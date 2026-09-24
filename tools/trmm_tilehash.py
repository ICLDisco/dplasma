#!/usr/bin/env python3
"""Audit the TILEHASH records ztrmm_RLT.jdf writes into parsec's debug history.

Every task hashes each tile it is handed, naming the tile of the original
matrix the copy came from, so the dump says what each task *saw* rather than
what the dataflow promised it. The checker adds the value the data collection
holds before and after each trmm, hashed the same way.

Within one iteration, three things must hold:

  A is read-only        every task that touches descA(i,j) must see one value.
  no write-after-read   a task reading descB(m,k) must see what read_B(m,k)
                        published, because the ztrmm that overwrites that tile
                        is supposed to be held back by ctl0.
  the C chain           whatever a task wrote into descB(m,c) must be what the
                        next task in the accumulation chain picks up, and the
                        last writer must be what the data collection ends up
                        holding.

Across iterations the input is rebuilt identically every time, so a fourth and
much sharper check applies: every task must see and produce the same bytes in
every iteration. The first task in time that does not is the culprit, and
whether its inputs or only its output moved says whether it was fed bad data
or produced it.

Usage: trmm_tilehash.py <dump> [...]
"""

import re
import sys
from collections import defaultdict

RECORD = re.compile(
    r"^ 0x[0-9a-f]+/\d+ \(\s*([0-9.e+-]+) s\) -- TILEHASH r(\d+) "
    r"(descA|descB)\((\d+),(\d+)\) (\S+) ([0-9a-f]{16})(?: @(0x[0-9a-f]+))?"
)
# The checker is precision-generated, so the marker names ztrmm/dtrmm/...
MARKER = re.compile(r"-- === \wtrmm iteration (\d+) on rank (\d+) begins")
TASK = re.compile(r"(read_A|read_B|trmm|gemm)\(([\d, ]+)\)(?:\.(\w+))?")
# The checker's own records name no task, just the stage they were taken at.
STAGE = ("lacpy", "final")


def load(paths):
    """Records grouped by rank then iteration, in the dump's own order.

    The dump is merge-sorted by timestamp within a rank, so every record
    following an iteration marker belongs to that iteration.
    """
    runs = defaultdict(lambda: defaultdict(list))
    current = {}
    for path in paths:
        for line in open(path):
            mark = MARKER.search(line)
            if mark:
                current[int(mark.group(2))] = int(mark.group(1))
                continue
            m = RECORD.match(line)
            if not m:
                continue
            ts, rank, mat, i, j, task, digest, ptr = m.groups()
            rank = int(rank)
            if task in STAGE:
                key = {"task": task, "locals": [], "role": task}
            else:
                t = TASK.match(task)
                if t is None:
                    continue
                key = {
                    "task": t.group(1),
                    "locals": [int(x) for x in t.group(2).split(",")],
                    "role": t.group(3) or "out",
                }
            key.update(time=float(ts), tile=(mat, int(i), int(j)), hash=digest,
                       ptr=ptr)
            runs[rank][current.get(rank, -1)].append(key)
    return {r: dict(sorted(v.items())) for r, v in sorted(runs.items())}


def describe(r):
    if r["task"] in STAGE:
        return r["task"]
    return "%s(%s).%s" % (r["task"], ",".join(str(x) for x in r["locals"]), r["role"])


def key_of(r):
    return (r["task"], tuple(r["locals"]), r["role"], r["tile"])


def audit_A(records, report):
    seen = defaultdict(list)
    for r in records:
        if r["tile"][0] == "descA":
            seen[r["tile"]].append(r)
    for tile, rs in sorted(seen.items()):
        digests = {r["hash"] for r in rs}
        if len(digests) == 1:
            continue
        majority = max(digests, key=lambda d: sum(r["hash"] == d for r in rs))
        report("descA(%d,%d) read %d times with %d different values (%s is the "
               "majority)" % (tile[1], tile[2], len(rs), len(digests), majority))
        for r in rs:
            if r["hash"] != majority:
                report("    %8.4fs  %-22s saw %s" % (r["time"], describe(r), r["hash"]))


def audit_B_reads(records, report):
    """Every read of descB(m,k) as an operand must match what read_B published."""
    published = {r["tile"]: r["hash"] for r in records if r["task"] == "read_B"}
    for r in records:
        if r["tile"][0] != "descB":
            continue
        if not ((r["task"] == "trmm" and r["role"] == "B")
                or (r["task"] == "gemm" and r["role"] == "A")):
            continue
        want = published.get(r["tile"])
        if want is not None and want != r["hash"]:
            report("descB(%d,%d) was %s at read_B, but %s saw %s at %.4fs"
                   % (r["tile"][1], r["tile"][2], want, describe(r), r["hash"],
                      r["time"]))


def audit_C_chain(records, report):
    """A tile handed on must arrive as it was left, and land in the collection."""
    last = {}
    for r in records:
        if r["tile"][0] != "descB":
            continue
        if r["role"] == "out":
            last[r["tile"]] = r
        elif (r["task"] == "gemm" and r["role"] == "C") or r["role"] == "final":
            prev = last.get(r["tile"])
            if prev is not None and prev["hash"] != r["hash"]:
                report("descB(%d,%d) left %s as %s but arrived at %s as %s (%.4fs)"
                       % (r["tile"][1], r["tile"][2], describe(prev), prev["hash"],
                          describe(r), r["hash"], r["time"]))


def audit_aliasing(records, report):
    """A buffer must not be two tiles at once.

    Every record names the tile its task believes it holds; the address says
    which buffer that actually is. If one address carries two tile names with
    overlapping lifetimes, a task is writing over somebody else's data while
    both of them report perfectly consistent hashes.
    """
    span = defaultdict(lambda: defaultdict(list))
    for r in records:
        if r["ptr"] is not None:
            span[r["ptr"]][r["tile"]].append(r["time"])
    for ptr, tiles in sorted(span.items()):
        if len(tiles) < 2:
            continue
        ranges = sorted((min(ts), max(ts), tile) for tile, ts in tiles.items())
        for (_, prev_end, prev), (nxt_start, _, nxt) in zip(ranges, ranges[1:]):
            # Timestamps land on a millisecond grid, so two uses that merely
            # touch cannot be told from an arena buffer legitimately recycled
            # after its last reader. Only a strict overlap means anything.
            if nxt_start < prev_end:
                report("%s is %s(%d,%d) until %.4fs and %s(%d,%d) from %.4fs"
                       % (ptr, prev[0], prev[1], prev[2], prev_end,
                          nxt[0], nxt[1], nxt[2], nxt_start))
                break


def audit_iterations(iterations, report):
    """Every task must see and produce the same bytes in every iteration.

    The baseline is what most iterations agreed on, not iteration 0. Taking
    the first iteration as truth makes anything it does differently -- BLAS
    dispatching on its first call, say -- look like nine failures instead of
    one oddity, and buries the real thing underneath.
    """
    if len(iterations) < 2:
        return

    readings = defaultdict(dict)
    for it, records in iterations.items():
        for r in records:
            readings[key_of(r)][it] = r

    diverged = []
    for key, per_it in readings.items():
        tally = defaultdict(int)
        for r in per_it.values():
            tally[r["hash"]] += 1
        if len(tally) == 1:
            continue
        agreed = max(tally, key=lambda h: tally[h])
        want = next(r for r in per_it.values() if r["hash"] == agreed)
        for it, r in per_it.items():
            if r["hash"] != agreed:
                diverged.append((it, r, want, tally[r["hash"]], len(per_it)))

    if not diverged:
        return

    by_iteration = defaultdict(list)
    for it, r, want, odd, total in diverged:
        by_iteration[it].append((r, want, odd, total))
    for items in by_iteration.values():
        items.sort(key=lambda d: d[0]["time"])

    for it, items in sorted(by_iteration.items()):
        first, want, odd, total = items[0]
        report("iteration %d: %d readings disagree with what the other "
               "iterations agreed on, earliest at %.4fs"
               % (it, len(items), first["time"]))
        report("    %-22s %s(%d,%d) %s, %d of %d iterations say %s"
               % (describe(first), first["tile"][0], first["tile"][1],
                  first["tile"][2], first["hash"], total - odd, total,
                  want["hash"]))
        if first["task"] in STAGE:
            continue
        siblings = [r for r in iterations[it]
                    if r["task"] == first["task"] and r["locals"] == first["locals"]]
        inputs = [r for r in siblings if r["role"] != "out"]
        bad_in = [r for r in inputs
                  if r["hash"] != readings[key_of(r)][
                      max(set(readings[key_of(r)]),
                          key=lambda i: sum(
                              readings[key_of(r)][i]["hash"] == readings[key_of(r)][j]["hash"]
                              for j in readings[key_of(r)]))]["hash"]]
        if bad_in:
            report("    it was fed bad data: " + ", ".join(
                "%s %s(%d,%d)" % (r["role"], r["tile"][0], r["tile"][1], r["tile"][2])
                for r in bad_in))
        elif inputs:
            report("    its %d inputs all agree with the other iterations, so "
                   "this task produced the divergence" % len(inputs))


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)

    runs = load(sys.argv[1:])
    if not runs:
        sys.exit("no TILEHASH records found; was parsec built with "
                 "-DPARSEC_DEBUG_HISTORY=ON?")

    # descA is read-only and replicated, so every rank must agree about every
    # tile of it, whichever iteration that rank happened to dump.
    seen = defaultdict(dict)
    for rank, iterations in runs.items():
        for records in iterations.values():
            for r in records:
                if r["tile"][0] == "descA":
                    seen[r["tile"]].setdefault(rank, set()).add(r["hash"])
    split = {t: p for t, p in seen.items() if len(set().union(*p.values())) > 1}
    print("=== descA across ranks: %d tiles, %d read inconsistently"
          % (len(seen), len(split)))
    for tile, per in sorted(split.items())[:20]:
        print("    descA(%d,%d): %s" % (tile[1], tile[2], ", ".join(
            "r%d=%s" % (r, "/".join(sorted(h))) for r, h in sorted(per.items()))))
    print()

    total, comparable = 0, 0
    total += len(split)
    for rank, iterations in runs.items():
        its = sorted(iterations)
        print("=== rank %d: iterations %s, %d records"
              % (rank, its if its != [-1] else "unmarked",
                 sum(len(v) for v in iterations.values())))

        found = []
        for name, audit in (("A is read-only", audit_A),
                            ("no write-after-read on B", audit_B_reads),
                            ("C chain and write-back", audit_C_chain),
                            ("no two tiles share a buffer", audit_aliasing)):
            before = len(found)
            for records in iterations.values():
                audit(records, found.append)
            print("  %-26s %s" % (name, "ok" if len(found) == before
                                  else "%d violations" % (len(found) - before)))

        before = len(found)
        audit_iterations(iterations, found.append)
        if len(iterations) < 2:
            print("  %-26s NOT CHECKED - only one iteration in this dump"
                  % "same bytes every iteration")
        else:
            comparable += 1
            print("  %-26s %s" % ("same bytes every iteration",
                                  "ok" if len(found) == before
                                  else "%d divergences" % (len(found) - before)))
        for line in found:
            print("    " + line)
        total += len(found)

    if not comparable:
        print("\nNo rank's dump held more than one iteration, so the strongest "
              "check never ran.\nA single iteration audits clean by "
              "construction when it is the one that behaved:\nrun with "
              "DPLASMA_CHECK_TRACE=1 so every iteration is kept.")
        return 2

    print("\n%d violations total" % total)
    return 1 if total else 0


if __name__ == "__main__":
    sys.exit(main())
