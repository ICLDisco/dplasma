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
    r"^ (0x[0-9a-f]+)/\d+ \(\s*([0-9.e+-]+) s\) -- TILEHASH r(\d+) "
    r"(descA|descB)\((\d+),(\d+)\) (\S+) ([0-9a-f]{16})(?: @(0x[0-9a-f]+))?"
)
# The checker is precision-generated, so the marker names ztrmm/dtrmm/...
MARKER = re.compile(r"-- === \wtrmm iteration (\d+) on rank (\d+) begins")
TASK = re.compile(r"(read_A|read_B|trmm|gemm)\(([\d, ]+)\)(?:\.([\w-]+))?")
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
        for lineno, line in enumerate(open(path)):
            mark = MARKER.search(line)
            if mark:
                current[int(mark.group(2))] = int(mark.group(1))
                continue
            m = RECORD.match(line)
            if not m:
                continue
            thread, ts, rank, mat, i, j, task, digest, ptr = m.groups()
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
                       ptr=ptr, seq=lineno, thread=thread)
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
                or (r["task"] == "gemm" and r["role"].startswith("A"))):
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
            span[r["ptr"]][r["tile"]].append(r["seq"])
    for ptr, tiles in sorted(span.items()):
        if len(tiles) < 2:
            continue
        # The printed timestamp carries three significant digits, so hundreds
        # of records share one value and two spans that merely touch say
        # nothing. Position in the dump is a total order over the rank -- it
        # is merge-sorted across that rank's threads -- so use that instead,
        # and require the two spans to genuinely interleave rather than abut.
        ranges = sorted((min(s), max(s), tile) for tile, s in tiles.items())
        for (_, prev_end, prev), (nxt_start, nxt_end, nxt) in zip(ranges, ranges[1:]):
            if nxt_start < prev_end:
                report("%s carries %s(%d,%d) over records %d-%d and %s(%d,%d) "
                       "over %d-%d; the two interleave"
                       % (ptr, prev[0], prev[1], prev[2],
                          min(tiles[prev]), prev_end,
                          nxt[0], nxt[1], nxt[2], nxt_start, nxt_end))
                break


def audit_kernel_window(records, report):
    """An operand must not move while the kernel is reading it.

    Each input is hashed on the way into the body and again once the kernel
    returns. The two readings bracket the only interval the other checks
    cannot see, so a difference here says the tile was written during the
    call rather than before or after it.
    """
    before = {}
    for r in records:
        if r["task"] in STAGE:
            continue
        ident = (r["task"], tuple(r["locals"]), r["role"].split("-")[0])
        if r["role"].endswith("-after"):
            was = before.get(ident)
            if was is not None and was["hash"] != r["hash"]:
                report("%s %s(%d,%d) was %s when %s started and %s when the "
                       "kernel returned"
                       % (ident[2], r["tile"][0], r["tile"][1], r["tile"][2],
                          was["hash"], describe(was).rsplit(".", 1)[0], r["hash"]))
        else:
            before[ident] = r


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
        inputs = [r for r in siblings
                  if r["role"] != "out" and not r["role"].endswith("-after")]
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
            for line in intruders(iterations[it], siblings, first):
                report("    " + line)
            for line in provenance(readings, iterations, it, siblings, first):
                report("    " + line)


def provenance(readings, iterations, it, siblings, out):
    """Where else in the run does the wrong value appear?

    Bytes that land in a buffer came from somewhere. If the wrong output is
    a value this run computed elsewhere, that names what overwrote it far
    more directly than any timing argument: the same tile at an earlier point
    in its own chain means a stale copy arrived, another tile's value means
    that tile's buffer was written here, and its own C input means the
    kernel's write never took.
    """
    lines = []
    for r in siblings:
        if r["role"] != "out" and r["hash"] == out["hash"]:
            lines.append("the wrong value is its own %s input, so the kernel's "
                         "write did not take" % r["role"])
    matches = []
    for key, per_it in readings.items():
        if key == key_of(out):
            continue
        for other_it, r in per_it.items():
            if r["hash"] == out["hash"]:
                matches.append((other_it, r))
    if matches:
        same_tile = [m for m in matches if m[1]["tile"] == out["tile"]]
        pick = same_tile or matches
        lines.append("the wrong value also appears as %d other reading(s), e.g. %s"
                     % (len(matches), ", ".join(
                         "%s %s(%d,%d) in iteration %d"
                         % (describe(r), r["tile"][0], r["tile"][1], r["tile"][2], i)
                         for i, r in sorted(pick, key=lambda m: m[1]["seq"])[:3])))
    elif not lines:
        lines.append("the wrong value appears nowhere else in the run, so it "
                     "was computed, not copied in")
    return lines


def intruders(records, siblings, out):
    """Who else held this task's output buffer while the task was running.

    A deterministic kernel handed correct inputs cannot produce a wrong
    answer on its own. If the buffer it wrote was clobbered between the hash
    taken on the way in and the one taken on the way out, everything
    downstream still agrees -- so the only way to name the culprit is to ask
    which other task was holding the same address at the time.
    """
    if out["ptr"] is None or not siblings:
        return []
    start, end = min(r["seq"] for r in siblings), out["seq"]
    mine = (out["task"], tuple(out["locals"]))
    others = [r for r in records
              if r["ptr"] == out["ptr"]
              and (r["task"], tuple(r["locals"])) != mine
              and start < r["seq"] < end]
    # How much of the task's run the window actually covers: the dump groups
    # records by thread when timestamps tie, so a window holding nothing from
    # any other thread has not observed the concurrency it claims to rule out.
    window = [r for r in records if start < r["seq"] < end]
    elsewhere = {r["thread"] for r in window} - {out["thread"]}
    if not others:
        return ["nothing else touched %s over records %d-%d (%d records there, "
                "%d from other threads)"
                % (out["ptr"], start, end, len(window), 
                   sum(r["thread"] in elsewhere for r in window))]
    lines = ["%s was also held by %d other records while this task ran:"
             % (out["ptr"], len(others))]
    for r in sorted(others, key=lambda r: r["seq"])[:8]:
        lines.append("      record %-7d %-22s %s(%d,%d) %s"
                     % (r["seq"], describe(r), r["tile"][0], r["tile"][1],
                        r["tile"][2], r["hash"]))
    return lines


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
                            ("no two tiles share a buffer", audit_aliasing),
                            ("operands hold still", audit_kernel_window)):
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
