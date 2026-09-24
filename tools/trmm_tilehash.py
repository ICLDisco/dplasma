#!/usr/bin/env python3
"""Audit the TILEHASH records ztrmm_RLT.jdf writes into parsec's debug history.

Every task hashes each tile it is handed, naming the tile of the original
matrix the copy came from, so the dump says what each task *saw* rather than
what the dataflow promised it. Three things must hold, and each failure points
somewhere different:

  A is read-only        every task that touches descA(i,j) must see one value.
                        A disagreement means a tile arrived wrong, was released
                        early, or was overwritten in a buffer still in use.

  no write-after-read   a task reading descB(m,k) must see what read_B(m,k)
                        published, because the ztrmm that overwrites that tile
                        is supposed to be held back by ctl0.

  the C chain           whatever a task wrote into descB(m,c) must be what the
                        next task in the accumulation chain picks up.

Usage: trmm_tilehash.py <dump> [...]
"""

import re
import sys
from collections import defaultdict

RECORD = re.compile(
    r"^ 0x[0-9a-f]+/\d+ \(\s*([0-9.e+-]+) s\) -- TILEHASH r(\d+) "
    r"(descA|descB)\((\d+),(\d+)\) (\S+) ([0-9a-f]{16})"
)
TASK = re.compile(r"(read_A|read_B|trmm|gemm)\(([\d, ]+)\)(?:\.(\w+))?")


def load(paths):
    """Records per rank, in the dump's own (timestamp-merged) order.

    Ranks write to one interleaved stream, so every record carries its own.
    """
    ranks = defaultdict(list)
    for path in paths:
        with open(path) as fh:
            for line in fh:
                m = RECORD.match(line)
                if not m:
                    continue
                ts, rank, mat, i, j, task, digest = m.groups()
                t = TASK.match(task)
                ranks[int(rank)].append(
                    {
                        "time": float(ts),
                        "tile": (mat, int(i), int(j)),
                        "task": t.group(1),
                        "locals": [int(x) for x in t.group(2).split(",")],
                        "role": t.group(3) or "out",
                        "hash": digest,
                    }
                )
    return {k: v for k, v in sorted(ranks.items())}


def describe(r):
    return "%s(%s).%s" % (r["task"], ",".join(str(x) for x in r["locals"]), r["role"])


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
        report(
            "descA(%d,%d) read %d times with %d different values (%s is the majority)"
            % (tile[1], tile[2], len(rs), len(digests), majority)
        )
        for r in rs:
            if r["hash"] != majority:
                report("    %8.4fs  %-22s saw %s" % (r["time"], describe(r), r["hash"]))


def audit_B_reads(records, report):
    """Every read of descB(m,k) as an operand must match what read_B published."""
    published = {}
    for r in records:
        if r["task"] == "read_B":
            published[r["tile"]] = r["hash"]
    for r in records:
        if r["tile"][0] != "descB":
            continue
        # trmm's B and gemm's A are reads of the entry-point value; C and out
        # are the accumulator, which is supposed to change.
        if not (
            (r["task"] == "trmm" and r["role"] == "B")
            or (r["task"] == "gemm" and r["role"] == "A")
        ):
            continue
        want = published.get(r["tile"])
        if want is None or want == r["hash"]:
            continue
        report(
            "descB(%d,%d) was %s at read_B, but %s saw %s at %.4fs"
            % (r["tile"][1], r["tile"][2], want, describe(r), r["hash"], r["time"])
        )


def audit_C_chain(records, report):
    """A tile handed down the accumulation chain must arrive as it was left."""
    last = {}
    for r in records:
        if r["tile"][0] != "descB":
            continue
        if r["role"] == "out":
            last[r["tile"]] = r
        elif r["task"] == "gemm" and r["role"] == "C":
            prev = last.get(r["tile"])
            if prev is not None and prev["hash"] != r["hash"]:
                report(
                    "descB(%d,%d) left %s as %s but arrived at %s as %s (%.4fs)"
                    % (
                        r["tile"][1],
                        r["tile"][2],
                        describe(prev),
                        prev["hash"],
                        describe(r),
                        r["hash"],
                        r["time"],
                    )
                )


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)

    ranks = load(sys.argv[1:])
    if not ranks:
        sys.exit("no TILEHASH records found; was parsec built with "
                 "-DPARSEC_DEBUG_HISTORY=ON?")

    total = 0
    for rank, records in ranks.items():
        rows = sorted({r["locals"][0] for r in records if r["task"] == "read_B"})
        print("=== rank %d: %d records, rows %s" % (rank, len(records), rows))
        found = []
        for name, audit in (
            ("A is read-only", audit_A),
            ("no write-after-read on B", audit_B_reads),
            ("C accumulation chain", audit_C_chain),
        ):
            before = len(found)
            audit(records, found.append)
            print("  %-26s %s" % (name, "ok" if len(found) == before
                                  else "%d violations" % (len(found) - before)))
        for line in found:
            print("    " + line)
        total += len(found)

    print("\n%d violations total" % total)
    return 1 if total else 0


if __name__ == "__main__":
    sys.exit(main())
