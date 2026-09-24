#!/usr/bin/env python3
#
# Copyright (c) 2026      NVIDIA Corporation.  All rights reserved.
#
"""Find tiles that enter a JDF's dataflow through more than one task.

A flow that reads a data collection directly, '<- ddescA(i,j)', is where a
tile enters the DAG. Two tasks doing that for the same tile each get their
own parsec_data_get_copy of the same device-0 copy, with no dependency
between them: the runtime has no way to know the two are related. While both
only read, that costs nothing. As soon as one of them holds the tile RW and
overwrites it in place, the ordering rests entirely on whatever control edge
the author remembered to add, and a control edge is not checked against the
flows by anything.

The same goes the other way: two tasks writing the same tile back to the
collection, '-> ddescA(i,j)', are two last writers for one tile.

Two entry points naming different tiles are not a problem, and nearly all
of them do: a factorization typically reads A(k,k), A(m,k) and A(k,n) out
of the collection at k == 0, and those never coincide. So the execution
space of each task is enumerated over a small concrete matrix and the index
expressions evaluated, and only entry points that genuinely land on the
same tile are reported. Anything whose space or condition cannot be
evaluated is listed separately rather than quietly dropped.

Usage: jdf_entrypoints.py [-v] [file.jdf ...]   (default: all of src/*.jdf)
"""

import glob
import os
import re
import sys
from collections import defaultdict

# A global declared as a data collection, e.g.
#   ddescA     [type = "dplasma_data_collection_t*"]
GLOBAL = re.compile(r'^(\w+)\s*\[\s*type\s*=\s*"?([^"\]]+)')
COLLECTION_TYPES = ("data_collection_t*", "tiled_matrix_t*", "matrix_t*")

# A task header sits at the left margin: name(params) possibly then '['.
TASK = re.compile(r"^(\w+)\s*\(([^)]*)\)\s*(\[|$)")
# A flow opens with its access mode; the deps may continue on later lines.
FLOW = re.compile(r"^\s*(READ|RW|WRITE|CTL)\s+(\w+)\s+(.*)$")
DEP = re.compile(r"^\s*(<-|->)\s*(.*)$")
# 'cond ? name(args)' is the collection; 'cond ? flow name(args)' is a task.
TARGET = re.compile(r"^(?:(.*?)\?\s*)?(\w+)\s*\(")


def arguments(text, open_paren):
    """The argument list, counting nested parens rather than stopping at one."""
    depth, i = 0, open_paren
    while i < len(text):
        if text[i] == "(":
            depth += 1
        elif text[i] == ")":
            depth -= 1
            if depth == 0:
                return text[open_paren + 1:i]
        i += 1
    return text[open_paren + 1:]


# A small matrix, and wide enough that a range like k+1 .. nt-1 is not empty.
SHAPE = {"mt": 5, "nt": 5, "lmt": 5, "lnt": 5, "lm": 20, "ln": 20,
         "mb": 4, "nb": 4, "m": 20, "n": 20, "i": 0, "j": 0, "bsiz": 16,
         "dtype": 0, "storage": 0, "rtile": 0, "llm": 20, "lln": 20}
ASSIGN = re.compile(r"^\s*(\w+)\s*=\s*(.+?)\s*$")
DEFAULT = re.compile(r'\bdefault\s*=\s*"([^"]*)"')
RANGE = re.compile(r"^(.*?)\s*\.\.\s*(.*)$")


def pythonise(expr):
    """A JDF integer expression as something Python can evaluate."""
    # An execution space often bounds itself with a one-line inline_c.
    expr = re.sub(r"%\{\s*return\s+(.*?)\s*;\s*%\}", r"(\1)", expr)
    expr = re.sub(r"\b\w+\s*->\s*", "", expr)          # descA->mt  ->  mt
    expr = expr.replace("&&", " and ").replace("||", " or ")
    expr = re.sub(r"!(?!=)", " not ", expr)
    expr = re.sub(r"/(?!/)", "//", expr)               # integer division
    # 'c ? a : b' is right-associative and nests, so rewrite innermost first.
    while True:
        m = re.search(r"\(([^()?:]+)\?([^()?:]+):([^()?:]+)\)", expr)
        if not m:
            return expr
        expr = expr[:m.start()] + "((%s) if (%s) else (%s))" % (
            m.group(2), m.group(1), m.group(3)) + expr[m.end():]


HELPERS = {"dplasma_imin": min, "dplasma_imax": max}


def value(expr, env):
    return eval(pythonise(expr), {"__builtins__": {}}, dict(env, **HELPERS))


def tiles(task, index, when, shape):
    """Every (i,j) a dependency can name, over the task's execution space.

    Returns None when the space, the condition or the index uses something
    that cannot be evaluated here -- an inline_c local, say -- so that an
    undecidable case is never mistaken for a safe one.
    """
    params = task["params"]
    if not params:
        return None
    found = set()

    def walk(left, env):
        if not left:
            if when and not value(when, env):
                return
            found.add(tuple(value(p, env) for p in index.split(",")))
            return
        # A range may be written in terms of another parameter, in whichever
        # order reads best, so take whichever one can be evaluated now.
        for name in left:
            if name not in task["ranges"]:
                raise ValueError("no range for '%s'" % name)
            lo, hi = task["ranges"][name]
            try:
                span = range(value(lo, env), value(hi, env) + 1)
            except Exception:
                continue
            rest = [p for p in left if p != name]
            for v in span:
                env[name] = v
                for local, rhs in task["locals"].items():
                    try:
                        env[local] = value(rhs, env)
                    except Exception:
                        pass
                walk(rest, env)
            return
        raise ValueError("circular execution space %s" % left)

    try:
        walk(list(params), dict(shape))
    except Exception:
        return None
    return found


def collections(lines):
    """The data collections, and the scalar globals with an evaluable default.

    Execution spaces lean on globals such as KT or minMN, declared hidden
    with a default expression; without them most spaces cannot be walked.
    """
    names, scalars = set(), {}
    for line in lines:
        m = GLOBAL.match(line)
        if not m:
            continue
        if any(t in m.group(2) for t in COLLECTION_TYPES):
            names.add(m.group(1))
            continue
        default = DEFAULT.search(line)
        if default and "int" in m.group(2):
            scalars[m.group(1)] = default.group(1)

    env = dict(SHAPE)
    for _ in range(len(scalars)):          # they may refer to one another
        pending = False
        for name, expr in scalars.items():
            if name in env:
                continue
            try:
                env[name] = value(expr, env)
            except Exception:
                pending = True
        if not pending:
            break
    return names, env


def parse(path):
    """Every direct reference to a data collection, by task and flow."""
    lines = open(path, errors="replace").read().splitlines()
    dcs, shape = collections(lines)
    refs = []

    tasks = {}
    task, mode, flow = None, None, None
    in_body = in_prologue = False
    for line in lines:
        stripped = line.strip()

        if in_prologue:
            if stripped.startswith("%}"):
                in_prologue = False
            continue
        if stripped.startswith('extern "C"') or stripped == "%{":
            in_prologue = True
            continue
        if in_body:
            if stripped == "END":
                in_body = False
            continue
        if stripped.startswith("BODY"):
            in_body = True
            continue

        head = TASK.match(line)
        if head:
            task, mode, flow = head.group(1), None, None
            tasks[task] = {
                "params": [p.strip() for p in head.group(2).split(",")
                           if p.strip()],
                "ranges": {}, "locals": {},
            }
            continue

        # Between the header and the first flow sit the execution space and
        # the locals; an inline_c local cannot be evaluated, so it is skipped
        # and any dependency needing it comes back undecidable.
        if (task and mode is None and "<-" not in line
                and line.count("%{") == line.count("%}")):
            assign = ASSIGN.match(line.split("/*")[0].split("//")[0])
            if assign and not assign.group(1) in ("loc", ):
                name, rhs = assign.group(1), assign.group(2)
                span = RANGE.match(rhs)
                if span and name in tasks[task]["params"]:
                    tasks[task]["ranges"][name] = (span.group(1), span.group(2))
                elif not span:
                    tasks[task]["locals"][name] = rhs
                continue

        opening = FLOW.match(line)
        if opening:
            mode, flow = opening.group(1), opening.group(2)
            rest = opening.group(3)
        else:
            rest = line

        # Annotations sit in brackets and mention the collection only as an
        # argument to a helper, never as 'name(' in dependency position.
        dep = DEP.match(rest.split("[")[0])
        if not dep or task is None or flow is None:
            continue
        body = dep.group(2).strip()
        target = TARGET.match(body)
        if target and target.group(2) in dcs:
            refs.append({
                "task": task, "flow": flow, "mode": mode,
                "dir": dep.group(1), "dc": target.group(2),
                "when": re.sub(r"\s+", " ", (target.group(1) or "").strip()),
                "index": re.sub(r"\s+", "", arguments(body, target.end() - 1)),
            })
    for r in refs:
        r["tiles"] = tiles(tasks[r["task"]], r["index"], r["when"], shape)
    return refs


def describe(r):
    return ("%-5s %-3s %-26s in %-20s %s"
            % (r["mode"], r["flow"], "%s(%s)" % (r["dc"], r["index"]),
               r["task"], "when " + r["when"] if r["when"] else "always"))


def report(path, refs, verbose):
    clashes, unknown = [], []
    for direction, word in (("<-", "enter"), ("->", "leave")):
        points = {}
        for r in refs:
            if r["dir"] == direction:
                # Several conditional deps on one flow are one entry point,
                # and between them they cover every tile it can name.
                key = (r["dc"], r["task"], r["flow"])
                if key in points and None not in (points[key]["tiles"],
                                                  r["tiles"]):
                    points[key]["tiles"] = points[key]["tiles"] | r["tiles"]
                else:
                    points.setdefault(key, dict(r))

        ordered = sorted(points.values(), key=lambda r: (r["dc"], r["task"]))
        for i, a in enumerate(ordered):
            for b in ordered[i+1:]:
                if a["dc"] != b["dc"]:
                    continue
                if a["tiles"] is None or b["tiles"] is None:
                    unknown.append((word, a, b))
                    continue
                shared = a["tiles"] & b["tiles"]
                if not shared:
                    continue
                # Two readers of a tile nobody writes is not a hazard.
                if word == "enter" and not any(
                        r["mode"] in ("RW", "WRITE") for r in (a, b)):
                    continue
                clashes.append((word, a, b, shared))

    if clashes or (verbose and unknown):
        print("\n%s" % os.path.basename(path))
    for word, a, b, shared in clashes:
        print("  !! %s: the same tile can %s the dataflow twice, e.g. (%s)"
              % (a["dc"], word, ",".join(str(x) for x in sorted(shared)[0])))
        print("       %s" % describe(a))
        print("       %s" % describe(b))
    if verbose:
        for word, a, b in unknown:
            print("  ?? %s: could not decide whether these %s together"
                  % (a["dc"], word))
            print("       %s" % describe(a))
            print("       %s" % describe(b))
    return len(clashes), len(unknown)


def main():
    args = sys.argv[1:]
    verbose = "-v" in args
    paths = [a for a in args if a != "-v"]
    if not paths:
        here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        paths = sorted(glob.glob(os.path.join(here, "src", "*.jdf")))

    bad, undecided, files, refs, blind = 0, 0, 0, 0, 0
    for path in paths:
        parsed = parse(path)
        refs += len(parsed)
        blind += sum(1 for r in parsed if r["tiles"] is None)
        n, u = report(path, parsed, verbose)
        bad, undecided, files = bad + n, undecided + u, files + (1 if n else 0)

    print("\n%d of %d files have a tile reachable from two places (%d in all)."
          % (files, len(paths), bad))
    print("%d of %d collection references were resolved; %d pair(s) rest on "
          "one that was not%s." % (refs - blind, refs, undecided,
                                   "" if verbose else ", and -v lists them"))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
