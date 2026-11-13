#!/usr/bin/env python3
"""Copy an as-run script with some literal strings replaced (authored for the repository).

Several as-run scripts write into fixed folders of the original runs. To re-run them without touching the stored
results, run_paper.sh runs a copy in which only those folder strings point to a new folder under OUT. This tool makes
the copy, checks that each string occurs in the script, and checks that every changed line is a line holding one of
the replaced strings, so that nothing else in the code changes. The list of changed lines is written next to the copy
(<dst>.redirect.txt).

    python3 redirect.py <src> <dst> <old> <new> [<old> <new> ...]
"""
import sys

src, dst, pairs = sys.argv[1], sys.argv[2], sys.argv[3:]
if not pairs or len(pairs) % 2:
    sys.exit("usage: redirect.py <src> <dst> <old> <new> [<old> <new> ...]")
pairs = list(zip(pairs[0::2], pairs[1::2]))
text = open(src, encoding="utf-8").read()
for old, _ in pairs:
    if old not in text:
        sys.exit(f"redirect.py: {old!r} not found in {src}")
lines = text.split("\n")
out, log = [], []
for i, line in enumerate(lines, 1):
    new = line
    for old, rep in pairs:
        new = new.replace(old, rep)
    if new != line:
        if not any(old in line for old, _ in pairs):
            sys.exit(f"redirect.py: line {i} changed without a replaced string")
        log.append(f"{i}\t{line.strip()}\t->\t{new.strip()}")
    out.append(new)
assert len(out) == len(lines)
open(dst, "w", encoding="utf-8").write("\n".join(out))
open(dst + ".redirect.txt", "w", encoding="utf-8").write(f"source\t{src}\n" + "\n".join(log) + "\n")
print(f"{dst}: {len(log)} line(s) redirected")
