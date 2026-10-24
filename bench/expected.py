"""List the datasets a SAS program creates, read mechanically from its code, for the structure check.

    python3 bench/expected.py corpus/phuse/X.sas      # prints JSON: {"datasets": [...], "skipped": [...]}

Counted: DATA step outputs, OUT= on any PROC, PROC SQL CREATE TABLE, ODS OUTPUT targets. Names
built from macro variables (containing &) can't be resolved without running SAS, so they are
listed under "skipped" and not required. Comments are removed first (same rules as scan.py).
"""

import json
import re
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from scan import strip_comments  # noqa: E402

NAME = r"(?:\w+\.)?([A-Za-z_&][\w&.]*)"


def clean(n):
    n = n.split(".")[-1] if "." in n and "&" not in n else n
    return n.upper()


def expected(text):
    code = strip_comments(text)
    found, skipped = [], []
    def add(n):
        n = n.strip()
        if not re.fullmatch(r"(?:\w+\.)?[A-Za-z_&][\w&.]*", n) or n.upper() in ("_NULL_", "_LAST_", "_DATA_"):
            return
        (skipped if "&" in n else found).append(clean(n))
    for m in re.finditer(r"(?im)^\s*data\s+([^;]+);", code):
        spec = re.sub(r"\([^)]*\)", " ", m.group(1))       # drop dataset options
        for tok in spec.split():
            if "=" in tok:
                break
            add(tok)
    for m in re.finditer(r"(?i)\bout\s*=\s*" + NAME, code):
        add(m.group(1))
    for m in re.finditer(r"(?i)\bcreate\s+table\s+" + NAME, code):
        add(m.group(1))
    for m in re.finditer(r"(?i)\bods\s+output\s+([^;]+);", code):
        for pair in re.findall(r"\w+\s*=\s*" + NAME, m.group(1)):
            add(pair)
    uniq = lambda xs: sorted(set(xs))
    return {"datasets": uniq(found), "skipped": uniq(skipped)}


if __name__ == "__main__":
    print(json.dumps(expected(open(sys.argv[1], encoding="utf-8", errors="replace").read()), indent=2))
