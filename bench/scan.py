"""Static scanner: flag SAS constructs that are known to break naive SAS-to-R translation.

Usage:
    python3 bench/scan.py FILE.sas [FILE.sas ...]      # print findings as JSON
    python3 bench/scan.py --summary FILE.sas ...        # counts per construct

Each construct is a regex over comment-stripped code, applied per DATA/PROC step where the rule
needs step context (for example MERGE without BY). Line numbers refer to the original file.
"""

import json
import re
import sys

# (id, description, pattern, scope) - scope "line" matches anywhere; "step" gets a whole step.
RULES = [
    ("retain", "RETAIN: values carry across rows; R has no implicit row loop",
     r"\bretain\b", "line"),
    ("first_last", "first./last. BY-group flags depend on sort order",
     r"\b(first|last)\.\w+", "line"),
    ("lag_dif", "LAG/DIF are queue-based, not 'previous row' when called conditionally",
     r"\b(lag|dif)\d*\s*\(", "line"),
    ("double_set", "more than one SET in a step, or FIRSTOBS look-ahead",
     None, "step"),
    ("merge_no_by", "MERGE without BY joins by row position",
     None, "step"),
    ("merge_in_flags", "MERGE with IN= flags (left/inner join semantics written by hand)",
     r"\bin\s*=\s*\w+", "line"),
    ("nodupkey", "PROC SORT NODUPKEY/NODUP keeps the first row in the incoming order",
     r"\bnodup(key)?\b", "line"),
    ("missing_compare", "comparisons with missing values (. or .z): SAS missing sorts below every number",
     r"(>\s*\.z?\b|<\s*\.z?\b|=\s*\.\s*(and|or|then|;)|\.z\b)", "line"),
    ("chained_compare", "chained comparisons (a >= b > .z): SAS means a >= b AND b > .z; R evaluates left to right",
     r"[\w.]+\s*(<=|>=|<|>|\^=)\s*[\w.]+\s*(<=|>=|<|>)\s*[\w.]+|\b(lt|le|gt|ge)\b\s*[\w.]+\s*\b(lt|le|gt|ge)\b", "line"),
    ("silent_informat", "?? / ? informat modifiers suppress errors and return missing",
     r"\?\?|,\s*\?\s*\w+\.", "line"),
    ("informat_format", "INPUT/PUT with formats (rounding, padding, date parsing rules)",
     r"\b(input|put)\s*\(", "line"),
    ("date_handling", "SAS dates: day counts since 1960, date literals, INTCK/INTNX/MDY",
     r"'\d{1,2}[a-z]{3}\d{2,4}'d|\b(intck|intnx|mdy|today|datepart|yymmdd\d*|date9)\b", "line"),
    ("proc_format_cntlin", "PROC FORMAT CNTLIN builds formats from data",
     r"\bcntlin\s*=", "line"),
    ("macro", "macro language: %macro, &var, CALL SYMPUT, CALL EXECUTE",
     r"%macro\b|&\w+|call\s+symputx?\b|call\s+execute\b", "line"),
    ("ods_output", "ODS OUTPUT: results captured from a procedure's printed tables",
     r"\bods\s+output\b", "line"),
    ("exact_tests", "exact tests (Fisher) in PROC FREQ: check two-sided definition and table layout",
     r"\b(fisher|exact)\b", "line"),
    ("sql_remerge", "PROC SQL summary functions next to non-grouped columns (SAS remerges silently)",
     r"\b(count|sum|avg|mean|min|max)\s*\(", "line"),
    ("string_semantics", "SAS strings: 1-based INDEX/SUBSTR, blank-padded comparisons, fixed LENGTH",
     r"\b(index|substr|compress|strip|trim|left|length)\s*\(|\blength\s+\w+\s+\$", "line"),
    ("proc_transpose", "PROC TRANSPOSE: column naming from ID values and prefixes",
     r"\bproc\s+transpose\b", "line"),
]


def strip_comments(text):
    """Blank out /* */ and *...; comments, keeping newlines so line numbers stay correct."""
    keep_lines = lambda m: re.sub(r"[^\n]", " ", m.group(0))
    text = re.sub(r"/\*.*?\*/", keep_lines, text, flags=re.S)
    # statement comments: '*' at the start of a statement, up to the next ';'
    text = re.sub(r"(^|;)(\s*)\*[^;]*;", lambda m: m.group(1) + m.group(2) + keep_lines(m)[len(m.group(1) + m.group(2)):], text, flags=re.M)
    return text


def steps(text):
    """Split into (start_line, step_text) chunks at DATA/PROC boundaries."""
    out, start, buf = [], 1, []
    for i, line in enumerate(text.split("\n"), 1):
        if re.match(r"\s*(data|proc)\b", line, re.I) and buf:
            out.append((start, "\n".join(buf)))
            start, buf = i, []
        buf.append(line)
    out.append((start, "\n".join(buf)))
    return out


def line_of(step_start, step_text, pos):
    return step_start + step_text.count("\n", 0, pos)


def scan(path):
    raw = open(path, encoding="utf-8", errors="replace").read()
    code = strip_comments(raw)
    found = {}
    for rid, _, pattern, scope in RULES:
        if scope != "line":
            continue
        rx = re.compile(pattern, re.I)
        hits = [i for i, line in enumerate(code.split("\n"), 1) if rx.search(line)]
        if hits:
            found[rid] = hits
    for start, step in steps(code):
        if not re.match(r"\s*data\b", step, re.I):
            continue
        sets = [m.start() for m in re.finditer(r"\bset\b", step, re.I)]
        if len(sets) > 1 or re.search(r"\bfirstobs\s*=", step, re.I):
            found.setdefault("double_set", []).append(line_of(start, step, sets[-1] if sets else 0))
        m = re.search(r"\bmerge\b", step, re.I)
        if m and not re.search(r"\bby\b", step, re.I):
            found.setdefault("merge_no_by", []).append(line_of(start, step, m.start()))
    return {"file": path, "lines": raw.count("\n") + 1, "constructs": found}


def main(argv):
    summary = "--summary" in argv
    files = [a for a in argv if not a.startswith("--")]
    results = [scan(f) for f in files]
    if summary:
        counts = {rid: sum(1 for r in results if rid in r["constructs"]) for rid, *_ in RULES}
        print(json.dumps({"programs": len(results), "programs_with_construct": counts}, indent=2))
    else:
        print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main(sys.argv[1:])
