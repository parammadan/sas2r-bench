"""Scanner golden file: what bench/scan.py finds on the PhUSE corpus programs and on one small case
per rule (plus comment and step edge cases). The R port in forward2028 (core/R/scan.R) is checked
against a copy of the same file, so CI fails in whichever repo changes its scanner alone.

    python3 -I bench/scanner_golden.py           check scan.py against tests/scanner_golden.json
    python3 -I bench/scanner_golden.py --write   rewrite it (only after a deliberate rule change;
                                                 then copy it to forward2028 tests/fixtures/)
"""
import hashlib, json, os, sys, tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "bench"))
from scan import scan  # noqa: E402

GOLDEN = os.path.join(ROOT, "tests", "scanner_golden.json")
PHUSE = ["ae_common", "ae_pref", "ae_serious", "demo_summary", "disposition_a", "mean_time",
         "medication_listing", "outliers_Shift_Table", "outliers_TEHiLo_Anytime", "pk_param_summary"]
CASES = [
    ("retain", "data a; set b;\n  retain n 0;\nrun;"),
    ("first_last", "data a; set b; by usubjid;\n  if first.usubjid then n = 0;\nrun;"),
    ("lag_dif", "data a; set b;\n  x = lag(y);\n  d = dif2(y);\nrun;"),
    ("double_set", "data a;\n  set b;\n  set c;\nrun;"),
    ("firstobs_lookahead", "data a;\n  merge b b(firstobs=2 keep=x rename=(x=nx));\nrun;"),
    ("merge_no_by", "data a;\n  merge b c;\nrun;"),
    ("merge_with_by", "data a;\n  merge b c;\n  by id;\nrun;"),
    ("merge_in_flags", "data a;\n  merge b(in=x) c(in=y);\n  by id;\n  if x;\nrun;"),
    ("nodupkey", "proc sort data=a out=b nodupkey;\n  by id;\nrun;"),
    ("missing_compare", "data a; set b;\n  if x > . then y = 1;\n  if z = . then w = 1;\n  if v < .z then u = 1;\nrun;"),
    ("decimal_not_missing", "data a; set b;\n  if x > .5 then y = 1;\nrun;"),
    ("chained_compare", "data a; set b;\n  if ASTDT>=TRTSDT>.z then TRTEMFL='Y';\n  if 1 lt x le 5 then k = 1;\nrun;"),
    ("silent_informat", "data a; set b;\n  x = input(c, ?? best.);\nrun;"),
    ("informat_format", "data a; set b;\n  c = put(x, 8.1);\nrun;"),
    ("date_handling", "data a; set b;\n  d = '01jan2014'd;\n  n = intck('day', s, e);\nrun;"),
    ("proc_format_cntlin", "proc format cntlin=fmt;\nrun;"),
    ("macro", "%macro m(v);\n  %put &v;\n%mend;\ndata _null_; call symputx('a', 1); run;"),
    ("ods_output", "ods output FishersExact = f1;\nproc freq data=a; tables t*e / fisher; run;"),
    ("sql_remerge", "proc sql;\n  select id, count(*) as n from a;\nquit;"),
    ("string_semantics", "data a; set b;\n  length name $ 20;\n  x = substr(y, 1, 3);\nrun;"),
    ("proc_transpose", "proc transpose data=a out=b prefix=v;\n  by id;\nrun;"),
    ("comments_ignored", "/* retain x; */\n* retain y;\nx = 1;"),
    ("comment_keeps_line_numbers", "/* a\n   b */\nretain z;"),
    ("nothing_risky", "x = 1;"),
]


def digest(path, algo):
    return hashlib.new(algo, open(path, "rb").read()).hexdigest()


def build():
    g = {"about": "bench/scan.py output; forward2028 core/R/scan.R must give the same", "programs": {}, "cases": []}
    for pid in PHUSE:
        p = os.path.join(ROOT, "corpus", "phuse", pid + ".sas")
        g["programs"][pid] = {"sha256": digest(p, "sha256"), "md5": digest(p, "md5"), "constructs": scan(p)["constructs"]}
    with tempfile.TemporaryDirectory() as d:
        for name, sas in CASES:
            p = os.path.join(d, name + ".sas")
            with open(p, "w") as f:
                f.write(sas)
            g["cases"].append({"name": name, "sas": sas, "constructs": scan(p)["constructs"]})
    return g


if __name__ == "__main__":
    g = build()
    if "--write" in sys.argv:
        with open(GOLDEN, "w") as f:
            json.dump(g, f, indent=1)
            f.write("\n")
        print(f"wrote {GOLDEN}")
    else:
        want = json.load(open(GOLDEN))
        bad = [k for k in want["programs"] if want["programs"][k] != g["programs"].get(k)]
        bad += [c["name"] for c, h in zip(want["cases"], g["cases"]) if c != h]
        bad += ["case count"] if len(want["cases"]) != len(g["cases"]) else []
        for b in bad:
            print("differs:", b)
        print(f"scanner golden: {len(g['programs'])} programs, {len(g['cases'])} cases, {len(bad)} differ")
        sys.exit(1 if bad else 0)
