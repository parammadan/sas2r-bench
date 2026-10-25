#!/usr/bin/env python3
"""Evaluation Set v1: the recorded runs, references and scores the project reports.

    python3 -I eval/v1.py            verify every listed file against eval/v1_manifest.json
    python3 -I eval/v1.py --create   write the manifest once; refuses if it already exists

Exit status 1 if any listed file is missing or has changed. Files in ../sas2r-bench are checked
when that checkout is present (locally) and reported as skipped when it isn't (CI); their copies
under core/snapshot/sas2r are always checked. New files are never flagged: v1 is the listed set,
and new work adds a v2 set instead of editing this one. See docs/EVALUATION_SET_V1.md.
"""
import glob, hashlib, json, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
MANIFEST = os.path.join(HERE, "v1_manifest.json")

# What belongs to v1, per evaluation. Patterns are relative to the repo (or to ../sas2r-bench).
SETS = {
    "sap_tables": ["tlf/generated/*/meta.json", "tlf/generated/*/program.R", "tlf/generated/*/shell.md",
                   "tlf/refs/*.csv", "tlf/results/*.json", "tlf/results/*_cells.csv", "tlf/FINDINGS.md",
                   "core/snapshot/tables/*/meta.json", "core/snapshot/tables/*/program.R",
                   "core/snapshot/tables/*/shell.md", "core/snapshot/tables/*/validation.json"],
    "sas2r": ["core/snapshot/sas2r/corpus.json", "core/snapshot/sas2r/FINDINGS.md",
              "core/snapshot/sas2r/*/translation.R", "core/snapshot/sas2r/*/translation_meta.json",
              "core/snapshot/sas2r/*/result.json"],
    "helpdesk": ["helpdesk/questions.yaml", "helpdesk/recorded.json", "helpdesk/tickets.json",
                 "core/snapshot/helpdesk/kb.json", "core/snapshot/helpdesk/recorded.json",
                 "core/snapshot/helpdesk/tickets.json"],
    "amendment_impact": ["loa/results.json", "loa/graded.json", "loa/gold.yaml",
                         "core/snapshot/loa/results.json", "core/snapshot/loa/graded.json", "core/snapshot/loa/gold.yaml"],
    "data_checks": ["core/snapshot/dq.json"],
}
BENCH_SETS = {"sas2r": ["translations/**/*", "results/**/*"]}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def tracked(repo):
    out = subprocess.run(["git", "-C", repo, "ls-files"], capture_output=True, text=True, check=True).stdout
    return set(out.splitlines())


def collect(repo, patterns):
    files = tracked(repo)
    hits = set()
    for p in patterns:
        hits.update(os.path.relpath(f, repo) for f in glob.glob(os.path.join(repo, p), recursive=True) if os.path.isfile(f))
    return sorted(h for h in hits if h in files)   # committed files only, never local leftovers


def run_info(repo, path):
    """Model, prompt hash and date recorded inside a run's own metadata file, when it has them."""
    if not path.endswith(".json"):
        return None
    try:
        d = json.load(open(os.path.join(repo, path)))
    except (ValueError, OSError):
        return None
    if not isinstance(d, dict):
        return None
    info = {k: d[k] for k in ("model", "translator", "prompt_sha", "when") if k in d and isinstance(d[k], str)}
    return info or None


def create(bench):
    if os.path.exists(MANIFEST):
        sys.exit(f"{MANIFEST} already exists; v1 is frozen. Add a v2 set instead.")
    head = subprocess.run(["git", "-C", ROOT, "rev-parse", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
    m = {"set": "Evaluation Set v1", "frozen_at_commit": head, "repos": {"forward2028": {}, "sas2r-bench": {}}}
    for repo_key, repo, sets in (("forward2028", ROOT, SETS), ("sas2r-bench", bench, BENCH_SETS)):
        if repo_key == "sas2r-bench":
            if not os.path.isdir(bench):
                sys.exit(f"{bench} not found; --create needs the sas2r-bench checkout")
            m["repos"][repo_key]["commit"] = subprocess.run(["git", "-C", bench, "rev-parse", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
        for name, pats in sets.items():
            entries = []
            for f in collect(repo, pats):
                e = {"path": f, "sha256": sha256(os.path.join(repo, f)), "bytes": os.path.getsize(os.path.join(repo, f))}
                info = run_info(repo, f)
                if info:
                    e["run"] = info
                entries.append(e)
            if not entries:
                sys.exit(f"no files for {repo_key}/{name}: {pats}")
            m["repos"][repo_key][name] = entries
    with open(MANIFEST, "w") as f:
        json.dump(m, f, indent=1)
        f.write("\n")
    n = sum(len(v) for r in m["repos"].values() for k, v in r.items() if isinstance(v, list))
    print(f"wrote {MANIFEST}: {n} files")


def verify(root, bench):
    m = json.load(open(MANIFEST))
    bad, checked, skipped = [], 0, {}
    for repo_key, repo in (("forward2028", root), ("sas2r-bench", bench)):
        present = os.path.isdir(repo)
        for name, entries in m["repos"][repo_key].items():
            if not isinstance(entries, list):
                continue
            if not present:
                skipped[repo_key] = skipped.get(repo_key, 0) + len(entries)
                continue
            for e in entries:
                p = os.path.join(repo, e["path"])
                if not os.path.isfile(p):
                    bad.append(f"missing  {repo_key}/{e['path']}")
                elif sha256(p) != e["sha256"]:
                    bad.append(f"changed  {repo_key}/{e['path']}")
                checked += 1
    for b in bad:
        print(b)
    note = "".join(f", {n} skipped ({k} not checked out)" for k, n in skipped.items())
    print(f"Evaluation Set v1: {checked} files checked, {len(bad)} changed or missing{note}")
    return 1 if bad else 0


if __name__ == "__main__":
    args = sys.argv[1:]
    def opt(flag, default):
        return args[args.index(flag) + 1] if flag in args else default
    root = os.path.abspath(opt("--root", ROOT))
    bench = os.path.abspath(opt("--bench", os.path.join(os.path.dirname(root), "sas2r-bench")))
    if "--create" in args:
        create(bench)
    else:
        sys.exit(verify(root, bench))
