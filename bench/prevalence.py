"""How common is each risky construct? Scan every .sas file in phuse-org/phuse-scripts plus the
CDISC Pilot programs, and write results/construct_prevalence.json.

    python3 bench/prevalence.py
"""

import json
import pathlib
import subprocess
import sys
import urllib.parse
import urllib.request

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from corpus import WORK, download, fetch_all  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent
TREE = "https://api.github.com/repos/phuse-org/phuse-scripts/git/trees/master?recursive=1"
RAW = "https://raw.githubusercontent.com/phuse-org/phuse-scripts/master/"


def main():
    fetch_all()
    with urllib.request.urlopen(TREE) as r:
        tree = json.load(r)["tree"]
    paths = [t["path"] for t in tree if t["type"] == "blob" and t["path"].lower().endswith(".sas")]
    files = []
    for p in paths:
        dest = WORK / "phuse_all" / p.replace("/", "__").replace(" ", "_")
        download(RAW + urllib.parse.quote(p), dest)
        files.append(str(dest))
    files += [str(p) for p in sorted((WORK / "corpus").glob("*.sas"))]

    out = subprocess.run([sys.executable, "-I", str(ROOT / "bench" / "scan.py"), "--summary", *files],
                         capture_output=True, text=True, check=True).stdout
    summary = json.loads(out)
    summary["sources"] = {"phuse-scripts (MIT)": len(paths), "CDISC Pilot": len(files) - len(paths)}
    (ROOT / "results" / "construct_prevalence.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
