"""Corpus loading and downloads shared by the bench scripts (stdlib only, no PyYAML)."""

import json
import pathlib
import re
import subprocess
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parent.parent
WORK = ROOT / "work"
DATA = ROOT / "data"


def manifest():
    """Parse corpus.yaml through R's yaml package so the manifest stays plain YAML."""
    out = subprocess.run(
        ["Rscript", "-e", "cat(jsonlite::toJSON(yaml::read_yaml('corpus.yaml'), auto_unbox=TRUE))"],
        cwd=ROOT, capture_output=True, text=True, check=True,
    ).stdout
    return json.loads(out)


def load_corpus():
    return manifest()["programs"]


def download(url, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists():
        urllib.request.urlretrieve(url, dest)
    return dest


def fetch_all():
    m = manifest()
    for name, path in m["data"].items():
        download(m["data_base"] + path, DATA / f"{name}.xpt")
    for prog in m["programs"]:
        if prog["source"] == "cdisc":
            download(prog["url"], WORK / "corpus" / f"{prog['id']}.sas")


def program_path(prog):
    if prog["source"] == "cdisc":
        return WORK / "corpus" / f"{prog['id']}.sas"
    return ROOT / prog["path"]


def program_source(prog):
    p = program_path(prog)
    if not p.exists():
        fetch_all()
    text = p.read_text(encoding="utf-8", errors="replace")
    return re.sub(r"\r\n?", "\n", text)
