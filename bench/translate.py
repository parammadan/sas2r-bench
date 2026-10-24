"""Translate each corpus program from SAS to R with a chosen translator.

Usage:
    python3 bench/translate.py --translator ollama:qwen2.5-coder:7b [--only ID ...]
    python3 bench/translate.py --translator claude-code:claude-opus-5-5
    python3 bench/translate.py --translator anthropic:claude-opus-5-5

Translators:
    ollama:<model>       local open-source model through Ollama (free; `ollama serve` must be running)
    claude-code:<model>  the Claude Code CLI in print mode with every tool disabled (uses your
                         Claude Code login, no API key)
    anthropic:<model>    the Anthropic API through the official SDK (needs ANTHROPIC_API_KEY;
                         real per-token cost)

Every translator gets the same prompt (prompts/translate.md + the SAS source + the input file
mapping) and nothing else: no data, no references, no tools. Output goes to
translations/<translator>/<id>.R with a .json sidecar recording the model, prompt hash,
time, tokens and reported cost.
"""

import argparse
import datetime as dt
import hashlib
import json
import pathlib
import re
import subprocess
import sys
import time
import urllib.request

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from corpus import load_corpus, program_source  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent
SYSTEM = (ROOT / "prompts" / "translate.md").read_text()


def user_message(prog, sas):
    lines = [f"  {ref} -> data/{name}.xpt" for ref, name in prog["inputs"].items()]
    return (
        f"SAS program `{prog['id']}.sas`:\n\n```sas\n{sas}\n```\n\n"
        "Input mapping (SAS reference -> local file, relative to the working directory):\n"
        + "\n".join(lines)
        + "\nAny URL in the program that points to adae.xpt or adsl.xpt means the local file of the same name.\n"
    )


def extract_r(text):
    blocks = re.findall(r"```(?:r|R)?\s*\n(.*?)```", text, flags=re.S)
    return max(blocks, key=len) if blocks else text


def via_ollama(model, system, user):
    body = json.dumps({
        "model": model,
        "stream": False,
        "options": {"temperature": 0, "seed": 1, "num_ctx": 16384},
        "messages": [{"role": "system", "content": system}, {"role": "user", "content": user}],
    }).encode()
    req = urllib.request.Request("http://localhost:11434/api/chat", data=body,
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=3600) as r:
        d = json.load(r)
    return d["message"]["content"], {
        "input_tokens": d.get("prompt_eval_count"),
        "output_tokens": d.get("eval_count"),
        "cost_usd": 0.0,
    }


def via_claude_code(model, system, user):
    cmd = ["claude", "-p", "--tools", "", "--model", model, "--output-format", "json",
           "--no-session-persistence", "--system-prompt", system]
    p = subprocess.run(cmd, input=user, capture_output=True, text=True, timeout=3600,
                       cwd="/tmp")  # run outside the repo so no project files are in reach
    if p.returncode != 0:
        raise RuntimeError(p.stderr.strip() or p.stdout.strip())
    d = json.loads(p.stdout)
    if d.get("is_error"):
        raise RuntimeError(d.get("result"))
    usage = d.get("usage", {})
    return d["result"], {
        "input_tokens": usage.get("input_tokens"),
        "output_tokens": usage.get("output_tokens"),
        "cost_usd": d.get("total_cost_usd"),
        "cost_note": "equivalent API cost reported by Claude Code; on a subscription it counts toward plan usage",
    }


def via_anthropic(model, system, user):
    import anthropic  # pip install anthropic

    client = anthropic.Anthropic()
    with client.messages.stream(
        model=model,
        max_tokens=64000,
        output_config={"effort": "high"},
        system=system,
        messages=[{"role": "user", "content": user}],
    ) as stream:
        msg = stream.get_final_message()
    if msg.stop_reason == "refusal":
        raise RuntimeError(f"refused: {msg.stop_details}")
    text = "".join(b.text for b in msg.content if b.type == "text")
    price_in, price_out = {"claude-opus-5-5": (4.0, 20.0)}.get(model, (None, None))
    cost = None
    if price_in is not None:
        cost = (msg.usage.input_tokens * price_in + msg.usage.output_tokens * price_out) / 1e6
    return text, {"input_tokens": msg.usage.input_tokens,
                  "output_tokens": msg.usage.output_tokens, "cost_usd": cost}


BACKENDS = {"ollama": via_ollama, "claude-code": via_claude_code, "anthropic": via_anthropic}


def slug(translator):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", translator)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--translator", default="ollama:qwen2.5-coder:7b")
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--force", action="store_true", help="redo translations that already exist")
    a = ap.parse_args()

    backend, _, model = a.translator.partition(":")
    call = BACKENDS[backend]
    outdir = ROOT / "translations" / slug(a.translator)
    outdir.mkdir(parents=True, exist_ok=True)
    prompt_hash = hashlib.sha256(SYSTEM.encode()).hexdigest()[:12]

    for prog in load_corpus():
        if a.only and prog["id"] not in a.only:
            continue
        dest = outdir / f"{prog['id']}.R"
        if dest.exists() and not a.force:
            print(f"{prog['id']}: exists, skipping")
            continue
        user = user_message(prog, program_source(prog))
        t0 = time.time()
        try:
            text, usage = call(model, SYSTEM, user)
            error = None
        except Exception as e:  # recorded, not raised: a failed translation is a result
            text, usage, error = "", {}, str(e)
        meta = {
            "id": prog["id"], "translator": a.translator, "prompt_sha": prompt_hash,
            "when": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            "seconds": round(time.time() - t0, 1), "error": error, **usage,
        }
        dest.write_text(extract_r(text) if text else "")
        (outdir / f"{prog['id']}.json").write_text(json.dumps(meta, indent=2) + "\n")
        (outdir / f"{prog['id']}.raw.md").write_text(text)
        print(f"{prog['id']}: {'ERROR ' + error if error else 'ok'} ({meta['seconds']}s)")


if __name__ == "__main__":
    main()
