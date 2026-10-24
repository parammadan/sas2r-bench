#!/bin/sh
# One command: fetch data, translate, run, score, report.
#   ./run_all.sh                                   # default: claude-code:claude-opus-5-5 (your Claude Code login)
#   ./run_all.sh ollama:qwen2.5-coder:7b           # local open-source model (needs 16 GB+ RAM)
#   ./run_all.sh anthropic:claude-opus-5-5         # Anthropic API key, per-token cost
set -e
cd "$(dirname "$0")"
translator="${1:-claude-code:claude-opus-5-5}"
slug=$(printf '%s' "$translator" | sed 's/[^A-Za-z0-9._-][^A-Za-z0-9._-]*/_/g')

python3 -I -c "import sys; sys.path.insert(0, 'bench'); import corpus; corpus.fetch_all()"
python3 -I bench/translate.py --translator "$translator"
Rscript bench/run.R "$slug"
python3 -I bench/report.py
