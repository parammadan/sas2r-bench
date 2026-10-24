# sas2r-bench

Can an AI translate a clinical-trial SAS program into R that gives the same answer? This bench
runs the translation, executes the R on real data, and checks the output against published
results, so the answer is a measurement instead of an impression.

Part of [Forward 2028](https://github.com/parammadan/forward2028): AI drafts, the harness
verifies, a programmer signs off.

## How it works

1. **Corpus** (`corpus.yaml`): SAS programs that run on the public CDISC Pilot Study data.
2. **Scan** (`bench/scan.py`): flags constructs that commonly break naive translation (RETAIN,
   `first.`/`last.`, MERGE without BY, missing-value comparisons, `??` informats, macros, ODS
   OUTPUT, ...), with line numbers.
3. **Translate** (`bench/translate.py`): one prompt (`prompts/translate.md`), the SAS source and
   the input file mapping. No data, no references, no tools.
4. **Run** (`bench/run.R`): executes each translation under macOS `sandbox-exec` (no network,
   writes only to its own folder).
5. **Compare** (`compare/*.R`): scores the output. Comparators were committed before any
   translation existed.
6. **Report** (`bench/report.py`): `results/README.md` and `results/summary.json`.

## Run it

```sh
./run_all.sh                                # free, local: ollama:qwen2.5-coder:7b
./run_all.sh claude-code:claude-opus-5-5    # uses your Claude Code login, no API key
./run_all.sh anthropic:claude-opus-5-5      # Anthropic API key, per-token cost
```

Needs R (haven, dplyr, jsonlite, yaml) and Python 3. The local option needs
[Ollama](https://ollama.com) with `ollama pull qwen2.5-coder:7b`. The API option needs
`pip install anthropic` and `ANTHROPIC_API_KEY`.

## Reference levels

| Level | Meaning |
|---|---|
| exact | The program's own output is published, so a faithful translation must match it |
| published | A published CSR table is the check; a mismatch can be the SAS program's own bug |
| none | Needs a real SAS run before it can be scored |

## Results

See [results/README.md](results/README.md) and [results/FINDINGS.md](results/FINDINGS.md).
Every number comes from a real run; nothing is filled in by hand except the CSR reference
values, which are copied from the report and cited in [refs/README.md](refs/README.md).

## Sources and licenses

- CDISC SDTM/ADaM Pilot Project data and programs (credit: CDISC). Downloaded at run time,
  never committed: <https://github.com/cdisc-org/sdtm-adam-pilot-project>
- PhUSE scripts (MIT), `corpus/phuse/`, license in `corpus/phuse/LICENSE-phuse-scripts.md`:
  <https://github.com/phuse-org/phuse-scripts>
- This repo's own code: MIT (`LICENSE`).

Independent project. Not affiliated with, endorsed by, or developed for Harvard University or
the Center for Biostatistics in AIDS Research (CBAR).
