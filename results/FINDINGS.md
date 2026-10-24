# Findings, run 1

Translator: `claude-code:claude-opus-5-5` (Claude Code print mode, all tools disabled, prompt in
`prompts/translate.md`). Each program was translated once, blind: no data, no references.
The comparators were committed before any translation existed (see git history).

| Program | Pre-registered check | What actually happened |
|---|---|---|
| `adae.sas` | **fail** | 1,191/1,191 rows matched; **53 of 55 released columns identical in every row**. The 2 "missing" columns are a reference issue, not a translation error (below). |
| `at14-5-02.sas` | pass | 7/7 non-zero count cells match CSR Table 14-5.02. **No p-values produced**, matching what the SAS program itself selects (below). |
| `ae_serious.sas` | pass | Subject counts match the CSR in all 3 rows. The program's own mistakes were reproduced, not fixed (below). |

## What the run found

### 1. The released ADAE was not produced by `adae.sas` alone
The program derives `ADUR` (line 15) and keeps it under that name (lines 120, 125, 131). The
released `adae.xpt` has no `ADUR`; it has `ADURN` (label "AE Duration (N)") and `ADURU` ("DAY"
in 714 rows, blank in 477). The translation's `ADUR` equals the released `ADURN` in every row.
So some later step, not in the CDISC package, renamed the variable and added the unit. A
"faithful translation" check against the released file can therefore never fully pass. This is
recorded rather than patched: the comparator stays as registered.

### 2. `at14-5-02.sas` cannot produce its p-values on this data
The Fisher tests select `TRTAN in (0,2)` and `TRTAN in (0,1)` (lines 110, 117, 140, 147, 167,
174). The data codes arms as 0, 54 and 81, so every test table holds placebo only. The
translation reproduced this and produced no p-values. The CSR's p-values (0.494, 0.243) came
from a different program (`saetable.sas`, named in the CSR footer). Confirming SAS's exact
behavior here (empty ODS output vs. a warning) still needs a real SAS run.

### 3. `ae_serious.sas` differs from the CSR by design and by mistake
- **Percent column shows the count.** The cell text is built with `put(count,8.1)` where a
  percentage was meant, so it prints `1(1.0)` instead of `1 (1.2%)`. The translation kept it.
- **One 3-arm test, not two pairwise tests.** `table trtan*flg / fisher` runs one Fisher exact
  test across all three arms. Checked independently in R: 0.3255 (any SAE), 0.5503 (syncope),
  0.6614 (partial seizures), exactly what the translation produced. The CSR reports
  placebo-vs-low and placebo-vs-high instead.

## What this means for the workflow
- On these three programs the translation was **faithful**: every difference from the published
  numbers traced back to the SAS program or to the reference, not to the R.
- That is exactly why "matches the published table" is not enough on its own. A faithful
  translation of a program with a bug reproduces the bug; a translation that "fixes" it silently
  is a different program. Both checks are needed, and a human decides which one applies.
- Three programs is a small corpus. These are observations, not a pass rate.

# Run 1 addition: `demo_summary.sas` (PhUSE, 303 lines)

Same translator and prompt, translated once, blind. The comparator was committed before the
translation existed.

**Pre-registered check: fail, "no REPORT dataset". That is the correct outcome**, because the SAS
program itself never creates REPORT:

1. **Wrong informat name.** The `%freqstats` call for age groups passes `ordfmt=agegrfmto`, and the
   macro runs `ord = input(agegr1, agegrfmto.)`. The program defines that informat as `agegrfmt`
   (line 43), not `agegrfmto`. In SAS that is a compile error, so the step stops and `AGEGRRES` is
   never created. The final `data report; set ageres weightres sexres agegrres raceres ethnicres;`
   then fails on the missing dataset, so `REPORT` is never created either. The translation
   reproduced both failures and said why in comments.
2. **Every arm gets the placebo percentage.** Inside the macro, `col&i = put(n&i,6.0-R) || pct1;`
   uses `pct1` for every arm, so each arm prints the first arm's percentage (Low Dose shows
   "50 (61.6%)"; the true figure is 59.5%). The translation kept this too.

**Post-hoc, labelled as such:** scoring the summary datasets the program does create (`AGERES`,
`WEIGHTRES`, `SEXRES`) with the same comparator gives 42 of 60: every age and weight statistic
(36) and every sex count (6) matches CSR Table 14-2.01. The other 18 are the age-group rows that the
informat bug prevents.
