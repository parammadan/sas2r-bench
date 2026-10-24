# demo_summary.sas (PhUSE scriptathon) -> published reference: CSR Table 14-2.01.
# The program's REPORT dataset has NAME (variable label, first row of each block), COL0 (the
# statistic or category label) and COL1-COL3 (text cells for arms in TRT01AN order: 0, 54, 81).
# Scored: age and weight n / mean / SD / median / min / max, sex counts, age-group counts and
# percentages, compared at the precision the program prints (it prints SD with one decimal and
# percentages with one decimal, the CSR with two and none). Where the two print at different
# precisions, they match when both could come from the same true value: the rounding interval of
# the program's printed number overlaps the rounding interval of the CSR's.
# Not scored: race (the 2006 CSR's groups differ from the released ADSL's RACE coding) and
# ethnicity (not in the CSR).

compare_demo_summary <- function(out, data_dir) {
  csr <- read.csv("refs/csr_t14_2_01.csv", stringsAsFactors = FALSE)
  d <- pick(out, "REPORT")
  if (is.null(d)) return(list(reference = "published", pass = FALSE, reason = "no REPORT dataset in out.rds"))
  for (v in c("NAME", "COL0", "COL1", "COL2", "COL3")) if (!v %in% names(d)) d[[v]] <- NA
  d$NAME <- trimws(as.character(d$NAME)); d$COL0 <- trimws(as.character(d$COL0))
  # carry the block name down to the rows under it
  blk <- d$NAME; for (i in seq_along(blk)) if ((is.na(blk[i]) || blk[i] == "") && i > 1) blk[i] <- blk[i - 1]
  d$BLOCK <- tolower(blk)
  arms <- c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")
  nums <- function(x) as.numeric(regmatches(x, gregexpr("-?\\d+\\.?\\d*", x))[[1]])
  r2 <- function(x, k) sign(x) * floor(abs(x) * 10^k + 0.5 + 1e-9) / 10^k

  checks <- list()
  # k_got: decimals the program prints; k_want: decimals the CSR prints
  add <- function(what, arm, want, got, k_got, k_want = k_got) {
    ok <- length(got) == 1 && !is.na(got) &&
      abs(got - want) <= 0.5 * 10^-k_got + 0.5 * 10^-k_want + 1e-9 &&
      (k_got != k_want || abs(r2(got, k_got) - r2(want, k_want)) < 1e-9)
    checks[[length(checks) + 1]] <<- list(what = what, arm = arm, csr = want, got = if (length(got)) got else NA, match = ok)
  }
  cell <- function(block_pat, label_pat, j) {
    rows <- d[grepl(block_pat, d$BLOCK) & grepl(label_pat, d$COL0, ignore.case = TRUE), , drop = FALSE]
    if (nrow(rows) != 1) return(NA_character_)
    as.character(rows[[paste0("COL", j)]])
  }
  ref <- function(row, arm, stat) as.numeric(csr$value[csr$row == row & csr$column == arm & csr$statistic == stat])

  for (blk in list(list(pat = "^age", csr = "Age (y)"), list(pat = "weight", csr = "Baseline weight(kg)"))) {
    for (j in 1:3) {
      a <- arms[j]
      add(paste(blk$csr, "n"), a, ref(blk$csr, a, "n"), nums(cell(blk$pat, "^n", j))[1], 0)
      add(paste(blk$csr, "mean"), a, ref(blk$csr, a, "mean"), nums(cell(blk$pat, "^mean", j))[1], 1)
      add(paste(blk$csr, "sd"), a, ref(blk$csr, a, "sd"), nums(cell(blk$pat, "^sd", j))[1], 1, 2)
      add(paste(blk$csr, "median"), a, ref(blk$csr, a, "median"), nums(cell(blk$pat, "^median", j))[1], 1)
      mm <- nums(cell(blk$pat, "^min", j))
      add(paste(blk$csr, "min"), a, ref(blk$csr, a, "min"), mm[1], 0)
      add(paste(blk$csr, "max"), a, ref(blk$csr, a, "max"), if (length(mm) >= 2) mm[2] else NA, 0)
    }
  }
  for (cat in list(c("^m$|^male", "Male"), c("^f$|^female", "Female"))) for (j in 1:3)
    add(paste("Sex", cat[2], "count"), arms[j], ref(cat[2], arms[j], "count"), nums(cell("^sex", cat[1], j))[1], 0)
  for (cat in list(c("<\\s*65", "<65 yrs"), c("65\\s*-\\s*80", "65-80 yrs"), c(">\\s*80", ">80 yrs"))) for (j in 1:3) {
    v <- nums(cell("^age cat", cat[1], j))
    add(paste("Age group", cat[2], "count"), arms[j], ref(cat[2], arms[j], "count"), v[1], 0)
    add(paste("Age group", cat[2], "pct"), arms[j], ref(cat[2], arms[j], "pct"), if (length(v) >= 2) v[2] else NA, 1, 0)
  }
  ok <- vapply(checks, `[[`, logical(1), "match")
  list(reference = "published", numbers = length(checks), matched = sum(ok),
       not_matched = lapply(checks[!ok], function(x) x[c("what", "arm", "csr", "got")]),
       pass = all(ok))
}
