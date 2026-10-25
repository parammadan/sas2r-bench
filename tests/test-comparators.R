# Comparator self-tests. For each comparator: an out.rds built straight from its reference must
# pass, and the same out.rds with one planted error must fail and name what is wrong. The outputs
# are synthetic (made from refs/), so no CDISC data or translation is needed.
#   Rscript tests/run.R

setwd(Sys.getenv("BENCH_ROOT"))   # comparators read refs/ relative to the repo root
csr501 <- read.csv("refs/csr_t14_5_01.csv", stringsAsFactors = FALSE)
csr502 <- read.csv("refs/csr_t14_5_02.csv", stringsAsFactors = FALSE)
csr101 <- read.csv("refs/csr_t14_1_01.csv", stringsAsFactors = FALSE)
csr201 <- read.csv("refs/csr_t14_2_01.csv", stringsAsFactors = FALSE)
ARMS <- c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")

plant <- function(out, ds, col, row = 1, by = 1) { out[[ds]][[col]][row] <- out[[ds]][[col]][row] + by; out }
expected_datasets <- function(id) jsonlite::fromJSON(file.path("refs", "expected", paste0(id, ".json")))$datasets
with_structure <- function(id, out) {        # add the other datasets the program creates, empty
  for (d in setdiff(expected_datasets(id), names(out))) out[[d]] <- data.frame()
  out
}

test_that("ae_pref: perfect counts pass; one wrong preferred-term count fails", {
  lv <- function(l, term_col) {
    r <- csr501[csr501$level == l & csr501$n > 0, ]
    d <- data.frame(TRTAN = r$trtan, `_FREQ_` = r$n, check.names = FALSE)
    if (!is.null(term_col)) d[[term_col]] <- r$term
    d
  }
  out <- with_structure("ae_pref", list(STAT1 = lv("any", NULL), STAT2 = lv("soc", "AESOC"), STAT3 = lv("pt", "AEDECOD")))
  r <- compare_ae_pref(out, NULL)
  expect_true(r$pass); expect_equal(r$matched, r$cells); expect_equal(r$cells, sum(csr501$n > 0))
  bad <- compare_ae_pref(plant(out, "STAT3", "_FREQ_"), NULL)
  expect_false(bad$pass); expect_equal(bad$matched, r$cells - 1)
  expect_equal(nrow(bad$levels$pt$not_matched), 1)
})

test_that("ae_common: perfect counts pass; one wrong count fails", {
  r0 <- csr501[csr501$level == "pt" & csr501$n > 0, ]
  out <- with_structure("ae_common", list(COUNTAE = data.frame(TRTAN = r0$trtan, AEDECOD = r0$term, NUMAE = r0$n)))
  r <- compare_ae_common(out, NULL)
  expect_true(r$pass); expect_true(r$structure$pass)
  expect_false(compare_ae_common(plant(out, "COUNTAE", "NUMAE", 5), NULL)$pass)
  expect_false(compare_ae_common(list(), NULL)$pass)            # missing dataset is a fail, not an error
})

test_that("at14_5_02: perfect counts and event counts pass; a wrong event count fails", {
  mk <- function(l, term_col) {
    r <- csr502[csr502$level == l & csr502$n > 0, ]
    d <- data.frame(TRTAN = r$trtan, COUNT = r$n, EVENTS = r$events)
    if (!is.null(term_col)) d[[term_col]] <- r$term
    d
  }
  out <- list(RESULTSUB = mk("any", NULL), RESULTSOC = mk("soc", "AEBODSYS"), RESULTPT = mk("pt", "AEDECOD"))
  r <- compare_at14_5_02(out, NULL)
  expect_true(r$pass); expect_equal(r$count_cells_matched, sum(csr502$n > 0)); expect_equal(r$p_values_found, 0)
  bad <- compare_at14_5_02(plant(out, "RESULTPT", "EVENTS"), NULL)
  expect_false(bad$pass); expect_equal(bad$count_cells_matched, r$count_cells_expected - 1)
})

test_that("ae_serious: perfect leading counts pass; one wrong arm count fails", {
  cell <- function(term, a) { n <- csr502$n[csr502$term == term & csr502$trtan == a]; sprintf("%d(%.1f)", n, n) }
  terms <- c("ANY BODY SYSTEM", "SYNCOPE", "PARTIAL SEIZURES WITH SECONDARY GENERALISATION")
  rows <- c("Number of subjects reporting serious adverse events", terms[2:3])
  fr <- data.frame(AEDECOD = rows, TRT0 = sapply(terms, cell, 0), TRT54 = sapply(terms, cell, 54),
                   TRT81 = sapply(terms, cell, 81), CVALUE1 = "0.5", row.names = NULL)
  r <- compare_ae_serious(list(FINAL_REPORT = fr), NULL)
  expect_true(r$pass); expect_equal(r$rows_matched, 3)
  fr$TRT81[2] <- "9(9.0)"
  bad <- compare_ae_serious(list(FINAL_REPORT = fr), NULL)
  expect_false(bad$pass); expect_false(bad$rows$SYNCOPE$counts_match); expect_true(bad$rows$any$counts_match)
})

test_that("disposition_a: completers and non-completers from the CSR pass; a wrong count fails", {
  v <- function(row, arm) csr101$value[csr101$row == row & csr101$column == arm & csr101$statistic == "count"]
  done <- sapply(ARMS, v, row = "Complete Study"); itt <- sapply(ARMS, v, row = "Intent-To-Treat (ITT)")
  tl <- data.frame(TRT01PN = rep(c(0, 54, 81), each = 2), DISCONFL = rep(c(NA, "Y"), 3),
                   COUNT = as.vector(rbind(done, itt - done)))
  out <- with_structure("disposition_a", list(TOPLINE = tl))
  r <- compare_disposition_a(out, NULL)
  expect_true(r$pass); expect_equal(r$matched, 3)
  bad <- compare_disposition_a(plant(out, "TOPLINE", "COUNT", 1), NULL)   # placebo completers off by one
  expect_false(bad$pass); expect_false(bad$arms[[1]]$match); expect_true(bad$arms[[2]]$match)
})

test_that("demo_summary: a REPORT printed from the CSR at the program's precision passes; a wrong mean fails", {
  ref <- function(row, arm, stat) csr201$value[csr201$row == row & csr201$column == arm & csr201$statistic == stat]
  f1 <- function(x) sprintf("%.1f", x)
  block <- function(name, csr_row) {
    st <- list(n = function(a) sprintf("%d", ref(csr_row, a, "n")), Mean = function(a) f1(ref(csr_row, a, "mean")),
               SD = function(a) f1(ref(csr_row, a, "sd")), Median = function(a) f1(ref(csr_row, a, "median")),
               `Min, Max` = function(a) sprintf("%g, %g", ref(csr_row, a, "min"), ref(csr_row, a, "max")))
    do.call(rbind, lapply(seq_along(st), function(i) data.frame(NAME = if (i == 1) name else "", COL0 = names(st)[i],
      COL1 = st[[i]](ARMS[1]), COL2 = st[[i]](ARMS[2]), COL3 = st[[i]](ARMS[3]))))
  }
  cats <- function(name, labels, csr_rows, pct) do.call(rbind, lapply(seq_along(labels), function(i) {
    cellv <- function(a) if (pct) sprintf("%d (%.1f%%)", ref(csr_rows[i], a, "count"), ref(csr_rows[i], a, "pct")) else sprintf("%d", ref(csr_rows[i], a, "count"))
    data.frame(NAME = if (i == 1) name else "", COL0 = labels[i], COL1 = cellv(ARMS[1]), COL2 = cellv(ARMS[2]), COL3 = cellv(ARMS[3]))
  }))
  rep_ <- rbind(block("Age (y)", "Age (y)"), block("Weight (kg)", "Baseline weight(kg)"),
                cats("Sex", c("M", "F"), c("Male", "Female"), FALSE),
                cats("Age Category", c("<65", "65-80", ">80"), c("<65 yrs", "65-80 yrs", ">80 yrs"), TRUE))
  r <- compare_demo_summary(list(REPORT = rep_), NULL)
  expect_equal(r$numbers, 60); expect_equal(r$matched, 60); expect_true(r$pass)
  rep_$COL2[2] <- f1(ref("Age (y)", ARMS[2], "mean") + 0.3)
  bad <- compare_demo_summary(list(REPORT = rep_), NULL)
  expect_false(bad$pass); expect_equal(bad$matched, 59)
  expect_equal(bad$not_matched[[1]]$what, "Age (y) mean")
})

test_that("structure: every expected dataset present passes; one missing fails and is named", {
  for (id in c("mean_time", "outliers_Shift_Table", "outliers_TEHiLo_Anytime", "pk_param_summary", "medication_listing")) {
    want <- expected_datasets(id)
    out <- setNames(lapply(want, function(x) data.frame()), paste0("WORK.", tolower(want)))   # libref and case ignored
    expect_true(get(paste0("compare_", id))(out, NULL)$pass, info = id)
    bad <- get(paste0("compare_", id))(out[-1], NULL)
    expect_false(bad$pass, info = id); expect_equal(bad$datasets_missing, toupper(want[1]), info = id)
  }
})

test_that("adae: an exact copy of the reference passes; one changed cell fails on that column only", {
  dir <- tempfile("adae"); dir.create(dir)
  ref <- data.frame(USUBJID = c("01-701-1015", "01-701-1015", "01-701-1023"), AESEQ = c(1, 2, 1),
                    AEDECOD = c("APPLICATION SITE ERYTHEMA", "DIARRHOEA", "ERYTHEMA"),
                    ASTDT = as.Date(c("2014-01-03", "2014-01-09", "2012-08-07")), stringsAsFactors = FALSE)
  haven::write_xpt(ref, file.path(dir, "adae.xpt"))
  got <- ref; got$ASTDT <- as.numeric(got$ASTDT - as.Date("1960-01-01"))     # SAS day counts are accepted
  r <- compare_adae(list(ADAE = got[c(3, 1, 2), ]), dir)                     # row order doesn't matter
  expect_true(r$pass); expect_equal(r$cell_match, 1)
  got$AEDECOD[2] <- "DIARRHEA"
  bad <- compare_adae(list(ADAE = got), dir)
  expect_false(bad$pass); expect_equal(names(bad$columns_partial), "AEDECOD")
  expect_false(compare_adae(list(ADAE = got[-1]), dir)$pass)                 # a missing row fails
})
