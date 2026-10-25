# Self-tests for the bench: comparators (planted errors) and the scanner (golden file).
#   Rscript tests/run.R        (from the repo root)
library(testthat)
Sys.setenv(BENCH_ROOT = normalizePath("."))
source("compare/common.R")
for (f in list.files("compare", pattern = "\\.R$", full.names = TRUE)) source(f)
r <- test_dir("tests", reporter = "summary", stop_on_failure = TRUE)
