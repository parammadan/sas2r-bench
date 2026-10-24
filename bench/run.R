# Run each translation and score it.
#
#   Rscript bench/run.R <translator-slug> [id ...]
#
# Each translation runs in work/run/<slug>/<id>/ with data/ linked in, under macOS sandbox-exec:
# no network, and writes only inside that folder and R's temp dir. The comparator for the program
# (compare/<id>.R) then scores out.rds, and the result goes to results/<slug>/<id>.json.

suppressPackageStartupMessages({
  library(jsonlite)
})
source("compare/common.R")
for (f in list.files("compare", pattern = "\\.R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
slug <- args[1]
corpus <- yaml::read_yaml("corpus.yaml")
ids <- vapply(corpus$programs, `[[`, character(1), "id")
if (length(args) > 1) ids <- intersect(ids, args[-1])

root <- normalizePath(".")
data_dir <- file.path(root, "data")
timeout_s <- 600

sandbox_profile <- function(run_dir) {
  tmp <- normalizePath(tempdir(), mustWork = FALSE)
  sprintf('(version 1)
(allow default)
(deny network*)
(deny file-write* (subpath "%s"))
(allow file-write* (subpath "%s") (subpath "%s") (subpath "/private/var/folders") (subpath "/private/tmp") (literal "/dev/null"))',
    Sys.getenv("HOME"), run_dir, tmp)
}

dir.create(file.path("results", slug), recursive = TRUE, showWarnings = FALSE)

for (id in ids) {
  src <- file.path(root, "translations", slug, paste0(id, ".R"))
  result <- list(id = id, translator = slug)
  if (!file.exists(src) || file.size(src) == 0) {
    result$stage <- "translate"
    result$pass <- FALSE
    result$reason <- "no translation"
  } else {
    run_dir <- file.path(root, "work", "run", slug, id)
    unlink(run_dir, recursive = TRUE)
    dir.create(run_dir, recursive = TRUE)
    file.symlink(data_dir, file.path(run_dir, "data"))
    file.copy(src, file.path(run_dir, "translation.R"))
    profile <- file.path(run_dir, ".sandbox.sb")
    writeLines(sandbox_profile(run_dir), profile)

    t0 <- Sys.time()
    old_wd <- setwd(run_dir)
    status <- system2(
      "sandbox-exec", c("-f", profile, "Rscript", "--vanilla", "translation.R"),
      stdout = file.path(run_dir, "stdout.txt"), stderr = file.path(run_dir, "stderr.txt"),
      timeout = timeout_s,
      wait = TRUE
    )
    setwd(old_wd)
    result$seconds <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    result$exit_status <- status
    err <- readLines(file.path(run_dir, "stderr.txt"), warn = FALSE)
    result$stderr_tail <- tail(err[nzchar(err)], 8)

    out_file <- file.path(run_dir, "out.rds")
    if (status != 0 || !file.exists(out_file)) {
      result$stage <- "run"
      result$pass <- FALSE
      result$reason <- if (status == 124) "timeout" else if (status != 0) "R error" else "no out.rds"
    } else {
      out <- tryCatch(readRDS(out_file), error = function(e) NULL)
      result$datasets <- names(out)
      scored <- tryCatch(get(paste0("compare_", id))(out, data_dir),
                         error = function(e) list(pass = FALSE, reason = paste("comparator error:", conditionMessage(e))))
      result$stage <- "compare"
      result <- c(result, scored)
    }
  }
  write(toJSON(result, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = NA),
        file.path("results", slug, paste0(id, ".json")))
  cat(sprintf("%-12s %-8s pass=%s %s\n", id, result$stage, isTRUE(result$pass),
              if (!is.null(result$reason)) result$reason else ""))
}
