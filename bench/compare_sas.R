# Compare each R translation with its program's real SAS output, dataset by dataset, cell by cell.
#
#   Rscript bench/compare_sas.R <translator-slug> [sas_out_dir]
#
# sas_out_dir (default work/sas_out) holds <id>/<DATASET>.csv files exported by the SAS kit
# (bench/sas_kit.py). R outputs come from work/run/<slug>/<id>/out.rds (bench/run.R).
# Matching rules: datasets by name; rows in order after sorting both sides by all shared
# columns; numbers equal to 8 decimals; text trimmed; SAS dates (exported as DATE9. text such as
# 07MAR2013) compared with R Dates. Results go to results/sas_vs_r/<slug>/<id>.json.

suppressPackageStartupMessages(library(jsonlite))

args <- commandArgs(trailingOnly = TRUE)
slug <- args[1]
sas_dir <- if (length(args) > 1) args[2] else file.path("work", "sas_out")
out_dir <- file.path("results", "sas_vs_r", slug)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

as_cmp <- function(x) {
  if (inherits(x, "Date")) return(format(x, "%Y-%m-%d"))
  if (is.factor(x)) x <- as.character(x)
  if (is.numeric(x)) return(ifelse(is.na(x), NA, trimws(formatC(round(x, 8), format = "fg", digits = 15))))
  x <- trimws(as.character(x))
  d9 <- grepl("^\\d{2}[A-Z]{3}\\d{4}$", x)
  if (any(d9, na.rm = TRUE) && mean(d9 | is.na(x) | x == "", na.rm = TRUE) == 1) {
    x <- format(as.Date(x, format = "%d%b%Y"), "%Y-%m-%d")
  } else {
    num <- suppressWarnings(as.numeric(x))
    if (all(is.na(num) == (is.na(x) | x == ""))) return(ifelse(is.na(num), NA, trimws(formatC(round(num, 8), format = "fg", digits = 15))))
  }
  x[x == ""] <- NA
  x
}

compare_one <- function(sas, r) {
  names(sas) <- toupper(names(sas)); names(r) <- toupper(names(r))
  cols <- intersect(names(sas), names(r))
  if (!length(cols)) return(list(rows_sas = nrow(sas), rows_r = nrow(r), columns_shared = 0, cell_match = 0))
  s <- as.data.frame(lapply(sas[cols], as_cmp), stringsAsFactors = FALSE, check.names = FALSE)
  t <- as.data.frame(lapply(as.data.frame(r, check.names = FALSE)[cols], as_cmp), stringsAsFactors = FALSE, check.names = FALSE)
  # sort by text columns first (identifiers, labels), then numbers, so one differing number
  # changes its own row's position as little as possible
  is_num <- vapply(s, function(x) all(is.na(x) | grepl("^-?[0-9.eE+-]+$", x)), logical(1))
  ord <- c(cols[!is_num], cols[is_num])
  key <- function(d) do.call(order, c(unname(as.list(d[ord])), list(na.last = TRUE)))
  s <- s[key(s), , drop = FALSE]; t <- t[key(t), , drop = FALSE]
  n <- min(nrow(s), nrow(t))
  per_col <- vapply(cols, function(v) if (n == 0) 0 else
    mean((is.na(s[[v]][1:n]) & is.na(t[[v]][1:n])) | (!is.na(s[[v]][1:n]) & !is.na(t[[v]][1:n]) & s[[v]][1:n] == t[[v]][1:n])), numeric(1))
  list(rows_sas = nrow(sas), rows_r = nrow(r),
       columns_sas_only = setdiff(names(sas), names(r)), columns_r_only = setdiff(names(r), names(sas)),
       columns_shared = length(cols), columns_exact = sum(per_col == 1),
       columns_partial = as.list(round(per_col[per_col < 1], 4)),
       cell_match = round(mean(per_col), 4),
       identical = nrow(sas) == nrow(r) && all(per_col == 1) && setequal(names(sas), names(r)))
}

ids <- list.dirs(sas_dir, full.names = FALSE, recursive = FALSE)
for (id in ids) {
  csvs <- list.files(file.path(sas_dir, id), pattern = "\\.csv$", full.names = TRUE)
  rds <- file.path("work", "run", slug, id, "out.rds")
  res <- list(id = id, translator = slug)
  if (!length(csvs)) { res$reason <- "no SAS output exported" }
  else if (!file.exists(rds)) { res$reason <- "no R output: run bench/run.R first" }
  else {
    r_out <- readRDS(rds)
    names(r_out) <- toupper(sub("^.*\\.", "", names(r_out)))
    per_ds <- lapply(csvs, function(f) {
      nm <- toupper(sub("\\.csv$", "", basename(f)))
      sas <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "."))
      if (!nm %in% names(r_out)) return(list(dataset = nm, in_r = FALSE, rows_sas = nrow(sas)))
      c(list(dataset = nm, in_r = TRUE), compare_one(sas, r_out[[nm]]))
    })
    res$datasets_sas <- length(per_ds)
    res$datasets_in_r <- sum(vapply(per_ds, `[[`, logical(1), "in_r"))
    res$datasets_identical <- sum(vapply(per_ds, function(x) isTRUE(x$identical), logical(1)))
    res$datasets_r_only <- setdiff(names(r_out), toupper(sub("\\.csv$", "", basename(csvs))))
    res$detail <- per_ds
  }
  write(toJSON(res, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = NA),
        file.path(out_dir, paste0(id, ".json")))
  cat(sprintf("%-24s %s\n", id, if (!is.null(res$reason)) res$reason else
    sprintf("SAS datasets %d, in R %d, identical %d", res$datasets_sas, res$datasets_in_r, res$datasets_identical)))
}
