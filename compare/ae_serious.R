# ae_serious.sas (PhUSE scriptathon) -> published reference: CSR Table 14-5.02.
# The program's FINAL_REPORT has AEDECOD, TRT0, TRT54, TRT81 (text "count(x)") and CVALUE1
# (Fisher p-value). Subject counts are the leading integer of each cell and are checked
# against the CSR. The text in brackets is reported, not scored: reading the source, the program
# formats the count where a percentage would go, so a faithful translation should copy that.

compare_ae_serious <- function(out, data_dir) {
  csr <- read.csv("refs/csr_t14_5_02.csv", stringsAsFactors = FALSE)
  d <- pick(out, "FINAL_REPORT")
  if (is.null(d)) {
    return(list(reference = "published", pass = FALSE, reason = "no FINAL_REPORT dataset in out.rds"))
  }
  arm_cols <- c(`0` = "TRT0", `54` = "TRT54", `81` = "TRT81")
  lead_int <- function(x) suppressWarnings(as.numeric(sub("^\\s*(\\d+).*$", "\\1", x)))

  rows <- list(
    any = "Number of subjects reporting serious adverse events",
    SYNCOPE = "SYNCOPE",
    PSG = "PARTIAL SEIZURES WITH SECONDARY GENERALISATION"
  )
  want_term <- c(any = "ANY BODY SYSTEM", SYNCOPE = "SYNCOPE",
                 PSG = "PARTIAL SEIZURES WITH SECONDARY GENERALISATION")

  per_row <- lapply(names(rows), function(k) {
    r <- d[toupper(trimws(d$AEDECOD)) == toupper(rows[[k]]), , drop = FALSE]
    if (nrow(r) != 1) return(list(found = FALSE))
    w <- csr[csr$term == want_term[[k]], ]
    got <- vapply(arm_cols, function(col) {
      x <- if (col %in% names(r)) r[[col]] else NA
      v <- lead_int(x)
      if (is.na(v)) 0 else v      # an absent arm cell means nobody had the event
    }, numeric(1))
    want <- w$n[match(as.numeric(names(arm_cols)), w$trtan)]
    list(found = TRUE, counts_got = unname(got), counts_csr = want,
         counts_match = all(got == want),
         cells_text = unname(vapply(arm_cols, function(col) if (col %in% names(r)) as.character(r[[col]]) else NA_character_, character(1))),
         p_value_text = if ("CVALUE1" %in% names(r)) trimws(as.character(r$CVALUE1)) else NA)
  })
  names(per_row) <- names(rows)

  ok <- vapply(per_row, function(x) isTRUE(x$counts_match), logical(1))
  list(reference = "published", rows = per_row, rows_expected = length(rows),
       rows_matched = sum(ok), pass = all(ok))
}
