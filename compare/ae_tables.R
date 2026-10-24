# ae_pref.sas and ae_common.sas -> published reference: CSR Table 14-5.01 (refs/csr_t14_5_01.csv).
# Both count subjects with treatment-emergent AEs through the first-occurrence flags, so their
# counts can be checked against the CSR's n per arm. A CSR term printed at its column width
# (about 38-40 characters) is matched as a prefix of the full term.
# Scored: every non-zero CSR cell at the levels the program produces. Each program's structure
# (refs/expected) is reported alongside.

csr_14_5_01 <- function() read.csv("refs/csr_t14_5_01.csv", stringsAsFactors = FALSE)

term_match <- function(data_terms, csr_term) {
  d <- toupper(trimws(data_terms)); t <- toupper(trimws(csr_term))
  d == t | (nchar(t) >= 36 & startsWith(d, t))
}

score_counts <- function(got, ref, term_col, count_col) {
  # got: data frame with TRTAN, the term column (or none for "any") and a count column
  ok <- vapply(seq_len(nrow(ref)), function(i) {
    rows <- got[as.numeric(as.character(got$TRTAN)) == ref$trtan[i], , drop = FALSE]
    if (!is.null(term_col)) rows <- rows[term_match(rows[[term_col]], ref$term[i]), , drop = FALSE]
    nrow(rows) >= 1 && sum(as.numeric(as.character(rows[[count_col]]))) == ref$n[i]
  }, logical(1))
  list(cells = nrow(ref), matched = sum(ok),
       not_matched = head(ref[!ok, c("level", "term", "trtan", "n")], 20))
}

compare_ae_pref <- function(out, data_dir) {
  csr <- csr_14_5_01(); csr <- csr[csr$n > 0, ]
  cnt <- function(d) if ("_FREQ_" %in% names(d)) "_FREQ_" else if ("N" %in% names(d)) "N" else NA
  levels <- list(any = list("STAT1", NULL), soc = list("STAT2", "AESOC"), pt = list("STAT3", "AEDECOD"))
  res <- lapply(names(levels), function(lv) {
    d <- pick(out, levels[[lv]][[1]])
    if (is.null(d) || is.na(cnt(d))) return(list(found = FALSE, cells = sum(csr$level == lv), matched = 0))
    c(list(found = TRUE, count_column = cnt(d)), score_counts(d, csr[csr$level == lv, ], levels[[lv]][[2]], cnt(d)))
  })
  names(res) <- names(levels)
  cells <- sum(vapply(res, `[[`, numeric(1), "cells")); matched <- sum(vapply(res, `[[`, numeric(1), "matched"))
  list(reference = "published", levels = res, cells = cells, matched = matched,
       structure = compare_structure("ae_pref", out), pass = matched == cells)
}

compare_ae_common <- function(out, data_dir) {
  csr <- csr_14_5_01(); csr <- csr[csr$level == "pt" & csr$n > 0, ]
  d <- pick(out, "COUNTAE")
  if (is.null(d) || !all(c("TRTAN", "AEDECOD", "NUMAE") %in% names(d))) {
    return(list(reference = "published", pass = FALSE, reason = "no COUNTAE with TRTAN, AEDECOD, NUMAE",
                structure = compare_structure("ae_common", out)))
  }
  s <- score_counts(d, csr, "AEDECOD", "NUMAE")
  c(list(reference = "published"), s, list(structure = compare_structure("ae_common", out), pass = s$matched == s$cells))
}
