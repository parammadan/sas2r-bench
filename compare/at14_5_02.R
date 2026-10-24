# at14-5-02.sas -> published reference: CSR Table 14-5.02.
# The program builds RESULTSUB (any body system), RESULTSOC and RESULTPT with COUNT (subjects),
# EVENTS (records) and DENOM, and runs Fisher tests into F1*/F2* datasets via ODS OUTPUT.
# Counts are checked against the CSR. P-values are reported separately: the program's WHERE
# clauses use TRTAN in (0,1) and (0,2), but the data codes arms as 0/54/81, so what SAS itself
# outputs for them can only be settled by a real SAS run.

compare_at14_5_02 <- function(out, data_dir) {
  csr <- read.csv("refs/csr_t14_5_02.csv", stringsAsFactors = FALSE)
  levels <- list(any = "RESULTSUB", soc = "RESULTSOC", pt = "RESULTPT")

  checks <- list()
  for (lv in names(levels)) {
    d <- pick(out, levels[[lv]])
    want <- csr[csr$level == lv & csr$n > 0, ]
    if (is.null(d)) {
      checks[[lv]] <- list(found = FALSE, expected_cells = nrow(want), matched = 0)
      next
    }
    term <- switch(lv, any = NULL, soc = "AEBODSYS", pt = "AEDECOD")
    ok <- vapply(seq_len(nrow(want)), function(i) {
      r <- d[as.numeric(d$TRTAN) == want$trtan[i], , drop = FALSE]
      if (!is.null(term)) r <- r[trimws(r[[term]]) == want$term[i], , drop = FALSE]
      nrow(r) == 1 && as.numeric(r$COUNT) == want$n[i] && as.numeric(r$EVENTS) == want$events[i]
    }, logical(1))
    checks[[lv]] <- list(found = TRUE, expected_cells = nrow(want), matched = sum(ok),
                         rows_got = nrow(d))
  }

  pnames <- c("F1PT", "F2PT", "F1SOC", "F2SOC", "F1SUB", "F2SUB")
  pvals <- lapply(pnames, function(n) {
    d <- pick(out, n)
    if (is.null(d) || nrow(d) == 0) return(NULL)
    v <- intersect(c("NVALUE1", "CVALUE1"), names(d))[1]
    if (is.na(v)) return(NULL)
    suppressWarnings(as.numeric(d[[v]]))
  })
  names(pvals) <- pnames
  found <- unlist(pvals)

  expected <- sum(vapply(checks, `[[`, numeric(1), "expected_cells"))
  matched <- sum(vapply(checks, `[[`, numeric(1), "matched"))
  list(
    reference = "published",
    count_cells_expected = expected,
    count_cells_matched = matched,
    levels = checks,
    p_values_found = length(found[!is.na(found)]),
    p_value_datasets = names(Filter(Negate(is.null), pvals)),
    p_values = lapply(Filter(Negate(is.null), pvals), function(x) round(x, 4)),
    pass = matched == expected
  )
}
