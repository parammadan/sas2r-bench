# disposition_a.sas -> published reference: CSR Table 14-1.01, "Complete Study" row
# (refs/csr_t14_1_01.csv): 58 / 25 / 27 completed. TOPLINE holds PROC FREQ counts of DISCONFL by
# TRT01PN; DISCONFL is "Y if DCREASCD ^= 'Completed', null otherwise" (define.xml), so the blank
# level is the completers and "Y" is everyone else: 86 - 58, 84 - 25, 84 - 27 = 28 / 59 / 57.

compare_disposition_a <- function(out, data_dir) {
  csr <- read.csv("refs/csr_t14_1_01.csv", stringsAsFactors = FALSE)
  arms <- c(`0` = "Placebo", `54` = "Xanomeline Low Dose", `81` = "Xanomeline High Dose")
  done <- csr$value[csr$row == "Complete Study" & csr$statistic == "count"][match(arms, csr$column[csr$row == "Complete Study" & csr$statistic == "count"])]
  itt <- csr$value[csr$row == "Intent-To-Treat (ITT)" & csr$statistic == "count"][match(arms, csr$column[csr$row == "Intent-To-Treat (ITT)" & csr$statistic == "count"])]
  d <- pick(out, "TOPLINE")
  st <- compare_structure("disposition_a", out)
  if (is.null(d) || !all(c("TRT01PN", "DISCONFL", "COUNT") %in% names(d))) {
    return(list(reference = "published", pass = FALSE, reason = "no TOPLINE with TRT01PN, DISCONFL, COUNT", structure = st))
  }
  fl <- trimws(as.character(d$DISCONFL)); fl[is.na(fl)] <- ""
  checks <- lapply(seq_along(arms), function(j) {
    a <- as.numeric(names(arms)[j])
    comp <- sum(as.numeric(as.character(d$COUNT)[as.numeric(as.character(d$TRT01PN)) == a & fl == ""]))
    disc <- sum(as.numeric(as.character(d$COUNT)[as.numeric(as.character(d$TRT01PN)) == a & fl == "Y"]))
    list(arm = arms[[j]], completed = comp, completed_csr = done[j], discontinued = disc, discontinued_csr = itt[j] - done[j],
         match = comp == done[j] && disc == itt[j] - done[j])
  })
  ok <- vapply(checks, `[[`, logical(1), "match")
  list(reference = "published", arms = checks, matched = sum(ok), cells = length(ok), structure = st, pass = all(ok))
}
