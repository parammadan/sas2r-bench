```r
library(haven)
library(dplyr)

adsl <- read_xpt("data/adsl.xpt")
adae <- read_xpt("data/adae.xpt")

sort_by_vars <- function(df, vars) {
  if (nrow(df) == 0) return(df)
  keys <- lapply(vars, function(v) df[[v]])
  o <- do.call(order, c(keys, list(na.last = FALSE)))
  df[o, , drop = FALSE]
}

up <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[is.na(x)] <- ""
  x
}

aestat <- function(adae, byvar, wfun) {
  adae <- sort_by_vars(adae, byvar)
  d <- adae[wfun(adae), , drop = FALSE]
  numv <- names(d)[vapply(d, is.numeric, logical(1))]
  numv <- setdiff(numv, byvar)
  first <- numv[1]
  key <- do.call(paste, c(lapply(byvar, function(v) ifelse(is.na(d[[v]]), "<NA>", as.character(d[[v]]))), sep = "\r"))
  ukey <- unique(key)
  rows <- lapply(ukey, function(k) {
    idx <- which(key == k)
    r <- d[idx[1], byvar, drop = FALSE]
    r$`_TYPE_` <- 0
    r$`_FREQ_` <- length(idx)
    r$n <- sum(!is.na(d[[first]][idx]))
    r
  })
  res <- bind_rows(rows)
  list(res = as.data.frame(res), adae = adae)
}

r1 <- aestat(adae, "TRTAN", function(d) up(d$SAFFL) == "Y" & up(d$AOCCFL) == "Y")
stat1 <- r1$res
r2 <- aestat(r1$adae, c("TRTAN", "AESOC"), function(d) up(d$SAFFL) == "Y" & up(d$AOCCSFL) == "Y")
stat2 <- r2$res
r3 <- aestat(r2$adae, c("TRTAN", "AESOC", "AEDECOD"), function(d) up(d$SAFFL) == "Y" & up(d$AOCCPFL) == "Y")
stat3 <- r3$res

saveRDS(list(ADSL = adsl, ADAE = r3$adae, STAT1 = stat1, STAT2 = stat2, STAT3 = stat3), "out.rds")
```