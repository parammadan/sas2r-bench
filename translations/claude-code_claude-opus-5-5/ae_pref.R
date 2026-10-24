library(haven)
library(dplyr)

adsl <- read_xpt("data/adsl.xpt")
adae <- read_xpt("data/adae.xpt")

ADSL <- adsl
ADAE <- adae

sas_sort <- function(df, byvars) {
  keys <- lapply(byvars, function(v) df[[v]])
  ord <- do.call(order, c(keys, list(na.last = FALSE, method = "radix")))
  df[ord, , drop = FALSE]
}

up_strip <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  toupper(trimws(x))
}

ae_stat <- function(df, byvars, flagvar) {
  # PROC SORT (stable) on the full dataset
  df <- sas_sort(df, byvars)
  # WHERE clause
  sub <- df[up_strip(df$SAFFL) == "Y" & up_strip(df[[flagvar]]) == "Y", , drop = FALSE]
  # first numeric analysis variable (no VAR statement: all numeric vars not in BY)
  numvars <- names(df)[vapply(df, is.numeric, logical(1))]
  numvars <- setdiff(numvars, byvars)
  firstvar <- numvars[1]
  if (nrow(sub) == 0) {
    out <- sub[0, byvars, drop = FALSE]
    out$`_TYPE_` <- numeric(0)
    out$`_FREQ_` <- numeric(0)
    out$n <- numeric(0)
  } else {
    # group consecutive rows by BY values (data already sorted)
    keydf <- sub[, byvars, drop = FALSE]
    keystr <- do.call(paste, c(lapply(keydf, function(x) ifelse(is.na(x), "\001NA", as.character(x))), sep = "\002"))
    newgrp <- c(TRUE, keystr[-1] != keystr[-length(keystr)])
    gid <- cumsum(newgrp)
    firstrows <- which(newgrp)
    out <- sub[firstrows, byvars, drop = FALSE]
    out$`_TYPE_` <- 0
    out$`_FREQ_` <- as.numeric(tabulate(gid))
    out$n <- as.numeric(tapply(!is.na(sub[[firstvar]]), gid, sum))
  }
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  list(data = df, stat = out)
}

r1 <- ae_stat(ADAE, c("TRTAN"), "AOCCFL")
ADAE <- r1$data
STAT1 <- r1$stat

r2 <- ae_stat(ADAE, c("TRTAN", "AESOC"), "AOCCSFL")
ADAE <- r2$data
STAT2 <- r2$stat

r3 <- ae_stat(ADAE, c("TRTAN", "AESOC", "AEDECOD"), "AOCCPFL")
ADAE <- r3$data
STAT3 <- r3$stat

saveRDS(list(ADSL = ADSL, ADAE = ADAE, STAT1 = STAT1, STAT2 = STAT2, STAT3 = STAT3), "out.rds")
