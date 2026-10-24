library(haven)
library(dplyr)
library(stringr)

SAS_ORIGIN <- as.Date("1960-01-01")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Strip labels; right-trim character values; blank character = missing (SAS)
clean_df <- function(df) {
  df <- zap_label(zap_labels(df))
  df %>% mutate(across(where(is.character), function(x) {
    x <- str_trim(x, side = "right")
    x[!is.na(x) & x == ""] <- NA_character_
    x
  }))
}

# Date variables read as Date class -> SAS date (days since 1960-01-01)
to_sas_num <- function(x) {
  if (inherits(x, "Date")) as.numeric(x - SAS_ORIGIN) else as.numeric(x)
}

# SAS input(x, ?? yymmdd10.) : reads first 10 columns of the stripped string
sas_input_iso <- function(x) {
  x <- substr(str_trim(x, side = "both"), 1, 10)
  ok <- !is.na(x) & str_detect(x, "^\\d{4}\\D\\d{2}\\D\\d{2}$")
  out <- rep(NA_real_, length(x))
  if (any(ok)) {
    d <- as.Date(str_replace_all(x[ok], "\\D", "-"), format = "%Y-%m-%d")
    out[ok] <- as.numeric(d - SAS_ORIGIN)
  }
  out
}

# SAS numeric comparisons: missing is smaller than any number; result is 1/0
sas_num <- function(x) ifelse(is.na(x), -Inf, as.numeric(x))
sas_ge  <- function(a, b) as.numeric(sas_num(a) >= sas_num(b))
sas_gt  <- function(a, b) as.numeric(sas_num(a) >  sas_num(b))
# "expr > .z" : true for any non-missing value (0/1 results always qualify)
gt_dotz <- function(v) v > -Inf

# SAS DATA step MERGE by BY variables (no FIRST./LAST.). For each BY group the
# output has max(nx, ny) rows; row k pairs with obs min(k, n) of each dataset,
# so the last obs is retained when a dataset runs out. Missing when group absent.
sas_merge <- function(x, y, by, drop_y_only = FALSE) {
  key_of <- function(d) {
    do.call(paste, c(lapply(by, function(v) as.character(d[[v]])), sep = "\r"))
  }
  xk <- key_of(x)
  yk <- key_of(y)
  xs <- ave(seq_along(xk), xk, FUN = seq_along)
  ys <- ave(seq_along(yk), yk, FUN = seq_along)

  kx <- !duplicated(xk)
  ky <- !duplicated(yk) & !(yk %in% xk)
  keys <- bind_rows(x[kx, by, drop = FALSE], y[ky, by, drop = FALSE])
  keyv <- c(xk[kx], yk[ky])

  nx <- tabulate(match(xk, keyv), nbins = length(keyv))
  ny <- tabulate(match(yk, keyv), nbins = length(keyv))

  ord <- do.call(order, c(unname(as.list(keys)), list(method = "radix", na.last = FALSE)))
  keys <- keys[ord, , drop = FALSE]
  keyv <- keyv[ord]
  nx <- nx[ord]
  ny <- ny[ord]

  keep_g <- if (drop_y_only) nx > 0 else rep(TRUE, length(nx))
  nmax <- ifelse(keep_g, pmax(nx, ny), 0L)
  gi <- rep(seq_along(nmax), nmax)
  k <- sequence(nmax)
  kk <- keyv[gi]

  xrow <- match(paste(kk, pmin(k, nx[gi]), sep = "\u0001"),
                paste(xk, xs, sep = "\u0001"))
  yrow <- match(paste(kk, pmin(k, ny[gi]), sep = "\u0001"),
                paste(yk, ys, sep = "\u0001"))

  out <- x[xrow, , drop = FALSE]
  yv <- setdiff(names(y), by)
  yd <- y[yrow, yv, drop = FALSE]
  for (v in yv) {
    if (v %in% names(out)) {
      out[[v]] <- coalesce(yd[[v]], out[[v]])
    } else {
      out[[v]] <- yd[[v]]
    }
  }
  for (v in by) out[[v]] <- keys[[v]][gi]
  out
}

# Macro %flagfrst
flagfrst <- function(inds, sortvars = character(0), invar = NULL, outvar,
                     ser = FALSE, derm = FALSE) {
  w <- inds$TRTEMFL %in% "Y"
  if (!is.null(invar)) {
    w <- w & !is.na(inds[[invar]]) & inds[[invar]] != ""
  }
  if (ser)  w <- w & (inds$AESER %in% "Y")
  if (derm) w <- w & (inds$CQ01NAM %in% "DERMATOLOGIC EVENTS")

  ds01 <- inds[w, ]
  ord <- do.call(order, c(unname(as.list(ds01[c("USUBJID", sortvars, "ASTDT")])),
                          list(method = "radix", na.last = FALSE)))
  ds01 <- ds01[ord, ]

  bycols <- if (is.null(invar)) {
    "USUBJID"
  } else {
    c("USUBJID", sortvars[seq_len(match(invar, sortvars))])
  }
  first <- !duplicated(ds01[bycols])

  ds02 <- ds01[first, c("USUBJID", "AESEQ")]
  ds02[[outvar]] <- rep("Y", nrow(ds02))

  out <- ds02[order(ds02$USUBJID, ds02$AESEQ, method = "radix", na.last = FALSE), ]
  list(ds01 = ds01, ds02 = ds02, out = out)
}

# ---------------------------------------------------------------------------
# Step 1: proc sort sdtmnew.AE -> AE by USUBJID AESEQ
# ---------------------------------------------------------------------------
AE <- clean_df(read_xpt("data/ae.xpt"))
AE <- AE[order(AE$USUBJID, AE$AESEQ, method = "radix", na.last = FALSE), ]

# ---------------------------------------------------------------------------
# Step 2: ASTDT, AENDT, ADUR
# ---------------------------------------------------------------------------
AE$ASTDT <- sas_input_iso(AE$AESTDTC)
AE$AENDT <- sas_input_iso(AE$AEENDTC)
both <- !is.na(AE$ASTDT) & !is.na(AE$AENDT)
AE$ADUR <- ifelse(both, AE$AENDT - AE$ASTDT + 1, NA_real_)

# ---------------------------------------------------------------------------
# proc sort adamnew.ADSL -> TEMP (keep/rename)
# ---------------------------------------------------------------------------
adsl <- clean_df(read_xpt("data/adsl.xpt"))
TEMP <- adsl %>%
  mutate(TRTSDT = to_sas_num(TRTSDT), TRTEDT = to_sas_num(TRTEDT)) %>%
  select(SITEID, USUBJID, TRTA = TRT01A, TRTAN = TRT01AN, AGE, AGEGR1, AGEGR1N,
         RACE, RACEN, SEX, SAFFL, TRTSDT, TRTEDT)
TEMP <- TEMP[order(TEMP$USUBJID, method = "radix", na.last = FALSE), ]

# ---------------------------------------------------------------------------
# Merge AE(in=inAE) TEMP; by USUBJID; if inAE; plus derivations
# ---------------------------------------------------------------------------
AE <- sas_merge(AE, TEMP, "USUBJID", drop_y_only = TRUE)

# ASTDT imputation: if ASTDT=. and length(AESTDTC)=7 then YYYY-MM-01
fix <- is.na(AE$ASTDT) & !is.na(AE$AESTDTC) &
  nchar(str_trim(AE$AESTDTC, side = "right")) == 7
AE$ASTDTF <- ifelse(fix, "D", NA_character_)
AE$ASTDT[fix] <- sas_input_iso(str_replace_all(paste0(AE$AESTDTC[fix], "-01"), " ", ""))

# ASTDY: SAS evaluates "ASTDT>=TRTSDT>.z" left to right as (ASTDT>=TRTSDT)>.z,
# which is always true, so the ELSE IF branch is never reached.
c1 <- gt_dotz(sas_ge(AE$ASTDT, AE$TRTSDT))
c2 <- !c1 & gt_dotz(sas_gt(AE$TRTSDT, AE$ASTDT))
AE$ASTDY <- ifelse(c1, AE$ASTDT - AE$TRTSDT + 1,
                   ifelse(c2, AE$ASTDT - AE$TRTSDT, NA_real_))

# AENDY: same left-to-right comparison quirk
d1 <- gt_dotz(sas_ge(AE$AENDT, AE$TRTSDT))
d2 <- !d1 & gt_dotz(sas_gt(AE$TRTSDT, AE$AENDT))
AE$AENDY <- ifelse(d1, AE$AENDT - AE$TRTSDT + 1,
                   ifelse(d2, AE$AENDT - AE$TRTSDT, NA_real_))

# TRTEMFL: IF (ASTDT>=TRTSDT>.z) -> always 'Y' under SAS left-to-right evaluation
AE$TRTEMFL <- ifelse(gt_dotz(sas_ge(AE$ASTDT, AE$TRTSDT)), "Y", "N")

# CQ01NAM (index() is case-sensitive; blank/missing gives 0)
derm_cond <- grepl("APPLICATION", AE$AEDECOD, fixed = TRUE) |
  grepl("DERMATITIS", AE$AEDECOD, fixed = TRUE) |
  grepl("ERYTHEMA", AE$AEDECOD, fixed = TRUE) |
  grepl("BLISTER", AE$AEDECOD, fixed = TRUE) |
  (AE$AEBODSYS %in% "SKIN AND SUBCUTANEOUS TISSUE DISORDERS" &
     !(AE$AEDECOD %in% c("COLD SWEAT", "HYPERHIDROSIS", "ALOPECIA")))
AE$CQ01NAM <- ifelse(derm_cond, "DERMATOLOGIC EVENTS", NA_character_)

# ---------------------------------------------------------------------------
# Flag macro calls (each call overwrites ds01/ds02 in WORK; last call persists)
# ---------------------------------------------------------------------------
AE <- AE[order(AE$USUBJID, AE$AESEQ, method = "radix", na.last = FALSE), ]

r <- flagfrst(AE, character(0), NULL, "AOCCFL");                       AEA  <- r$out
r <- flagfrst(AE, "AEBODSYS", "AEBODSYS", "AOCCSFL");                  AEB  <- r$out
r <- flagfrst(AE, c("AEBODSYS", "AEDECOD"), "AEDECOD", "AOCCPFL");     AEC  <- r$out
r <- flagfrst(AE, character(0), NULL, "AOCC02FL", ser = TRUE);        SAEA <- r$out
r <- flagfrst(AE, "AEBODSYS", "AEBODSYS", "AOCC03FL", ser = TRUE);    SAEB <- r$out
r <- flagfrst(AE, c("AEBODSYS", "AEDECOD"), "AEDECOD", "AOCC04FL", ser = TRUE); SAEC <- r$out
r <- flagfrst(AE, character(0), NULL, "AOCC01FL", derm = TRUE);       DERM <- r$out
DS01 <- r$ds01
DS02 <- r$ds02

# ---------------------------------------------------------------------------
# proc sort AE by usubjid aeseq; data AE; merge AE AEA AEB AEC SAEA SAEB SAEC DERM
# ---------------------------------------------------------------------------
AE <- AE[order(AE$USUBJID, AE$AESEQ, method = "radix", na.last = FALSE), ]
bykey <- c("USUBJID", "AESEQ")
AE <- sas_merge(AE, AEA,  bykey)
AE <- sas_merge(AE, AEB,  bykey)
AE <- sas_merge(AE, AEC,  bykey)
AE <- sas_merge(AE, SAEA, bykey)
AE <- sas_merge(AE, SAEB, bykey)
AE <- sas_merge(AE, SAEC, bykey)
AE <- sas_merge(AE, DERM, bykey)

# ---------------------------------------------------------------------------
# data AE; retain ...; set AE(keep=...); by usubjid aeseq;
# ---------------------------------------------------------------------------
retain_vars <- c("STUDYID", "SITEID", "USUBJID", "TRTA", "TRTAN", "AGE", "AGEGR1",
                 "AGEGR1N", "RACE", "RACEN", "SEX", "SAFFL", "TRTSDT", "TRTEDT",
                 "ASTDT", "ASTDTF", "ASTDY", "AENDT", "AENDY", "ADUR", "AETERM",
                 "AELLT", "AELLTCD", "AEDECOD", "AEPTCD", "AEHLT", "AEHLTCD",
                 "AEHLGT", "AEHLGTCD", "AEBODSYS", "AESOC", "AESOCCD", "AESEV",
                 "AESER", "AESCAN", "AESCONG", "AESDISAB", "AESDTH", "AESHOSP",
                 "AESLIFE", "AESOD", "AEREL", "AEACN", "AEOUT", "AESEQ", "TRTEMFL",
                 "AOCCFL", "AOCCSFL", "AOCCPFL", "AOCC02FL", "AOCC03FL", "AOCC04FL",
                 "CQ01NAM", "AOCC01FL")

# RETAIN creates any missing variable as numeric missing
for (v in setdiff(retain_vars, names(AE))) {
  AE[[v]] <- rep(NA_real_, nrow(AE))
}
AE <- AE[, retain_vars]

# ---------------------------------------------------------------------------
# data adamnew.ADAE
# ---------------------------------------------------------------------------
ADAE <- AE

# endsas; -> nothing after this point is executed

saveRDS(list(
  AE   = AE,
  TEMP = TEMP,
  AEA  = AEA,
  AEB  = AEB,
  AEC  = AEC,
  SAEA = SAEA,
  SAEB = SAEB,
  SAEC = SAEC,
  DERM = DERM,
  DS01 = DS01,
  DS02 = DS02,
  ADAE = ADAE
), "out.rds")
