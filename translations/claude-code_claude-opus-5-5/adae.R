suppressPackageStartupMessages({
  library(dplyr)
  library(haven)
})

ae_raw <- haven::read_xpt("data/ae.xpt")
adsl_raw <- haven::read_xpt("data/adsl.xpt")
names(ae_raw) <- toupper(names(ae_raw))
names(adsl_raw) <- toupper(names(adsl_raw))

rtrim <- function(x) sub(" +$", "", x)
to_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  as.Date(as.numeric(x), origin = "1960-01-01")
}
parse_ymd10 <- function(x) {
  x <- ifelse(is.na(x), "", trimws(x))
  s <- substr(x, 1, 10)
  ok <- grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", s)
  out <- rep(as.Date(NA), length(x))
  out[ok] <- as.Date(s[ok], format = "%Y-%m-%d")
  out
}
sas_order <- function(df, vars) {
  if (nrow(df) == 0) return(integer(0))
  do.call(order, c(lapply(vars, function(v) df[[v]]), list(na.last = FALSE, method = "radix")))
}
nz <- function(x) ifelse(is.na(x), "", x)

# ---- Sort AE and derive dates ----
AE <- ae_raw[sas_order(ae_raw, c("USUBJID", "AESEQ")), , drop = FALSE]
AE$ASTDT <- parse_ymd10(AE$AESTDTC)
AE$AENDT <- parse_ymd10(AE$AEENDTC)
AE$ADUR <- ifelse(!is.na(AE$ASTDT) & !is.na(AE$AENDT),
                  as.numeric(AE$AENDT - AE$ASTDT) + 1, NA_real_)

# ---- TEMP from ADSL ----
adsl_s <- adsl_raw[sas_order(adsl_raw, "USUBJID"), , drop = FALSE]
keepv <- c("SITEID","USUBJID","TRT01A","TRT01AN","AGE","AGEGR1","AGEGR1N","RACE","RACEN","SEX","SAFFL","TRTSDT","TRTEDT")
TEMP <- adsl_s[, names(adsl_s) %in% keepv, drop = FALSE]
names(TEMP)[names(TEMP) == "TRT01A"] <- "TRTA"
names(TEMP)[names(TEMP) == "TRT01AN"] <- "TRTAN"
if ("TRTSDT" %in% names(TEMP)) TEMP$TRTSDT <- to_date(TEMP$TRTSDT)
if ("TRTEDT" %in% names(TEMP)) TEMP$TRTEDT <- to_date(TEMP$TRTEDT)

# ---- Merge AE with TEMP (if inAE) ----
idx <- match(AE$USUBJID, TEMP$USUBJID)
tvars <- setdiff(names(TEMP), "USUBJID")
for (v in tvars) {
  tv <- TEMP[[v]][idx]
  if (v %in% names(AE)) {
    AE[[v]] <- ifelse(is.na(idx), AE[[v]], tv)
  } else {
    AE[[v]] <- tv
  }
}

aestdtc <- nz(AE$AESTDTC)
imp <- is.na(AE$ASTDT) & nchar(rtrim(aestdtc)) == 7
ASTDTF <- rep(NA_character_, nrow(AE))
if (any(imp)) {
  s <- gsub(" ", "", paste0(aestdtc[imp], "-01"))
  AE$ASTDT[imp] <- parse_ymd10(s)
  ASTDTF[imp] <- "D"
}
AE$ASTDTF <- ifelse(is.na(ASTDTF), "", ASTDTF)

TRTSDT <- AE$TRTSDT
ASTDT <- AE$ASTDT
AENDT <- AE$AENDT
c1 <- !is.na(TRTSDT) & !is.na(ASTDT) & ASTDT >= TRTSDT
c2 <- !is.na(TRTSDT) & !is.na(ASTDT) & TRTSDT > ASTDT
dd <- as.numeric(ASTDT - TRTSDT)
ASTDY <- ifelse(c1, dd + 1, ifelse(c2, dd, NA_real_))
e1 <- !is.na(TRTSDT) & !is.na(AENDT) & AENDT >= TRTSDT
e2 <- !is.na(TRTSDT) & !is.na(AENDT) & TRTSDT > AENDT
de <- as.numeric(AENDT - TRTSDT)
AENDY <- ifelse(e1, de + 1, ifelse(e2, de, NA_real_))
ANLSTDY <- ifelse(!is.na(ASTDT) & !is.na(TRTSDT), dd + 1, NA_real_)
TRTEMFL <- ifelse(c1, "Y", "N")

aedecod <- nz(AE$AEDECOD)
aebodsys <- nz(AE$AEBODSYS)
dermc <- grepl("APPLICATION", aedecod, fixed = TRUE) |
  grepl("DERMATITIS", aedecod, fixed = TRUE) |
  grepl("ERYTHEMA", aedecod, fixed = TRUE) |
  grepl("BLISTER", aedecod, fixed = TRUE) |
  (rtrim(aebodsys) == "SKIN AND SUBCUTANEOUS TISSUE DISORDERS" &
     !(rtrim(aedecod) %in% c("COLD SWEAT", "HYPERHIDROSIS", "ALOPECIA")))
AE$CQ01NAM <- ifelse(dermc, "DERMATOLOGIC EVENTS", "")
AE$ASTDY <- ASTDY
AE$AENDY <- AENDY
AE$ANLSTDY <- ANLSTDY
AE$TRTEMFL <- TRTEMFL

# ---- flagfrst macro ----
DS01 <- NULL
DS02 <- NULL
flagfrst <- function(inds, sortvars, invar, outvar, ser = FALSE, derm = FALSE) {
  keep <- nz(inds$TRTEMFL) == "Y"
  if (!is.null(invar)) keep <- keep & nz(inds[[invar]]) != ""
  if (ser) keep <- keep & rtrim(nz(inds$AESER)) == "Y"
  if (derm) keep <- keep & rtrim(nz(inds$CQ01NAM)) == "DERMATOLOGIC EVENTS"
  d <- inds[keep, , drop = FALSE]
  byv <- c("USUBJID", sortvars, "ASTDT")
  ds01 <- d[sas_order(d, byv), , drop = FALSE]
  grp <- c("USUBJID", sortvars)
  n <- nrow(ds01)
  if (n > 0) {
    key <- do.call(paste, c(lapply(grp, function(v) as.character(ds01[[v]])), sep = "\r"))
    first <- c(TRUE, key[-1] != key[-n])
  } else first <- logical(0)
  ds02 <- ds01[first, , drop = FALSE]
  ds02[[outvar]] <- rep("Y", nrow(ds02))
  kv <- c(names(ds01)[names(ds01) %in% c("USUBJID", "AESEQ")], outvar)
  ds02 <- ds02[, kv, drop = FALSE]
  out <- ds02[sas_order(ds02, c("USUBJID", "AESEQ")), , drop = FALSE]
  DS01 <<- ds01
  DS02 <<- ds02
  out
}

AEA  <- flagfrst(AE, NULL, NULL, "AOCCFL")
AEB  <- flagfrst(AE, "AEBODSYS", "AEBODSYS", "AOCCSFL")
AEC  <- flagfrst(AE, c("AEBODSYS", "AEDECOD"), "AEDECOD", "AOCCPFL")
SAEA <- flagfrst(AE, NULL, NULL, "AOCC02FL", ser = TRUE)
SAEB <- flagfrst(AE, "AEBODSYS", "AEBODSYS", "AOCC03FL", ser = TRUE)
SAEC <- flagfrst(AE, c("AEBODSYS", "AEDECOD"), "AEDECOD", "AOCC04FL", ser = TRUE)
DERM <- flagfrst(AE, NULL, NULL, "AOCC01FL", derm = TRUE)

# ---- Merge flags ----
AE <- AE[sas_order(AE, c("USUBJID", "AESEQ")), , drop = FALSE]
for (fd in list(AEA, AEB, AEC, SAEA, SAEB, SAEC, DERM)) {
  fv <- setdiff(names(fd), c("USUBJID", "AESEQ"))
  k1 <- paste(AE$USUBJID, AE$AESEQ, sep = "\r")
  k2 <- paste(fd$USUBJID, fd$AESEQ, sep = "\r")
  m <- match(k1, k2)
  for (v in fv) {
    val <- fd[[v]][m]
    AE[[v]] <- ifelse(is.na(val), "", val)
  }
}

# ---- Reorder / keep ----
retain_list <- c("STUDYID","SITEID","USUBJID","TRTA","TRTAN","AGE","AGEGR1","AGEGR1N","RACE","RACEN","SEX","SAFFL","TRTSDT","TRTEDT",
  "ASTDT","ASTDTF","ASTDY","AENDT","AENDY","ADUR","AETERM",
  "AELLT","AELLTCD","AEDECOD","AEPTCD","AEHLT","AEHLTCD","AEHLGT","AEHLGTCD","AEBODSYS","AESOC","AESOCCD","AESEV","AESER","AESCAN","AESCONG","AESDISAB","AESDTH","AESHOSP","AESLIFE","AESOD","AEREL","AEACN",
  "AEOUT","AESEQ","TRTEMFL","AOCCFL","AOCCSFL","AOCCPFL","AOCC02FL","AOCC03FL","AOCC04FL","CQ01NAM","AOCC01FL")
keep_exact <- c("STUDYID","SITEID","USUBJID","TRTAN","TRTA","AESEQ","AETERM","AEDECOD","AEPTCD","AEBODSYS",
  "AESEV","AESER","AEACN","AEREL","AEOUT","AESCAN","AESCONG","AESDISAB","AESDTH","AESHOSP","AESLIFE","AESOD",
  "TRTEMFL","ASTDT","AENDT","ASTDTF","CQ01NAM","ASTDY","AENDY","AOCCFL","AOCCSFL","AOCCPFL",
  "AOCC01FL","AOCC02FL","AOCC03FL","AOCC04FL","AGE","AGEGR1","AGEGR1N","RACE","RACEN","SEX","SAFFL","TRTSDT","TRTEDT","ADUR")
keep_pref <- c("AELLT", "AEHL", "AESOC")
nm <- names(AE)
kept <- nm[nm %in% keep_exact | Reduce(`|`, lapply(keep_pref, function(p) startsWith(nm, p)))]
AEk <- AE[, kept, drop = FALSE]
for (v in setdiff(retain_list, kept)) AEk[[v]] <- rep(NA_real_, nrow(AEk))
final_order <- c(retain_list, setdiff(kept, retain_list))
AE <- AEk[, final_order, drop = FALSE]

ADAE <- AE

saveRDS(list(
  AE = AE,
  TEMP = TEMP,
  DS01 = DS01,
  DS02 = DS02,
  AEA = AEA,
  AEB = AEB,
  AEC = AEC,
  SAEA = SAEA,
  SAEB = SAEB,
  SAEC = SAEC,
  DERM = DERM,
  ADAE = ADAE
), "out.rds")
