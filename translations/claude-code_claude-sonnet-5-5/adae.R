library(haven)
library(dplyr)
library(stringr)

ae <- read_xpt("data/ae.xpt")
adsl <- read_xpt("data/adsl.xpt")

pdate <- function(x) {
  x <- trimws(ifelse(is.na(x), "", x))
  s <- substr(x, 1, 10)
  ok <- grepl("^\\d{4}-\\d{2}-\\d{2}$", s)
  out <- rep(as.Date(NA), length(s))
  out[ok] <- suppressWarnings(as.Date(s[ok], format = "%Y-%m-%d"))
  out
}

ae <- ae[order(ae$USUBJID, ae$AESEQ, method = "radix"), ]
ae$ASTDT <- pdate(ae$AESTDTC)
ae$AENDT <- pdate(ae$AEENDTC)
ae$ADUR <- ifelse(!is.na(ae$ASTDT) & !is.na(ae$AENDT),
                  as.numeric(ae$AENDT) - as.numeric(ae$ASTDT) + 1, NA_real_)

temp <- adsl %>%
  select(SITEID, USUBJID, TRT01A, TRT01AN, AGE, AGEGR1, AGEGR1N, RACE, RACEN, SEX,
         SAFFL, TRTSDT, TRTEDT) %>%
  rename(TRTA = TRT01A, TRTAN = TRT01AN) %>%
  arrange(USUBJID)
temp <- temp[!duplicated(temp$USUBJID), ]

aek <- ae[, setdiff(names(ae), setdiff(names(temp), "USUBJID")), drop = FALSE]
ae <- left_join(aek, temp, by = "USUBJID")
ae$ASTDTF <- rep("", nrow(ae))

s7 <- nchar(trimws(ifelse(is.na(ae$AESTDTC), "", ae$AESTDTC), "right")) == 7
idx <- is.na(ae$ASTDT) & s7
if (any(idx)) {
  d <- as.Date(paste0(trimws(ae$AESTDTC[idx]), "-01"), format = "%Y-%m-%d")
  ae$ASTDT[idx] <- d
  ae$ASTDTF[idx] <- ifelse(is.na(d), "", "D")
}

a <- as.numeric(ae$ASTDT); e <- as.numeric(ae$AENDT); t <- as.numeric(ae$TRTSDT)
ae$ASTDY <- dplyr::case_when(
  !is.na(a) & !is.na(t) & a >= t ~ a - t + 1,
  !is.na(a) & !is.na(t) & t > a ~ a - t,
  TRUE ~ NA_real_)
ae$AENDY <- dplyr::case_when(
  !is.na(e) & !is.na(t) & e >= t ~ e - t + 1,
  !is.na(e) & !is.na(t) & t > e ~ e - t,
  TRUE ~ NA_real_)
ae$ANLSTDY <- ifelse(!is.na(a) & !is.na(t), a - t + 1, NA_real_)
ae$TRTEMFL <- ifelse(!is.na(a) & !is.na(t) & a >= t, "Y", "N")

dec <- ifelse(is.na(ae$AEDECOD), "", ae$AEDECOD)
bod <- ifelse(is.na(ae$AEBODSYS), "", ae$AEBODSYS)
cond <- grepl("APPLICATION", dec, fixed = TRUE) | grepl("DERMATITIS", dec, fixed = TRUE) |
  grepl("ERYTHEMA", dec, fixed = TRUE) | grepl("BLISTER", dec, fixed = TRUE) |
  (bod == "SKIN AND SUBCUTANEOUS TISSUE DISORDERS" &
     !(dec %in% c("COLD SWEAT", "HYPERHIDROSIS", "ALOPECIA")))
ae$CQ01NAM <- ifelse(cond, "DERMATOLOGIC EVENTS", "")

flagfrst <- function(inds, sortvars, invar, outvar, ser = FALSE, derm = FALSE) {
  d <- inds[inds$TRTEMFL == "Y", ]
  if (invar != "") d <- d[!is.na(d[[invar]]) & trimws(d[[invar]]) != "", ]
  if (ser) d <- d[!is.na(d$AESER) & d$AESER == "Y", ]
  if (derm) d <- d[d$CQ01NAM == "DERMATOLOGIC EVENTS", ]
  keys <- c("USUBJID", sortvars)
  o <- do.call(order, c(lapply(keys, function(k) { x <- d[[k]]; x[is.na(x)] <- ""; x }),
                        list(as.numeric(d$ASTDT)), list(method = "radix")))
  d <- d[o, ]
  if (nrow(d) == 0) {
    r <- data.frame(USUBJID = character(), AESEQ = numeric(), stringsAsFactors = FALSE)
    r[[outvar]] <- character()
    return(r)
  }
  grp <- apply(as.data.frame(lapply(keys, function(k) { x <- d[[k]]; x[is.na(x)] <- ""; x })),
               1, paste, collapse = "\r")
  first <- !duplicated(grp)
  r <- d[first, c("USUBJID", "AESEQ")]
  r[[outvar]] <- "Y"
  r <- r[order(r$USUBJID, r$AESEQ, method = "radix"), ]
  as.data.frame(r)
}

aea <- flagfrst(ae, character(0), "", "AOCCFL")
aeb <- flagfrst(ae, "AEBODSYS", "AEBODSYS", "AOCCSFL")
aec <- flagfrst(ae, c("AEBODSYS", "AEDECOD"), "AEDECOD", "AOCCPFL")
saea <- flagfrst(ae, character(0), "", "AOCC02FL", ser = TRUE)
saeb <- flagfrst(ae, "AEBODSYS", "AEBODSYS", "AOCC03FL", ser = TRUE)
saec <- flagfrst(ae, c("AEBODSYS", "AEDECOD"), "AEDECOD", "AOCC04FL", ser = TRUE)
derm <- flagfrst(ae, character(0), "", "AOCC01FL", derm = TRUE)

ae <- ae[order(ae$USUBJID, ae$AESEQ, method = "radix"), ]
for (f in list(aea, aeb, aec, saea, saeb, saec, derm)) {
  ae <- left_join(ae, f, by = c("USUBJID", "AESEQ"))
}
for (v in c("AOCCFL", "AOCCSFL", "AOCCPFL", "AOCC02FL", "AOCC03FL", "AOCC04FL", "AOCC01FL")) {
  ae[[v]][is.na(ae[[v]])] <- ""
}

ord <- c("STUDYID","SITEID","USUBJID","TRTA","TRTAN","AGE","AGEGR1","AGEGR1N","RACE","RACEN","SEX","SAFFL","TRTSDT","TRTEDT",
  "ASTDT","ASTDTF","ASTDY","AENDT","AENDY","ADUR","AETERM","AELLT","AELLTCD","AEDECOD","AEPTCD","AEHLT","AEHLTCD",
  "AEHLGT","AEHLGTCD","AEBODSYS","AESOC","AESOCCD","AESEV","AESER","AESCAN","AESCONG","AESDISAB","AESDTH","AESHOSP",
  "AESLIFE","AESOD","AEREL","AEACN","AEOUT","AESEQ","TRTEMFL","AOCCFL","AOCCSFL","AOCCPFL","AOCC02FL","AOCC03FL",
  "AOCC04FL","CQ01NAM","AOCC01FL")
for (v in ord) if (!(v %in% names(ae))) ae[[v]] <- NA_real_
extra <- setdiff(names(ae)[grepl("^(AELLT|AEHL|AESOC)", names(ae))], ord)
ADAE <- as.data.frame(ae[, c(ord, extra)])
ADAE <- ADAE[order(ADAE$USUBJID, ADAE$AESEQ, method = "radix"), ]
AE <- ADAE

saveRDS(list(AE = AE, ADAE = ADAE, TEMP = temp, AEA = aea, AEB = aeb, AEC = aec,
             SAEA = saea, SAEB = saeb, SAEC = saec, DERM = derm), "out.rds")
