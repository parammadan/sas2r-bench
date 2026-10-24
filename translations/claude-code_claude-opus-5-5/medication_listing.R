library(haven)
library(dplyr)

adcm_src <- haven::read_xpt("data/adcm.xpt")

# PROC SORT with WHERE SAFFL='Y'
ADCM <- adcm_src %>%
  filter(!is.na(SAFFL) & sub(" +$", "", SAFFL) == "Y") %>%
  arrange(TRTA, STUDYID, USUBJID, CMSTDTC, CMENDTC, CMDECOD)

# helpers
char_width <- function(x) {
  w <- attr(x, "width")
  if (is.null(w) || is.na(w)) {
    w <- suppressWarnings(max(nchar(x), na.rm = TRUE))
    if (!is.finite(w)) w <- 1
  }
  as.integer(w)
}
pad_to <- function(x, w) {
  x[is.na(x)] <- ""
  x <- substr(x, 1, w)
  n <- nchar(x)
  paste0(x, strrep(" ", pmax(0, w - n)))
}
compbl <- function(x) gsub(" +", " ", x)
rtrim <- function(x) sub(" +$", "", x)
best12 <- function(v) {
  sapply(v, function(x) {
    if (is.na(x)) return(".")
    if (x == round(x) && abs(x) < 1e12) return(format(x, scientific = FALSE, trim = TRUE))
    for (d in 11:0) {
      s <- formatC(x, format = "f", digits = d)
      if (grepl("\\.", s)) s <- sub("0+$", "", s)
      s <- sub("\\.$", "", s)
      if (nchar(s) <= 12) return(s)
    }
    formatC(x, format = "E", digits = 5)
  }, USE.NAMES = FALSE)
}
num_to_char12 <- function(v) {
  s <- best12(v)
  paste0(strrep(" ", pmax(0, 12 - nchar(s))), s)
}

w_st <- char_width(ADCM$CMSTDTC)
w_en <- char_width(ADCM$CMENDTC)

ADCM2 <- ADCM
cmst <- ADCM2$CMSTDTC; cmst[is.na(cmst)] <- ""
cmen <- ADCM2$CMENDTC; cmen[is.na(cmen)] <- ""
cmst[rtrim(cmst) == ""] <- substr("NK", 1, w_st)
cmen[rtrim(cmen) == ""] <- substr("NK", 1, w_en)
ADCM2$CMSTDTC <- cmst
ADCM2$CMENDTC <- cmen

ss <- compbl(paste0(pad_to(cmst, w_st), "/", pad_to(cmen, w_en)))
ss <- substr(rtrim(ss), 1, 25)
ADCM2$StartStop <- ss

# CMDOSE as character per SAS automatic conversion
if (is.numeric(ADCM2$CMDOSE)) {
  dose_chr <- num_to_char12(ADCM2$CMDOSE)
} else {
  dose_chr <- pad_to(ADCM2$CMDOSE, char_width(ADCM2$CMDOSE))
}

has_cmunit <- "cmunit" %in% tolower(names(ADCM2))
if (has_cmunit) {
  unm <- names(ADCM2)[tolower(names(ADCM2)) == "cmunit"][1]
  u <- ADCM2[[unm]]
  if (is.numeric(u)) {
    cond <- !is.na(u)
    u_chr <- num_to_char12(u)
  } else {
    u2 <- u; u2[is.na(u2)] <- ""
    cond <- rtrim(u2) != "" & rtrim(u2) != "."
    u_chr <- pad_to(u2, char_width(u))
  }
  du <- ifelse(cond,
               compbl(paste0(dose_chr, "(", u_chr, ")")),
               compbl(dose_chr))
} else {
  du <- compbl(dose_chr)
}
du <- substr(rtrim(du), 1, 50)
ADCM2$Doseunit <- du
ADCM2$Period <- "Study Period"
if (!has_cmunit) ADCM2$cmunit <- NA_real_

aenrf <- ADCM2$AENRF; aenrf[is.na(aenrf)] <- ""
ADCM2$Ongoing <- ifelse(rtrim(aenrf) == "Ongoing", "Y", "N")

saveRDS(list(ADCM = ADCM, ADCM2 = ADCM2), "out.rds")
