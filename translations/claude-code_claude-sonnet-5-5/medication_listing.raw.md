```r
library(haven)
library(dplyr)

src <- read_xpt("data/adcm.xpt")
for (n in names(src)) if (is.character(src[[n]])) src[[n]][is.na(src[[n]])] <- ""

adcm <- src[!is.na(src$SAFFL) & trimws(src$SAFFL) == "Y", ]
ord <- do.call(order, c(lapply(c("TRTA","STUDYID","USUBJID","CMSTDTC","CMENDTC","CMDECOD"),
                               function(v) {
                                 x <- adcm[[v]]
                                 if (is.character(x)) x <- trimws(x, "right")
                                 x
                               }),
                        list(method = "radix")))
adcm <- adcm[ord, ]
adcm <- as.data.frame(adcm)

compbl <- function(x) gsub(" +", " ", x)
best12 <- function(x) {
  s <- ifelse(is.na(x), ".", as.character(signif(x, 12)))
  formatC(s, width = 12, flag = " ")
}

adcm2 <- adcm
adcm2$Period <- "Study Period"
adcm2$CMSTDTC[trimws(adcm2$CMSTDTC) == ""] <- "NK"
adcm2$CMENDTC[trimws(adcm2$CMENDTC) == ""] <- "NK"
adcm2$StartStop <- substr(compbl(paste0(adcm2$CMSTDTC, "/", adcm2$CMENDTC)), 1, 25)

if (!"CMUNIT" %in% names(adcm2)) adcm2$CMUNIT <- NA_real_
dose <- if (is.character(adcm2$CMDOSE)) adcm2$CMDOSE else best12(adcm2$CMDOSE)
if (is.numeric(adcm2$CMUNIT)) {
  adcm2$DoseUnit <- ifelse(!is.na(adcm2$CMUNIT),
                           compbl(paste0(dose, "(", adcm2$CMUNIT, ")")),
                           compbl(dose))
} else {
  adcm2$DoseUnit <- ifelse(adcm2$CMUNIT != "",
                           compbl(paste0(dose, "(", adcm2$CMUNIT, ")")),
                           compbl(dose))
}
adcm2$DoseUnit <- substr(adcm2$DoseUnit, 1, 50)
adcm2$Ongoing <- ifelse(!is.na(adcm2$AENRF) & trimws(adcm2$AENRF) == "Ongoing", "Y", "N")

saveRDS(list(ADCM = adcm, ADCM2 = adcm2), "out.rds")
```