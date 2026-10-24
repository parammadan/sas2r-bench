```r
library(haven)
library(dplyr)

ord <- function(df, ...) {
  keys <- list(...)
  df[do.call(order, c(keys, list(method = "radix"))), , drop = FALSE]
}

adsl <- read_xpt("data/adsl.xpt")
names(adsl) <- toupper(names(adsl))
adae <- read_xpt("data/adae.xpt")
names(adae) <- toupper(names(adae))
adae <- adae[, c("USUBJID", "TRTAN", "AEDECOD", "TRTEMFL", "AOCCPFL", "SAFFL")]
adae <- ord(adae, adae$AEDECOD)

adae1 <- adae %>% filter(TRTEMFL == "Y" & AOCCPFL == "Y" & SAFFL == "Y")

countae <- adae1 %>%
  group_by(TRTAN, AEDECOD) %>%
  summarise(NUMAE = sum(!is.na(USUBJID)), .groups = "drop") %>%
  select(NUMAE, AEDECOD, TRTAN) %>%
  as.data.frame()
countae <- ord(countae, countae$TRTAN, countae$AEDECOD)

nsubj <- adsl %>%
  group_by(TRT01AN) %>%
  summarise(NUMUSUBJID = n_distinct(USUBJID[!is.na(USUBJID)]), .groups = "drop") %>%
  select(NUMUSUBJID, TRT01AN) %>%
  as.data.frame()
nsubj <- ord(nsubj, nsubj$TRT01AN)

count <- full_join(countae, rename(nsubj, TRTAN = TRT01AN), by = "TRTAN")
count <- ord(count, count$TRTAN)
count$AEPERCENT <- ifelse(!is.na(count$NUMAE) & !is.na(count$NUMUSUBJID),
                          count$NUMAE / count$NUMUSUBJID, NA_real_)
count <- ord(count, count$AEDECOD)
rownames(count) <- NULL

plc <- count %>% filter(TRTAN == 0) %>%
  rename(AEPERCENT_1 = AEPERCENT, TRTAN_1 = TRTAN, NUMAE_P = NUMAE, NUMUSUBJID_P = NUMUSUBJID)
tcount <- left_join(count, plc, by = "AEDECOD")
first <- !duplicated(tcount$AEDECOD) & !is.na(tcount$TRTAN_1)
tcount$NUMAE[first] <- tcount$NUMAE_P[first]
tcount$NUMUSUBJID[first] <- tcount$NUMUSUBJID_P[first]
tcount$NUMAE_P <- NULL
tcount$NUMUSUBJID_P <- NULL
tcount$AEPERCENT[is.na(tcount$AEPERCENT)] <- 0
tcount$AEPERCENT_1[is.na(tcount$AEPERCENT_1)] <- 0
tcount$RELRISK <- tcount$AEPERCENT - tcount$AEPERCENT_1
tcount$AESEQNO <- seq_len(nrow(tcount))
tcount <- tcount[, c("NUMAE", "AEDECOD", "TRTAN", "NUMUSUBJID", "AEPERCENT",
                     "TRTAN_1", "AEPERCENT_1", "RELRISK", "AESEQNO")]

maxpercen2 <- tcount %>%
  group_by(AEDECOD) %>%
  summarise(MINUSMAXAEPERCENT = -max(AEPERCENT), .groups = "drop") %>%
  select(MINUSMAXAEPERCENT, AEDECOD) %>%
  as.data.frame()
maxpercen2 <- ord(maxpercen2, maxpercen2$MINUSMAXAEPERCENT)
maxpercen2$DISPSEQNO <- seq_len(nrow(maxpercen2))
maxpercen2 <- ord(maxpercen2, maxpercen2$AEDECOD)
rownames(maxpercen2) <- NULL

tcount2 <- left_join(tcount, maxpercen2, by = "AEDECOD")

AEgrid <- merge(data.frame(TRTAN = sort(unique(tcount2$TRTAN), na.last = FALSE)),
                unique(tcount2[, c("AEDECOD", "DISPSEQNO")]), by = NULL)
AEgrid <- AEgrid[, c("TRTAN", "AEDECOD", "DISPSEQNO")]
AEgrid <- ord(AEgrid, AEgrid$AEDECOD, AEgrid$TRTAN)
rownames(AEgrid) <- NULL

tcount3 <- full_join(tcount2, AEgrid, by = c("AEDECOD", "TRTAN"), suffix = c("", ".g"))
tcount3$DISPSEQNO <- ifelse(is.na(tcount3$DISPSEQNO.g), tcount3$DISPSEQNO, tcount3$DISPSEQNO.g)
tcount3$DISPSEQNO.g <- NULL
tcount3 <- ord(tcount3, tcount3$AEDECOD, tcount3$TRTAN)
tcount3$AEPERCENT[is.na(tcount3$AEPERCENT)] <- 0
tcount3 <- ord(tcount3, tcount3$DISPSEQNO)
rownames(tcount3) <- NULL

saveRDS(list(ADSL = adsl, ADAE = adae, ADAE1 = adae1, COUNTAE = countae,
             NSUBJ = nsubj, COUNT = count, TCOUNT = tcount,
             MAXPERCEN2 = maxpercen2, TCOUNT2 = tcount2,
             AEGRID = AEgrid, TCOUNT3 = tcount3), "out.rds")
```