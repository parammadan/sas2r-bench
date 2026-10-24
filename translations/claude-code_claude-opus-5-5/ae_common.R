suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(haven)
})

# SAS sort order: missing numeric values come before all numbers
na_first <- function(x) !is.na(x)

# Read inputs
adsl <- read_xpt("data/adsl.xpt")
adsl <- as.data.frame(adsl)

adae_src <- as.data.frame(read_xpt("data/adae.xpt"))
keepv <- c("USUBJID", "TRTAN", "AEDECOD", "TRTEMFL", "AOCCPFL", "SAFFL")
adae <- adae_src[, names(adae_src)[toupper(names(adae_src)) %in% keepv], drop = FALSE]
for (v in names(adae)) if (is.character(adae[[v]])) adae[[v]][is.na(adae[[v]])] <- ""
for (v in names(adsl)) if (is.character(adsl[[v]])) adsl[[v]][is.na(adsl[[v]])] <- ""

# proc sort by aedecod (stable sort)
adae <- adae %>% arrange(AEDECOD)

# adae1
adae1 <- adae %>% filter(TRTEMFL %in% "Y" & AOCCPFL %in% "Y" & SAFFL %in% "Y")

# countae
countae <- adae1 %>%
  group_by(TRTAN, AEDECOD) %>%
  summarise(numAE = sum(USUBJID != ""), .groups = "drop") %>%
  arrange(na_first(TRTAN), TRTAN, AEDECOD) %>%
  select(numAE, AEDECOD, TRTAN) %>%
  as.data.frame()

# nsubj
nsubj <- adsl %>%
  group_by(TRT01AN) %>%
  summarise(numusubjid = n_distinct(USUBJID[USUBJID != ""]), .groups = "drop") %>%
  arrange(na_first(TRT01AN), TRT01AN) %>%
  select(numusubjid, TRT01AN) %>%
  as.data.frame()

# count: merge by trtan
ns_r <- nsubj %>% rename(TRTAN = TRT01AN)
if (any(!(countae$TRTAN %in% ns_r$TRTAN))) stop("ABORT CANCEL: trtan not in nsubj")
count <- countae %>%
  mutate(.src = 1L, .ord = row_number()) %>%
  full_join(ns_r, by = "TRTAN") %>%
  mutate(AEDECOD = ifelse(is.na(AEDECOD), "", AEDECOD),
         .ord = ifelse(is.na(.ord), 0L, .ord)) %>%
  arrange(na_first(TRTAN), TRTAN, .ord) %>%
  mutate(AEpercent = ifelse(!is.na(numAE) & !is.na(numusubjid), numAE / numusubjid, NA_real_)) %>%
  select(numAE, AEDECOD, TRTAN, numusubjid, AEpercent) %>%
  as.data.frame()

# sort count by aedecod (stable)
count <- count %>% arrange(AEDECOD)

# tcount: merge count with count(where=(trtan_1=0) rename=...)
a <- count %>% group_by(AEDECOD) %>% mutate(.i = row_number()) %>% ungroup()
b <- count %>%
  rename(AEpercent_1 = AEpercent, trtan_1 = TRTAN) %>%
  filter(!is.na(trtan_1) & trtan_1 == 0) %>%
  group_by(AEDECOD) %>% mutate(.i = row_number()) %>% ungroup() %>%
  rename(numAE_b = numAE, numusubjid_b = numusubjid)
na_n <- a %>% count(AEDECOD, name = ".na")
nb_n <- b %>% count(AEDECOD, name = ".nb")

tc <- full_join(a %>% rename(numAE_a = numAE, numusubjid_a = numusubjid),
                b, by = c("AEDECOD", ".i")) %>%
  left_join(na_n, by = "AEDECOD") %>%
  left_join(nb_n, by = "AEDECOD") %>%
  mutate(.na = ifelse(is.na(.na), 0L, .na), .nb = ifelse(is.na(.nb), 0L, .nb)) %>%
  arrange(AEDECOD, .i) %>%
  group_by(AEDECOD) %>%
  fill(TRTAN, AEpercent, numAE_a, numusubjid_a, .direction = "down") %>%
  fill(AEpercent_1, trtan_1, .direction = "down") %>%
  ungroup() %>%
  mutate(numAE = ifelse(.i <= .nb, numAE_b, numAE_a),
         numusubjid = ifelse(.i <= .nb, numusubjid_b, numusubjid_a))
# columns from b when b exhausted: retained via fill only if group has b rows
tc <- tc %>%
  mutate(AEpercent_1 = ifelse(.nb == 0, NA_real_, AEpercent_1),
         trtan_1 = ifelse(.nb == 0, NA_real_, trtan_1))

tcount <- tc %>%
  mutate(AEpercent = ifelse(is.na(AEpercent), 0, AEpercent),
         AEpercent_1 = ifelse(is.na(AEpercent_1), 0, AEpercent_1),
         RelRisk = AEpercent - AEpercent_1,
         aeseqno = as.numeric(row_number())) %>%
  select(numAE, AEDECOD, TRTAN, numusubjid, AEpercent, AEpercent_1, trtan_1, RelRisk, aeseqno) %>%
  as.data.frame()

# maxpercen2
maxpercen2 <- tcount %>%
  group_by(AEDECOD) %>%
  summarise(minusmaxAEpercent = -max(AEpercent), .groups = "drop") %>%
  arrange(minusmaxAEpercent, AEDECOD) %>%
  select(minusmaxAEpercent, AEDECOD) %>%
  mutate(dispseqno = as.numeric(row_number())) %>%
  arrange(AEDECOD) %>%
  as.data.frame()

# tcount2
tcount2 <- tcount %>%
  mutate(.o = row_number()) %>%
  left_join(maxpercen2, by = "AEDECOD") %>%
  arrange(.o) %>% select(-.o) %>%
  as.data.frame()

# AEgrid
trts <- tcount2 %>% distinct(TRTAN)
aeds <- tcount2 %>% distinct(AEDECOD, dispseqno)
AEgrid <- tidyr::crossing(trts, aeds) %>%
  arrange(AEDECOD, na_first(TRTAN), TRTAN) %>%
  select(TRTAN, AEDECOD, dispseqno) %>%
  as.data.frame()

# tcount3
tcount3 <- full_join(tcount2, AEgrid %>% rename(dispseqno_g = dispseqno),
                     by = c("AEDECOD", "TRTAN")) %>%
  mutate(dispseqno = ifelse(!is.na(dispseqno_g), dispseqno_g, dispseqno)) %>%
  select(-dispseqno_g) %>%
  arrange(AEDECOD, na_first(TRTAN), TRTAN) %>%
  mutate(AEpercent = ifelse(is.na(AEpercent), 0, AEpercent)) %>%
  select(names(tcount2)) %>%
  arrange(na_first(dispseqno), dispseqno) %>%
  as.data.frame()

saveRDS(list(
  ADSL = adsl,
  ADAE = adae,
  ADAE1 = adae1,
  COUNTAE = countae,
  NSUBJ = nsubj,
  COUNT = count,
  TCOUNT = tcount,
  MAXPERCEN2 = maxpercen2,
  TCOUNT2 = tcount2,
  AEGRID = AEgrid,
  TCOUNT3 = tcount3
), "out.rds")
