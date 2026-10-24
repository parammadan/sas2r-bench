```r
library(dplyr)
library(tidyr)
library(haven)
library(stringr)

chr_fix <- function(d) {
  d[] <- lapply(d, function(x) if (is.character(x)) { x[is.na(x)] <- ""; x } else x)
  d
}
adsl_raw <- chr_fix(as.data.frame(read_xpt("data/adsl.xpt")))
names(adsl_raw) <- toupper(names(adsl_raw))

ADSL <- adsl_raw %>%
  filter(TRT01P != "") %>%
  mutate(TRT01PN = case_when(TRT01PN == 0 ~ 3,
                             TRT01PN == 54 ~ 1,
                             TRT01PN == 81 ~ 2,
                             TRUE ~ TRT01PN))

freq <- function(d, vars) {
  ok <- rep(TRUE, nrow(d))
  for (v in vars) {
    x <- d[[v]]
    ok <- ok & !is.na(x) & !(is.character(x) & x == "")
  }
  d <- d[ok, , drop = FALSE]
  d %>% group_by(across(all_of(vars))) %>%
    summarise(COUNT = as.numeric(n()), .groups = "drop") %>%
    arrange(across(all_of(vars))) %>% as.data.frame()
}

TOPLINE <- freq(ADSL, c("STUDYID", "TRT01PN", "TRT01P", "DISCONFL"))

d1 <- ADSL %>% filter(DSRAEFL == "Y" | DTHFL == "Y")
MIDLINE1 <- freq(d1, c("STUDYID", "TRT01PN", "TRT01P", "DCREASCD"))
SORTORD1 <- freq(d1, c("STUDYID", "DCREASCD"))

d2 <- ADSL %>% filter(DCDECOD == "LACK OF EFFICACY")
MIDLINE2 <- freq(d2, c("STUDYID", "TRT01PN", "TRT01P", "DCREASCD"))
SORTORD2 <- freq(d2, c("STUDYID", "DCREASCD"))

d3 <- ADSL %>% filter(DCDECOD != "LACK OF EFFICACY", DSRAEFL != "Y", DTHFL != "Y", DISCONFL == "Y")
MIDLINE3 <- freq(d3, c("STUDYID", "TRT01PN", "TRT01P", "DCREASCD"))
SORTORD3 <- freq(d3, c("STUDYID", "DCREASCD"))

ALL <- bind_rows(
  TOPLINE %>% mutate(ORD1 = 1, ORD2 = 1),
  MIDLINE1 %>% mutate(ORD1 = 2, ORD2 = 1),
  MIDLINE2 %>% mutate(ORD1 = 3, ORD2 = 2),
  MIDLINE3 %>% mutate(ORD1 = 4, ORD2 = 3)
) %>% arrange(STUDYID, TRT01PN, TRT01P) %>% as.data.frame()

ALLS <- bind_rows(
  SORTORD1 %>% mutate(ORD1 = 2, ORD2 = 1),
  SORTORD2 %>% mutate(ORD1 = 3, ORD2 = 2),
  SORTORD3 %>% mutate(ORD1 = 4, ORD2 = 3)
) %>% arrange(STUDYID) %>%
  arrange(ORD1, ORD2, desc(COUNT)) %>% as.data.frame()

ALLS2 <- ALLS %>% group_by(ORD1, ORD2) %>%
  mutate(ORD3 = row_number()) %>% ungroup() %>% as.data.frame()

POP <- freq(ADSL, c("STUDYID", "TRT01PN", "TRT01P")) %>% rename(TOT = COUNT)

rnd <- function(x) floor(x * 10 + 0.5 + 1e-9) / 10
ALLP <- full_join(ALL, POP, by = c("STUDYID", "TRT01PN", "TRT01P")) %>%
  arrange(STUDYID, TRT01PN, TRT01P)
ALLP$DCREASCD <- as.character(ALLP$DCREASCD)
ALLP$DCREASCD[is.na(ALLP$DCREASCD)] <- ""
ALLP$COUNTC <- ifelse(is.na(ALLP$COUNT), NA_character_,
                      paste0(sprintf("%4d", as.integer(ALLP$COUNT)), " (",
                             sprintf("%5.1f", rnd(ALLP$COUNT / ALLP$TOT)), ")"))
ALLP$DCREASCD <- ifelse(!is.na(ALLP$ORD1) & ALLP$ORD1 == 1,
                        ifelse(!is.na(ALLP$DISCONFL) & ALLP$DISCONFL == "Y",
                               "Discontinued", "Completed the Study"),
                        ALLP$DCREASCD)
ALLP <- as.data.frame(ALLP)

ALLP <- ALLP %>% arrange(STUDYID, ORD1, ORD2, DCREASCD, TRT01PN)
ALLPT <- ALLP %>%
  pivot_wider(id_cols = c(STUDYID, ORD1, ORD2, DCREASCD),
              names_from = TRT01PN, names_prefix = "TRT", values_from = COUNTC) %>%
  mutate(`_NAME_` = "COUNTC", .after = DCREASCD) %>% as.data.frame()

maxt <- nrow(POP)
ALLPT2 <- ALLPT
for (i in seq_len(maxt)) {
  v <- paste0("TRT", i)
  if (!v %in% names(ALLPT2)) ALLPT2[[v]] <- NA_character_
  x <- ALLPT2[[v]]
  ALLPT2[[v]] <- ifelse(is.na(x) | x == "", "###0", gsub(" ", "#", x, fixed = TRUE))
}

# positional (SAS-style) merge by STUDYID ORD1 ORD2
A <- ALLPT2 %>% group_by(STUDYID, ORD1, ORD2) %>% mutate(RN = row_number()) %>% ungroup()
B <- ALLS2 %>% select(STUDYID, ORD1, ORD2, DCREASCD_B = DCREASCD, ORD3) %>%
  group_by(STUDYID, ORD1, ORD2) %>% mutate(RN = row_number()) %>% ungroup()
M <- full_join(A, B, by = c("STUDYID", "ORD1", "ORD2", "RN")) %>%
  arrange(STUDYID, ORD1, ORD2, RN)
M$DCREASCD <- ifelse(!is.na(M$DCREASCD_B), M$DCREASCD_B, M$DCREASCD)
M$ORD3 <- ifelse(M$ORD1 == 1 & toupper(M$DCREASCD) == "DISCONTINUED", 2,
          ifelse(M$ORD1 == 1, 1, M$ORD3))
ALLPT3 <- M %>%
  select(STUDYID, DCREASCD, ORD1, ORD2, ORD3, starts_with("TRT")) %>%
  arrange(STUDYID, ORD1, ORD2, ORD3) %>% as.data.frame()

FINAL <- ALLPT3
FINAL$PAGEX <- 1
FINAL$COL1 <- ifelse(FINAL$ORD1 == 1, trimws(FINAL$DCREASCD),
              ifelse(FINAL$ORD1 == 3, paste0("####", trimws(FINAL$DCREASCD)), NA_character_))

saveRDS(list(ADSL = ADSL, TOPLINE = TOPLINE, MIDLINE1 = MIDLINE1, SORTORD1 = SORTORD1,
             MIDLINE2 = MIDLINE2, SORTORD2 = SORTORD2, MIDLINE3 = MIDLINE3,
             SORTORD3 = SORTORD3, ALL = ALL, ALLS = ALLS, ALLS2 = ALLS2, POP = POP,
             ALLP = ALLP, ALLPT = ALLPT, ALLPT2 = ALLPT2, ALLPT3 = ALLPT3, FINAL = FINAL),
        "out.rds")
```