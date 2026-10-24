library(haven)
library(dplyr)

src <- read_xpt("data/advs.xpt")
names(src) <- toupper(names(src))

advs <- src %>%
  filter(PARAMCD == "DIABP", ATPTN == 815, ANL01FL == "Y",
         !is.na(AVISITN), AVISITN < 99, SAFFL == "Y")
keepv <- c("USUBJID","TRTPN","ANL01FL","PARAM","PARAMCD","AVISIT","AVISITN","ATPTN","SAFFL","AVAL")
advs <- advs[, intersect(names(src), keepv)]
advs <- as.data.frame(advs)

# sort by trtpn avisitn (stable)
advs <- advs[order(advs$TRTPN, advs$AVISITN), ]
rownames(advs) <- NULL

temp <- advs %>%
  group_by(TRTPN, AVISITN, AVISIT) %>%
  summarise(
    N = sum(!is.na(AVAL)),
    MEAN = if (sum(!is.na(AVAL)) > 0) mean(AVAL, na.rm = TRUE) else NA_real_,
    STDERR = if (sum(!is.na(AVAL)) > 1) sd(AVAL, na.rm = TRUE) / sqrt(sum(!is.na(AVAL))) else NA_real_,
    .groups = "drop"
  ) %>%
  arrange(TRTPN, AVISITN) %>%
  as.data.frame()

temp1 <- temp %>%
  mutate(LO = MEAN - qt(0.95, N - 1) * STDERR,
         HI = MEAN + qt(0.95, N - 1) * STDERR) %>%
  select(-STDERR)

A <- temp1 %>% filter(TRTPN == 0)
B <- temp1 %>% filter(TRTPN == 54)
C <- temp1 %>% filter(TRTPN == 81)

A1 <- A %>% rename(A_N = N, A_MEAN = MEAN, A_LCL = LO, A_UCL = HI)
B1 <- B %>% rename(B_N = N, B_MEAN = MEAN, B_LCL = LO, B_UCL = HI)
C1 <- C %>% rename(C_N = N, C_MEAN = MEAN, C_LCL = LO, C_UCL = HI)

ALL <- A1 %>% select(-TRTPN) %>%
  full_join(B1 %>% select(-TRTPN), by = c("AVISITN", "AVISIT")) %>%
  full_join(C1 %>% select(-TRTPN), by = c("AVISITN", "AVISIT")) %>%
  arrange(AVISITN, AVISIT) %>%
  rename(XC = AVISITN, X = AVISIT) %>%
  as.data.frame()

saveRDS(list(ADVS = advs, TEMP = temp, TEMP1 = as.data.frame(temp1),
             A = as.data.frame(A), B = as.data.frame(B), C = as.data.frame(C),
             A1 = as.data.frame(A1), B1 = as.data.frame(B1), C1 = as.data.frame(C1),
             ALL = ALL), "out.rds")
