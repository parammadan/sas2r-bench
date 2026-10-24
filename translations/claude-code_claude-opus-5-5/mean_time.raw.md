```r
library(haven)
library(dplyr)

src <- read_xpt("data/advs.xpt")
names(src) <- toupper(names(src))

keepv <- c("USUBJID","TRTPN","ANL01FL","PARAM","PARAMCD","AVISIT","AVISITN","ATPTN","SAFFL","AVAL")
ADVS <- src %>%
  filter(PARAMCD == "DIABP" & !is.na(ATPTN) & ATPTN == 815 & ANL01FL == "Y" &
           !is.na(AVISITN) & AVISITN < 99 & SAFFL == "Y")
ADVS <- ADVS[, names(ADVS)[names(ADVS) %in% keepv]]

# PROC SORT by trtpn avisitn (stable, missing first)
ADVS <- ADVS %>% arrange(!is.na(TRTPN), TRTPN, !is.na(AVISITN), AVISITN)

# PROC MEANS BY trtpn avisitn avisit (consecutive groups)
key <- paste(ADVS$TRTPN, ADVS$AVISITN, ADVS$AVISIT, sep = "\r")
grp <- cumsum(c(TRUE, key[-1] != key[-length(key)]))
if (nrow(ADVS) == 0) grp <- integer(0)
ADVS$.grp <- grp
TEMP <- ADVS %>%
  group_by(.grp) %>%
  summarise(TRTPN = first(TRTPN), AVISITN = first(AVISITN), AVISIT = first(AVISIT),
            `_TYPE_` = 0, `_FREQ_` = n(),
            N = sum(!is.na(AVAL)),
            MEAN = if (N > 0) mean(AVAL, na.rm = TRUE) else NA_real_,
            SD = if (N > 1) sd(AVAL, na.rm = TRUE) else NA_real_,
            .groups = "drop") %>%
  mutate(STDERR = ifelse(N > 1, SD / sqrt(N), NA_real_),
         tq = ifelse(N > 1, suppressWarnings(qt(0.975, pmax(N - 1, 1))), NA_real_),
         LCLM = MEAN - tq * STDERR,
         UCLM = MEAN + tq * STDERR) %>%
  select(TRTPN, AVISITN, AVISIT, `_TYPE_`, `_FREQ_`, N, MEAN, STDERR, LCLM, UCLM) %>%
  as.data.frame()
ADVS$.grp <- NULL
ADVS <- as.data.frame(ADVS)

tinv95 <- function(df) ifelse(!is.na(df) & df > 0, suppressWarnings(qt(0.95, pmax(df, 1))), NA_real_)
TEMP1 <- TEMP %>%
  mutate(LO = MEAN - (tinv95(N - 1) * STDERR),
         HI = MEAN + (tinv95(N - 1) * STDERR)) %>%
  select(-`_TYPE_`, -`_FREQ_`, -LCLM, -UCLM, -STDERR) %>%
  as.data.frame()

A <- TEMP1[!is.na(TEMP1$TRTPN) & TEMP1$TRTPN == 0, , drop = FALSE]
B <- TEMP1[!is.na(TEMP1$TRTPN) & TEMP1$TRTPN == 54, , drop = FALSE]
C <- TEMP1[!is.na(TEMP1$TRTPN) & TEMP1$TRTPN == 81, , drop = FALSE]
rownames(A) <- NULL; rownames(B) <- NULL; rownames(C) <- NULL

rn <- function(d, p) {
  d %>% rename(!!paste0(p, "_N") := N, !!paste0(p, "_MEAN") := MEAN,
               !!paste0(p, "_LCL") := LO, !!paste0(p, "_UCL") := HI)
}
A1 <- rn(A, "A"); B1 <- rn(B, "B"); C1 <- rn(C, "C")

# SAS match-merge by AVISITN AVISIT
sas_merge <- function(dfs, by) {
  dfs <- lapply(dfs, function(d) d %>% group_by(across(all_of(by))) %>%
                  mutate(.seq = row_number(), .cnt = n()) %>% ungroup())
  keys <- bind_rows(lapply(dfs, function(d) d %>% distinct(across(all_of(c(by, ".cnt")))))) %>%
    group_by(across(all_of(by))) %>% summarise(.max = max(.cnt), .groups = "drop")
  keys <- keys[rep(seq_len(nrow(keys)), keys$.max), ] %>%
    group_by(across(all_of(by))) %>% mutate(.seq = row_number()) %>% ungroup() %>% select(-.max)
  out <- keys
  allv <- by
  for (d in dfs) {
    cnts <- d %>% distinct(across(all_of(c(by, ".cnt"))))
    o2 <- out %>% left_join(cnts, by = by) %>%
      mutate(.s2 = ifelse(is.na(.cnt), NA_integer_, pmin(.seq, .cnt))) %>% select(-.cnt)
    dd <- d %>% select(-.cnt) %>% rename(.s2 = .seq)
    newv <- setdiff(names(dd), c(by, ".s2"))
    o2 <- o2 %>% left_join(dd, by = c(by, ".s2"), suffix = c("", ".new"))
    for (v in intersect(newv, names(out))) {
      nv <- paste0(v, ".new")
      has <- !is.na(o2$.s2)
      o2[[v]][has] <- o2[[nv]][has]
      o2[[nv]] <- NULL
    }
    o2$.s2 <- NULL
    out <- o2
    allv <- c(allv, setdiff(newv, allv))
  }
  out <- out %>% arrange(!is.na(AVISITN), AVISITN, AVISIT, .seq)
  first_order <- unique(c(names(dfs[[1]]), unlist(lapply(dfs, names))))
  first_order <- setdiff(first_order, c(".seq", ".cnt"))
  as.data.frame(out[, first_order])
}
ALL <- sas_merge(list(A1, B1, C1), c("AVISITN", "AVISIT"))
ALL <- ALL %>% rename(XC = AVISITN, X = AVISIT) %>% select(-TRTPN) %>% as.data.frame()

saveRDS(list(ADVS = ADVS, TEMP = TEMP, TEMP1 = TEMP1, A = A, B = B, C = C,
             A1 = as.data.frame(A1), B1 = as.data.frame(B1), C1 = as.data.frame(C1),
             ALL = ALL), "out.rds")
```