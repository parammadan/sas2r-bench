library(haven)
library(dplyr)

advs0 <- read_xpt("data/advs.xpt")
names(advs0) <- toupper(names(advs0))
miss <- function(x) is.na(x) | trimws(x) == ""
advs <- advs0 %>%
  filter(PARAMCD == "DIABP", !is.na(ATPTN), ATPTN == 815, !is.na(ADY), ADY > 1,
         SAFFL == "Y", (!miss(CRIT1FL) | !miss(CRIT2FL))) %>%
  select(USUBJID, TRTPN, PARAM, PARAMCD, AVAL, CRIT1, CRIT1FL, CRIT2, CRIT2FL)

cnt <- function(d, nm, tefl) {
  if (nrow(d) == 0) return(tibble(TRTPN = numeric(), PARAM = character(), !!nm := integer(), TEFL = character()))
  d %>% group_by(TRTPN, PARAM) %>%
    summarise(!!nm := n_distinct(USUBJID), .groups = "drop") %>%
    mutate(TEFL = tefl)
}
TOTN_LOW <- cnt(advs %>% filter(!miss(CRIT1FL)), "TOTN", "Low")
N_LOW <- cnt(advs %>% filter(!is.na(CRIT1FL) & CRIT1FL == "Y"), "N", "Low")
TOTN_HIGH <- cnt(advs %>% filter(!miss(CRIT2FL)), "TOTN", "High")
N_HIGH <- cnt(advs %>% filter(!is.na(CRIT2FL) & CRIT2FL == "Y"), "N", "High")

TOTN <- bind_rows(TOTN_HIGH, TOTN_LOW)
TE_N <- bind_rows(N_HIGH, N_LOW)
TOTN <- TOTN[order(TOTN$TRTPN, TOTN$PARAM), ]
TE_N <- TE_N[order(TE_N$TRTPN, TE_N$PARAM), ]

keys <- distinct(bind_rows(TOTN[, c("TRTPN", "PARAM")], TE_N[, c("TRTPN", "PARAM")]))
keys <- keys[order(keys$TRTPN, keys$PARAM), ]
rows <- list()
for (i in seq_len(nrow(keys))) {
  a <- TOTN[TOTN$TRTPN == keys$TRTPN[i] & TOTN$PARAM == keys$PARAM[i], ]
  b <- TE_N[TE_N$TRTPN == keys$TRTPN[i] & TE_N$PARAM == keys$PARAM[i], ]
  k <- max(nrow(a), nrow(b))
  for (j in seq_len(k)) {
    ia <- if (nrow(a) > 0) min(j, nrow(a)) else NA
    ib <- if (nrow(b) > 0) min(j, nrow(b)) else NA
    rows[[length(rows) + 1]] <- tibble(
      TRTPN = keys$TRTPN[i], PARAM = keys$PARAM[i],
      TOTN = if (is.na(ia)) NA_integer_ else a$TOTN[ia],
      TEFL = if (!is.na(ib)) b$TEFL[ib] else a$TEFL[ia],
      N = if (is.na(ib)) NA_integer_ else b$N[ib])
  }
}
ALL_N <- bind_rows(rows)
ALL_N$PCT <- ALL_N$N / ALL_N$TOTN * 100

fisher_for <- function(flg, tefl) {
  d <- advs
  d <- d[order(d$USUBJID, d[[flg]], na.last = FALSE, method = "radix"), ]
  d <- d[!duplicated(d$USUBJID, fromLast = TRUE), ]
  d <- d[!miss(d[[flg]]), ]
  res <- list()
  for (p in sort(unique(d$PARAM))) {
    dd <- d[d$PARAM == p, ]
    tb <- table(dd$TRTPN, dd[[flg]])
    if (nrow(tb) >= 2 && ncol(tb) >= 2) {
      pv <- tryCatch(fisher.test(tb)$p.value, error = function(e) NA_real_)
      res[[length(res) + 1]] <- tibble(TEFL = tefl, NVALUE1 = pv)
    }
  }
  if (length(res) == 0) tibble(TEFL = character(), NVALUE1 = numeric()) else bind_rows(res)
}
FISHER_CRIT1FL <- fisher_for("CRIT1FL", "Low")
FISHER_CRIT2FL <- fisher_for("CRIT2FL", "High")
ALL_FISHER <- bind_rows(FISHER_CRIT1FL, FISHER_CRIT2FL)
ALL_FISHER <- ALL_FISHER[order(ALL_FISHER$TEFL), ]

ALL_N <- ALL_N[order(ALL_N$TEFL, ALL_N$TRTPN), ]
TABLE_VALS <- left_join(ALL_N, ALL_FISHER, by = "TEFL", relationship = "many-to-many")
TABLE_VALS$PCTC <- paste0("(", ifelse(is.na(TABLE_VALS$PCT), ".", sprintf("%.1f", TABLE_VALS$PCT)), ")")
TABLE_VALS$NC <- ifelse(is.na(TABLE_VALS$N), ".", as.character(TABLE_VALS$N))

o <- order(TABLE_VALS$PARAM, -xtfrm(TABLE_VALS$TEFL), TABLE_VALS$NVALUE1, na.last = FALSE)
TABLE_VALS <- TABLE_VALS[o, ]

pv5 <- function(x) {
  ifelse(is.na(x), ".",
         ifelse(x < 0.001, "<.001", sub("^0\\.", ".", sprintf("%.3f", x))))
}
n <- nrow(TABLE_VALS)
first <- rep(TRUE, n)
if (n > 1) {
  same <- TABLE_VALS$PARAM[-1] == TABLE_VALS$PARAM[-n] &
    TABLE_VALS$TEFL[-1] == TABLE_VALS$TEFL[-n] &
    ((is.na(TABLE_VALS$NVALUE1[-1]) & is.na(TABLE_VALS$NVALUE1[-n])) |
       (!is.na(TABLE_VALS$NVALUE1[-1]) & !is.na(TABLE_VALS$NVALUE1[-n]) &
          TABLE_VALS$NVALUE1[-1] == TABLE_VALS$NVALUE1[-n]))
  first[-1] <- !same
}
TABLE_VALS$NVALUE1C <- ifelse(first, pv5(TABLE_VALS$NVALUE1), "")

saveRDS(list(ADVS = advs, TOTN_LOW = TOTN_LOW, N_LOW = N_LOW, TOTN_HIGH = TOTN_HIGH,
             N_HIGH = N_HIGH, TOTN = TOTN, TE_N = TE_N, ALL_N = ALL_N,
             FISHER_CRIT1FL = FISHER_CRIT1FL, FISHER_CRIT2FL = FISHER_CRIT2FL,
             ALL_FISHER = ALL_FISHER, TABLE_VALS = TABLE_VALS), "out.rds")
