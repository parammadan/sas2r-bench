library(haven)
library(dplyr)
library(tidyr)
library(stringr)

advs_raw <- read_xpt("data/advs.xpt")
names(advs_raw) <- toupper(names(advs_raw))

propcase <- function(x) {
  x <- ifelse(is.na(x), "", x)
  vapply(x, function(s) {
    w <- strsplit(tolower(s), " ", fixed = TRUE)[[1]]
    if (length(w) == 0) return("")
    w <- paste0(toupper(substr(w, 1, 1)), substr(w, 2, nchar(w)))
    paste(w, collapse = " ")
  }, character(1), USE.NAMES = FALSE)
}
catx2 <- function(a, b) {
  a <- trimws(ifelse(is.na(a), "", a)); b <- trimws(ifelse(is.na(b), "", b))
  ifelse(a != "" & b != "", paste(a, b), ifelse(a != "", a, b))
}

ADVS <- advs_raw
ADVS$PARAM <- catx2(ADVS$PARAM, propcase(ADVS$ATPT))
ADVS$PARAMN <- ADVS$ATPTN

invar <- c("USUBJID","TRTPN","TRTP","PARAM","PARAMCD","PARAMN","ATPT","ATPTN","AVAL",
           "AVISITN","ANL01FL","BNRIND","ANRIND","SHIFT1","CHGCAT1","SAFFL")
INDS <- ADVS[, invar]
INDS <- INDS[which(INDS$PARAMCD == "DIABP" & INDS$ANL01FL == "Y" & INDS$AVISITN == 99 &
                     INDS$SAFFL == "Y" & INDS$BNRIND != "Missing" & INDS$ANRIND != "Missing"), ]

TRTCOMB <- INDS %>% distinct(START = TRTPN, LABEL = TRTP) %>%
  mutate(FMTNAME = "TRTFT", TYPE = "n") %>% arrange(START)
TRTCOMB <- as.data.frame(TRTCOMB)

# proc summary
d <- INDS %>% filter(!is.na(PARAM) & PARAM != "" & !is.na(BNRIND) & BNRIND != "" &
                       !is.na(ANRIND) & ANRIND != "" & !is.na(PARAMN) & !is.na(TRTPN))
mk <- function(g, type) {
  g %>% summarise(N = sum(!is.na(AVAL)), .groups = "drop") %>% mutate(`_TYPE_` = type)
}
t28 <- mk(d %>% group_by(PARAM, PARAMN, TRTPN), 28) %>% mutate(BNRIND = "", ANRIND = "")
t29 <- mk(d %>% group_by(PARAM, PARAMN, TRTPN, ANRIND), 29) %>% mutate(BNRIND = "")
t30 <- mk(d %>% group_by(PARAM, PARAMN, TRTPN, BNRIND), 30) %>% mutate(ANRIND = "")
t31 <- mk(d %>% group_by(PARAM, PARAMN, TRTPN, BNRIND, ANRIND), 31)
COUNTS <- bind_rows(t28, t29, t30, t31) %>%
  arrange(`_TYPE_`, PARAM, PARAMN, TRTPN, BNRIND, ANRIND) %>%
  arrange(PARAM, PARAMN, TRTPN)
COUNTS <- as.data.frame(COUNTS)

shiftord <- function(x) unname(c(Low = 1, Normal = 2, High = 3, Total = 4)[x])
shifttx <- function(i) c("Low", "Normal", "High", "Total")[i]

tmp <- COUNTS
tmp$BNRIND[is.na(tmp$BNRIND) | tmp$BNRIND == ""] <- "Total"
tmp$ANRIND[is.na(tmp$ANRIND) | tmp$ANRIND == ""] <- "Total"
tmp$BNRINDN <- shiftord(tmp$BNRIND)
tmp$ANRINDN <- shiftord(tmp$ANRIND)
isden <- tmp$BNRIND == "Total" & tmp$ANRIND == "Total"
DENOM <- tmp[isden, c("PARAM","PARAMN","TRTPN","N")]
names(DENOM)[4] <- "DENOM"
FIXTOT <- tmp[!isden, c("PARAM","PARAMN","TRTPN","ANRIND","BNRIND","ANRINDN","BNRINDN","N")]
FIXTOT <- FIXTOT %>% arrange(PARAM, PARAMN, TRTPN, BNRINDN, ANRINDN)
FIXTOT <- as.data.frame(FIXTOT)

ALLCAT <- as.data.frame(FIXTOT %>% distinct(PARAM, PARAMN, TRTPN) %>% arrange(PARAM, PARAMN, TRTPN))

FULL <- ALLCAT %>% crossing_helper <- NULL
FULL <- do.call(rbind, lapply(seq_len(nrow(ALLCAT)), function(i) {
  g <- expand.grid(ANRINDN = 1:4, BNRINDN = 1:4)
  cbind(ALLCAT[rep(i, 16), ], BNRINDN = g$BNRINDN, ANRINDN = g$ANRINDN, N = 0)
}))
rownames(FULL) <- NULL

ALLRECS <- FULL %>%
  left_join(FIXTOT %>% rename(N_F = N), by = c("PARAM","PARAMN","TRTPN","BNRINDN","ANRINDN")) %>%
  mutate(N = ifelse(is.na(N_F), N, N_F)) %>% select(-N_F) %>%
  select(PARAM, PARAMN, TRTPN, BNRINDN, ANRINDN, N, ANRIND, BNRIND) %>%
  arrange(PARAM, PARAMN, TRTPN, BNRINDN, ANRINDN)
ALLRECS <- as.data.frame(ALLRECS)

ALLDENOM <- ALLRECS %>% left_join(DENOM, by = c("PARAM","PARAMN","TRTPN"))
ALLDENOM$BNRIND <- shifttx(ALLDENOM$BNRINDN)
ALLDENOM$ANRIND <- shifttx(ALLDENOM$ANRINDN)
ALLDENOM <- as.data.frame(ALLDENOM[!is.na(ALLDENOM$DENOM), ])

rnd1 <- function(x) floor(x * 10 + 0.5 + 1e-9) / 10
trtlab <- function(v) TRTCOMB$LABEL[match(v, TRTCOMB$START)]

grp <- ALLDENOM %>% group_by(PARAM, PARAMN, TRTPN, BNRINDN) %>% group_keys()
rows <- lapply(seq_len(nrow(grp)), function(i) {
  g <- ALLDENOM[ALLDENOM$PARAM == grp$PARAM[i] & ALLDENOM$PARAMN == grp$PARAMN[i] &
                  ALLDENOM$TRTPN == grp$TRTPN[i] & ALLDENOM$BNRINDN == grp$BNRINDN[i], ]
  ns <- rep(NA_real_, 4); pc <- rep(NA_character_, 4)
  for (k in seq_len(nrow(g))) {
    a <- g$ANRINDN[k]
    ns[a] <- g$N[k]
    pc[a] <- if (g$DENOM[k] > 0) paste0("(", sprintf("%5.1f", rnd1(100 * g$N[k] / g$DENOM[k])), ")") else "(  0.0)"
  }
  last <- g[nrow(g), ]
  out <- data.frame(N1 = ns[1], N2 = ns[2], N3 = ns[3], N4 = ns[4],
                    PCT1 = pc[1], PCT2 = pc[2], PCT3 = pc[3], PCT4 = pc[4],
                    PARAM = last$PARAM, PARAMN = last$PARAMN, TRTPN = last$TRTPN,
                    BNRINDN = last$BNRINDN, BNRIND = last$BNRIND, DENOM = last$DENOM,
                    stringsAsFactors = FALSE)
  out$TRTTX <- catx2(trtlab(last$TRTPN), paste0("(N = ", format(last$DENOM, scientific = FALSE, trim = TRUE), ")"))
  out
})
TRANS <- do.call(rbind, rows)
rownames(TRANS) <- NULL

saveRDS(list(ADVS = ADVS, INDS = INDS, TRTCOMB = TRTCOMB, COUNTS = COUNTS,
             FIXTOT = FIXTOT, DENOM = DENOM, ALLCAT = ALLCAT, FULL = FULL,
             ALLRECS = ALLRECS, ALLDENOM = ALLDENOM, TRANS = TRANS), "out.rds")
