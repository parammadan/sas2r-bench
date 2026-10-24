suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
})

# ---------- helpers ----------
trunc_w <- function(x, w) if (!is.null(w)) substr(x, 1, w) else x

propcase <- function(x) {
  x <- ifelse(is.na(x), "", x)
  x <- tolower(x)
  gsub("(^|[ /\\-(.\t])([a-z])", "\\1\\U\\2", x, perl = TRUE)
}

catx_sp <- function(...) {
  m <- cbind(...)
  m[is.na(m)] <- ""
  m <- trimws(m)
  apply(m, 1, function(r) paste(r[r != ""], collapse = " "))
}

sas_round <- function(x, u) {
  y <- x / u
  sign(y) * floor(abs(y) + 0.5 + 1e-9) / (1 / u)
}

shiftord_in <- function(x) {
  v <- trimws(substr(ifelse(is.na(x), "", x), 1, 6))
  map <- c(Low = 1, Normal = 2, High = 3, Total = 4)
  out <- unname(map[v])
  nomatch <- is.na(out) & v != ""
  out[nomatch] <- suppressWarnings(as.numeric(v[nomatch]))
  out
}

shifttx_put <- function(x) {
  lab <- c("Low", "Normal", "High", "Total")
  out <- ifelse(!is.na(x) & x %in% 1:4, lab[match(x, 1:4)],
                ifelse(is.na(x), ".", as.character(x)))
  out
}

# ---------- ADVS ----------
advs <- read_xpt("data/advs.xpt")
w_param <- attr(advs$PARAM, "width")
advs$PARAM <- trunc_w(catx_sp(advs$PARAM, propcase(advs$ATPT)), w_param)
advs <- advs %>% mutate(PARAMN = ATPTN)
ADVS <- advs

# ---------- INDS ----------
invar <- c("USUBJID","TRTPN","TRTP","PARAM","PARAMCD","PARAMN","ATPT","ATPTN","AVAL",
           "AVISITN","ANL01FL","BNRIND","ANRIND","SHIFT1","CHGCAT1","SAFFL")
inds <- advs[, names(advs) %in% invar, drop = FALSE]
eqc <- function(x, v) !is.na(x) & trimws(x, "right") == v
inds <- inds[eqc(inds$PARAMCD, "DIABP") & eqc(inds$ANL01FL, "Y") &
               !is.na(inds$AVISITN) & inds$AVISITN == 99 & eqc(inds$SAFFL, "Y") &
               !eqc(inds$BNRIND, "Missing") & !eqc(inds$ANRIND, "Missing"), , drop = FALSE]
INDS <- inds
w_bnr <- attr(inds$BNRIND, "width")
w_anr <- attr(inds$ANRIND, "width")

# ---------- TRTCOMB ----------
TRTCOMB <- inds %>%
  transmute(START = TRTPN, LABEL = TRTP, FMTNAME = "TRTFT", TYPE = "n") %>%
  distinct() %>%
  arrange(!is.na(START), START, LABEL) %>%
  as.data.frame()

fmt_start <- TRTCOMB$START
fmt_label <- trimws(TRTCOMB$LABEL, "right")
trtft_put <- function(x) {
  idx <- match(x, fmt_start)
  ifelse(is.na(idx), ifelse(is.na(x), ".", as.character(x)), fmt_label[idx])
}

# ---------- COUNTS (proc summary) ----------
chr_miss <- function(x) is.na(x) | trimws(x) == ""
base <- inds[!chr_miss(inds$PARAM) & !is.na(inds$PARAMN) & !is.na(inds$TRTPN) &
               !chr_miss(inds$BNRIND) & !chr_miss(inds$ANRIND), , drop = FALSE]

mk_type <- function(type, use_b, use_a) {
  gv <- c("PARAM", "PARAMN", "TRTPN", if (use_b) "BNRIND", if (use_a) "ANRIND")
  r <- base %>% group_by(across(all_of(gv))) %>%
    summarise(`_FREQ_` = n(), N = sum(!is.na(AVAL)), .groups = "drop")
  if (!use_b) r$BNRIND <- ""
  if (!use_a) r$ANRIND <- ""
  r$`_TYPE_` <- type
  r %>% select(PARAM, PARAMN, TRTPN, BNRIND, ANRIND, `_TYPE_`, `_FREQ_`, N) %>%
    arrange(PARAM, PARAMN, TRTPN, BNRIND, ANRIND)
}
counts <- bind_rows(mk_type(28, FALSE, FALSE), mk_type(29, FALSE, TRUE),
                    mk_type(30, TRUE, FALSE), mk_type(31, TRUE, TRUE))
counts <- counts %>% arrange(PARAM, PARAMN, TRTPN) %>% as.data.frame()
COUNTS <- counts

# ---------- FIXTOT / DENOM ----------
ft <- counts
ft$BNRIND[chr_miss(ft$BNRIND)] <- trunc_w("Total", w_bnr)
ft$ANRIND[chr_miss(ft$ANRIND)] <- trunc_w("Total", w_anr)
ft$BNRINDN <- shiftord_in(ft$BNRIND)
ft$ANRINDN <- shiftord_in(ft$ANRIND)
isden <- eqc(ft$BNRIND, "Total") & eqc(ft$ANRIND, "Total")
DENOM <- ft[isden, c("PARAM", "PARAMN", "TRTPN", "N")] %>% rename(DENOM = N) %>% as.data.frame()
fixtot <- ft[!isden, c("PARAM", "PARAMN", "TRTPN", "BNRIND", "ANRIND", "N", "BNRINDN", "ANRINDN")]
fixtot <- fixtot %>%
  arrange(PARAM, !is.na(PARAMN), PARAMN, !is.na(TRTPN), TRTPN,
          !is.na(BNRINDN), BNRINDN, !is.na(ANRINDN), ANRINDN) %>% as.data.frame()
FIXTOT <- fixtot

# ---------- ALLCAT ----------
ALLCAT <- fixtot %>% distinct(PARAM, PARAMN, TRTPN) %>%
  arrange(PARAM, !is.na(PARAMN), PARAMN, !is.na(TRTPN), TRTPN) %>% as.data.frame()

# ---------- FULL ----------
grid <- expand.grid(ANRINDN = 1:4, BNRINDN = 1:4)[, c("BNRINDN", "ANRINDN")]
FULL <- ALLCAT[rep(seq_len(nrow(ALLCAT)), each = 16), , drop = FALSE]
FULL$BNRINDN <- rep(as.numeric(grid$BNRINDN), nrow(ALLCAT))
FULL$ANRINDN <- rep(as.numeric(grid$ANRINDN), nrow(ALLCAT))
FULL$N <- 0
rownames(FULL) <- NULL
FULL <- as.data.frame(FULL)

# ---------- ALLRECS ----------
keys5 <- c("PARAM", "PARAMN", "TRTPN", "BNRINDN", "ANRINDN")
ft2 <- fixtot %>% rename(N_FT = N) %>% mutate(.m = TRUE)
allrecs <- FULL %>% left_join(ft2, by = keys5, multiple = "all", relationship = "many-to-many")
allrecs$N <- ifelse(!is.na(allrecs$.m), allrecs$N_FT, allrecs$N)
allrecs$BNRIND[is.na(allrecs$BNRIND)] <- ""
allrecs$ANRIND[is.na(allrecs$ANRIND)] <- ""
allrecs <- allrecs %>% select(PARAM, PARAMN, TRTPN, BNRINDN, ANRINDN, N, BNRIND, ANRIND) %>%
  as.data.frame()
ALLRECS <- allrecs

# ---------- ALLDENOM ----------
alldenom <- allrecs %>% left_join(DENOM, by = c("PARAM", "PARAMN", "TRTPN"),
                                  multiple = "all", relationship = "many-to-many")
alldenom$BNRIND <- trunc_w(shifttx_put(alldenom$BNRINDN), w_bnr)
alldenom$ANRIND <- trunc_w(shifttx_put(alldenom$ANRINDN), w_anr)
alldenom <- alldenom[!is.na(alldenom$DENOM), , drop = FALSE]
rownames(alldenom) <- NULL
ALLDENOM <- as.data.frame(alldenom)

# ---------- TRANS ----------
ad <- ALLDENOM
if (nrow(ad) > 0) {
  kstr <- paste(ad$PARAM, ad$PARAMN, ad$TRTPN, ad$BNRINDN, sep = "\r")
  gid <- cumsum(c(TRUE, kstr[-1] != kstr[-length(kstr)]))
} else gid <- integer(0)
res <- lapply(split(seq_len(nrow(ad)), gid), function(ix) {
  ns <- rep(NA_real_, 4)
  pcts <- rep("", 4)
  for (i in ix) {
    k <- ad$ANRINDN[i]
    ns[k] <- ad$N[i]
    d <- ad$DENOM[i]
    if (!is.na(d) && d > 0) {
      pcts[k] <- paste0("(", sprintf("%5.1f", sas_round((100 * ad$N[i]) / d, 0.1)), ")")
    } else pcts[k] <- "(  0.0)"
  }
  l <- ix[length(ix)]
  dn <- ad$DENOM[l]
  dtxt <- if (is.na(dn)) "." else format(dn, scientific = FALSE, trim = TRUE)
  trttx <- substr(catx_sp(trtft_put(ad$TRTPN[l]), paste0("(N = ", dtxt, ")")), 1, 30)
  data.frame(N1 = ns[1], N2 = ns[2], N3 = ns[3], N4 = ns[4],
             PCT1 = pcts[1], PCT2 = pcts[2], PCT3 = pcts[3], PCT4 = pcts[4],
             PARAM = ad$PARAM[l], PARAMN = ad$PARAMN[l], TRTPN = ad$TRTPN[l],
             BNRINDN = ad$BNRINDN[l], BNRIND = ad$BNRIND[l], DENOM = dn,
             TRTTX = trttx, stringsAsFactors = FALSE)
})
TRANS <- if (length(res)) bind_rows(res) else
  data.frame(N1 = numeric(0), N2 = numeric(0), N3 = numeric(0), N4 = numeric(0),
             PCT1 = character(0), PCT2 = character(0), PCT3 = character(0), PCT4 = character(0),
             PARAM = character(0), PARAMN = numeric(0), TRTPN = numeric(0),
             BNRINDN = numeric(0), BNRIND = character(0), DENOM = numeric(0),
             TRTTX = character(0))
TRANS <- as.data.frame(TRANS)

saveRDS(list(ADVS = ADVS, INDS = INDS, TRTCOMB = TRTCOMB, COUNTS = COUNTS,
             FIXTOT = FIXTOT, DENOM = DENOM, ALLCAT = ALLCAT, FULL = FULL,
             ALLRECS = ALLRECS, ALLDENOM = ALLDENOM, TRANS = TRANS), "out.rds")
