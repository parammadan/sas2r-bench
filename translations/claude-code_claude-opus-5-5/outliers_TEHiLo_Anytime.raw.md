```r
suppressPackageStartupMessages({
  library(dplyr)
  library(haven)
})

# ---------- helpers ----------
sas_order <- function(df, keys, desc = rep(FALSE, length(keys))) {
  cols <- lapply(keys, function(k) {
    x <- df[[k]]
    if (is.numeric(x)) { x[is.na(x)] <- -Inf } else { x <- as.character(x); x[is.na(x)] <- "" ; x <- sub(" +$", "", x) }
    x
  })
  do.call(order, c(unname(cols), list(decreasing = desc, method = "radix")))
}
sas_sort <- function(df, keys, desc = rep(FALSE, length(keys))) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  df[sas_order(df, keys, desc), , drop = FALSE] |> `rownames<-`(NULL)
}
keystr <- function(df, by) {
  do.call(paste, c(lapply(by, function(b) { x <- df[[b]]; ifelse(is.na(x), "\001NA", as.character(x)) }), sep = "\r"))
}
merge_by <- function(A, B, by) {
  A <- as.data.frame(A, stringsAsFactors = FALSE); B <- as.data.frame(B, stringsAsFactors = FALSE)
  allcols <- c(names(A), setdiff(names(B), names(A)))
  typ <- sapply(allcols, function(cc) {
    x <- if (cc %in% names(A)) A[[cc]] else B[[cc]]
    if (is.numeric(x)) "num" else "chr"
  })
  miss <- lapply(allcols, function(cc) if (typ[[cc]] == "num") NA_real_ else "")
  names(miss) <- allcols
  K <- rbind(A[by], B[by])
  K <- K[!duplicated(keystr(K, by)), , drop = FALSE]
  K <- K[sas_order(K, by), , drop = FALSE]
  kA <- keystr(A, by); kB <- keystr(B, by); kK <- keystr(K, by)
  rows <- list()
  for (g in seq_len(nrow(K))) {
    ia <- which(kA == kK[g]); ib <- which(kB == kK[g])
    pdv <- miss
    for (b in by) pdv[[b]] <- K[[b]][g]
    for (i in seq_len(max(length(ia), length(ib)))) {
      if (i <= length(ia)) for (cc in names(A)) pdv[[cc]] <- A[[cc]][ia[i]]
      if (i <= length(ib)) for (cc in names(B)) pdv[[cc]] <- B[[cc]][ib[i]]
      rows[[length(rows) + 1]] <- pdv
    }
  }
  out <- lapply(allcols, function(cc) {
    if (length(rows) == 0) return(if (typ[[cc]] == "num") numeric(0) else character(0))
    v <- lapply(rows, `[[`, cc)
    if (typ[[cc]] == "num") as.numeric(unlist(v)) else as.character(unlist(v))
  })
  names(out) <- allcols
  as.data.frame(out, stringsAsFactors = FALSE)
}
rnd <- function(x, d) sign(x) * floor(abs(x) * 10^d + 0.5 + 1e-9) / 10^d
strip_attrs <- function(df) { df <- as.data.frame(df, stringsAsFactors = FALSE); for (n in names(df)) attributes(df[[n]]) <- NULL; df }

# ---------- data advs ----------
raw <- strip_attrs(read_xpt("data/advs.xpt"))
names(raw) <- toupper(names(raw))
for (n in names(raw)) if (is.character(raw[[n]])) raw[[n]][is.na(raw[[n]])] <- ""

keep <- with(raw,
  sub(" +$", "", PARAMCD) == "DIABP" & !is.na(ATPTN) & ATPTN == 815 &
  !is.na(ADY) & ADY > 1 & sub(" +$", "", SAFFL) == "Y" &
  (CRIT1FL != "" | CRIT2FL != ""))
keep[is.na(keep)] <- FALSE
keepvars <- c("USUBJID", "TRTPN", "PARAM", "PARAMCD", "AVAL", "CRIT1", "CRIT1FL", "CRIT2", "CRIT2FL")
advs <- raw[keep, names(raw)[names(raw) %in% keepvars], drop = FALSE]
rownames(advs) <- NULL

# ---------- proc sql ----------
sqlcount <- function(d, cntname, tefl) {
  out <- d |>
    group_by(TRTPN, PARAM) |>
    summarise(V = n_distinct(USUBJID[USUBJID != ""]), .groups = "drop") |>
    as.data.frame(stringsAsFactors = FALSE)
  out$TEFL <- rep(tefl, nrow(out))
  names(out)[names(out) == "V"] <- cntname
  out <- out[, c("TRTPN", "PARAM", cntname, "TEFL"), drop = FALSE]
  sas_sort(out, c("TRTPN", "PARAM", cntname, "TEFL"))
}
TOTN_LOW  <- sqlcount(advs[advs$CRIT1FL != "", , drop = FALSE], "TOTN", "Low")
N_LOW     <- sqlcount(advs[advs$CRIT1FL == "Y", , drop = FALSE], "N", "Low")
TOTN_HIGH <- sqlcount(advs[advs$CRIT2FL != "", , drop = FALSE], "TOTN", "High")
N_HIGH    <- sqlcount(advs[advs$CRIT2FL == "Y", , drop = FALSE], "N", "High")

TOTN <- rbind(TOTN_HIGH, TOTN_LOW)
TE_N <- rbind(N_HIGH, N_LOW)
TOTN <- sas_sort(TOTN, c("TRTPN", "PARAM"))
TE_N <- sas_sort(TE_N, c("TRTPN", "PARAM"))

ALL_N <- merge_by(TOTN, TE_N, c("TRTPN", "PARAM"))
ALL_N$PCT <- ALL_N$N / ALL_N$TOTN * 100

# ---------- fisher macro ----------
fisherpval <- function(advs, flg, tefl) {
  advs <- sas_sort(advs, c("USUBJID", flg))
  last <- c(advs$USUBJID[-1] != advs$USUBJID[-nrow(advs)], TRUE)
  if (nrow(advs) == 0) last <- logical(0)
  sub_ds <- advs[last, , drop = FALSE]; rownames(sub_ds) <- NULL
  f <- sub_ds[sub_ds$PARAM != "" & !is.na(sub_ds$TRTPN) & sub_ds[[flg]] != "", , drop = FALSE]
  pvals <- numeric(0)
  for (p in sort(unique(f$PARAM), method = "radix")) {
    s <- f[f$PARAM == p, , drop = FALSE]
    tab <- table(s$TRTPN, s[[flg]])
    tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
    if (nrow(tab) >= 2 && ncol(tab) >= 2) {
      pv <- tryCatch(fisher.test(tab, workspace = 2e8)$p.value,
                     error = function(e) fisher.test(tab, simulate.p.value = TRUE, B = 1e6)$p.value)
      pvals <- c(pvals, pv)
    }
  }
  fish <- data.frame(NVALUE1 = pvals, TEFL = rep(tefl, length(pvals)), stringsAsFactors = FALSE)
  list(advs = advs, sub = sub_ds, fish = fish)
}

r1 <- fisherpval(advs, "CRIT1FL", "Low")
advs <- r1$advs; CRIT1FL <- r1$sub; FISHER_CRIT1FL <- r1$fish
r2 <- fisherpval(advs, "CRIT2FL", "High")
advs <- r2$advs; CRIT2FL <- r2$sub; FISHER_CRIT2FL <- r2$fish

ALL_FISHER <- rbind(FISHER_CRIT1FL, FISHER_CRIT2FL)
ALL_FISHER <- sas_sort(ALL_FISHER, "TEFL")
ALL_N <- sas_sort(ALL_N, c("TEFL", "TRTPN"))

# ---------- table_vals ----------
TABLE_VALS <- merge_by(ALL_N, ALL_FISHER, "TEFL")
fmt51 <- function(x) ifelse(is.na(x), ".", sprintf("%.1f", rnd(x, 1)))
TABLE_VALS$PCTC <- paste0("(", fmt51(TABLE_VALS$PCT), ")")
TABLE_VALS$NC <- ifelse(is.na(TABLE_VALS$N), ".", format(TABLE_VALS$N, scientific = FALSE, trim = TRUE))
TABLE_VALS$NC <- ifelse(is.na(TABLE_VALS$N), ".", as.character(TABLE_VALS$N))

TABLE_VALS <- sas_sort(TABLE_VALS, c("PARAM", "TEFL", "NVALUE1"), c(FALSE, TRUE, FALSE))
pvfmt <- function(x) {
  ifelse(is.na(x), "    .",
         ifelse(x < 0.001, "<.001", sprintf("%5.3f", rnd(x, 3))))
}
nr <- nrow(TABLE_VALS)
if (nr > 0) {
  k <- keystr(TABLE_VALS, c("PARAM", "TEFL", "NVALUE1"))
  first <- c(TRUE, k[-1] != k[-nr])
  TABLE_VALS$NVALUE1C <- ifelse(first, pvfmt(TABLE_VALS$NVALUE1), "")
} else {
  TABLE_VALS$NVALUE1C <- character(0)
}

saveRDS(list(
  ADVS = advs,
  TOTN_LOW = TOTN_LOW,
  N_LOW = N_LOW,
  TOTN_HIGH = TOTN_HIGH,
  N_HIGH = N_HIGH,
  TOTN = TOTN,
  TE_N = TE_N,
  ALL_N = ALL_N,
  CRIT1FL = CRIT1FL,
  FISHER_CRIT1FL = FISHER_CRIT1FL,
  CRIT2FL = CRIT2FL,
  FISHER_CRIT2FL = FISHER_CRIT2FL,
  ALL_FISHER = ALL_FISHER,
  TABLE_VALS = TABLE_VALS
), "out.rds")
```