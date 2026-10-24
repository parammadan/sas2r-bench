# Translation of ae_serious.sas (target26) to base R + haven

# ---------------------------------------------------------------------------
# SAS-compatible helpers
# ---------------------------------------------------------------------------
read_sas_xpt <- function(path) {
  d <- as.data.frame(haven::read_xpt(path))
  names(d) <- toupper(names(d))
  for (v in names(d)) {
    x <- d[[v]]
    if (inherits(x, "haven_labelled")) x <- as.vector(haven::zap_labels(x))
    if (is.character(x)) x <- trimws(x, which = "right")
    d[[v]] <- x
  }
  d
}

# SAS missing: NA or blank character
is_miss <- function(x) {
  m <- is.na(x)
  if (is.character(x)) m <- m | (trimws(x) %in% "")
  m
}

# SAS "x = 'v'" (trailing blanks ignored; missing -> FALSE inside where)
eqs <- function(x, v) !is.na(x) & trimws(x, which = "right") == v

# Stable SAS-like ordering: character by bytes (C locale), numeric missing first
# (ascending) or last (descending)
sas_order_idx <- function(d, by, desc = character(0)) {
  keys <- list()
  add <- function(k) keys[[length(keys) + 1]] <<- k
  for (v in by) {
    x <- d[[v]]
    if (is.character(x)) {
      add(ifelse(is.na(x), "", x))
    } else if (v %in% desc) {
      add(is.na(x))
      add(-x)
    } else {
      add(!is.na(x))
      add(x)
    }
  }
  if (length(keys) == 0) return(seq_len(nrow(d)))
  do.call(order, c(keys, list(method = "radix")))
}

sas_sort <- function(d, by, desc = character(0)) {
  d <- d[sas_order_idx(d, by, desc), , drop = FALSE]
  rownames(d) <- NULL
  d
}

# SAS DATA step MERGE (no-BY allowed). Within each BY group, output max(nx, ny)
# rows; the shorter side is paired with its last observation (retained).
# Returns all columns plus logical flags .A (in x) and .B (in y).
sas_merge <- function(x, y, by = character(0)) {
  mk <- function(d) {
    if (length(by) == 0) return(rep("", nrow(d)))
    do.call(paste, c(unname(lapply(by, function(v) as.character(d[[v]]))), sep = "\r"))
  }
  kx <- mk(x)
  ky <- mk(y)
  if (length(by) == 0) {
    keys <- ""
  } else {
    kd <- rbind(x[, by, drop = FALSE], y[, by, drop = FALSE])
    kd <- kd[!duplicated(kd), , drop = FALSE]
    kd <- kd[sas_order_idx(kd, by), , drop = FALSE]
    keys <- mk(kd)
  }
  lx <- split(seq_len(nrow(x)), factor(kx, levels = keys))
  ly <- split(seq_len(nrow(y)), factor(ky, levels = keys))
  xi <- integer(0); yi <- integer(0)
  fa <- logical(0); fb <- logical(0)
  for (i in seq_along(keys)) {
    a <- lx[[i]]
    b <- ly[[i]]
    n <- max(length(a), length(b))
    if (n == 0) next
    ia <- if (length(a) > 0) a[pmin(seq_len(n), length(a))] else rep(NA_integer_, n)
    ib <- if (length(b) > 0) b[pmin(seq_len(n), length(b))] else rep(NA_integer_, n)
    xi <- c(xi, ia)
    yi <- c(yi, ib)
    fa <- c(fa, rep(length(a) > 0, n))
    fb <- c(fb, rep(length(b) > 0, n))
  }
  xs <- x[xi, , drop = FALSE]
  ys <- y[yi, , drop = FALSE]
  out <- xs
  for (v in names(ys)) {
    if (v %in% by) next
    out[[v]] <- ys[[v]]
  }
  for (v in by) out[[v]] <- ifelse(fa, xs[[v]], ys[[v]])
  out$.A <- fa
  out$.B <- fb
  rownames(out) <- NULL
  out
}

drop_flags <- function(d) d[, setdiff(names(d), c(".A", ".B")), drop = FALSE]

merge_all <- function(x, y, by) drop_flags(sas_merge(x, y, by))

merge_in <- function(x, y, by, both = FALSE) {
  m <- sas_merge(x, y, by)
  m <- m[if (both) m$.A & m$.B else m$.A, , drop = FALSE]
  drop_flags(m)
}

# PROC FREQ table output (COUNT = frequency or weighted frequency, PERCENT of table).
# sparse = TRUE: all combinations of the observed levels are output (zeros included).
freqN <- function(d, vars, weight = NULL, sparse = FALSE) {
  ok <- rep(TRUE, nrow(d))
  for (v in vars) ok <- ok & !is_miss(d[[v]])
  if (!is.null(weight)) ok <- ok & !is.na(d[[weight]])
  d <- d[ok, , drop = FALSE]
  w <- if (is.null(weight)) rep(1, nrow(d)) else as.numeric(d[[weight]])
  if (sparse) {
    lv <- lapply(vars, function(v) sort(unique(d[[v]]), method = "radix"))
    names(lv) <- vars
    grid <- expand.grid(lv, stringsAsFactors = FALSE)
  } else {
    grid <- d[vars]
    grid <- grid[sas_order_idx(grid, vars), , drop = FALSE]
    grid <- grid[!duplicated(grid), , drop = FALSE]
  }
  keyd <- do.call(paste, c(unname(as.list(d[vars])), sep = "\r"))
  keyg <- do.call(paste, c(unname(as.list(grid)), sep = "\r"))
  cnt <- vapply(keyg, function(k) sum(w[keyd == k]), numeric(1), USE.NAMES = FALSE)
  out <- data.frame(grid, COUNT = cnt, stringsAsFactors = FALSE)
  out$PERCENT <- 100 * out$COUNT / sum(out$COUNT)
  rownames(out) <- NULL
  out
}

# Fisher's exact p-value (Pr <= P) for an r x c table of COUNT
fisher_p <- function(f, r, c) {
  if (nrow(f) == 0) return(NA_real_)
  rl <- sort(unique(f[[r]]), method = "radix")
  cl <- sort(unique(f[[c]]), method = "radix")
  if (length(rl) < 2 || length(cl) < 2) return(NA_real_)
  M <- matrix(0, length(rl), length(cl))
  ri <- match(f[[r]], rl)
  ci <- match(f[[c]], cl)
  for (k in seq_along(ri)) M[ri[k], ci[k]] <- M[ri[k], ci[k]] + f$COUNT[k]
  tryCatch(fisher.test(M, workspace = 2e8)$p.value, error = function(e) NA_real_)
}

fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 0.0001, "<.0001", sprintf("%.4f", p)))
}

# PROC TRANSPOSE (var CP, id TRTAN, prefix trt, by = by vars)
transpose_cp <- function(d, by, trt_levels) {
  dk <- do.call(paste, c(unname(as.list(d[by])), sep = "\r"))
  g <- d[!duplicated(dk), by, drop = FALSE]
  gk <- do.call(paste, c(unname(as.list(g)), sep = "\r"))
  ik <- paste(dk, d$TRTAN, sep = "\r")
  out <- g
  out[["_NAME_"]] <- "cp"
  for (t in trt_levels) {
    m <- match(paste(gk, t, sep = "\r"), ik)
    out[[paste0("trt", t)]] <- d$CP[m]
  }
  rownames(out) <- NULL
  out
}

# SAS SET of two datasets (column union)
bind_rows_base <- function(a, b) {
  cols <- union(names(a), names(b))
  for (v in setdiff(cols, names(a))) a[[v]] <- NA
  for (v in setdiff(cols, names(b))) b[[v]] <- NA
  out <- rbind(a[cols], b[cols])
  rownames(out) <- NULL
  out
}

# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
adae_raw <- read_sas_xpt("data/adae.xpt")
adsl_raw <- read_sas_xpt("data/adsl.xpt")

# ---------------------------------------------------------------------------
# SAE data: work.adae
# ---------------------------------------------------------------------------
adae <- adae_raw[, intersect(names(adae_raw),
                             c("USUBJID", "TRTAN", "AEDECOD", "AESER",
                               "AOCC02FL", "AOCC04FL", "SAFFL")),
                 drop = FALSE]
keep_rows <- (eqs(adae$AOCC02FL, "Y") | eqs(adae$AOCC04FL, "Y")) &
  eqs(adae$AESER, "Y") & eqs(adae$SAFFL, "Y")
adae <- adae[keep_rows, , drop = FALSE]
adae <- sas_sort(adae, "USUBJID")

# ---------------------------------------------------------------------------
# Subject info: work.adsl (nodupkey by usubjid trt01an)
# ---------------------------------------------------------------------------
adsl <- adsl_raw[, intersect(names(adsl_raw), c("USUBJID", "TRT01AN")), drop = FALSE]
adsl <- adsl[!is.na(adsl$TRT01AN), , drop = FALSE]
adsl <- sas_sort(adsl, c("USUBJID", "TRT01AN"))
adsl <- adsl[!duplicated(adsl[c("USUBJID", "TRT01AN")]), , drop = FALSE]
rownames(adsl) <- NULL

adae <- sas_sort(adae, "USUBJID")

# merge adsl(in=a) adae(in=b); by usubjid; if a and b;
adae_sl <- merge_in(adsl, adae, "USUBJID", both = TRUE)

# ---------------------------------------------------------------------------
# Denominator: sum_1 and macro variables (trtnum, numtrt)
# ---------------------------------------------------------------------------
sum_1 <- freqN(adsl, "TRT01AN")[, c("TRT01AN", "COUNT")]
trtnum <- sum_1$TRT01AN
numtrt <- length(trtnum)

# proc sort adae_sl out=adae_sl2 nodupkey by usubjid aedecod
adae_sl2 <- sas_sort(adae_sl, c("USUBJID", "AEDECOD"))
adae_sl2 <- adae_sl2[!duplicated(adae_sl2[c("USUBJID", "AEDECOD")]), , drop = FALSE]
rownames(adae_sl2) <- NULL

# proc sort adae_sl out=adae_dummy(keep=aedecod) nodupkey; by aedecod; where aocc02fl='Y'
adae_dummy <- adae_sl[eqs(adae_sl$AOCC02FL, "Y"), , drop = FALSE]
adae_dummy <- adae_dummy["AEDECOD"]
adae_dummy <- sas_sort(adae_dummy, "AEDECOD")
adae_dummy <- adae_dummy[!duplicated(adae_dummy$AEDECOD), , drop = FALSE]
rownames(adae_dummy) <- NULL

adae_sl2 <- sas_sort(adae_sl2, "AEDECOD")

# merge adae_sl2(in=a) adae_dummy(in=b); by aedecod; if a and b;
adae_sl3 <- merge_in(adae_sl2, adae_dummy, "AEDECOD", both = TRUE)

# proc freq adae_sl3: table aedecod*trtan / out=sum_2(keep=aedecod trtan count)
sum_2 <- freqN(adae_sl3, c("AEDECOD", "TRTAN"))[, c("AEDECOD", "TRTAN", "COUNT")]

# %dummy macro: one row per aedecod for each treatment level (ascending)
dummy_2 <- data.frame(
  AEDECOD = rep(adae_dummy$AEDECOD, each = numtrt),
  TRTAN   = rep(trtnum, times = nrow(adae_dummy)),
  stringsAsFactors = FALSE
)

# merge dummy_2(in=a) sum_2(in=b); by aedecod trtan;
adae_sl4 <- merge_all(dummy_2, sum_2, c("AEDECOD", "TRTAN"))

# proc sort adae_sl4; by trtan;
adae_sl4 <- sas_sort(adae_sl4, "TRTAN")

# merge adae_sl4(in=a) sum_1(rename=(trt01an=trtan count=bign)); by trtan; if a;
bign_df <- data.frame(TRTAN = sum_1$TRT01AN, BIGN = sum_1$COUNT)
adae_sl5 <- merge_in(adae_sl4, bign_df, "TRTAN", both = FALSE)

# data adae_sl6: set adae_sl5; by trtan aedecod; do j=0 to 1; ... output; end;
n5 <- nrow(adae_sl5)
k5 <- paste(adae_sl5$TRTAN, adae_sl5$AEDECOD, sep = "\r")
first5 <- if (n5 > 0) c(TRUE, k5[-1] != k5[-n5]) else logical(0)
ri  <- rep(seq_len(n5), each = 2)
jj  <- rep(c(0, 1), times = n5)
fl  <- first5[ri]
cnt5 <- adae_sl5$COUNT[ri]
bg5  <- adae_sl5$BIGN[ri]
adae_sl6 <- adae_sl5[ri, , drop = FALSE]
# resp is only assigned on the first observation of each aedecod group;
# otherwise it is missing (not retained), and so is count_final
adae_sl6$RESP <- ifelse(fl, jj, NA_real_)
adae_sl6$COUNT_FINAL <- ifelse(
  !fl, NA_real_,
  ifelse(jj == 1,
         ifelse(is.na(cnt5), 0, cnt5),
         ifelse(is.na(cnt5), bg5, bg5 - cnt5))
)
rownames(adae_sl6) <- NULL

adae_sl6 <- sas_sort(adae_sl6, "AEDECOD")

# ---------------------------------------------------------------------------
# Fisher / frequency per aedecod (by-group PROC FREQ, weight count_final)
# ---------------------------------------------------------------------------
gl <- sort(unique(adae_sl6$AEDECOD[!is_miss(adae_sl6$AEDECOD)]), method = "radix")
gidx <- split(seq_len(nrow(adae_sl6)), factor(adae_sl6$AEDECOD, levels = gl))

sum_final_parts <- list()
fish1_parts <- list()
for (i in seq_along(gl)) {
  g <- gl[i]
  f <- freqN(adae_sl6[gidx[[i]], , drop = FALSE], c("TRTAN", "RESP"),
             weight = "COUNT_FINAL", sparse = TRUE)
  f <- data.frame(AEDECOD = g, f, stringsAsFactors = FALSE)
  sum_final_parts[[i]] <- f
  p <- fisher_p(f, "TRTAN", "RESP")
  if (!is.na(p)) {
    fish1_parts[[length(fish1_parts) + 1]] <- data.frame(
      AEDECOD = g, NAME1 = "XP2_FISH", LABEL1 = "Pr <= P",
      CVALUE1 = fmt_p(p), NVALUE1 = p, stringsAsFactors = FALSE)
  }
}
sum_final <- do.call(rbind, sum_final_parts)
rownames(sum_final) <- NULL

fish1 <- if (length(fish1_parts) > 0) {
  do.call(rbind, fish1_parts)
} else {
  data.frame(AEDECOD = character(0), NAME1 = character(0), LABEL1 = character(0),
             CVALUE1 = character(0), NVALUE1 = numeric(0), stringsAsFactors = FALSE)
}
rownames(fish1) <- NULL

# data fish2(keep=aedecod cvalue1); where name1='XP2_FISH';
fish2 <- fish1[, c("AEDECOD", "CVALUE1"), drop = FALSE]

# data sum_final: where resp=1; cp; ord
sf <- sum_final[!is.na(sum_final$RESP) & sum_final$RESP == 1, , drop = FALSE]
sf$CP <- paste0(sprintf("%.0f", sf$COUNT), "(", sprintf("%.1f", sf$COUNT), ")")
sf$ORD <- ifelse(sf$TRTAN == 81, 1, 2)
sf <- sas_sort(sf, c("AEDECOD", "ORD"))

# data sum_final; retain count_max; if ord=1 then count_max=count;
cm <- NA_real_
cmv <- numeric(nrow(sf))
for (i in seq_len(nrow(sf))) {
  if (sf$ORD[i] == 1) cm <- sf$COUNT[i]
  cmv[i] <- cm
}
sf$COUNT_MAX <- cmv
sf <- sas_sort(sf, c("AEDECOD", "COUNT_MAX"))

# proc transpose by aedecod count_max; id trtan; var cp; prefix trt
trt_levels <- sort(unique(sf$TRTAN), method = "radix")
sum_final2 <- transpose_cp(sf, c("AEDECOD", "COUNT_MAX"), trt_levels)

# data sum_final2; merge sum_final2 fish2; by aedecod;
sum_final2 <- merge_all(sum_final2, fish2, "AEDECOD")

# proc sort; by descending count_max aedecod;
sum_final2 <- sas_sort(sum_final2, c("COUNT_MAX", "AEDECOD"), desc = "COUNT_MAX")

# ---------------------------------------------------------------------------
# At least one SAE (subject level)
# ---------------------------------------------------------------------------
adae_sl_sub <- adae_sl[eqs(adae_sl$AOCC02FL, "Y"), "USUBJID", drop = FALSE]
adae_sl_sub <- sas_sort(adae_sl_sub, "USUBJID")
adae_sl_sub <- adae_sl_sub[!duplicated(adae_sl_sub$USUBJID), , drop = FALSE]
rownames(adae_sl_sub) <- NULL

# data adae_sl_sub2; merge adsl(in=a) adae_sl_sub(in=b); by usubjid; if a;
tmpm <- sas_merge(adsl, adae_sl_sub, "USUBJID")
tmpm <- tmpm[tmpm$.A, , drop = FALSE]
tmpm$FLG <- ifelse(tmpm$.B, 1, 0)
adae_sl_sub2 <- drop_flags(tmpm)
rownames(adae_sl_sub2) <- NULL

# Fisher test on subject-level table (trtan renamed from trt01an)
tmp_sub <- adae_sl_sub2
names(tmp_sub)[names(tmp_sub) == "TRT01AN"] <- "TRTAN"
f_sub <- freqN(tmp_sub, c("TRTAN", "FLG"), sparse = TRUE)
p_sub <- fisher_p(f_sub, "TRTAN", "FLG")

fish_sub <- if (!is.na(p_sub)) {
  data.frame(NAME1 = "XP2_FISH", LABEL1 = "Pr <= P", CVALUE1 = fmt_p(p_sub),
             NVALUE1 = p_sub, stringsAsFactors = FALSE)
} else {
  data.frame(NAME1 = character(0), LABEL1 = character(0), CVALUE1 = character(0),
             NVALUE1 = numeric(0), stringsAsFactors = FALSE)
}

# data fish_sub2(keep=cvalue1); where name1='XP2_FISH';
fish_sub2 <- fish_sub[, "CVALUE1", drop = FALSE]

# data sum_final_sub; length aedecod $100; where flg=1; cp; aedecod=...
sum_final_sub <- f_sub[f_sub$FLG == 1, , drop = FALSE]
sum_final_sub$CP <- paste0(sprintf("%.0f", sum_final_sub$COUNT), "(",
                           sprintf("%.1f", sum_final_sub$COUNT), ")")
sum_final_sub$AEDECOD <- "Number of subjects reporting serious adverse events"
sum_final_sub <- sas_sort(sum_final_sub, "AEDECOD")

# proc transpose by aedecod; id trtan; var cp; prefix trt
sum_final_sub2 <- transpose_cp(sum_final_sub, "AEDECOD",
                               sort(unique(sum_final_sub$TRTAN), method = "radix"))

# data sum_final_sub2; merge sum_final_sub2 fish_sub2;  (no BY: sequential)
sum_final_sub2 <- merge_all(sum_final_sub2, fish_sub2, character(0))

# data final_report; set sum_final_sub2 sum_final2;
final_report <- bind_rows_base(sum_final_sub2, sum_final2)

# ---------------------------------------------------------------------------
# Save all datasets
# ---------------------------------------------------------------------------
saveRDS(
  list(
    ADAE           = adae,
    ADSL           = adsl,
    ADAE_SL        = adae_sl,
    SUM_1          = sum_1,
    ADAE_SL2       = adae_sl2,
    ADAE_DUMMY     = adae_dummy,
    ADAE_SL3       = adae_sl3,
    SUM_2          = sum_2,
    DUMMY_2        = dummy_2,
    ADAE_SL4       = adae_sl4,
    ADAE_SL5       = adae_sl5,
    ADAE_SL6       = adae_sl6,
    FISH1          = fish1,
    SUM_FINAL      = sf,
    FISH2          = fish2,
    SUM_FINAL2     = sum_final2,
    ADAE_SL_SUB    = adae_sl_sub,
    ADAE_SL_SUB2   = adae_sl_sub2,
    FISH_SUB       = fish_sub,
    SUM_FINAL_SUB  = sum_final_sub,
    FISH_SUB2      = fish_sub2,
    SUM_FINAL_SUB2 = sum_final_sub2,
    FINAL_REPORT   = final_report
  ),
  "out.rds"
)
