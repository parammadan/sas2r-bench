```r
suppressPackageStartupMessages({
  library(dplyr)
  library(haven)
})

# ---------- helpers ----------
sas_sort <- function(df, vars) {
  keys <- lapply(vars, function(v) {
    x <- df[[v]]
    if (is.character(x)) sub(" +$", "", x) else x
  })
  ord <- do.call(order, c(keys, list(na.last = FALSE, method = "radix")))
  df[ord, , drop = FALSE]
}

put_num1 <- function(x, w, d) {
  if (is.na(x)) return(formatC(".", width = w))
  for (dd in d:0) {
    r <- sign(x) * floor(abs(x) * 10^dd + 0.5 + 1e-9) / 10^dd
    r <- r + 0
    s <- sprintf(paste0("%.", dd, "f"), r)
    if (dd > 0 && nchar(s) > w && grepl("^-?0\\.", s)) s <- sub("0\\.", ".", s)
    if (nchar(s) <= w) return(formatC(s, width = w))
  }
  strrep("*", w)
}
put_num <- function(x, w, d = 0) vapply(x, put_num1, character(1), w = w, d = d)

compbl <- function(x) gsub(" {2,}", " ", x)
rtrim <- function(x) sub(" +$", "", x)

pk_fmt <- c(AUCINFO = "AUC", AUCTAU = "AUCtau", AUCLAST = "AUC(0-tlast)",
            CMAX = "Cmax", CMAXD = "Cmax/Dose", MRTLAST = "MRT", TLAST = "tlast",
            TMAX = "tmax", CLRFO = "CL/F", VZFO = "VzF", LAMZ = "lambdaZ",
            LAMZHL = "t1/2", CLAST = "Clast")
put_pk <- function(x) {
  k <- rtrim(ifelse(is.na(x), "", x))
  out <- ifelse(k %in% names(pk_fmt), pk_fmt[k], k)
  out <- substr(out, 1, 12)
  formatC(out, width = -12)
}

# ---------- ADPP ----------
src <- read_xpt("data/adpp.xpt")
src <- as.data.frame(src)

keep <- (!is.na(src$PKFN) & src$PKFN == 1) &
  (is.na(src$ANL01FN) | src$ANL01FN != 1) &
  (rtrim(src$PARAMCD) %in% c("AUCINFO", "AUCLAST", "CMAX", "TMAX", "MRTLAST", "LAMZHL")) &
  (!is.na(src$PPSTRESN)) &
  (gsub(" ", "", ifelse(is.na(src$PPCATL), "", src$PPCATL)) == "ANAL2")
adpp <- src[keep, , drop = FALSE]

uni <- ifelse(is.na(adpp$PPORRESU), "", trimws(adpp$PPORRESU))
adpp$`_PARAM` <- rtrim(compbl(paste0(put_pk(adpp$PPTESTCD), " (", uni, ")")))
adpp$AVAL_log <- ifelse(!is.na(adpp$AVAL) & adpp$AVAL > 0, log(adpp$AVAL), NA_real_)

byv <- c("ANLPRNT", "PPSPECL", "PPCATL", "TRTAN", "PPSDY", "_PARAM")
adpp <- sas_sort(adpp, byv)
rownames(adpp) <- NULL

# ---------- PROC UNIVARIATE ----------
sstat <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  mn <- if (n > 0) mean(x) else NA_real_
  sd_ <- if (n > 1) sd(x) else NA_real_
  cv <- if (!is.na(sd_) && !is.na(mn) && mn != 0) sd_ / mn * 100 else NA_real_
  md <- if (n > 0) median(x) else NA_real_
  mi <- if (n > 0) min(x) else NA_real_
  ma <- if (n > 0) max(x) else NA_real_
  c(n = n, mean = mn, std = sd_, cv = cv, med = md, min = mi, max = ma)
}

grp_key <- do.call(paste, c(lapply(byv, function(v) {
  x <- adpp[[v]]; ifelse(is.na(x), "\001NA", as.character(x))
}), sep = "\002"))
first_idx <- which(!duplicated(grp_key))
adpp_sum <- adpp[first_idx, byv, drop = FALSE]
st <- lapply(first_idx, function(i) {
  idx <- which(grp_key == grp_key[i])
  s1 <- sstat(adpp$AVAL[idx]); s2 <- sstat(adpp$AVAL_log[idx])
  c(n1 = s1[["n"]], n2 = s2[["n"]], mn1 = s1[["mean"]], mn2 = s2[["mean"]],
    std1 = s1[["std"]], std2 = s2[["std"]], cv1 = s1[["cv"]], cv2 = s2[["cv"]],
    med1 = s1[["med"]], med2 = s2[["med"]], min1 = s1[["min"]], min2 = s2[["min"]],
    max1 = s1[["max"]], max2 = s2[["max"]])
})
st <- as.data.frame(do.call(rbind, st))
adpp_sum <- cbind(adpp_sum, st)
rownames(adpp_sum) <- NULL

# ---------- ADPP_SUM_T ----------
P <- adpp_sum$`_PARAM`
has <- function(s) grepl(s, P, fixed = TRUE)
isAUC <- has("AUC"); isC <- has("Cmax"); isM <- has("MRT"); isT <- has("t1/2"); isTm <- has("tmax")
nr <- nrow(adpp_sum)

v1 <- put_num(adpp_sum$n1, 2, 0)

v2 <- rep(" ", nr)
v2 <- ifelse(isAUC, paste0(put_num(adpp_sum$mn1, 6), " (", put_num(adpp_sum$std1, 5), ")"),
      ifelse(isC, paste0(put_num(adpp_sum$mn1, 5), " (", put_num(adpp_sum$std1, 4), ")"),
      ifelse(isM, paste0(put_num(adpp_sum$mn1, 4, 1), " (", put_num(adpp_sum$std1, 4, 2), ")"),
      ifelse(isT, paste0(put_num(adpp_sum$mn1, 4, 2), " (", put_num(adpp_sum$std1, 4, 2), ")"),
             " "))))

v3 <- ifelse(!isTm, put_num(adpp_sum$cv1, 4, 1), " ")

gm <- exp(adpp_sum$mn2)
v4 <- ifelse(isAUC, put_num(gm, 6),
      ifelse(isC, put_num(gm, 5),
      ifelse(isM, put_num(gm, 4, 1),
      ifelse(isT, put_num(gm, 4, 2), " "))))

v5 <- ifelse(!isTm, put_num(sqrt(exp(adpp_sum$std2^2) - 1) * 100, 4, 1), " ")

v6 <- ifelse(isAUC, put_num(adpp_sum$med1, 6),
      ifelse(isC, put_num(adpp_sum$med1, 5),
      ifelse(isM, put_num(adpp_sum$med1, 4, 1),
      ifelse(isT, put_num(adpp_sum$med1, 4, 2),
      ifelse(isTm, put_num(adpp_sum$med1, 4, 2), v5)))))

v7 <- ifelse(isAUC, paste0("[", put_num(adpp_sum$min1, 6), ";", put_num(adpp_sum$max1, 6), "]"),
      ifelse(isC, paste0("[", put_num(adpp_sum$min1, 4), ";", put_num(adpp_sum$max1, 5), "]"),
      ifelse(isM, paste0("[", put_num(adpp_sum$min1, 4, 2), ";", put_num(adpp_sum$max1, 4, 1), "]"),
      ifelse(isT, paste0("[", put_num(adpp_sum$min1, 4, 2), ";", put_num(adpp_sum$max1, 4, 1), "]"),
      ifelse(isTm, paste0("[", put_num(adpp_sum$min1, 4, 2), ";", put_num(adpp_sum$max1, 4, 2), "]"),
             v6)))))

vals <- list(v1, v2, v3, v4, v5, v6, v7)
parts <- lapply(1:7, function(k) {
  d <- adpp_sum
  d$value <- rtrim(substr(vals[[k]], 1, 30))
  d$s_order <- k
  d$.row <- seq_len(nr)
  d
})
adpp_sum_t <- do.call(rbind, parts)
adpp_sum_t <- adpp_sum_t[order(adpp_sum_t$.row, adpp_sum_t$s_order), ]
adpp_sum_t$.row <- NULL
rownames(adpp_sum_t) <- NULL

adpp_sum_t <- sas_sort(adpp_sum_t, c("ANLPRNT", "PPSPECL", "PPCATL", "TRTAN", "PPSDY", "s_order", "_PARAM"))
rownames(adpp_sum_t) <- NULL

# ---------- PROC REPORT OUT=CLXCHECKNOBS (approximation) ----------
plev <- unique(rtrim(adpp_sum_t$`_PARAM`))
plev <- plev[order(plev, method = "radix")]
ccols <- paste0("_C", 6 + seq_along(plev), "_")

rowv <- c("ANLPRNT", "PPSPECL", "PPCATL", "TRTAN", "PPSDY", "s_order")
mk_key <- function(d, vars) do.call(paste, c(lapply(vars, function(v) {
  x <- d[[v]]; ifelse(is.na(x), "\001NA", as.character(x))
}), sep = "\002"))

out_rows <- list()
if (nrow(adpp_sum_t) > 0) {
  bykey <- mk_key(adpp_sum_t, c("ANLPRNT", "PPSPECL", "PPCATL"))
  for (bk in unique(bykey)) {
    sub <- adpp_sum_t[bykey == bk, , drop = FALSE]
    rk <- mk_key(sub, rowv)
    first <- which(!duplicated(rk))
    det <- sub[first, rowv, drop = FALSE]
    for (j in seq_along(plev)) {
      det[[ccols[j]]] <- vapply(first, function(i) {
        m <- which(rk == rk[i] & rtrim(sub$`_PARAM`) == plev[j])
        if (length(m)) sub$value[m[1]] else ""
      }, character(1))
    }
    det$`_BREAK_` <- ""
    # page row
    pg <- det[1, , drop = FALSE]
    pg$TRTAN <- NA; pg$PPSDY <- NA; pg$s_order <- NA
    for (cc in ccols) pg[[cc]] <- ""
    pg$`_BREAK_` <- "_PAGE_"
    # break after PPSDY rows
    bkey <- mk_key(det, c("TRTAN", "PPSDY"))
    pieces <- list(pg)
    for (b in unique(bkey)) {
      dd <- det[bkey == b, , drop = FALSE]
      br <- dd[nrow(dd), , drop = FALSE]
      br$s_order <- NA
      for (cc in ccols) br[[cc]] <- ""
      br$`_BREAK_` <- "PPSDY"
      pieces[[length(pieces) + 1]] <- dd
      pieces[[length(pieces) + 1]] <- br
    }
    out_rows[[length(out_rows) + 1]] <- do.call(rbind, pieces)
  }
  clx <- do.call(rbind, out_rows)
} else {
  clx <- adpp_sum_t[0, rowv, drop = FALSE]
  for (cc in ccols) clx[[cc]] <- character(0)
  clx$`_BREAK_` <- character(0)
}
rownames(clx) <- NULL

saveRDS(list(ADPP = adpp, ADPP_SUM = adpp_sum, ADPP_SUM_T = adpp_sum_t,
             CLXCHECKNOBS = clx), "out.rds")
```