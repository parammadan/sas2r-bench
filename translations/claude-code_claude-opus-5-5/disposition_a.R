suppressPackageStartupMessages({
  library(dplyr)
  library(haven)
})

rtrim <- function(x) { x <- as.character(x); x[is.na(x)] <- ""; sub(" +$", "", x) }

sas_order <- function(df, keys, desc = character()) {
  if (nrow(df) == 0) return(df)
  args <- lapply(keys, function(k) {
    x <- df[[k]]
    if (is.character(x)) {
      x <- rtrim(x)
      r <- match(x, sort(unique(x), method = "radix"))
    } else {
      x <- as.numeric(x); x[is.na(x)] <- -Inf; r <- x
    }
    if (k %in% desc) r <- -r
    r
  })
  df[do.call(order, c(args, list(method = "radix"))), , drop = FALSE]
}

sas_freq <- function(d, vars) {
  keep <- rep(TRUE, nrow(d))
  for (v in vars) {
    x <- d[[v]]
    if (is.character(x)) keep <- keep & rtrim(x) != "" else keep <- keep & !is.na(x)
  }
  d <- d[keep, , drop = FALSE]
  out <- d %>% group_by(across(all_of(vars))) %>% summarise(COUNT = as.numeric(n()), .groups = "drop")
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  sas_order(out, vars)
}

interleave <- function(lst, by) {
  d <- bind_rows(lapply(seq_along(lst), function(i) {
    x <- lst[[i]]; x$.src <- i; x$.rn <- seq_len(nrow(x)); x
  }))
  d <- sas_order(d, c(by, ".src", ".rn"))
  d
}

sas_merge <- function(A, B, by) {
  A <- as.data.frame(A, stringsAsFactors = FALSE)
  B <- as.data.frame(B, stringsAsFactors = FALSE)
  keyf <- function(d) do.call(paste, c(lapply(by, function(b) {
    x <- d[[b]]; if (is.character(x)) rtrim(x) else as.character(x)
  }), sep = "\r"))
  kA <- keyf(A); kB <- keyf(B)
  keysdf <- bind_rows(A[by], B[by])
  keysdf$.k <- keyf(keysdf)
  keysdf <- keysdf[!duplicated(keysdf$.k), , drop = FALSE]
  keysdf <- sas_order(keysdf, by)
  idxA <- c(); idxB <- c(); useB <- c()
  for (k in keysdf$.k) {
    ia <- which(kA == k); ib <- which(kB == k)
    na <- length(ia); nb <- length(ib); n <- max(na, nb)
    for (i in seq_len(n)) {
      idxA <- c(idxA, if (na > 0) ia[min(i, na)] else NA)
      idxB <- c(idxB, if (nb > 0) ib[min(i, nb)] else NA)
      useB <- c(useB, nb > 0 && i <= nb)
    }
  }
  cols <- c(names(A), setdiff(names(B), names(A)))
  out <- list()
  for (cc in cols) {
    inA <- cc %in% names(A); inB <- cc %in% names(B)
    if (inA && inB) {
      va <- A[[cc]][idxA]; vb <- B[[cc]][idxB]
      if (is.character(va) || is.character(vb)) { va <- as.character(va); vb <- as.character(vb) }
      v <- ifelse(useB, vb, va)
    } else if (inA) {
      v <- A[[cc]][idxA]
    } else {
      v <- B[[cc]][idxB]
    }
    if (is.character(v)) v[is.na(v)] <- ""
    out[[cc]] <- v
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}

## ADSL
src <- read_xpt("data/adsl.xpt")
src <- zap_labels(zap_formats(src))
names(src) <- toupper(names(src))
src <- as.data.frame(src, stringsAsFactors = FALSE)
ADSL <- src[rtrim(src$TRT01P) != "", , drop = FALSE]
ADSL$TRT01PN <- ifelse(!is.na(ADSL$TRT01PN) & ADSL$TRT01PN == 0, 3,
                ifelse(!is.na(ADSL$TRT01PN) & ADSL$TRT01PN == 54, 1,
                ifelse(!is.na(ADSL$TRT01PN) & ADSL$TRT01PN == 81, 2, ADSL$TRT01PN)))

TOPLINE <- sas_freq(ADSL, c("STUDYID", "TRT01PN", "TRT01P", "DISCONFL"))

w1 <- ADSL[rtrim(ADSL$DSRAEFL) == "Y" | rtrim(ADSL$DTHFL) == "Y", , drop = FALSE]
MIDLINE1 <- sas_freq(w1, c("STUDYID", "TRT01PN", "TRT01P", "DCREASCD"))
SORTORD1 <- sas_freq(w1, c("STUDYID", "DCREASCD"))

w2 <- ADSL[rtrim(ADSL$DCDECOD) == "LACK OF EFFICACY", , drop = FALSE]
MIDLINE2 <- sas_freq(w2, c("STUDYID", "TRT01PN", "TRT01P", "DCREASCD"))
SORTORD2 <- sas_freq(w2, c("STUDYID", "DCREASCD"))

w3 <- ADSL[rtrim(ADSL$DCDECOD) != "LACK OF EFFICACY" & rtrim(ADSL$DSRAEFL) != "Y" &
           rtrim(ADSL$DTHFL) != "Y" & rtrim(ADSL$DISCONFL) == "Y", , drop = FALSE]
MIDLINE3 <- sas_freq(w3, c("STUDYID", "TRT01PN", "TRT01P", "DCREASCD"))
SORTORD3 <- sas_freq(w3, c("STUDYID", "DCREASCD"))

## ALL
d <- interleave(list(TOPLINE, MIDLINE1, MIDLINE2, MIDLINE3), c("STUDYID", "TRT01PN", "TRT01P"))
d$ORD1 <- c(1, 2, 3, 4)[d$.src]
d$ORD2 <- c(1, 1, 2, 3)[d$.src]
d$DISCONFL[is.na(d$DISCONFL)] <- ""
d$DCREASCD[is.na(d$DCREASCD)] <- ""
ALL <- d[, c("STUDYID", "TRT01PN", "TRT01P", "DISCONFL", "COUNT", "DCREASCD", "ORD1", "ORD2")]
rownames(ALL) <- NULL

## ALLS
d <- interleave(list(SORTORD1, SORTORD2, SORTORD3), c("STUDYID"))
d$ORD1 <- c(2, 3, 4)[d$.src]
d$ORD2 <- c(1, 2, 3)[d$.src]
ALLS <- d[, c("STUDYID", "DCREASCD", "COUNT", "ORD1", "ORD2")]
ALLS <- sas_order(ALLS, c("ORD1", "ORD2", "COUNT"), desc = "COUNT")
rownames(ALLS) <- NULL

ALLS2 <- ALLS
ord3 <- numeric(nrow(ALLS2))
for (i in seq_len(nrow(ALLS2))) {
  first <- i == 1 || !identical(ALLS2$ORD1[i], ALLS2$ORD1[i - 1]) || !identical(ALLS2$ORD2[i], ALLS2$ORD2[i - 1])
  ord3[i] <- if (first) 1 else ord3[i - 1] + 1
}
ALLS2$ORD3 <- ord3

## POP
POP <- sas_freq(ADSL, c("STUDYID", "TRT01PN", "TRT01P"))
names(POP)[names(POP) == "COUNT"] <- "TOT"
rownames(POP) <- NULL

## ALLP
m <- sas_merge(ALL, POP, c("STUDYID", "TRT01PN", "TRT01P"))
m$COUNTC <- ""
ok <- !is.na(m$COUNT)
r <- floor(m$COUNT[ok] / m$TOT[ok] * 10 + 0.5) / 10
m$COUNTC[ok] <- paste0(sprintf("%4d", as.integer(m$COUNT[ok])), " (", sprintf("%5.1f", r), ")")
o1 <- !is.na(m$ORD1) & m$ORD1 == 1
m$DCREASCD[o1] <- ifelse(rtrim(m$DISCONFL[o1]) == "Y", "Discontinued", "Completed the Study")
ALLP <- m[, c("DCREASCD", setdiff(names(m), "DCREASCD"))]
ALLP <- sas_order(ALLP, c("STUDYID", "ORD1", "ORD2", "DCREASCD", "TRT01PN"))
rownames(ALLP) <- NULL

## ALLPT
byv <- c("STUDYID", "ORD1", "ORD2", "DCREASCD")
gk <- do.call(paste, c(lapply(byv, function(b) as.character(ALLP[[b]])), sep = "\r"))
ug <- unique(gk)
ids <- unique(as.character(ALLP$TRT01PN))
ALLPT <- ALLP[match(ug, gk), byv, drop = FALSE]
ALLPT$`_NAME_` <- rep("countc", nrow(ALLPT))
for (id in ids) {
  v <- rep("", length(ug))
  sel <- which(as.character(ALLP$TRT01PN) == id)
  v[match(gk[sel], ug)] <- ALLP$COUNTC[sel]
  ALLPT[[paste0("TRT", id)]] <- v
}
rownames(ALLPT) <- NULL

## ALLPT2
maxt <- nrow(POP)
ALLPT2 <- ALLPT
for (i in seq_len(maxt)) {
  cn <- paste0("TRT", i)
  if (!cn %in% names(ALLPT2)) ALLPT2[[cn]] <- rep("", nrow(ALLPT2))
  x <- rtrim(ALLPT2[[cn]])
  ALLPT2[[cn]] <- ifelse(x == "", "###0", gsub(" ", "#", x))
}

## ALLPT3
m <- sas_merge(ALLPT2, ALLS2, c("STUDYID", "ORD1", "ORD2"))
o1 <- !is.na(m$ORD1) & m$ORD1 == 1
m$ORD3[o1 & toupper(rtrim(m$DCREASCD)) == "DISCONTINUED"] <- 2
m$ORD3[o1 & toupper(rtrim(m$DCREASCD)) != "DISCONTINUED"] <- 1
keepv <- names(m)[names(m) %in% c("STUDYID", "DCREASCD", "ORD1", "ORD2", "ORD3") | grepl("^TRT", toupper(names(m)))]
ALLPT3 <- m[, keepv, drop = FALSE]
ALLPT3 <- sas_order(ALLPT3, c("STUDYID", "ORD1", "ORD2", "ORD3"))
rownames(ALLPT3) <- NULL

## FINAL
FINAL <- ALLPT3
col1 <- rep("", nrow(FINAL))
o1 <- !is.na(FINAL$ORD1) & FINAL$ORD1 == 1
o3 <- !is.na(FINAL$ORD1) & FINAL$ORD1 == 3
col1[o1] <- trimws(FINAL$DCREASCD[o1])
col1[o3] <- paste0("####", trimws(FINAL$DCREASCD[o3]))
FINAL <- cbind(data.frame(COL1 = col1, stringsAsFactors = FALSE), FINAL)
FINAL$PAGEX <- 1

saveRDS(list(
  ADSL = ADSL, TOPLINE = TOPLINE,
  MIDLINE1 = MIDLINE1, SORTORD1 = SORTORD1,
  MIDLINE2 = MIDLINE2, SORTORD2 = SORTORD2,
  MIDLINE3 = MIDLINE3, SORTORD3 = SORTORD3,
  ALL = ALL, ALLS = ALLS, ALLS2 = ALLS2, POP = POP,
  ALLP = ALLP, ALLPT = ALLPT, ALLPT2 = ALLPT2, ALLPT3 = ALLPT3,
  FINAL = FINAL
), "out.rds")
