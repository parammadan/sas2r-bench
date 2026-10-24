suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
})

# ---------- helpers ----------
sas_round <- function(x, d) sign(x) * floor(abs(x) * 10^d + 0.5 + 1e-9) / 10^d

put_wd <- function(x, w, d) {
  out <- character(length(x))
  na <- is.na(x)
  out[na] <- formatC(".", width = w)
  if (any(!na)) {
    s <- formatC(sas_round(x[!na], d), format = "f", digits = d)
    out[!na] <- formatC(s, width = w)
  }
  out
}

rtrim <- function(x) sub("\\s+$", "", x)

sas_order_num <- function(x) order(!is.na(x), x, method = "radix")

q_def5 <- function(x, p) {
  x <- sort(x[!is.na(x)])
  n <- length(x)
  if (n == 0) return(NA_real_)
  np <- n * p
  j <- floor(np + 1e-10)
  g <- np - j
  if (abs(g) < 1e-10) {
    if (j == 0) return(x[1])
    if (j >= n) return(x[n])
    (x[j] + x[j + 1]) / 2
  } else {
    x[min(j + 1, n)]
  }
}

# ---------- data ----------
src <- read_xpt("data/adsl.xpt")
names(src) <- toupper(names(src))

ADSL <- src[!is.na(src$ITTFL) & rtrim(src$ITTFL) == "Y", ]
ADSL <- as.data.frame(ADSL)

# proc sort nodupkey -> treatments
o <- sas_order_num(ADSL$TRT01AN)
tmp <- ADSL[o, ]
tmp <- tmp[!duplicated(tmp$TRT01AN), ]
keepv <- names(ADSL)[names(ADSL) %in% c("TRT01AN", "TRT01A")]
TREATMENTS <- tmp[, keepv, drop = FALSE]
TREATMENTS$`__DUMMYTRT` <- seq_len(nrow(TREATMENTS))
rownames(TREATMENTS) <- NULL

# sort adsl
ADSL <- ADSL[sas_order_num(ADSL$TRT01AN), ]
rownames(ADSL) <- NULL

# merge
ADSL_WORK <- ADSL
m <- match(ADSL_WORK$TRT01AN, TREATMENTS$TRT01AN)
ADSL_WORK$TRT01A <- TREATMENTS$TRT01A[m]
ADSL_WORK$`__DUMMYTRT` <- TREATMENTS$`__DUMMYTRT`[m]

trtnum <- nrow(TREATMENTS)

# ---------- descrstats ----------
destats_lbl <- c(NC = "n[a]", MEANC = "Mean", STDC = "SD", MEDIANC = "Median",
                 Q1Q3C = "Q1, Q3", MINMAXC = "Min, Max", MISSINGC = "Missing")
destats_ord <- c(NC = 0, MEANC = 1, STDC = 2, MEDIANC = 3, Q1Q3C = 4, MINMAXC = 5, MISSINGC = 6)

descrstats <- function(indata, varn, dec, name, ord) {
  ind <- indata[sas_order_num(indata$`__DUMMYTRT`), ]
  rownames(ind) <- NULL
  grps <- unique(ind$`__DUMMYTRT`)
  mr <- do.call(rbind, lapply(grps, function(g) {
    x <- ind[[varn]][ind$`__DUMMYTRT` == g]
    xx <- x[!is.na(x)]
    data.frame(
      `__DUMMYTRT` = g, `_TYPE_` = 0, `_FREQ_` = length(x),
      N = length(xx),
      MEAN = if (length(xx)) mean(xx) else NA_real_,
      STD = if (length(xx) > 1) sd(xx) else NA_real_,
      MIN = if (length(xx)) min(xx) else NA_real_,
      MAX = if (length(xx)) max(xx) else NA_real_,
      Q1 = q_def5(xx, 0.25), Q3 = q_def5(xx, 0.75),
      MEDIAN = q_def5(xx, 0.5),
      check.names = FALSE)
  }))
  d1 <- dec + 1
  mr$NC <- put_wd(mr$N, 6, 0)
  mr$MEANC <- put_wd(mr$MEAN, 6, d1)
  mr$MEDIANC <- put_wd(mr$MEDIAN, 6, d1)
  mr$STDC <- put_wd(mr$STD, 6, d1)
  mr$Q1Q3C <- rtrim(paste0(rtrim(put_wd(mr$Q1, 6, d1)), ", ", rtrim(put_wd(mr$Q3, 6, d1))))
  mr$MINMAXC <- rtrim(paste0(rtrim(put_wd(mr$MIN, 6, dec)), ", ", put_wd(mr$MAX, 6, dec)))

  vars <- c("NC", "MEANC", "STDC", "MINMAXC", "Q1Q3C", "MEDIANC")
  tr <- data.frame(`_NAME_` = tolower(vars), check.names = FALSE, stringsAsFactors = FALSE)
  for (k in seq_len(nrow(mr))) {
    tr[[paste0("COL", mr$`__DUMMYTRT`[k])]] <- rtrim(unlist(mr[k, vars]))
  }
  rownames(tr) <- NULL

  dr <- data.frame(row.names = seq_len(nrow(tr)))
  for (i in seq_len(trtnum)) {
    cn <- paste0("COL", i)
    dr[[cn]] <- if (cn %in% names(tr)) tr[[cn]] else ""
  }
  up <- toupper(tr$`_NAME_`)
  dr$NAME <- name
  dr$COL0 <- unname(destats_lbl[up])
  dr$ORD <- unname(destats_ord[up])
  dr$MAINORD <- ord
  rownames(dr) <- NULL
  list(indata = ind, meansres = mr, trout = tr, descrres = dr)
}

r_age <- descrstats(ADSL_WORK, "AGE", 0, "Age (years)", 2)
AGERES <- r_age$descrres
r_wt <- descrstats(ADSL_WORK, "WEIGHTBL", 0, "Weight (kg)", 5)
WEIGHTRES <- r_wt$descrres
`__INDATA` <- r_wt$indata
`__MEANSRES` <- r_wt$meansres
`__TROUT` <- r_wt$trout
`__DESCRRES` <- r_wt$descrres

# ---------- freqstats ----------
fmt_apply <- function(x, map) {
  w <- max(nchar(map))
  out <- unname(map[rtrim(x)])
  um <- is.na(out)
  out[um] <- substr(x[um], 1, w)
  out
}
infmt_apply <- function(x, map) unname(map[rtrim(x)])

sexfmt <- c(N = "n[a]", F = "F", M = "M")
sexfmto <- c(N = 1, F = 2, M = 3)
agegrfmt <- c(N = "n[a]", "<65" = "<65", "65-80" = "65-80", ">80" = ">80")
racefmt <- c(N = "n[a]", WHITE = "White",
             "BLACK OR AFRICAN AMERICAN" = "Black or African American",
             "AMERICAN INDIAN OR ALASKA NATIVE" = "American Indian or Alaska Native")
racefmto <- c(N = 1, WHITE = 2, "BLACK OR AFRICAN AMERICAN" = 3,
              "AMERICAN INDIAN OR ALASKA NATIVE" = 4)
ethfmt <- c(N = "n[a]", "NOT HISPANIC OR LATINO" = "Not Hispanic or Latino",
            "HISPANIC OR LATINO" = "Hispanic or Latino")
ethfmto <- c(N = 1, "NOT HISPANIC OR LATINO" = 2, "HISPANIC OR LATINO" = 3)

freqstats <- function(indata, var, varfmt, ordfmt, name, ord) {
  v <- indata[[var]]
  nm <- !(is.na(v) | rtrim(v) == "")
  d <- indata[nm, ]
  fr <- d %>% group_by(.data[[var]], `__DUMMYTRT`) %>%
    summarise(COUNT = n(), .groups = "drop") %>% as.data.frame()
  fr <- fr[order(fr[[var]], fr$`__DUMMYTRT`, method = "radix"), ]
  rownames(fr) <- NULL

  denoms <- sapply(seq_len(trtnum), function(i)
    sum(nm & indata$`__DUMMYTRT` == i & !is.na(indata$USUBJID) & indata$USUBJID != ""))

  ids <- unique(fr$`__DUMMYTRT`)
  lv <- unique(fr[[var]])
  tf <- data.frame(x = lv, stringsAsFactors = FALSE)
  names(tf) <- var
  for (id in ids) {
    sub <- fr[fr$`__DUMMYTRT` == id, ]
    tf[[paste0("N", id)]] <- as.numeric(sub$COUNT[match(lv, sub[[var]])])
  }

  den <- data.frame(x = "N", stringsAsFactors = FALSE)
  names(den) <- var
  for (i in seq_len(trtnum)) den[[paste0("N", i)]] <- as.numeric(denoms[i])

  out <- bind_rows(tf, den)
  for (i in seq_len(trtnum)) {
    cn <- paste0("N", i)
    if (!cn %in% names(out)) out[[cn]] <- NA_real_
  }
  out$NAME <- name
  out$COL0 <- NA_character_
  for (i in seq_len(trtnum)) out[[paste0("COL", i)]] <- NA_character_
  for (i in seq_len(trtnum)) {
    cn <- paste0("N", i)
    out[[cn]][is.na(out[[cn]])] <- 0
    n <- out[[cn]]
    dn <- denoms[i]
    pct <- ifelse(n == 0, "",
                  formatC(paste0("(", gsub(" ", "", put_wd(100 * n / dn, 5, 1)), "%)"), width = 9))
    out[[paste0("PCT", i)]] <- pct
    out[[paste0("COL", i)]] <- rtrim(paste0(put_wd(n, 6, 0), formatC(out$PCT1, width = 9)))
  }
  out$MAINORD <- ord
  out$COL0 <- fmt_apply(out[[var]], varfmt)
  out$ORD <- infmt_apply(out[[var]], ordfmt)
  out <- out[, c(var, paste0("N", seq_len(trtnum)), "NAME", "COL0",
                 paste0("COL", seq_len(trtnum)), paste0("PCT", seq_len(trtnum)),
                 "MAINORD", "ORD")]
  rownames(out) <- NULL
  list(frout = fr, trout_fr = tf, denom = den, out = out)
}

r_sex <- freqstats(ADSL_WORK, "SEX", sexfmt, sexfmto, "Sex n(%)", 1)
SEXRES <- r_sex$out

# agegr: informat AGEGRFMTO does not exist -> data step fails, AGEGRRES not created
r_agegr <- freqstats(ADSL_WORK, "AGEGR1", agegrfmt, c(N = 1), "Age categories n(%)", 3)

r_race <- freqstats(ADSL_WORK, "RACE", racefmt, racefmto, "Race n(%)", 4)
RACERES <- r_race$out

r_eth <- freqstats(ADSL_WORK, "ETHNIC", ethfmt, ethfmto, "Ethnicity n(%)", 5)
ETHNICRES <- r_eth$out

`__FROUT` <- r_eth$frout
`__TROUT_FR` <- r_eth$trout_fr
DENOM <- r_eth$denom

# REPORT not created: SET references nonexistent WORK.AGEGRRES

saveRDS(list(
  ADSL = ADSL,
  TREATMENTS = TREATMENTS,
  ADSL_WORK = ADSL_WORK,
  `__INDATA` = `__INDATA`,
  `__MEANSRES` = `__MEANSRES`,
  `__TROUT` = `__TROUT`,
  `__DESCRRES` = `__DESCRRES`,
  AGERES = AGERES,
  WEIGHTRES = WEIGHTRES,
  `__FROUT` = `__FROUT`,
  `__TROUT_FR` = `__TROUT_FR`,
  DENOM = DENOM,
  SEXRES = SEXRES,
  RACERES = RACERES,
  ETHNICRES = ETHNICRES
), "out.rds")
