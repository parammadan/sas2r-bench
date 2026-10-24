```r
library(haven)
library(dplyr)
library(tidyr)
library(stringr)

adsl0 <- read_xpt("data/adsl.xpt")
adsl <- adsl0 %>% filter(ITTFL == "Y")
names(adsl) <- toupper(names(adsl))

# proc sort by trt01an (missing first, stable)
adsl <- adsl[order(adsl$TRT01AN, na.last = FALSE), ]
ADSL <- adsl

treatments <- adsl[!duplicated(adsl$TRT01AN), c("TRT01AN", "TRT01A")]
treatments$`__dummytrt` <- seq_len(nrow(treatments))
TREATMENTS <- treatments

adsl_work <- adsl %>%
  select(-any_of("TRT01A")) %>%
  left_join(treatments, by = "TRT01AN")
# restore TRT01A in original position-ish (values identical)
adsl_work$TRT01A <- treatments$TRT01A[match(adsl_work$TRT01AN, treatments$TRT01AN)]
adsl_work <- adsl_work[, c(names(adsl), "__dummytrt")]
ADSL_WORK <- adsl_work

trtnum <- nrow(treatments)

fmtnum <- function(x, w, d) {
  r <- sign(x) * floor(abs(x) * 10^d + 0.5 + 1e-9) / 10^d
  s <- formatC(r, format = "f", digits = d, width = w)
  s[is.na(x)] <- formatC(".", width = w)
  s
}

destats <- c(NC = "n[a]", MEANC = "Mean", STDC = "SD", MEDIANC = "Median",
             Q1Q3C = "Q1, Q3", MINMAXC = "Min, Max", MISSINGC = "Missing")
destatso <- c(NC = 0, MEANC = 1, STDC = 2, MEDIANC = 3, Q1Q3C = 4, MINMAXC = 5, MISSINGC = 6)

descrstats <- function(indata, varn, dec, name, ord) {
  indata_s <- indata[order(indata$`__dummytrt`), ]
  groups <- sort(unique(indata_s$`__dummytrt`))
  rows <- lapply(groups, function(g) {
    x <- indata_s[[varn]][indata_s$`__dummytrt` == g]
    fr <- length(x)
    xn <- x[!is.na(x)]
    n <- length(xn)
    q <- if (n > 0) quantile(xn, c(0.25, 0.75), type = 2, names = FALSE) else c(NA_real_, NA_real_)
    data.frame(`__dummytrt` = g, `_TYPE_` = 0, `_FREQ_` = fr, n = n,
               mean = if (n > 0) mean(xn) else NA_real_,
               std = if (n > 1) sd(xn) else NA_real_,
               min = if (n > 0) min(xn) else NA_real_,
               max = if (n > 0) max(xn) else NA_real_,
               q1 = q[1], q3 = q[2],
               median = if (n > 0) median(xn) else NA_real_,
               check.names = FALSE)
  })
  m <- bind_rows(rows)
  m$nc <- fmtnum(m$n, 6, 0)
  m$meanc <- fmtnum(m$mean, 6, dec + 1)
  m$medianc <- fmtnum(m$median, 6, dec + 1)
  m$stdc <- fmtnum(m$std, 6, dec + 1)
  m$q1q3c <- paste0(str_trim(fmtnum(m$q1, 6, dec + 1), "right"), ", ",
                    str_trim(fmtnum(m$q3, 6, dec + 1), "right"))
  m$minmaxc <- paste0(str_trim(fmtnum(m$min, 6, dec), "right"), ", ", fmtnum(m$max, 6, dec))

  vars <- c("nc", "meanc", "stdc", "minmaxc", "q1q3c", "medianc")
  tr <- data.frame(`_NAME_` = vars, check.names = FALSE, stringsAsFactors = FALSE)
  for (k in seq_along(groups)) {
    tr[[paste0("col", groups[k])]] <- vapply(vars, function(v) m[[v]][k], character(1), USE.NAMES = FALSE)
  }
  descr <- data.frame(row.names = seq_len(nrow(tr)))
  for (i in seq_len(trtnum)) descr[[paste0("col", i)]] <- tr[[paste0("col", i)]]
  descr$name <- name
  descr$col0 <- unname(destats[toupper(tr$`_NAME_`)])
  descr$ord <- unname(destatso[toupper(tr$`_NAME_`)])
  descr$mainord <- ord
  list(indata = indata_s, meansres = m, trout = tr, descrres = descr)
}

age <- descrstats(adsl_work, "AGE", 0, "Age (years)", 2)
AGERES <- age$descrres
wt <- descrstats(adsl_work, "WEIGHTBL", 0, "Weight (kg)", 5)
WEIGHTRES <- wt$descrres
`__INDATA` <- wt$indata
`__MEANSRES` <- wt$meansres
`__TROUT` <- wt$trout
`__DESCRRES` <- wt$descrres

freqstats <- function(indata, var, fmtmap, ordmap, name, ord, ordfmt_exists = TRUE) {
  v <- indata[[var]]
  ok <- !is.na(v) & v != ""
  d <- indata[ok, ]
  fr <- d %>% count(.data[[var]], `__dummytrt`, name = "COUNT")
  fr <- fr[order(fr[[var]], fr$`__dummytrt`, method = "radix"), ]
  frout <- as.data.frame(fr)
  denoms <- vapply(seq_len(trtnum), function(i) sum(d$`__dummytrt` == i), numeric(1))
  vals <- unique(frout[[var]])
  ids <- unique(frout$`__dummytrt`)
  tr <- data.frame(tmp = vals, stringsAsFactors = FALSE)
  names(tr) <- var
  for (id in ids) {
    sub <- frout[frout$`__dummytrt` == id, ]
    tr[[paste0("n", id)]] <- sub$COUNT[match(vals, sub[[var]])]
  }
  denom <- data.frame(tmp = "N", stringsAsFactors = FALSE)
  names(denom) <- var
  for (i in seq_len(trtnum)) denom[[paste0("n", i)]] <- denoms[i]
  res <- bind_rows(tr, denom)
  res$name <- name
  res$col0 <- ifelse(res[[var]] %in% names(fmtmap), unname(fmtmap[res[[var]]]), res[[var]])
  pcts <- list()
  for (i in seq_len(trtnum)) {
    nn <- res[[paste0("n", i)]]
    nn[is.na(nn)] <- 0
    res[[paste0("n", i)]] <- nn
    p <- ifelse(nn == 0, "",
                paste0("(", str_trim(fmtnum(100 * nn / denoms[i], 5, 1)), "%)"))
    pcts[[i]] <- sprintf("%9s", p)
  }
  for (i in seq_len(trtnum)) {
    res[[paste0("col", i)]] <- str_trim(paste0(fmtnum(res[[paste0("n", i)]], 6, 0), pcts[[1]]), "right")
  }
  for (i in seq_len(trtnum)) res[[paste0("pct", i)]] <- pcts[[i]]
  res$mainord <- ord
  if (ordfmt_exists) {
    res$ord <- unname(ordmap[res[[var]]])
  }
  front <- c(var, grep("^n[0-9]+$", names(res), value = TRUE), "name", "col0",
             paste0("col", seq_len(trtnum)), paste0("pct", seq_len(trtnum)), "mainord")
  if (ordfmt_exists) front <- c(front, "ord")
  res <- res[, front]
  list(frout = frout, trout_fr = tr, denom = denom, res = if (ordfmt_exists) res else NULL)
}

sexf <- freqstats(adsl_work, "SEX",
                  c(N = "n[a]", F = "F", M = "M"),
                  c(N = 1, F = 2, M = 3), "Sex n(%)", 1)
SEXRES <- sexf$res

# agegrfmto informat does not exist in SAS -> data step fails, AGEGRRES not created
agef <- freqstats(adsl_work, "AGEGR1",
                  c(N = "n[a]", "<65" = "<65", "65-80" = "65-80", ">80" = ">80"),
                  NULL, "Age categories n(%)", 3, ordfmt_exists = FALSE)

racef <- freqstats(adsl_work, "RACE",
                   c(N = "n[a]", WHITE = "White", `BLACK OR AFRICAN AMERICAN` = "Black or African American",
                     `AMERICAN INDIAN OR ALASKA NATIVE` = "American Indian or Alaska Native"),
                   c(N = 1, WHITE = 2, `BLACK OR AFRICAN AMERICAN` = 3, `AMERICAN INDIAN OR ALASKA NATIVE` = 4),
                   "Race n(%)", 4)
RACERES <- racef$res

ethf <- freqstats(adsl_work, "ETHNIC",
                  c(N = "n[a]", `NOT HISPANIC OR LATINO` = "Not Hispanic or Latino",
                    `HISPANIC OR LATINO` = "Hispanic or Latino"),
                  c(N = 1, `NOT HISPANIC OR LATINO` = 2, `HISPANIC OR LATINO` = 3),
                  "Ethnicity n(%)", 5)
ETHNICRES <- ethf$res
`__FROUT` <- ethf$frout
`__TROUT_FR` <- ethf$trout_fr
DENOM <- ethf$denom

# REPORT is not created: AGEGRRES does not exist, so the SET statement fails

saveRDS(list(ADSL = ADSL, TREATMENTS = TREATMENTS, ADSL_WORK = ADSL_WORK,
             `__INDATA` = `__INDATA`, `__MEANSRES` = `__MEANSRES`,
             `__TROUT` = `__TROUT`, `__DESCRRES` = `__DESCRRES`,
             AGERES = AGERES, WEIGHTRES = WEIGHTRES,
             `__FROUT` = `__FROUT`, `__TROUT_FR` = `__TROUT_FR`, DENOM = DENOM,
             SEXRES = SEXRES, RACERES = RACERES, ETHNICRES = ETHNICRES),
        "out.rds")
```