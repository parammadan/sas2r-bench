```r
library(haven)
library(dplyr)
library(tidyr)
library(stringr)

pkmap <- c(AUCINFO="AUC", AUCTAU="AUCtau", AUCLAST="AUC(0-tlast)", CMAX="Cmax",
           CMAXD="Cmax/Dose", MRTLAST="MRT", TLAST="tlast", TMAX="tmax",
           CLRFO="CL/F", VZFO="VzF", LAMZ="lambdaZ", LAMZHL="t1/2", CLAST="Clast")

sas_order <- function(df, keys) {
  o <- do.call(order, c(unname(as.list(df[keys])), list(na.last = FALSE, method = "radix")))
  df[o, , drop = FALSE]
}

# SAS w.d numeric format
fmtn <- function(x, w, d) {
  vapply(x, function(v) {
    if (is.na(v)) return(formatC(".", width = w))
    for (dd in d:0) {
      s <- sprintf("%.*f", dd, v)
      if (nchar(s) <= w) return(formatC(s, width = w))
    }
    formatC(substr(format(v, digits = max(w - 2, 1)), 1, w), width = w)
  }, character(1), USE.NAMES = FALSE)
}

raw <- read_xpt("data/adpp.xpt")
for (nm in names(raw)) if (is.character(raw[[nm]])) raw[[nm]][is.na(raw[[nm]])] <- ""

adpp <- raw %>%
  filter(PKFN == 1 & (is.na(ANL01FN) | ANL01FN != 1) &
           str_trim(PARAMCD) %in% c("AUCINFO","AUCLAST","CMAX","TMAX","MRTLAST","LAMZHL") &
           !is.na(PPSTRESN) &
           gsub("\\s", "", PPCATL) == "ANAL2")

lab <- ifelse(adpp$PPTESTCD %in% names(pkmap), pkmap[adpp$PPTESTCD], adpp$PPTESTCD)
adpp$`_PARAM` <- gsub(" +", " ", paste0(lab, " (", str_trim(adpp$PPORRESU), ")"))
adpp$AVAL_log <- suppressWarnings(ifelse(!is.na(adpp$AVAL) & adpp$AVAL > 0, log(adpp$AVAL), NA_real_))

bykeys <- c("ANLPRNT","PPSPECL","PPCATL","TRTAN","PPSDY","_PARAM")
adpp <- sas_order(adpp, bykeys)

st <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  m <- if (n > 0) mean(x) else NA_real_
  s <- if (n > 1) sd(x) else NA_real_
  cv <- if (n > 1 && !is.na(m) && m != 0) s / m * 100 else NA_real_
  list(n = n, mn = m, sd = s, cv = cv,
       med = if (n > 0) median(x) else NA_real_,
       min = if (n > 0) min(x) else NA_real_,
       max = if (n > 0) max(x) else NA_real_)
}

grp <- adpp %>% group_by(across(all_of(bykeys)), .drop = TRUE) %>% group_keys()
idx <- adpp %>% group_by(across(all_of(bykeys))) %>% group_rows()
rows <- lapply(seq_along(idx), function(i) {
  a <- st(adpp$AVAL[idx[[i]]]); b <- st(adpp$AVAL_log[idx[[i]]])
  data.frame(n1 = a$n, n2 = b$n, mn1 = a$mn, mn2 = b$mn, std1 = a$sd, std2 = b$sd,
             cv1 = a$cv, cv2 = b$cv, med1 = a$med, med2 = b$med,
             min1 = a$min, min2 = b$min, max1 = a$max, max2 = b$max)
})
adpp_sum <- bind_cols(as.data.frame(grp), bind_rows(rows))
adpp_sum <- sas_order(adpp_sum, bykeys)
rownames(adpp_sum) <- NULL

D <- adpp_sum
P <- D$`_PARAM`
hasAUC <- grepl("AUC", P, fixed = TRUE)
hasCmax <- grepl("Cmax", P, fixed = TRUE)
hasMRT <- grepl("MRT", P, fixed = TRUE)
hasT12 <- grepl("t1/2", P, fixed = TRUE)
hasTmax <- grepl("tmax", P, fixed = TRUE)

pick <- function(auc, cmax, mrt, t12, tmax, other) {
  case_when(hasAUC ~ auc, hasCmax ~ cmax, hasMRT ~ mrt, hasT12 ~ t12,
            hasTmax ~ tmax, TRUE ~ other)
}

v1 <- fmtn(D$n1, 2, 0)
v2 <- pick(paste0(fmtn(D$mn1,6,0)," (",fmtn(D$std1,5,0),")"),
           paste0(fmtn(D$mn1,5,0)," (",fmtn(D$std1,4,0),")"),
           paste0(fmtn(D$mn1,4,1)," (",fmtn(D$std1,4,2),")"),
           paste0(fmtn(D$mn1,4,2)," (",fmtn(D$std1,4,2),")"),
           " ", " ")
v3 <- ifelse(!hasTmax, fmtn(D$cv1, 4, 1), " ")
g <- exp(D$mn2)
v4 <- pick(fmtn(g,6,0), fmtn(g,5,0), fmtn(g,4,1), fmtn(g,4,2), " ", " ")
v5 <- ifelse(!hasTmax, fmtn(sqrt(exp(D$std2^2) - 1) * 100, 4, 1), " ")
v6 <- pick(fmtn(D$med1,6,0), fmtn(D$med1,5,0), fmtn(D$med1,4,1), fmtn(D$med1,4,2),
           fmtn(D$med1,4,2), v5)
v7 <- pick(paste0("[",fmtn(D$min1,6,0),";",fmtn(D$max1,6,0),"]"),
           paste0("[",fmtn(D$min1,4,0),";",fmtn(D$max1,5,0),"]"),
           paste0("[",fmtn(D$min1,4,2),";",fmtn(D$max1,4,1),"]"),
           paste0("[",fmtn(D$min1,4,2),";",fmtn(D$max1,4,1),"]"),
           paste0("[",fmtn(D$min1,4,2),";",fmtn(D$max1,4,2),"]"), v6)

vals <- list(v1, v2, v3, v4, v5, v6, v7)
adpp_sum_t <- bind_rows(lapply(1:7, function(k) {
  x <- D; x$s_order <- k; x$value <- vals[[k]]; x$.i <- seq_len(nrow(D)); x
}))
adpp_sum_t <- adpp_sum_t[order(adpp_sum_t$.i, adpp_sum_t$s_order), ]
adpp_sum_t$.i <- NULL
rownames(adpp_sum_t) <- NULL

adpp_sum_t <- sas_order(adpp_sum_t, c("ANLPRNT","PPSPECL","PPCATL","TRTAN","PPSDY","s_order","_PARAM"))
rownames(adpp_sum_t) <- NULL

# Approximation of PROC REPORT OUT= dataset
params <- sort(unique(adpp_sum_t$`_PARAM`), method = "radix")
wide <- adpp_sum_t %>%
  select(ANLPRNT, PPSPECL, PPCATL, TRTAN, PPSDY, s_order, `_PARAM`, value) %>%
  mutate(`_PARAM` = factor(`_PARAM`, levels = params, labels = paste0("_C", 6 + seq_along(params), "_"))) %>%
  pivot_wider(names_from = `_PARAM`, values_from = value, values_fill = " ")
wide <- sas_order(as.data.frame(wide), c("ANLPRNT","PPSPECL","PPCATL","TRTAN","PPSDY","s_order"))
wide$`_BREAK_` <- ""
CLXCHECKNOBS <- wide
rownames(CLXCHECKNOBS) <- NULL

ADPP <- adpp
ADPP_SUM <- adpp_sum
ADPP_SUM_T <- adpp_sum_t

saveRDS(list(ADPP = ADPP, ADPP_SUM = ADPP_SUM, ADPP_SUM_T = ADPP_SUM_T,
             CLXCHECKNOBS = CLXCHECKNOBS), "out.rds")
```