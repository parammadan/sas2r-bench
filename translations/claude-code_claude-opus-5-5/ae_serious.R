suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(haven)
})

# ---------- helpers ----------
sas_sort <- function(df, cols, desc = rep(FALSE, length(cols))) {
  if (nrow(df) == 0) return(df)
  keys <- lapply(cols, function(cn) {
    x <- df[[cn]]
    if (is.character(x)) { x[is.na(x)] <- ""; x } else { x <- as.numeric(x); x[is.na(x)] <- -Inf; x }
  })
  o <- do.call(order, c(keys, list(method = "radix", decreasing = desc)))
  out <- df[o, , drop = FALSE]
  rownames(out) <- NULL
  out
}
nodupkey <- function(df, cols) {
  df <- sas_sort(df, cols)
  k <- df[, cols, drop = FALSE]
  df <- df[!duplicated(k), , drop = FALSE]
  rownames(df) <- NULL
  df
}
rtrim <- function(x) sub(" +$", "", x)
fmt_p <- function(p) if (is.na(p)) "" else if (p < 1e-4) "<.0001" else sprintf("%.4f", p)
fmt_tp <- function(p) if (is.na(p)) "" else if (p >= 1e-4) sprintf("%.4f", p) else formatC(p, format = "E", digits = 3)

# weighted 2-way freq (sparse), dropping zero/missing weights
wfreq <- function(d, rv, cv, w = NULL) {
  wt <- if (is.null(w)) rep(1, nrow(d)) else d[[w]]
  keep <- !is.na(wt) & wt > 0
  d <- d[keep, , drop = FALSE]; wt <- wt[keep]
  rl <- sort(unique(d[[rv]])); cl <- sort(unique(d[[cv]]))
  tab <- tapply(wt, list(factor(d[[rv]], levels = rl), factor(d[[cv]], levels = cl)), sum)
  tab[is.na(tab)] <- 0
  out <- data.frame(R = rep(rl, each = length(cl)), C = rep(cl, times = length(rl)),
                    COUNT = as.vector(t(tab)))
  names(out)[1:2] <- c(rv, cv)
  out$PERCENT <- out$COUNT / sum(tab) * 100
  list(out = out, tab = tab)
}
fisher_rows <- function(tab, tablab) {
  n <- sum(tab)
  ptab <- exp(sum(lfactorial(rowSums(tab))) + sum(lfactorial(colSums(tab))) -
                lfactorial(n) - sum(lfactorial(tab)))
  p2 <- tryCatch(fisher.test(tab, workspace = 2e7)$p.value, error = function(e) NA_real_)
  data.frame(Table = tablab,
             NAME1 = c("P_TABLE", "XP2_FISH"),
             LABEL1 = c("Table Probability (P)", "Pr <= P"),
             CVALUE1 = c(fmt_tp(ptab), fmt_p(p2)),
             NVALUE1 = c(ptab, p2),
             stringsAsFactors = FALSE)
}

# ---------- inputs ----------
src_adae <- read_xpt("data/adae.xpt"); names(src_adae) <- toupper(names(src_adae))
src_adsl <- read_xpt("data/adsl.xpt"); names(src_adsl) <- toupper(names(src_adsl))
src_adae <- as.data.frame(zap_formats(zap_labels(src_adae)), stringsAsFactors = FALSE)
src_adsl <- as.data.frame(zap_formats(zap_labels(src_adsl)), stringsAsFactors = FALSE)

# ---------- ADAE ----------
keepv <- c("USUBJID", "TRTAN", "AEDECOD", "AESER", "AOCC02FL", "AOCC04FL", "SAFFL")
ADAE <- src_adae[, names(src_adae)[names(src_adae) %in% keepv], drop = FALSE]
ADAE <- ADAE[(rtrim(ADAE$AOCC02FL) %in% "Y" | rtrim(ADAE$AOCC04FL) %in% "Y") &
               rtrim(ADAE$AESER) %in% "Y" & rtrim(ADAE$SAFFL) %in% "Y", , drop = FALSE]
ADAE <- sas_sort(ADAE, "USUBJID")

# ---------- ADSL ----------
ADSL <- src_adsl[, names(src_adsl)[names(src_adsl) %in% c("USUBJID", "TRT01AN")], drop = FALSE]
ADSL <- ADSL[!is.na(ADSL$TRT01AN), , drop = FALSE]
ADSL <- nodupkey(ADSL, c("USUBJID", "TRT01AN"))

# ---------- ADAE_SL ----------
ADAE_SL <- inner_join(ADSL, ADAE, by = "USUBJID", multiple = "all")
ADAE_SL <- sas_sort(ADAE_SL, "USUBJID")

# ---------- SUM_1 ----------
SUM_1 <- ADSL %>% count(TRT01AN, name = "COUNT") %>% arrange(TRT01AN) %>% as.data.frame()
trtnums <- SUM_1$TRT01AN
trtcn <- setNames(SUM_1$COUNT, as.character(SUM_1$TRT01AN))

# ---------- ADAE_SL2, ADAE_DUMMY ----------
ADAE_SL2 <- nodupkey(ADAE_SL, c("USUBJID", "AEDECOD"))
tmp <- ADAE_SL[rtrim(ADAE_SL$AOCC02FL) %in% "Y", , drop = FALSE]
ADAE_DUMMY <- nodupkey(tmp, "AEDECOD")[, "AEDECOD", drop = FALSE]
ADAE_SL2 <- sas_sort(ADAE_SL2, "AEDECOD")

ADAE_SL3 <- ADAE_SL2[ADAE_SL2$AEDECOD %in% ADAE_DUMMY$AEDECOD, , drop = FALSE]
rownames(ADAE_SL3) <- NULL

# ---------- SUM_2 ----------
SUM_2 <- ADAE_SL3 %>% count(AEDECOD, TRTAN, name = "COUNT") %>% as.data.frame()
SUM_2 <- sas_sort(SUM_2, c("AEDECOD", "TRTAN"))

# ---------- DUMMY_2 ----------
DUMMY_2 <- data.frame(AEDECOD = rep(ADAE_DUMMY$AEDECOD, each = length(trtnums)),
                      TRTAN = rep(trtnums, times = nrow(ADAE_DUMMY)),
                      stringsAsFactors = FALSE)

# ---------- ADAE_SL4, ADAE_SL5 ----------
ADAE_SL4 <- full_join(DUMMY_2, SUM_2, by = c("AEDECOD", "TRTAN"))
ADAE_SL4 <- sas_sort(ADAE_SL4, c("AEDECOD", "TRTAN"))
ADAE_SL4 <- sas_sort(ADAE_SL4, "TRTAN")

s1r <- SUM_1; names(s1r) <- c("TRTAN", "BIGN")
ADAE_SL5 <- left_join(ADAE_SL4, s1r, by = "TRTAN")

# ---------- ADAE_SL6 ----------
idx <- rep(seq_len(nrow(ADAE_SL5)), each = 2)
ADAE_SL6 <- ADAE_SL5[idx, , drop = FALSE]
ADAE_SL6$J <- rep(c(0, 1), times = nrow(ADAE_SL5))
ADAE_SL6$RESP <- ADAE_SL6$J
ADAE_SL6$COUNT_FINAL <- with(ADAE_SL6, ifelse(RESP == 1, ifelse(is.na(COUNT), 0, COUNT),
                                              ifelse(is.na(COUNT), BIGN, BIGN - COUNT)))
rownames(ADAE_SL6) <- NULL
ADAE_SL6 <- sas_sort(ADAE_SL6, "AEDECOD")

# ---------- FREQ by aedecod ----------
fl <- list(); sl <- list()
for (g in unique(ADAE_SL6$AEDECOD)) {
  d <- ADAE_SL6[ADAE_SL6$AEDECOD == g, , drop = FALSE]
  r <- wfreq(d, "TRTAN", "RESP", "COUNT_FINAL")
  sl[[length(sl) + 1]] <- cbind(data.frame(AEDECOD = g, stringsAsFactors = FALSE), r$out)
  fr <- fisher_rows(r$tab, "Table trtan * resp")
  fl[[length(fl) + 1]] <- cbind(data.frame(AEDECOD = g, stringsAsFactors = FALSE), fr)
}
FISH1 <- bind_rows(fl)
SUM_FINAL <- bind_rows(sl)

FISH2 <- FISH1[FISH1$NAME1 == "XP2_FISH", c("AEDECOD", "CVALUE1"), drop = FALSE]
rownames(FISH2) <- NULL

SUM_FINAL <- SUM_FINAL[SUM_FINAL$RESP %in% 1, , drop = FALSE]
SUM_FINAL$CP <- paste0(sprintf("%.0f", SUM_FINAL$COUNT), "(", sprintf("%.1f", SUM_FINAL$COUNT), ")")
SUM_FINAL$ORD <- ifelse(SUM_FINAL$TRTAN %in% 81, 1, 2)
SUM_FINAL <- sas_sort(SUM_FINAL, c("AEDECOD", "ORD"))
cm <- NA_real_; cmv <- numeric(nrow(SUM_FINAL))
for (i in seq_len(nrow(SUM_FINAL))) {
  if (SUM_FINAL$ORD[i] == 1) cm <- SUM_FINAL$COUNT[i]
  cmv[i] <- cm
}
SUM_FINAL$COUNT_MAX <- cmv
SUM_FINAL <- sas_sort(SUM_FINAL, c("AEDECOD", "COUNT_MAX"))

SUM_FINAL2 <- SUM_FINAL %>%
  select(AEDECOD, COUNT_MAX, TRTAN, CP) %>%
  pivot_wider(id_cols = c(AEDECOD, COUNT_MAX), names_from = TRTAN, values_from = CP,
              names_prefix = "TRT") %>%
  as.data.frame(stringsAsFactors = FALSE)
SUM_FINAL2 <- cbind(SUM_FINAL2[, c("AEDECOD", "COUNT_MAX")], `_NAME_` = "cp",
                    SUM_FINAL2[, setdiff(names(SUM_FINAL2), c("AEDECOD", "COUNT_MAX")), drop = FALSE],
                    stringsAsFactors = FALSE)
SUM_FINAL2 <- full_join(SUM_FINAL2, FISH2, by = "AEDECOD")
SUM_FINAL2 <- sas_sort(SUM_FINAL2, "AEDECOD")
SUM_FINAL2 <- sas_sort(SUM_FINAL2, c("COUNT_MAX", "AEDECOD"), c(TRUE, FALSE))

# ---------- subject-level ----------
tmp <- ADAE_SL[rtrim(ADAE_SL$AOCC02FL) %in% "Y", , drop = FALSE]
ADAE_SL_SUB <- nodupkey(tmp, "USUBJID")[, "USUBJID", drop = FALSE]

ADAE_SL_SUB2 <- ADSL
ADAE_SL_SUB2$FLG <- ifelse(ADAE_SL_SUB2$USUBJID %in% ADAE_SL_SUB$USUBJID, 1, 0)

d <- ADAE_SL_SUB2; names(d)[names(d) == "TRT01AN"] <- "TRTAN"
r <- wfreq(d, "TRTAN", "FLG")
FISH_SUB <- fisher_rows(r$tab, "Table trtan * flg")
SUM_FINAL_SUB <- r$out

FISH_SUB2 <- FISH_SUB[FISH_SUB$NAME1 == "XP2_FISH", "CVALUE1", drop = FALSE]
rownames(FISH_SUB2) <- NULL

SUM_FINAL_SUB <- SUM_FINAL_SUB[SUM_FINAL_SUB$FLG %in% 1, , drop = FALSE]
SUM_FINAL_SUB$CP <- paste0(sprintf("%.0f", SUM_FINAL_SUB$COUNT), "(", sprintf("%.1f", SUM_FINAL_SUB$COUNT), ")")
SUM_FINAL_SUB <- cbind(data.frame(AEDECOD = "Number of subjects reporting serious adverse events",
                                  stringsAsFactors = FALSE), SUM_FINAL_SUB)
rownames(SUM_FINAL_SUB) <- NULL
SUM_FINAL_SUB <- sas_sort(SUM_FINAL_SUB, "AEDECOD")

SUM_FINAL_SUB2 <- SUM_FINAL_SUB %>%
  select(AEDECOD, TRTAN, CP) %>%
  pivot_wider(id_cols = AEDECOD, names_from = TRTAN, values_from = CP, names_prefix = "TRT") %>%
  as.data.frame(stringsAsFactors = FALSE)
SUM_FINAL_SUB2 <- cbind(SUM_FINAL_SUB2[, "AEDECOD", drop = FALSE], `_NAME_` = "cp",
                        SUM_FINAL_SUB2[, setdiff(names(SUM_FINAL_SUB2), "AEDECOD"), drop = FALSE],
                        stringsAsFactors = FALSE)
# one-to-one merge
n1 <- max(nrow(SUM_FINAL_SUB2), nrow(FISH_SUB2))
a <- SUM_FINAL_SUB2[seq_len(n1), , drop = FALSE]
b <- FISH_SUB2[seq_len(n1), , drop = FALSE]
SUM_FINAL_SUB2 <- cbind(a, CVALUE1 = b$CVALUE1, stringsAsFactors = FALSE)
rownames(SUM_FINAL_SUB2) <- NULL

FINAL_REPORT <- bind_rows(SUM_FINAL_SUB2, SUM_FINAL2)

saveRDS(list(
  ADAE = ADAE, ADSL = ADSL, ADAE_SL = ADAE_SL, SUM_1 = SUM_1,
  ADAE_SL2 = ADAE_SL2, ADAE_DUMMY = ADAE_DUMMY, ADAE_SL3 = ADAE_SL3, SUM_2 = SUM_2,
  DUMMY_2 = DUMMY_2, ADAE_SL4 = ADAE_SL4, ADAE_SL5 = ADAE_SL5, ADAE_SL6 = ADAE_SL6,
  FISH1 = FISH1, FISH2 = FISH2, SUM_FINAL = SUM_FINAL, SUM_FINAL2 = SUM_FINAL2,
  ADAE_SL_SUB = ADAE_SL_SUB, ADAE_SL_SUB2 = ADAE_SL_SUB2, FISH_SUB = FISH_SUB,
  SUM_FINAL_SUB = SUM_FINAL_SUB, FISH_SUB2 = FISH_SUB2, SUM_FINAL_SUB2 = SUM_FINAL_SUB2,
  FINAL_REPORT = FINAL_REPORT
), "out.rds")
