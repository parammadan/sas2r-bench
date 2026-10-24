library(dplyr)
library(tidyr)

ordr <- function(df, ...) df[order(..., method = "radix"), , drop = FALSE]
pfmt <- function(p) if (p < 0.0001) "<.0001" else sprintf("%.4f", p)
fish_p <- function(m) {
  if (any(rowSums(m) == 0) || any(colSums(m) == 0)) return(NA_character_)
  r <- tryCatch(fisher.test(m)$p.value, error = function(e) NA_real_)
  if (is.na(r)) NA_character_ else pfmt(r)
}

adae_raw <- haven::read_xpt("data/adae.xpt")
names(adae_raw) <- tolower(names(adae_raw))
adsl_raw <- haven::read_xpt("data/adsl.xpt")
names(adsl_raw) <- tolower(names(adsl_raw))

adae <- adae_raw %>%
  filter((aocc02fl %in% "Y" | aocc04fl %in% "Y") & aeser %in% "Y" & saffl %in% "Y") %>%
  select(usubjid, trtan, aedecod, aeser, aocc02fl, aocc04fl, saffl)
adae <- adae[, c("usubjid","trtan","aedecod","aeser","aocc02fl","aocc04fl","saffl")]
adae <- ordr(adae, adae$usubjid)

adsl <- adsl_raw %>% filter(!is.na(trt01an)) %>% select(usubjid, trt01an)
adsl <- ordr(adsl, adsl$usubjid, adsl$trt01an)
adsl <- adsl[!duplicated(adsl[, c("usubjid", "trt01an")]), ]

adae_sl <- inner_join(adsl, adae, by = "usubjid", relationship = "many-to-many")
adae_sl <- ordr(adae_sl, adae_sl$usubjid)

sum_1 <- adsl %>% count(trt01an, name = "count") %>% mutate(count = as.numeric(count))
sum_1 <- ordr(sum_1, sum_1$trt01an)
trtnum <- sum_1$trt01an

adae_sl2 <- ordr(adae_sl, adae_sl$usubjid, adae_sl$aedecod)
adae_sl2 <- adae_sl2[!duplicated(adae_sl2[, c("usubjid", "aedecod")]), ]

d <- adae_sl[adae_sl$aocc02fl %in% "Y", ]
d <- ordr(d, d$aedecod)
d <- d[!duplicated(d$aedecod), ]
adae_dummy <- d[, "aedecod", drop = FALSE]

adae_sl2 <- ordr(adae_sl2, adae_sl2$aedecod)
adae_sl3 <- adae_sl2[adae_sl2$aedecod %in% adae_dummy$aedecod, ]

sum_2 <- adae_sl3 %>% group_by(aedecod, trtan) %>% summarise(count = as.numeric(n()), .groups = "drop") %>%
  as.data.frame()
sum_2 <- ordr(sum_2, sum_2$aedecod, sum_2$trtan)

dummy_2 <- data.frame(aedecod = rep(adae_dummy$aedecod, each = length(trtnum)),
                      trtan = rep(trtnum, times = nrow(adae_dummy)), stringsAsFactors = FALSE)

adae_sl4 <- full_join(dummy_2, sum_2, by = c("aedecod", "trtan"))
adae_sl4 <- ordr(adae_sl4, adae_sl4$aedecod, adae_sl4$trtan)
adae_sl4 <- ordr(adae_sl4, adae_sl4$trtan)

adae_sl5 <- left_join(adae_sl4, rename(sum_1, trtan = trt01an, bign = count), by = "trtan")

x <- adae_sl5[rep(seq_len(nrow(adae_sl5)), each = 2), ]
x$j <- rep(0:1, times = nrow(adae_sl5))
x$resp <- x$j
x$count_final <- ifelse(x$resp == 1,
                        ifelse(is.na(x$count), 0, x$count),
                        ifelse(is.na(x$count), x$bign, x$bign - x$count))
rownames(x) <- NULL
adae_sl6 <- ordr(x, x$aedecod)

# PROC FREQ with weight, by aedecod
sf <- adae_sl6 %>% group_by(aedecod, trtan, resp) %>%
  summarise(COUNT = sum(count_final), .groups = "drop") %>% as.data.frame()
sf <- ordr(sf, sf$aedecod, sf$trtan, sf$resp)
sf <- sf %>% group_by(aedecod) %>% mutate(PERCENT = COUNT / sum(COUNT) * 100) %>% ungroup() %>% as.data.frame()

fish2 <- do.call(rbind, lapply(split(adae_sl6, adae_sl6$aedecod), function(g) {
  tb <- tapply(g$count_final, list(g$trtan, g$resp), sum)
  data.frame(aedecod = g$aedecod[1], cvalue1 = fish_p(tb), stringsAsFactors = FALSE)
}))
fish2 <- fish2[!is.na(fish2$cvalue1), ]
rownames(fish2) <- NULL
fish2 <- ordr(fish2, fish2$aedecod)

sum_final <- sf[sf$resp == 1, ]
sum_final$cp <- paste0(trimws(sprintf("%.0f", sum_final$COUNT)), "(", trimws(sprintf("%.1f", sum_final$COUNT)), ")")
sum_final$ord <- ifelse(sum_final$trtan == 81, 1, 2)
sum_final <- ordr(sum_final, sum_final$aedecod, sum_final$ord)
sum_final$count_max <- ifelse(sum_final$ord == 1, sum_final$COUNT, NA_real_)
sum_final <- fill(sum_final, count_max, .direction = "down")
sum_final <- ordr(sum_final, sum_final$aedecod, sum_final$count_max)
rownames(sum_final) <- NULL

sum_final2 <- sum_final %>%
  pivot_wider(id_cols = c(aedecod, count_max), names_from = trtan, names_prefix = "trt", values_from = cp) %>%
  as.data.frame()
sum_final2$"_NAME_" <- "cp"
sum_final2 <- sum_final2[, c("aedecod", "count_max", "_NAME_", setdiff(names(sum_final2), c("aedecod", "count_max", "_NAME_")))]
sum_final2 <- left_join(sum_final2, fish2, by = "aedecod")
sum_final2 <- sum_final2[order(-sum_final2$count_max, sum_final2$aedecod, method = "radix"), ]
rownames(sum_final2) <- NULL

adae_sl_sub <- adae_sl[adae_sl$aocc02fl %in% "Y", ]
adae_sl_sub <- ordr(adae_sl_sub, adae_sl_sub$usubjid)
adae_sl_sub <- adae_sl_sub[!duplicated(adae_sl_sub$usubjid), "usubjid", drop = FALSE]

adae_sl_sub2 <- adsl
adae_sl_sub2$flg <- ifelse(adae_sl_sub2$usubjid %in% adae_sl_sub$usubjid, 1, 0)

s <- adae_sl_sub2 %>% rename(trtan = trt01an)
sfs <- s %>% group_by(trtan, flg) %>% summarise(COUNT = as.numeric(n()), .groups = "drop") %>% as.data.frame()
sfs <- tidyr::complete(sfs, trtan, flg, fill = list(COUNT = 0)) %>% as.data.frame()
sfs <- ordr(sfs, sfs$trtan, sfs$flg)
sfs$PERCENT <- sfs$COUNT / sum(sfs$COUNT) * 100
sum_final_sub <- sfs
tb <- tapply(sfs$COUNT, list(sfs$trtan, sfs$flg), sum)
pv <- fish_p(tb)
fish_sub2 <- data.frame(cvalue1 = pv, stringsAsFactors = FALSE)
if (is.na(pv)) fish_sub2 <- fish_sub2[0, , drop = FALSE]

sum_final_sub <- sum_final_sub[sum_final_sub$flg == 1, ]
sum_final_sub$cp <- paste0(trimws(sprintf("%.0f", sum_final_sub$COUNT)), "(", trimws(sprintf("%.1f", sum_final_sub$COUNT)), ")")
sum_final_sub$aedecod <- "Number of subjects reporting serious adverse events"
sum_final_sub <- sum_final_sub[, c("aedecod", "trtan", "flg", "COUNT", "PERCENT", "cp")]
rownames(sum_final_sub) <- NULL

sum_final_sub2 <- sum_final_sub %>%
  pivot_wider(id_cols = aedecod, names_from = trtan, names_prefix = "trt", values_from = cp) %>%
  as.data.frame()
sum_final_sub2$"_NAME_" <- "cp"
sum_final_sub2 <- sum_final_sub2[, c("aedecod", "_NAME_", setdiff(names(sum_final_sub2), c("aedecod", "_NAME_")))]
sum_final_sub2$cvalue1 <- if (nrow(fish_sub2) > 0) fish_sub2$cvalue1[1] else NA_character_

final_report <- bind_rows(sum_final_sub2, sum_final2)

up <- function(df) { df <- as.data.frame(df); names(df) <- toupper(names(df)); rownames(df) <- NULL; df }

saveRDS(list(
  ADAE = up(adae), ADSL = up(adsl), ADAE_SL = up(adae_sl), SUM_1 = up(sum_1),
  ADAE_SL2 = up(adae_sl2), ADAE_DUMMY = up(adae_dummy), ADAE_SL3 = up(adae_sl3),
  SUM_2 = up(sum_2), DUMMY_2 = up(dummy_2), ADAE_SL4 = up(adae_sl4),
  ADAE_SL5 = up(adae_sl5), ADAE_SL6 = up(adae_sl6),
  FISH2 = up(fish2), SUM_FINAL = up(sum_final), SUM_FINAL2 = up(sum_final2),
  ADAE_SL_SUB = up(adae_sl_sub), ADAE_SL_SUB2 = up(adae_sl_sub2),
  FISH_SUB2 = up(fish_sub2), SUM_FINAL_SUB = up(sum_final_sub),
  SUM_FINAL_SUB2 = up(sum_final_sub2), FINAL_REPORT = up(final_report)
), "out.rds")
