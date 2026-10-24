suppressPackageStartupMessages({
  library(dplyr)
  library(haven)
})

ADSL <- haven::read_xpt("data/adsl.xpt")
ADAE <- haven::read_xpt("data/adae.xpt")

# ---- helpers ----
sas_arrange <- function(df, ...) {
  keys <- c(...)
  ord_args <- list()
  for (k in keys) {
    v <- df[[k]]
    ord_args[[length(ord_args) + 1]] <- !is.na(v)
    if (is.character(v)) {
      ord_args[[length(ord_args) + 1]] <- rank(v, ties.method = "min", na.last = TRUE)
    } else {
      ord_args[[length(ord_args) + 1]] <- v
    }
  }
  if (is.character(df[[keys[1]]])) {
    # use C-locale ordering for character keys
    ord_args <- list()
    for (k in keys) {
      v <- df[[k]]
      ord_args[[length(ord_args) + 1]] <- !is.na(v)
      if (is.character(v)) {
        u <- sort(unique(v[!is.na(v)]), method = "radix")
        ord_args[[length(ord_args) + 1]] <- match(v, u)
      } else {
        ord_args[[length(ord_args) + 1]] <- v
      }
    }
  } else {
    ord_args <- list()
    for (k in keys) {
      v <- df[[k]]
      ord_args[[length(ord_args) + 1]] <- !is.na(v)
      if (is.character(v)) {
        u <- sort(unique(v[!is.na(v)]), method = "radix")
        ord_args[[length(ord_args) + 1]] <- match(v, u)
      } else {
        ord_args[[length(ord_args) + 1]] <- v
      }
    }
  }
  o <- do.call(order, c(ord_args, list(method = "radix")))
  df[o, , drop = FALSE]
}

cnt_distinct <- function(x) {
  x <- x[!is.na(x) & x != ""]
  length(unique(x))
}
cnt_nonmiss <- function(x) sum(!is.na(x) & x != "")

fmt_best <- function(x) ifelse(is.na(x), ".", format(x, scientific = FALSE, trim = TRUE))
fmt_81 <- function(x) {
  ifelse(is.na(x), ".",
         sprintf("%.1f", sign(x) * floor(abs(x) * 10 + 0.5 + 1e-9) / 10))
}
make_result <- function(count, denom, events) {
  r <- paste0(fmt_best(count), " (", fmt_81((count / denom) * 100), "%) [",
              fmt_best(events), "]")
  substr(r, 1, 50)
}
fmt_pvalue <- function(p) {
  ifelse(is.na(p), ".",
         ifelse(p < 0.0001, "<.0001", sprintf("%.4f", p)))
}

# ---- descriptive ----
ADSL_S <- ADSL %>% filter(SAFFL %in% "Y")
ADAE_S <- ADAE %>% filter(TRTEMFL %in% "Y" & AESER %in% "Y")

DENOM <- ADSL_S %>% mutate(TRTAN = TRT01AN) %>%
  group_by(TRTAN) %>% summarise(DENOM = cnt_distinct(USUBJID), .groups = "drop") %>%
  select(DENOM, TRTAN) %>% as.data.frame()
DENOM <- sas_arrange(DENOM, "TRTAN")

AEPTD <- ADAE_S %>% group_by(TRTAN, AEBODSYS, AEDECOD) %>%
  summarise(COUNT = cnt_distinct(USUBJID), .groups = "drop") %>%
  select(COUNT, AEBODSYS, AEDECOD, TRTAN) %>% as.data.frame()
AEPTD <- sas_arrange(AEPTD, c("TRTAN", "AEBODSYS", "AEDECOD"))

AEPTA <- ADAE_S %>% group_by(TRTAN, AEBODSYS, AEDECOD) %>%
  summarise(EVENTS = cnt_nonmiss(USUBJID), .groups = "drop") %>%
  select(EVENTS, AEBODSYS, AEDECOD, TRTAN) %>% as.data.frame()
AEPTA <- sas_arrange(AEPTA, c("TRTAN", "AEBODSYS", "AEDECOD"))

RESULTPT <- AEPTD %>% mutate(SUMMARYLEVEL = 3) %>%
  left_join(DENOM, by = "TRTAN") %>%
  left_join(AEPTA, by = c("TRTAN", "AEBODSYS", "AEDECOD")) %>%
  mutate(RESULT = make_result(COUNT, DENOM, EVENTS)) %>%
  select(COUNT, AEBODSYS, AEDECOD, TRTAN, SUMMARYLEVEL, DENOM, EVENTS, RESULT) %>%
  as.data.frame()
RESULTPT <- sas_arrange(RESULTPT, c("AEBODSYS", "AEDECOD", "TRTAN"))

AESOCD <- ADAE_S %>% group_by(TRTAN, AEBODSYS) %>%
  summarise(COUNT = cnt_distinct(USUBJID), .groups = "drop") %>%
  select(COUNT, AEBODSYS, TRTAN) %>% as.data.frame()
AESOCD <- sas_arrange(AESOCD, c("TRTAN", "AEBODSYS"))

AESOCA <- ADAE_S %>% group_by(TRTAN, AEBODSYS) %>%
  summarise(EVENTS = cnt_nonmiss(USUBJID), .groups = "drop") %>%
  select(EVENTS, AEBODSYS, TRTAN) %>% as.data.frame()
AESOCA <- sas_arrange(AESOCA, c("TRTAN", "AEBODSYS"))

RESULTSOC <- AESOCD %>% mutate(SUMMARYLEVEL = 2) %>%
  left_join(DENOM, by = "TRTAN") %>%
  left_join(AESOCA, by = c("TRTAN", "AEBODSYS")) %>%
  mutate(RESULT = make_result(COUNT, DENOM, EVENTS)) %>%
  select(COUNT, AEBODSYS, TRTAN, SUMMARYLEVEL, DENOM, EVENTS, RESULT) %>%
  as.data.frame()
RESULTSOC <- sas_arrange(RESULTSOC, c("AEBODSYS", "TRTAN"))

AESUBD <- ADAE_S %>% group_by(TRTAN) %>%
  summarise(COUNT = cnt_distinct(USUBJID), .groups = "drop") %>%
  select(COUNT, TRTAN) %>% as.data.frame()
AESUBD <- sas_arrange(AESUBD, "TRTAN")

AESUBA <- ADAE_S %>% group_by(TRTAN) %>%
  summarise(EVENTS = cnt_nonmiss(USUBJID), .groups = "drop") %>%
  select(EVENTS, TRTAN) %>% as.data.frame()
AESUBA <- sas_arrange(AESUBA, "TRTAN")

RESULTSUB <- AESUBD %>% mutate(SUMMARYLEVEL = 1) %>%
  left_join(DENOM, by = "TRTAN") %>%
  left_join(AESUBA, by = "TRTAN") %>%
  mutate(RESULT = make_result(COUNT, DENOM, EVENTS)) %>%
  select(COUNT, TRTAN, SUMMARYLEVEL, DENOM, EVENTS, RESULT) %>%
  as.data.frame()
RESULTSUB <- sas_arrange(RESULTSUB, "TRTAN")

# ---- inferential ----
subj <- ADSL_S %>% transmute(USUBJID, TRTAN = TRT01AN) %>% distinct() %>% as.data.frame()
subj <- sas_arrange(subj, c("USUBJID", "TRTAN"))

cross_join_ordered <- function(a, b) {
  ia <- rep(seq_len(nrow(a)), each = nrow(b))
  ib <- rep(seq_len(nrow(b)), times = nrow(a))
  out <- cbind(a[ia, , drop = FALSE], b[ib, , drop = FALSE])
  rownames(out) <- NULL
  out
}

run_fisher <- function(df, byvar, trts) {
  d <- df %>% filter(TRTAN %in% trts)
  do_one <- function(sub) {
    tab <- table(sub$TRTAN, sub$EVENT)
    if (nrow(tab) < 2 || ncol(tab) < 2) return(NULL)
    p <- fisher.test(tab)$p.value
    data.frame(Table = "Table TRTAN * EVENT", Name1 = "XP2_FISH",
               Label1 = "Two-sided Pr <= P", cValue1 = fmt_pvalue(p),
               nValue1 = p, stringsAsFactors = FALSE)
  }
  if (is.null(byvar)) {
    res <- do_one(d)
    if (is.null(res)) res <- data.frame(Table = character(), Name1 = character(),
                                        Label1 = character(), cValue1 = character(),
                                        nValue1 = numeric(), stringsAsFactors = FALSE)
    return(res)
  }
  lv <- unique(d[[byvar]])
  out <- list()
  for (v in lv) {
    sub <- d[d[[byvar]] == v, , drop = FALSE]
    r <- do_one(sub)
    if (!is.null(r)) {
      r <- cbind(setNames(data.frame(v, stringsAsFactors = FALSE), byvar), r)
      out[[length(out) + 1]] <- r
    }
  }
  if (length(out) == 0) {
    res <- data.frame(x = character(), Table = character(), Name1 = character(),
                      Label1 = character(), cValue1 = character(),
                      nValue1 = numeric(), stringsAsFactors = FALSE)
    names(res)[1] <- byvar
    return(res)
  }
  do.call(rbind, out)
}

# PT level
pts <- ADAE_S %>% distinct(AEDECOD) %>% as.data.frame()
pts <- sas_arrange(pts, "AEDECOD")
BACKBONE <- cross_join_ordered(subj, pts)

HADEVENT <- ADAE_S %>% distinct(STUDYID, USUBJID, TRTAN, AEDECOD) %>% as.data.frame()
HADEVENT <- sas_arrange(HADEVENT, c("STUDYID", "USUBJID", "TRTAN", "AEDECOD"))

FISHER <- BACKBONE %>%
  left_join(HADEVENT, by = c("USUBJID", "TRTAN", "AEDECOD")) %>%
  mutate(EVENT = ifelse(!is.na(STUDYID) & STUDYID == "CDISCPILOT01", 1, 0)) %>%
  select(USUBJID, TRTAN, AEDECOD, EVENT) %>% as.data.frame()
FISHER <- sas_arrange(FISHER, c("AEDECOD", "TRTAN"))

F1PT <- run_fisher(FISHER, "AEDECOD", c(0, 2))
F2PT <- run_fisher(FISHER, "AEDECOD", c(0, 1))

# SOC level
socs <- ADAE_S %>% distinct(AEBODSYS) %>% as.data.frame()
socs <- sas_arrange(socs, "AEBODSYS")
BACKBONE1 <- cross_join_ordered(subj, socs)

HADEVENT1 <- ADAE_S %>% distinct(STUDYID, USUBJID, TRTAN, AEBODSYS) %>% as.data.frame()
HADEVENT1 <- sas_arrange(HADEVENT1, c("STUDYID", "USUBJID", "TRTAN", "AEBODSYS"))

FISHER1 <- BACKBONE1 %>%
  left_join(HADEVENT1, by = c("USUBJID", "TRTAN", "AEBODSYS")) %>%
  mutate(EVENT = ifelse(!is.na(STUDYID) & STUDYID == "CDISCPILOT01", 1, 0)) %>%
  select(USUBJID, TRTAN, AEBODSYS, EVENT) %>% as.data.frame()
FISHER1 <- sas_arrange(FISHER1, c("AEBODSYS", "TRTAN"))

F1SOC <- run_fisher(FISHER1, "AEBODSYS", c(0, 2))
F2SOC <- run_fisher(FISHER1, "AEBODSYS", c(0, 1))

# Subject level
BACKBONE2 <- subj

HADEVENT2 <- ADAE_S %>% distinct(STUDYID, USUBJID, TRTAN) %>% as.data.frame()
HADEVENT2 <- sas_arrange(HADEVENT2, c("STUDYID", "USUBJID", "TRTAN"))

FISHER2 <- BACKBONE2 %>%
  left_join(HADEVENT2, by = c("USUBJID", "TRTAN")) %>%
  mutate(EVENT = ifelse(!is.na(STUDYID) & STUDYID == "CDISCPILOT01", 1, 0)) %>%
  select(USUBJID, TRTAN, EVENT) %>% as.data.frame()
FISHER2 <- sas_arrange(FISHER2, "TRTAN")

F1SUB <- run_fisher(FISHER2, NULL, c(0, 2))
F2SUB <- run_fisher(FISHER2, NULL, c(0, 1))

saveRDS(list(
  ADSL = ADSL, ADAE = ADAE,
  DENOM = DENOM, AEPTD = AEPTD, AEPTA = AEPTA, RESULTPT = RESULTPT,
  AESOCD = AESOCD, AESOCA = AESOCA, RESULTSOC = RESULTSOC,
  AESUBD = AESUBD, AESUBA = AESUBA, RESULTSUB = RESULTSUB,
  BACKBONE = BACKBONE, HADEVENT = HADEVENT, FISHER = FISHER,
  F1PT = F1PT, F2PT = F2PT,
  BACKBONE1 = BACKBONE1, HADEVENT1 = HADEVENT1, FISHER1 = FISHER1,
  F1SOC = F1SOC, F2SOC = F2SOC,
  BACKBONE2 = BACKBONE2, HADEVENT2 = HADEVENT2, FISHER2 = FISHER2,
  F1SUB = F1SUB, F2SUB = F2SUB
), "out.rds")
