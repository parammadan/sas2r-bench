# adae.sas -> exact reference: the released ADaM adae.xpt.
# Rows are matched on USUBJID + AESEQ, then every released column is compared value by value.
# Labels, formats and column order are not scored here.

compare_adae <- function(out, data_dir) {
  ref <- as.data.frame(haven::read_xpt(file.path(data_dir, "adae.xpt")))
  got <- pick(out, "ADAE")
  if (is.null(got)) {
    return(list(reference = "exact", pass = FALSE, reason = "no ADAE dataset in out.rds"))
  }

  date_cols <- names(ref)[vapply(ref, inherits, logical(1), "Date")]
  missing_cols <- setdiff(names(ref), names(got))
  common <- intersect(names(ref), names(got))

  if (!all(c("USUBJID", "AESEQ") %in% names(got))) {
    return(list(reference = "exact", pass = FALSE, reason = "ADAE has no USUBJID/AESEQ keys",
                columns_missing = missing_cols))
  }
  key <- function(d) paste(trimws(d$USUBJID), as.numeric(d$AESEQ))
  ref_k <- key(ref)
  got_k <- key(got)
  joined <- match(ref_k, got_k)

  per_col <- sapply(common, function(v) {
    a <- if (v %in% date_cols) as_iso_date(ref[[v]]) else norm(ref[[v]])
    b <- got[[v]][joined]
    b <- if (v %in% date_cols) as_iso_date(b) else norm(b)
    mean((is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b))
  })

  cells_total <- nrow(ref) * ncol(ref)
  cells_ok <- sum(per_col * nrow(ref))
  list(
    reference = "exact",
    rows_ref = nrow(ref),
    rows_got = nrow(got),
    rows_matched = sum(!is.na(joined)),
    columns_ref = ncol(ref),
    columns_missing = missing_cols,
    columns_extra = setdiff(names(got), names(ref)),
    columns_exact = names(per_col)[per_col == 1],
    columns_partial = as.list(round(per_col[per_col < 1], 4)),
    cell_match = round(cells_ok / cells_total, 4),
    pass = nrow(got) == nrow(ref) && length(missing_cols) == 0 && all(per_col == 1)
  )
}
