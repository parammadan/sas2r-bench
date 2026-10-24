# Helpers shared by the per-program comparators. A comparator reads the translation's out.rds
# (a named list of every dataset the program has at the end) and returns a list that
# bench/compare.R writes to results/<translator>/<id>.json.

suppressPackageStartupMessages(library(dplyr))

# Find a dataset by SAS name, ignoring case and any "libref." prefix.
pick <- function(out, name) {
  if (is.null(out)) return(NULL)
  keys <- sub("^.*\\.", "", toupper(names(out)))
  i <- match(toupper(name), keys)
  if (is.na(i)) return(NULL)
  d <- out[[i]]
  names(d) <- toupper(names(d))
  as.data.frame(d)
}

# Put one column into a comparable form: dates as ISO text, numbers rounded, text trimmed,
# empty text as NA (SAS stores missing character values as blanks).
norm <- function(x) {
  if (inherits(x, "Date")) return(format(x, "%Y-%m-%d"))
  if (inherits(x, "POSIXt")) return(format(x, "%Y-%m-%d %H:%M:%S"))
  if (is.factor(x)) x <- as.character(x)
  if (is.numeric(x)) return(as.character(round(x, 8)))
  x <- trimws(as.character(x))
  x[x == ""] <- NA
  x
}

# SAS stores dates as days since 1960-01-01; accept either form for date columns.
as_iso_date <- function(x) {
  if (inherits(x, "Date")) return(format(x, "%Y-%m-%d"))
  if (is.numeric(x)) return(format(as.Date(x, origin = "1960-01-01"), "%Y-%m-%d"))
  norm(x)
}
