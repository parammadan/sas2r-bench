# Structure check for programs with no published reference yet: did the translation run, and does
# out.rds hold every dataset the SAS code creates (refs/expected/<id>.json, read mechanically from
# the program by bench/expected.py)? Datasets named with macro variables can't be listed without
# SAS and are not required. Full scoring waits for real SAS output (SAS OnDemand).

compare_structure <- function(id, out) {
  exp <- jsonlite::fromJSON(file.path("refs", "expected", paste0(id, ".json")))
  got <- unique(toupper(sub("^.*\\.", "", names(out))))
  want <- toupper(exp$datasets)
  missing <- setdiff(want, got)
  list(reference = "structure", datasets_expected = length(want), datasets_present = length(want) - length(missing),
       datasets_missing = missing, datasets_extra = setdiff(got, want), not_required = exp$skipped,
       pass = length(missing) == 0)
}

compare_mean_time <- function(out, data_dir) compare_structure("mean_time", out)
compare_outliers_Shift_Table <- function(out, data_dir) compare_structure("outliers_Shift_Table", out)
compare_outliers_TEHiLo_Anytime <- function(out, data_dir) compare_structure("outliers_TEHiLo_Anytime", out)
compare_pk_param_summary <- function(out, data_dir) compare_structure("pk_param_summary", out)
compare_medication_listing <- function(out, data_dir) compare_structure("medication_listing", out)
