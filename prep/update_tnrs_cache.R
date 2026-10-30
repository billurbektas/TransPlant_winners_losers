# P1: add new names to the TNRS cache in prep/tnrs/, outside the pipeline.
# Run from the repository root when tar_make stops on names missing from the
# cache (a changed SQLite, community ledger or TRY file):
#   Rscript prep/update_tnrs_cache.R
# TNRS answers can change between runs (P6 found a TRY name that did), so a
# name is resolved once and kept: names already in the cache are never sent
# again. Each new row carries its date and the tidyTRY and TNRS versions; the
# warnings of the call go to prep/tnrs/warnings.csv.
# The cache was seeded on 30 Oct 2026 from the pipeline store, with the
# 28 Oct resolutions that P1-P4 and step 01 were approved on.

library(targets)
source("prep/R/taxonomy.R")

script <- "prep/_targets.R"
store <- "prep/_targets"
warnings_file <- "prep/tnrs/warnings.csv"

add_to_cache <- function(names, cache_file) {
  cache <- read_tnrs_cache(cache_file)
  new <- setdiff(names, cache$submitted_name)
  if (!length(new)) {
    message(cache_file, ": cache complete (", nrow(cache), " names)")
    return(invisible(0L))
  }
  message(cache_file, ": resolving ", length(new), " new names with TNRS")
  res <- tnrs_resolve(new)
  stopifnot(setequal(res$submitted_name, new))
  today <- as.character(Sys.Date())
  res$resolved_date <- today
  res$tidyTRY <- as.character(utils::packageVersion("tidyTRY"))
  res$TNRS <- as.character(utils::packageVersion("TNRS"))
  readr::write_csv(dplyr::bind_rows(cache, res), cache_file, na = "")
  warns <- attr(res, "warnings")
  if (length(warns))
    readr::write_csv(tibble::tibble(cache_file = basename(cache_file), resolved_date = today, message = warns),
                     warnings_file, append = TRUE)
  invisible(length(new))
}

# study names first: the TRY candidates are drawn from the genera of their
# resolved names
tar_make(names = study_names, script = script, store = store)
add_to_cache(tar_read(study_names, store = store)$name, "prep/tnrs/name_resolution_study.csv")
tar_make(names = try_candidate_names, script = script, store = store)
add_to_cache(tar_read(try_candidate_names, store = store)$name, "prep/tnrs/name_resolution_try.csv.gz")
message("Now rerun the pipeline: Rscript -e 'targets::tar_make(script = \"prep/_targets.R\", store = \"prep/_targets\")'")
