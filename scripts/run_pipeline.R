#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
selected_targets <- if (length(args)) args else NULL
crew_setting <- tolower(trimws(Sys.getenv("TRANSPLANT_USE_CREW", "true")))
if (!crew_setting %in% c("true", "false")) {
  stop("TRANSPLANT_USE_CREW must be true or false.", call. = FALSE)
}
use_crew <- identical(crew_setting, "true")

if (is.null(selected_targets)) {
  targets::tar_make(use_crew = use_crew, callr_function = NULL)
} else {
  rlang::inject(
    targets::tar_make(
      names = tidyselect::all_of(!!selected_targets),
      use_crew = use_crew,
      callr_function = NULL
    )
  )
}
