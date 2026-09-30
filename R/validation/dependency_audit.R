# Module: R/validation/dependency_audit
# Responsibility: Confirm that the active target graph is disconnected from
# historical code, fixtures, reports, and reference results.
# Inputs: The active `_targets.R` manifest.
# Outputs: One review row for target names and one for target commands.
# Side effects: Reads target definitions but never executes a target.

#' Audit disconnection from historical workflow material
#'
#' @param root Repository root containing `_targets.R`.
#' @return A two-row table checking active target names and commands.
audit_historical_dependency_disconnection <- function(root = NULL) {
  if (is.null(root)) {
    root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
    repeat {
      if (file.exists(file.path(root, "_targets.R"))) break
      parent <- dirname(root)
      if (identical(parent, root)) {
        stop("Could not locate _targets.R from the current directory.",
             call. = FALSE)
      }
      root <- parent
    }
  }
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  previous_directory <- setwd(root)
  on.exit(setwd(previous_directory), add = TRUE)
  manifest <- targets::tar_manifest(
    fields = c("name", "command"),
    callr_function = NULL,
    script = file.path(root, "_targets.R")
  )
  historical_name_pattern <- paste(
    "(^legacy_|equivalence|step(19|22)|transition_reference|",
    "frozen_release_reference|resident_consistency_checks)",
    sep = ""
  )
  historical_command_pattern <- paste(
    "(archive/|reports/reference_text|tests/fixtures|",
    "reference_artifacts[.]csv|legacy_|step(19|22)|transition_reference)",
    sep = ""
  )
  name_hits <- manifest$name[
    grepl(historical_name_pattern, manifest$name, ignore.case = TRUE)
  ]
  command_hits <- manifest$name[
    grepl(historical_command_pattern, manifest$command, ignore.case = TRUE)
  ]
  result <- tibble::tibble(
    check = c(
      "No historical targets in active manifest",
      "No historical paths or objects in active commands"
    ),
    pass = c(!length(name_hits), !length(command_hits)),
    detail = c(
      if (length(name_hits)) paste(name_hits, collapse = "; ") else
        paste(nrow(manifest), "active target names checked"),
      if (length(command_hits)) paste(command_hits, collapse = "; ") else
        paste(nrow(manifest), "active target commands checked")
    )
  )
  if (any(!result$pass)) {
    stop(
      "Historical material remains connected to the active target graph: ",
      paste(result$detail[!result$pass], collapse = " | "),
      call. = FALSE
    )
  }
  result
}
