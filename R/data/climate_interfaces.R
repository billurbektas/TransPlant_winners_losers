# Module: R/data/climate_interfaces
# Responsibility: Check completeness and relationships in experiment climate outputs.
# Inputs: Site and experiment climate tables.
# Outputs: Validated climate interface.
# Side effects: None.

# Narrow validation interface for experiment-level climate targets.

#' Validate Experiment Climate Interface.
#'
#' @param climate A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
validate_experiment_climate_interface <- function(climate) {
  required <- c("experiment_id", "T", "P", "W", "T_z", "P_z", "W_z")
  missing <- setdiff(required, names(climate))
  if (length(missing)) stop("Experiment climate is missing: ", paste(missing, collapse = ", "))
  tibble::tibble(
    experiments = dplyr::n_distinct(climate$experiment_id),
    complete = all(stats::complete.cases(climate[required])),
    pass = !anyDuplicated(climate$experiment_id) && all(stats::complete.cases(climate[required]))
  )
}
