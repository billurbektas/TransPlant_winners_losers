# Narrow validation interface for experiment-level climate targets.

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
