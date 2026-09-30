# Module: config/analysis
# Responsibility: Implement the named scientific workflow stage or shared infrastructure.
# Inputs: Function arguments and explicitly declared target dependencies.
# Outputs: Validated in-memory objects and declared workflow artifacts.
# Side effects: Limited to functions whose contracts explicitly write, copy, or render files.

#' Transplant Config.
#'
#' @param root A project-relative or absolute filesystem path.
#' @return An R object described by the function name; its structure is validated before later stages use it.
transplant_config <- function(root = here::here()) {
  list(
    root = root,
    input_rdata = file.path(root, "data", "TransPlantNetwork_060824.RData"),
    try_file = file.path(root, "data", "try_47163.txt"),
    chelsa_cache = file.path(root, "data", "cache", "chelsa_monthly"),
    output_dir = file.path(root, "output"),
    plot_dir = file.path(root, "plot"), # reserved output location; currently unused
    taxonomy_min_score = 0.90,
    error_risk_max = 4,
    excluded_try_climates = c("Af", "Am", "Aw", "BWh", "BWk", "BSh", "BSk"),
    horizons = 2L:9L,
    climate_variables = c("tas", "pr"),
    climate_required_years = 2008:2023,
    climate_download_years = 2008:2021,
    chelsa_workers = as.integer(Sys.getenv("TRANSPLANT_CHELSA_WORKERS", "2")),
    chelsa_retries = as.integer(Sys.getenv("TRANSPLANT_CHELSA_RETRIES", "4")),
    chelsa_backoff_seconds = as.numeric(Sys.getenv("TRANSPLANT_CHELSA_BACKOFF", "2")),
    seed = 20260824L
  )
}

# These taxa remain in the community data and winner/loser contrast
# calculations. They are excluded at the resident-model boundary for every
# ecological model and phylogenetic analysis. Do not delete their calculated
# winner/loser records: the exclusion is an analysis-scope decision.
#' Resident Model Exclusions.
#'
#' @return An R object described by the function name; its structure is validated before later stages use it.
resident_model_exclusions <- function() {
  c(
    "Gymnocarpium dryopteris",
    "Athyrium spinulosum",
    "Botrychium lunaria",
    "Equisetum arvense",
    "Selaginella selaginoides",
    "Juniperus communis"
  )
}

#' Filter Resident Model Species.
#'
#' @param data A data frame whose required columns are checked before use.
#' @param species_column A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
filter_resident_model_species <- function(data, species_column) {
  if (length(species_column) != 1L || !species_column %in% names(data)) {
    stop("A valid species column is required for the resident-model exclusion.",
         call. = FALSE)
  }
  keep <- !as.character(data[[species_column]]) %in%
    resident_model_exclusions()
  result <- data[keep, , drop = FALSE]
  attr(result, "resident_model_exclusions") <-
    resident_model_exclusions()
  result
}
