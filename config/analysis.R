# Module: config/analysis
# Responsibility: Declare input paths, fixed scientific settings, and resident-model exclusions.
# Inputs: Project root and the optional CHELSA runtime environment settings.
# Outputs: Analysis configuration and explicit excluded-species names.
# Side effects: None.

#' Paths and settings for the TransPlant analysis.
#'
#' @param root Project directory containing data/, config/, and output/.
#' @return A list of input paths, ecological thresholds, climate years, and
#'   runtime settings. `chelsa_offline` is TRUE only when the environment
#'   variable `TRANSPLANT_CHELSA_OFFLINE` is set to `true` or `1`.
transplant_config <- function(root = here::here()) {
  offline_setting <- tolower(Sys.getenv("TRANSPLANT_CHELSA_OFFLINE", "false"))
  if (!offline_setting %in% c("true", "false", "1", "0")) {
    stop("TRANSPLANT_CHELSA_OFFLINE must be true, false, 1, or 0.", call. = FALSE)
  }
  list(
    root = root,
    input_rdata = file.path(root, "data", "TransPlantNetwork_060824.RData"),
    try_file = file.path(root, "data", "try_47163.txt"),
    chelsa_cache = file.path(root, "data", "cache", "chelsa_monthly"),
    output_dir = file.path(root, "output"),
    try_raster_sha256 = "59b749fe6711e72c720ff903044eb0efb4373524ee434740ed7384a0284c155b",
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
    chelsa_offline = offline_setting %in% c("true", "1"),
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
#' @param excluded_species Species excluded at the resident-model boundary.
#' @return An R object described by the function name; its structure is validated before later stages use it.
filter_resident_model_species <- function(data, species_column, excluded_species) {
  if (length(species_column) != 1L || !species_column %in% names(data)) {
    stop("A valid species column is required for the resident-model exclusion.",
         call. = FALSE)
  }
  if (!is.character(excluded_species) || anyNA(excluded_species)) {
    stop("Resident-model exclusions must be an explicit character vector.",
         call. = FALSE)
  }
  keep <- !as.character(data[[species_column]]) %in% excluded_species
  result <- data[keep, , drop = FALSE]
  attr(result, "resident_model_exclusions") <- excluded_species
  result
}
