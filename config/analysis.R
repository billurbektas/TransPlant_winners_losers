transplant_config <- function(root = here::here()) {
  list(
    root = root,
    input_rdata = file.path(root, "data", "TransPlantNetwork_060824.RData"),
    try_file = file.path(root, "data", "try_47163.txt"),
    chelsa_cache = file.path(root, "data", "cache", "chelsa_monthly"),
    output_dir = file.path(root, "output"),
    plot_dir = file.path(root, "plot"),
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

# These taxa remain in the upstream community data and winner/loser contrast
# calculations. They are excluded immediately upstream of every downstream
# ecological model and phylogenetic analysis. Do not delete their calculated
# winner/loser records: the exclusion is an analysis-scope decision.
resident_downstream_model_exclusions <- function() {
  c(
    "Gymnocarpium dryopteris",
    "Athyrium spinulosum",
    "Botrychium lunaria",
    "Equisetum arvense",
    "Selaginella selaginoides",
    "Juniperus communis"
  )
}

filter_downstream_model_species <- function(data, species_column) {
  if (length(species_column) != 1L || !species_column %in% names(data)) {
    stop("A valid species column is required for the downstream exclusion.",
         call. = FALSE)
  }
  keep <- !as.character(data[[species_column]]) %in%
    resident_downstream_model_exclusions()
  result <- data[keep, , drop = FALSE]
  attr(result, "downstream_model_exclusions") <-
    resident_downstream_model_exclusions()
  result
}
