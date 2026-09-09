# ==============================================================================
# Script: 09_climate_significant_bars_review.R
#
# Purpose:
#   Manually reconstruct review-only summaries of significant loser, winner,
#   and colonizer species x experiment responses across destination-climate
#   groups. The generated review table is intentionally absent from the clean
#   checkpoint and is not produced by _targets.R.
#
# Percentage denominators:
#   - Loser and winner: all resident species x experiment records.
#   - Colonizer: all nonresident species x experiment records, retaining the
#     reversible colonizer/outsider calculation pools.
#
# Reconstructed output:
#   output/review/climate_significant_counts_pool_percentages.csv
# ============================================================================== 

source("config/analysis.R")
source("R/pipeline/utils.R")
source("R/pipeline/figures.R")

config <- transplant_config()
analysis_path <- file.path(config$output_dir, "tables", "analysis_dataset.csv")
review_output_dir <- file.path(config$output_dir, "review")

if (!file.exists(analysis_path)) {
  stop("The saved analysis dataset is missing.", call. = FALSE)
}
dir.create(review_output_dir, recursive = TRUE, showWarnings = FALSE)

climate_levels <- c(
  "Cold/Dry", "Intermediate/Dry", "Warm/Dry",
  "Cold/Intermediate", "Cold/Wet"
)
outcome_levels <- c("Loser", "Winner", "Colonizer")

analysis_data <- readr::read_csv(analysis_path, show_col_types = FALSE) |>
  dplyr::mutate(
    category = as.character(.data$analysis_category),
    species = dplyr::coalesce(.data$accepted_name, .data$species),
    experiment_id = as.character(.data$experiment_id),
    family = dplyr::if_else(
      is.na(.data$family) | trimws(.data$family) == "",
      "Unresolved", trimws(.data$family)
    )
  )

if (anyDuplicated(dplyr::select(analysis_data, "species", "experiment_id"))) {
  stop("The analysis dataset must contain one row per species x experiment.", call. = FALSE)
}

# Reproduce the reviewed equal-width display assignments -----------------------
experiment_climate <- analysis_data |>
  dplyr::distinct(
    .data$experiment_id, .data$destination_T_JJA, .data$destination_P_JJA
  ) |>
  dplyr::filter(
    is.finite(.data$destination_T_JJA),
    is.finite(.data$destination_P_JJA)
  )
if (nrow(experiment_climate) != 42L) {
  stop("Expected 42 distinct experiment climate profiles.", call. = FALSE)
}

temperature_edges <- seq(
  min(experiment_climate$destination_T_JJA),
  max(experiment_climate$destination_T_JJA),
  length.out = 4L
)
precipitation_edges <- seq(
  min(experiment_climate$destination_P_JJA),
  max(experiment_climate$destination_P_JJA),
  length.out = 4L
)

experiment_climate <- experiment_climate |>
  dplyr::mutate(
    temperature_context = cut(
      .data$destination_T_JJA,
      breaks = .env$temperature_edges,
      labels = c("cold", "intermediate", "warm"),
      include.lowest = TRUE
    ),
    precipitation_context = cut(
      .data$destination_P_JJA,
      breaks = .env$precipitation_edges,
      labels = c("dry", "intermediate", "wet"),
      include.lowest = TRUE
    ),
    grid_context = paste0(
      as.character(.data$temperature_context), "_",
      as.character(.data$precipitation_context)
    ),
    display_context = destination_climate_display_context(
      .data$experiment_id, .data$grid_context
    ),
    climate_context = factor(
      dplyr::recode(
        .data$display_context,
        cold_dry = "Cold/Dry",
        intermediate_dry = "Intermediate/Dry",
        warm_dry = "Warm/Dry",
        cold_intermediate = "Cold/Intermediate",
        cold_wet = "Cold/Wet"
      ),
      levels = climate_levels
    )
  )

expected_context_counts <- c(
  `Cold/Dry` = 14L,
  `Intermediate/Dry` = 9L,
  `Warm/Dry` = 11L,
  `Cold/Intermediate` = 3L,
  `Cold/Wet` = 5L
)
observed_context_counts <- table(experiment_climate$climate_context)
if (anyNA(experiment_climate$climate_context) ||
    !identical(
      as.integer(observed_context_counts[names(expected_context_counts)]),
      as.integer(expected_context_counts)
    )) {
  stop("The reviewed destination-climate assignments changed.", call. = FALSE)
}

analysis_data <- analysis_data |>
  dplyr::left_join(
    dplyr::select(experiment_climate, "experiment_id", "climate_context"),
    by = "experiment_id"
  )
if (anyNA(analysis_data$climate_context)) {
  stop("Every analysis row needs one reviewed climate assignment.", call. = FALSE)
}

# Event counts and pool-specific denominators ---------------------------------
make_event_rows <- function(data) {
  dplyr::bind_rows(
    data |>
      dplyr::filter(.data$pool == "resident") |>
      dplyr::mutate(outcome = "Loser", event = .data$category == "loser"),
    data |>
      dplyr::filter(.data$pool == "resident") |>
      dplyr::mutate(outcome = "Winner", event = .data$category == "winner"),
    data |>
      dplyr::filter(.data$pool %in% c("colonizer", "outsider")) |>
      dplyr::mutate(
        outcome = "Colonizer", event = .data$category == "colonizer"
      )
  ) |>
    dplyr::mutate(outcome = factor(.data$outcome, levels = outcome_levels))
}

summarise_entity <- function(data, entity_column, entity_type, entity_order) {
  event_rows <- make_event_rows(data) |>
    dplyr::mutate(entity = as.character(.data[[entity_column]]))
  observed <- event_rows |>
    dplyr::summarise(
      eligible_pool_count = dplyr::n(),
      significant_count = sum(.data$event),
      .by = c("entity", "climate_context", "outcome")
    )
  tidyr::crossing(
    entity = entity_order,
    climate_context = factor(climate_levels, levels = climate_levels),
    outcome = factor(outcome_levels, levels = outcome_levels)
  ) |>
    dplyr::left_join(
      observed, by = c("entity", "climate_context", "outcome")
    ) |>
    dplyr::mutate(
      entity_type = entity_type,
      entity = factor(.data$entity, levels = entity_order),
      significant_count = tidyr::replace_na(.data$significant_count, 0L),
      eligible_pool_count = tidyr::replace_na(.data$eligible_pool_count, 0L),
      pool_rate = dplyr::if_else(
        .data$eligible_pool_count > 0L,
        .data$significant_count / .data$eligible_pool_count,
        NA_real_
      ),
      eligible_pool = dplyr::if_else(
        .data$outcome %in% c("Loser", "Winner"),
        "resident", "nonresident"
      )
    )
}

species_data <- analysis_data |>
  dplyr::mutate(entity = "All species")
species_summary <- summarise_entity(
  species_data, "entity", "all_species", "All species"
)

family_order <- analysis_data |>
  dplyr::filter(.data$category %in% c("loser", "winner", "colonizer")) |>
  dplyr::count(.data$family, name = "significant_total") |>
  dplyr::arrange(dplyr::desc(.data$significant_total), .data$family) |>
  dplyr::pull(.data$family)
if (length(family_order) != 28L) {
  stop("Expected 28 families with at least one significant response.", call. = FALSE)
}
family_summary <- analysis_data |>
  dplyr::filter(.data$family %in% family_order) |>
  summarise_entity("family", "family", family_order)

functional_order <- c("Forb", "Graminoid")
functional_summary <- analysis_data |>
  dplyr::filter(.data$functional_group %in% c("forb", "graminoid")) |>
  dplyr::mutate(
    functional_group = stringr::str_to_title(.data$functional_group)
  ) |>
  summarise_entity("functional_group", "functional_group", functional_order)

all_summaries <- dplyr::bind_rows(
  species_summary, family_summary, functional_summary
) |>
  dplyr::mutate(
    climate_context = as.character(.data$climate_context),
    outcome = as.character(.data$outcome),
    entity = as.character(.data$entity),
    pool_percentage = 100 * .data$pool_rate
  )
write_csv_atomic(
  all_summaries,
  file.path(
    review_output_dir, "climate_significant_counts_pool_percentages.csv"
  )
)

# Validation -------------------------------------------------------------------
expected_events <- c(Loser = 132L, Winner = 69L, Colonizer = 34L)
observed_events <- species_summary |>
  dplyr::summarise(
    events = sum(.data$significant_count), .by = "outcome"
  )
if (!identical(
  as.integer(observed_events$events[match(names(expected_events), observed_events$outcome)]),
  as.integer(expected_events)
)) {
  stop("The significant-response totals changed.", call. = FALSE)
}
family_events <- family_summary |>
  dplyr::summarise(events = sum(.data$significant_count), .by = "outcome")
if (!identical(
  as.integer(family_events$events[match(names(expected_events), family_events$outcome)]),
  as.integer(expected_events)
)) {
  stop("Family panels do not conserve significant-response counts.", call. = FALSE)
}
if (any(all_summaries$significant_count > all_summaries$eligible_pool_count) ||
    any(all_summaries$pool_rate < 0 | all_summaries$pool_rate > 1, na.rm = TRUE)) {
  stop("A pool-specific response rate failed validation.", call. = FALSE)
}

cat("Climate significant-response summaries completed.\n")
cat("Percentages use resident denominators for loser/winner and nonresident denominators for colonizer.\n")
cat("Climate counts: Cold/Dry 14; Intermediate/Dry 9; Warm/Dry 11; Cold/Intermediate 3; Cold/Wet 5.\n")
cat("No maintained target or figure manifest was changed.\n")
