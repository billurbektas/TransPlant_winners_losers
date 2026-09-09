# ==============================================================================
# Script: 06_biological_context_models_review.R
#
# Purpose:
#   Manually reconstruct the review-only family, functional-group, and trait
#   context model trials. Figure rendering has been removed. Generated models,
#   tables, diagnostics, and prediction grids are intentionally absent from the
#   clean checkpoint and are not produced by _targets.R.
#
# Input:
#   output/tables/analysis_dataset.csv
#
# Reconstructed output:
#   output/review/biological_context_models/
# ==============================================================================

source("config/analysis.R")
source("R/pipeline/utils.R")
source("R/pipeline/figures.R")
source("R/pipeline/climate_response_model.R")
source("R/pipeline/biological_context_review.R")

config <- transplant_config()
input_path <- file.path(config$output_dir, "tables", "analysis_dataset.csv")
review_dir <- file.path(
  config$output_dir, "review", "biological_context_models"
)

if (!file.exists(input_path)) {
  stop("Missing biological-context review input: ", input_path, call. = FALSE)
}

dir.create(review_dir, recursive = TRUE, showWarnings = FALSE)

cat("Preparing review data and denominator-selected families...\n")
analysis_data <- readr::read_csv(input_path, show_col_types = FALSE)
prepared <- prepare_biological_context_review(
  analysis_data, family_minimum = 50L
)

cat("Fitting mixed multinomial trials and prescribed binomial fallbacks...\n")
fits <- fit_biological_context_review(
  prepared, mixed_timeout = 35, binary_timeout = 90
)

cat("Collecting coefficients, nested comparisons, and diagnostics...\n")
coefficients <- collect_biological_context_coefficients(fits)
likelihood_ratios <- collect_biological_context_likelihood_ratios(fits)
diagnostics <- collect_biological_context_diagnostics(fits)

cat("Generating population-level review predictions (2,000 draws)...\n")
predictions <- review_all_predictions(
  fits, prepared, n_draws = 2000L, seed = 47163L, grid_size = 60L
)
validation <- validate_biological_context_review(
  prepared, fits, predictions
)

trait_complete_cases <- dplyr::bind_rows(
  purrr::map_dfr(names(fits$traits), function(trait_name) {
    trait_fit <- fits$traits[[trait_name]]
    tibble::tibble(
      trait = trait_name,
      pool = c("resident", "nonresident"),
      n_species_experiment = c(
        nrow(trait_fit$resident_data),
        nrow(trait_fit$nonresident_data)
      ),
      n_species = c(
        dplyr::n_distinct(trait_fit$resident_data$accepted_species),
        dplyr::n_distinct(trait_fit$nonresident_data$accepted_species)
      ),
      n_experiments = c(
        dplyr::n_distinct(trait_fit$resident_data$experiment_id),
        dplyr::n_distinct(trait_fit$nonresident_data$experiment_id)
      )
    )
  })
)

cat("Writing review-only model objects and tables...\n")
write_rds_atomic(fits, file.path(review_dir, "model_objects.rds"))
write_rds_atomic(prepared, file.path(review_dir, "prepared_data.rds"))
write_rds_atomic(predictions, file.path(review_dir, "prediction_objects.rds"))
write_csv_atomic(prepared$snapshot, file.path(review_dir, "data_snapshot.csv"))
write_csv_atomic(
  prepared$family_selection, file.path(review_dir, "family_selection.csv")
)
write_csv_atomic(
  prepared$family_event_support,
  file.path(review_dir, "family_outcome_support.csv")
)
write_csv_atomic(
  prepared$functional_group_counts,
  file.path(review_dir, "functional_group_descriptive_counts.csv")
)
write_csv_atomic(
  prepared$trait_scaling, file.path(review_dir, "trait_scaling.csv")
)
write_csv_atomic(
  trait_complete_cases, file.path(review_dir, "trait_complete_cases.csv")
)
write_csv_atomic(
  diagnostics, file.path(review_dir, "model_diagnostics.csv")
)
write_csv_atomic(
  dplyr::filter(diagnostics, !.data$clean),
  file.path(review_dir, "failed_model_diagnostics.csv")
)
write_csv_atomic(
  coefficients, file.path(review_dir, "coefficients_odds_ratios.csv")
)
write_csv_atomic(
  likelihood_ratios, file.path(review_dir, "likelihood_ratio_tests.csv")
)
write_csv_atomic(
  predictions$family, file.path(review_dir, "family_context_changes.csv")
)
write_csv_atomic(
  predictions$functional_group,
  file.path(review_dir, "functional_group_context_probabilities.csv")
)
write_csv_atomic(
  predictions$traits, file.path(review_dir, "trait_context_probabilities.csv")
)
write_csv_atomic(validation, file.path(review_dir, "validation.csv"))

cat("\nBiological-context review completed.\n")
cat("Outputs: ", review_dir, "\n", sep = "")
