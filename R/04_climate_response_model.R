# ==============================================================================
# Script: 04_climate_response_model.R
#
# Purpose:
#   Manually reconstruct the provisional pool-specific climate-response trial.
#   The candidate combines a resident multinomial with a nonresident binomial;
#   three species-random binomials are also fitted as sensitivity diagnostics.
#   Model outputs are reconstructable and intentionally absent from the clean
#   checkpoint. This script is not called by _targets.R. Only the probability
#   prototype is retained in the curated figure collection.
#
# Input:
#   output/tables/analysis_dataset.csv
#   output/tables/climate_context_experiments.csv
#
# Reconstructed outputs:
#   output/models/ and output/tables/climate_response_*
#   plot/climate_response_data_model_probability_prototype.pdf
#   plot/climate_response_binary_sensitivity_odds_ratios.pdf (not curated)
# ==============================================================================

source("config/analysis.R")
source("R/pipeline/utils.R")
source("R/pipeline/figures.R")
source("R/pipeline/climate_response_model.R")

config <- transplant_config()
analysis_path <- file.path(config$output_dir, "tables", "analysis_dataset.csv")
context_path <- file.path(
  config$output_dir, "tables", "climate_context_experiments.csv"
)
if (!file.exists(analysis_path)) {
  stop("Missing pool-model input: ", analysis_path, call. = FALSE)
}
if (!file.exists(context_path)) {
  stop("Missing climate-context input: ", context_path, call. = FALSE)
}

analysis_data <- readr::read_csv(analysis_path, show_col_types = FALSE)
climate_contexts <- readr::read_csv(context_path, show_col_types = FALSE)
prepared <- prepare_pool_climate_response_data(analysis_data)
models <- fit_pool_climate_response_models(prepared)
coefficients <- tidy_pool_climate_coefficients(models)
likelihood_ratio_tests <- pool_climate_likelihood_ratio_tests(models)
performance <- pool_climate_model_performance(models)
leave_one_experiment_out <- pool_climate_leave_one_experiment_out(models)
sensitivity <- fit_pool_binary_sensitivity(prepared)
sensitivity_coefficients <- tidy_pool_binary_sensitivity(sensitivity)
predictions <- pool_climate_prediction_surfaces(
  models, climate_contexts,
  grid_size = 60L, n_draws = 2000L, seed = 47163L
)

write_rds_atomic(
  models,
  file.path(config$output_dir, "models", "climate_response_pool_models.rds")
)
write_rds_atomic(
  sensitivity,
  file.path(
    config$output_dir, "models", "climate_response_binary_sensitivity.rds"
  )
)
write_csv_atomic(
  prepared$scaling,
  file.path(config$output_dir, "tables", "climate_response_pool_scaling.csv")
)
write_csv_atomic(
  prepared$category_counts,
  file.path(
    config$output_dir, "tables", "climate_response_pool_category_counts.csv"
  )
)
write_csv_atomic(
  prepared$pool_by_species,
  file.path(
    config$output_dir, "tables", "climate_response_pool_membership_by_species.csv"
  )
)
write_csv_atomic(
  models$diagnostics,
  file.path(
    config$output_dir, "tables", "climate_response_pool_diagnostics.csv"
  )
)
write_csv_atomic(
  coefficients,
  file.path(
    config$output_dir, "tables", "climate_response_pool_coefficients.csv"
  )
)
write_csv_atomic(
  likelihood_ratio_tests,
  file.path(
    config$output_dir, "tables", "climate_response_pool_likelihood_ratio_tests.csv"
  )
)
write_csv_atomic(
  performance$metrics,
  file.path(
    config$output_dir, "tables", "climate_response_pool_performance.csv"
  )
)
write_csv_atomic(
  performance$calibration,
  file.path(
    config$output_dir, "tables", "climate_response_pool_calibration.csv"
  )
)
write_csv_atomic(
  leave_one_experiment_out,
  file.path(
    config$output_dir, "tables",
    "climate_response_pool_leave_one_experiment_out.csv"
  )
)
write_csv_atomic(
  predictions$warming_levels,
  file.path(
    config$output_dir, "tables", "climate_response_pool_warming_levels.csv"
  )
)
write_csv_atomic(
  predictions$grid,
  file.path(
    config$output_dir, "tables", "climate_response_pool_prediction_grid.csv"
  )
)
write_csv_atomic(
  predictions$simulation,
  file.path(
    config$output_dir, "tables", "climate_response_pool_prediction_simulation.csv"
  )
)
write_csv_atomic(
  sensitivity$diagnostics,
  file.path(
    config$output_dir, "tables",
    "climate_response_binary_sensitivity_diagnostics.csv"
  )
)
write_csv_atomic(
  sensitivity$replication,
  file.path(
    config$output_dir, "tables",
    "climate_response_binary_sensitivity_replication.csv"
  )
)
write_csv_atomic(
  sensitivity_coefficients,
  file.path(
    config$output_dir, "tables",
    "climate_response_binary_sensitivity_coefficients.csv"
  )
)

plot_pool_binary_sensitivity_odds_ratios(
  sensitivity_coefficients,
  file.path(
    config$plot_dir, "climate_response_binary_sensitivity_odds_ratios.pdf"
  )
)
plot_pool_climate_response_review_prototype(
  analysis_data, climate_contexts, coefficients, predictions,
  file.path(
    config$plot_dir, "climate_response_data_model_probability_prototype.pdf"
  )
)
write_csv_atomic(
  figure_manifest(config$root),
  file.path(config$output_dir, "manifests", "figures.csv")
)

cat("\nPool-specific climate-response analysis completed.\n")
print(prepared$snapshot)
print(models$diagnostics)
print(coefficients)
print(likelihood_ratio_tests)
print(sensitivity$diagnostics)
print(performance$metrics)
