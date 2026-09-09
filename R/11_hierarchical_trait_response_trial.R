# ==============================================================================
# Script: 11_hierarchical_trait_response_trial.R
#
# Purpose:
#   Manually reconstruct the feasibility trial of the paper-inspired
#   hierarchical LDMC structure. This exploratory trial is intentionally
#   outside _targets.R. Generated diagnostics and model objects are deliberately
#   absent from the clean checkpoint.
#
# Reference model idea:
#   Martinez-Almoyna et al., "Climate, soil resources and microbial activity
#   shape the distributions of mountain plants based on their functional
#   traits" (Ecography, article ecog.05269).
#
# Input:
#   output/tables/analysis_dataset.csv
#
# Reconstructed output:
#   output/review/hierarchical_ldmc_trial/
# ==============================================================================

source("config/analysis.R")
source("R/pipeline/utils.R")
source("R/pipeline/climate_response_model.R")
source("R/pipeline/hierarchical_trait_trial.R")

config <- transplant_config()
input_path <- file.path(config$output_dir, "tables", "analysis_dataset.csv")
output_dir <- file.path(config$output_dir, "review", "hierarchical_ldmc_trial")

if (!file.exists(input_path)) {
  stop("Missing hierarchical-trial input: ", input_path, call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Preparing raw-cover, maximum-duration LDMC pool data...\n")
analysis_data <- readr::read_csv(input_path, show_col_types = FALSE)
prepared <- prepare_hierarchical_ldmc_trial(analysis_data)
print(prepared$snapshot)
print(prepared$counts)

cat("\nFitting resident PQL multinomial with six random components...\n")
cat("Category covariance is diagonal: no loser/winner random correlations.\n")
cat("If that hierarchy is not clean, fitting loser and winner binomial GLMMs.\n")
cat("The nonresident colonizer binomial GLMM is always fitted.\n")
fits <- fit_hierarchical_ldmc_trial(prepared)
event_support <- summarize_hierarchical_trial_event_support(prepared)
coefficients <- tidy_hierarchical_trial_coefficients(fits)
random_effects <- collect_hierarchical_trial_random_effects(fits)
validation <- validate_hierarchical_ldmc_trial(prepared, fits)

write_rds_atomic(
  list(
    resident_multinomial = fits$resident_multinomial$fit,
    binary = purrr::map(fits$binary, "fit")
  ),
  file.path(output_dir, "model_objects.rds")
)
write_rds_atomic(prepared, file.path(output_dir, "prepared_data.rds"))
write_csv_atomic(prepared$snapshot, file.path(output_dir, "data_snapshot.csv"))
write_csv_atomic(prepared$counts, file.path(output_dir, "outcome_counts.csv"))
write_csv_atomic(
  prepared$climate_scaling, file.path(output_dir, "climate_scaling.csv")
)
write_csv_atomic(
  prepared$trait_scaling, file.path(output_dir, "LDMC_scaling.csv")
)
write_csv_atomic(
  prepared$experiment_gradient_map,
  file.path(output_dir, "experiment_gradient_map.csv")
)
write_csv_atomic(
  prepared$replication, file.path(output_dir, "species_replication.csv")
)
write_csv_atomic(
  event_support, file.path(output_dir, "event_replication_support.csv")
)
write_csv_atomic(
  fits$diagnostics, file.path(output_dir, "model_diagnostics.csv")
)
write_csv_atomic(
  coefficients, file.path(output_dir, "coefficients_odds_ratios.csv")
)
write_csv_atomic(
  random_effects, file.path(output_dir, "random_effect_covariances.csv")
)
write_csv_atomic(validation, file.path(output_dir, "validation.csv"))

if (!all(validation$passed)) {
  stop(
    "Hierarchical LDMC trial validation failed: ",
    paste(validation$check[!validation$passed], collapse = ", "),
    call. = FALSE
  )
}

cat("\nHierarchical LDMC trial completed.\n")
print(fits$diagnostics)
cat("Outputs: ", output_dir, "\n", sep = "")
