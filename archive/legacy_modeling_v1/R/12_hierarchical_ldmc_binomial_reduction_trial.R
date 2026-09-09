# ==============================================================================
# Script: 12_hierarchical_ldmc_binomial_reduction_trial.R
#
# Purpose:
#   Manually reconstruct progressively simpler species-hierarchical binomial
#   LDMC diagnostics after removing experiment and gradient random intercepts.
#   This is a review-only feasibility grid outside _targets.R; generated model
#   objects, tables, and figures are intentionally absent from the checkpoint.
#
# Input:
#   output/tables/analysis_dataset.csv
#
# Reconstructed output:
#   output/review/hierarchical_ldmc_binomial_reductions/
# ==============================================================================

source("config/analysis.R")
source("R/pipeline/utils.R")
source("R/pipeline/figures.R")
source("R/pipeline/climate_response_model.R")
source("R/pipeline/hierarchical_trait_trial.R")

config <- transplant_config()
input_path <- file.path(config$output_dir, "tables", "analysis_dataset.csv")
output_dir <- file.path(
  config$output_dir, "review", "hierarchical_ldmc_binomial_reductions"
)
if (!file.exists(input_path)) {
  stop("Missing binomial-reduction input: ", input_path, call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Preparing LDMC-complete resident and nonresident pools...\n")
analysis_data <- readr::read_csv(input_path, show_col_types = FALSE)
prepared <- prepare_hierarchical_ldmc_trial(analysis_data)
plan <- hierarchical_ldmc_reduction_plan()

cat("Fitting 9 species-only binomial specifications for each outcome...\n")
fits <- fit_hierarchical_ldmc_binomial_reductions(prepared)
coefficients <- tidy_hierarchical_ldmc_reduction_coefficients(fits)
random_effects <- collect_hierarchical_ldmc_reduction_random_effects(fits)
validation <- validate_hierarchical_ldmc_binomial_reductions(fits)
signed_species <- signed_species_consistency(analysis_data)
repeated_species_effects <-
  extract_hierarchical_ldmc_repeated_species_effects(
    fits$models, prepared, signed_species
  )
repeated_species_validation <-
  validate_hierarchical_ldmc_repeated_species_effects(
    repeated_species_effects, signed_species
  )

write_rds_atomic(fits$models, file.path(output_dir, "model_objects.rds"))
write_csv_atomic(plan, file.path(output_dir, "model_plan.csv"))
write_csv_atomic(
  fits$diagnostics, file.path(output_dir, "model_diagnostics.csv")
)
write_csv_atomic(
  coefficients, file.path(output_dir, "coefficients_odds_ratios.csv")
)
write_csv_atomic(
  random_effects, file.path(output_dir, "random_effect_variances.csv")
)
write_csv_atomic(validation, file.path(output_dir, "validation.csv"))
write_csv_atomic(
  repeated_species_effects,
  file.path(output_dir, "repeated_species_random_intercepts.csv")
)
write_csv_atomic(
  repeated_species_validation,
  file.path(output_dir, "repeated_species_random_intercepts_validation.csv")
)
plot_hierarchical_ldmc_resident_odds_ratios(
  coefficients, random_effects,
  file.path(output_dir, paste0("loser_winner_odds_ratios", ".pdf"))
)
plot_hierarchical_ldmc_repeated_species_effects(
  repeated_species_effects,
  file.path(
    output_dir, paste0("loser_winner_species_random_intercepts", ".pdf")
  )
)

if (!all(validation$passed) || !all(repeated_species_validation$passed)) {
  stop(
    "Binomial-reduction validation failed: ",
    paste(c(
      validation$check[!validation$passed],
      repeated_species_validation$check[!repeated_species_validation$passed]
    ), collapse = ", "),
    call. = FALSE
  )
}

cat("\nBinomial-reduction trial completed.\n")
print(
  fits$diagnostics |>
    dplyr::select(
      "outcome", "specification", "converged", "singular",
      "clean_hierarchical_fit", "minimum_random_variance", "warning"
    ),
  n = Inf
)
cat("Outputs: ", output_dir, "\n", sep = "")
