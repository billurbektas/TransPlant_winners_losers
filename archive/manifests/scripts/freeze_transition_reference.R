#!/usr/bin/env Rscript

# One-time, reviewable migration utility for the canonical refactor.
#
# It condenses the validated Step 19/20/22 transition artifacts into a single
# regression fixture. The fixture contains only the numerical interfaces used
# by the canonical target graph; it does not contain fitted-model environments,
# caches, figures, reports, or raw source data.

args <- commandArgs(trailingOnly = TRUE)
if (!identical(args, "--write")) {
  stop(
    "Refusing to replace the transition fixture without the explicit --write flag.",
    call. = FALSE
  )
}

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(root, "R", "functions", "model_core.R"))

read_csv <- function(...) {
  readr::read_csv(..., show_col_types = FALSE, name_repair = "minimal")
}
choose_tree <- function(active, archived, marker) {
  if (file.exists(file.path(active, marker))) return(active)
  if (file.exists(file.path(archived, marker))) return(archived)
  active
}
downstream <- choose_tree(
  file.path(root, "output", "step_by_step_modeling", "downstream_526"),
  file.path(root, "archive", "canonical_source_snapshots", "output", "downstream_526"),
  file.path("tables", "auc_data_526.csv")
)
phylogeny <- choose_tree(
  file.path(root, "output", "step_by_step_modeling", "resident_phylogeny_pcoa"),
  file.path(root, "archive", "canonical_source_snapshots", "output", "resident_phylogeny_pcoa"),
  "resident_phylogeny.rds"
)
trait_audit <- choose_tree(
  file.path(root, "output", "step_by_step_modeling", "trait_data_audit"),
  file.path(root, "archive", "canonical_source_snapshots", "output", "trait_data_audit"),
  "trait_species_level_audit.csv"
)

required <- c(
  file.path(downstream, "tables", "auc_data_526.csv"),
  file.path(downstream, "models", "V_auc_526.rds"),
  file.path(downstream, "models", "M_both_526.rds"),
  file.path(downstream, "models", "M_family_aggregated_526.rds"),
  file.path(phylogeny, "resident_phylogeny.rds"),
  file.path(phylogeny, "resident_phylogeny.newick"),
  file.path(phylogeny, "species_pev_scores.csv")
)
missing <- required[!file.exists(required)]
if (length(missing)) {
  stop("Missing transition inputs: ", paste(missing, collapse = "; "),
       call. = FALSE)
}

fixture <- list(
  schema_version = 1L,
  authority = paste(
    "Validated Step 19/20/22 transition interfaces; the Step 22 numerical",
    "authority remains downstream_526_reference_results.rds"
  ),
  downstream_auc_data = read_csv(file.path(downstream, "tables", "auc_data_526.csv")),
  downstream_V_auc = readRDS(file.path(downstream, "models", "V_auc_526.rds")),
  downstream_consistency_checks = read_csv(
    file.path(downstream, "tables", "required_consistency_checks_526.csv")
  ),
  mboth_model_signature = canonical_model_signature(
    readRDS(file.path(downstream, "models", "M_both_526.rds"))
  ),
  mboth_bundle = list(
    fixed = read_csv(file.path(downstream, "tables", "M_both_526_fixed_coefficients.csv")),
    CR1 = read_csv(file.path(downstream, "tables", "M_both_526_fixed_effects_model_vs_CR1.csv")),
    random = read_csv(file.path(downstream, "tables", "M_both_526_random_diagnostics.csv")),
    variance_profiles = read_csv(file.path(downstream, "tables", "M_both_526_variance_profile_CIs.csv")),
    fit_statistics = read_csv(file.path(downstream, "tables", "M_both_526_fit_statistics.csv"))
  ),
  species_classification = read_csv(
    file.path(downstream, "tables", "species_minmax_classification_526.csv")
  ),
  family = list(
    classification = read_csv(file.path(downstream, "tables", "family_minmax_classification_526.csv")),
    model_summary = read_csv(file.path(downstream, "tables", "family_model_diagnostics_526.csv")),
    auc_data = read_csv(file.path(downstream, "tables", "family_auc_data_526.csv")),
    V = readRDS(file.path(downstream, "models", "V_family_auc_526.rds")),
    model_signature = canonical_model_signature(
      readRDS(file.path(downstream, "models", "M_family_aggregated_526.rds"))
    ),
    support = read_csv(file.path(downstream, "tables", "family_support_526.csv")),
    random = read_csv(file.path(downstream, "tables", "family_random_coefficients_526.csv")),
    cover_audit = read_csv(file.path(downstream, "tables", "family_plot_cover_reconstruction_audit_526.csv"))
  ),
  trait_scaling = read_csv(file.path(downstream, "tables", "trait_scaling_526.csv")),
  trait_species = read_csv(file.path(trait_audit, "trait_species_level_audit.csv")),
  climate_scaling = read_csv(file.path(downstream, "tables", "climate_scaling_526.csv")),
  pev_scaling = read_csv(file.path(downstream, "tables", "PEV1_PEV3_scaling_526.csv")),
  functional_species = read_csv(
    file.path(downstream, "tables", "species_graminoid_woodiness_526.csv")
  ),
  predictor = list(
    fit = read_csv(file.path(downstream, "tables", "all_predictor_ML_fit_statistics.csv")),
    comparisons = read_csv(file.path(downstream, "tables", "all_predictor_ML_nested_comparisons.csv")),
    CR1 = read_csv(file.path(downstream, "tables", "all_predictor_full_REML_CR1.csv"))
  ),
  heterogeneity = list(
    components = read_csv(file.path(downstream, "tables", "response_scale_heterogeneity_components_526.csv")),
    synthesis = read_csv(file.path(downstream, "tables", "heterogeneity_resolution_synthesis_526.csv"))
  ),
  phylogeny = list(
    eigen_summary = read_csv(file.path(phylogeny, "pcoa_eigenvalue_summary.csv")),
    tree_rds = readRDS(file.path(phylogeny, "resident_phylogeny.rds")),
    tree_newick = ape::read.tree(file.path(phylogeny, "resident_phylogeny.newick")),
    species_order = read_csv(file.path(phylogeny, "phylogeny_species_order.csv")),
    pev_scores = read_csv(file.path(phylogeny, "species_pev_scores.csv")) |>
      dplyr::select(dplyr::any_of(c("species_id", "PEV1", "PEV2", "PEV3"))),
    distance = readRDS(file.path(phylogeny, "phylogenetic_distance_matrix.rds")),
    zero_distance_audit = read_csv(file.path(phylogeny, "zero_distance_pair_audit.csv"))
  )
)

destination <- file.path(root, "tests", "fixtures", "canonical_transition_reference.rds")
temporary <- paste0(destination, ".tmp")
saveRDS(fixture, temporary, version = 3, compress = "xz")
if (!file.rename(temporary, destination)) {
  unlink(temporary)
  stop("Could not install the transition fixture atomically.", call. = FALSE)
}

hash <- system2("shasum", c("-a", "256", destination), stdout = TRUE)
cat(hash, "\n")
