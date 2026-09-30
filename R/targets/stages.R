# Module: R/targets/stages
# Responsibility: Implement the named scientific workflow stage or shared infrastructure.
# Inputs: Function arguments and explicitly declared target dependencies.
# Outputs: Validated in-memory objects and declared workflow artifacts.
# Side effects: Limited to functions whose contracts explicitly write, copy, or render files.

# Stage registry for the active targets workflow.
#
# Target definitions remain ordinary targets objects. This registry assigns
# every definition to one scientific stage so `_targets.R` exposes a readable
# top-level workflow without changing dependency inference or cache semantics.

#' Return every target definition exactly once.
#'
#' @return A list of `targets` target-definition objects.
pipeline_target_definitions <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      cache <<- c(
        data_target_definitions(),
        model_target_definitions(),
        publication_target_definitions()
      )
    }
    cache
  }
})

# Target names are read directly from `tar_tangle()` rather than copied into a
# second workflow registry. The cache only avoids reparsing the QMD once per
# stage factory while the old factory files remain during migration review.
tangled_scientific_target_names <- local({
  cache <- NULL
  function(path = "analysis_workflow.qmd") {
    if (is.null(cache)) {
      targets <- tarchetypes::tar_tangle(path)
      cache <<- vapply(targets, function(target) target$name, character(1))
    }
    cache
  }
})

# These replaced reports, pass-through inputs, hidden formula registry, and
# approved performance outputs remain on disk until the migration map is
# approved, but they are no longer part of the active dependency graph.
#' Factory targets disconnected while migration is reviewed
#'
#' @return Character target names excluded from every active stage factory.
disconnected_factory_target_names <- function() {
  c(
    "report_data_foundation", "report_response_estimation",
    "report_resident_climate_model", "report_family_context_model",
    "report_predictor_models", "report_model_synthesis",
    "report_release_validation", "published_report_data_foundation",
    "published_report_response_estimation",
    "published_report_resident_climate_model",
    "published_report_family_context_model", "published_report_predictor_models",
    "published_report_model_synthesis", "published_report_release_validation",
    "winner_loser_input", "analysis_data_input", "climate_context_input",
    "predictor_formula_registry", "model_cross_validation_predictions",
    "model_predictive_skill", "model_cross_validation_predictions_file",
    "model_predictive_skill_file"
  )
}

#' Assign a target name to its workflow stage.
#'
#' Historical comparisons are deliberately assigned to `legacy` before any
#' scientific-name matching. The temporary legacy stage is removed only after
#' the clean rebuild is explicitly accepted.
#'
#' @param name A single target name.
#' @return One stage identifier.
pipeline_target_stage <- function(name) {
  stopifnot(length(name) == 1L, is.character(name), nzchar(name))

  if (
    startsWith(name, "legacy_") ||
      grepl("equivalence|migration", name) ||
      name %in% c(
        "resident_consistency_checks", "legacy_certification_validation",
        "legacy_certification_validation_file",
        "refactor_certification",
        "phylogeny_reproduction_validation",
        "phylogeny_reproduction_validation_file",
        "phylogeny_reference_artifact_validation",
        "phylogeny_reference_artifact_validation_file",
        "frozen_release_reference_file",
        "frozen_release_reference_validation"
      )
  ) {
    return("legacy")
  }

  if (
    grepl(
      "^(release_artifact|run_provenance|output_manifest|current_scientific_release|release_validation)",
      name
    )
  ) {
    return("release")
  }

  if (
    grepl(
      "^(figure|appendix|report_|published_report_|retained_figure|figures_manifest|figure_provenance|current_model_report|current_model_results)",
      name
    )
  ) {
    return("publication")
  }

  if (
    grepl(
      "^(model_performance|model_cross_validation|model_predictive|model_fit_r2|model_fixed_effects|model_random_|model_effects|heterogeneity)",
      name
    )
  ) {
    return("evaluation")
  }

  if (
    grepl(
      "^(trait_|prepared_|functional_|standardized_|pev_|phylogen|resident_phylogeny|resident_pcoa|pcoa_|zero_distance|predictor_|ldmc_|lma_|height_|graminoid_|woodiness_|full_)",
      name
    )
  ) {
    return("predictor_model")
  }

  if (
    grepl(
      "^(analysis_settings|resident_|species_|family_|climate_scaling|resident_climate)",
      name
    )
  ) {
    return("resident_model")
  }

  if (grepl("^(auc_|winner_loser|temporal_|flow_validation)", name)) {
    return("response")
  }

  if (
    grepl(
      "^(try_|chelsa_|site_year_climate|experiment_climate|taxonomy_family|functional_groups|analysis_data|descriptive_analysis_data|climate_context|transformed_data)",
      name
    )
  ) {
    return("covariate")
  }

  if (
    grepl(
      "^(current_input|input_inventory|site_metadata|species_mapping|community|sites$|sites_|experiments$|experiments_|taxonomy_overrides|vegetation_|omitted_|included_|analysis_community)",
      name
    )
  ) {
    return("community")
  }

  if (
    name == "config" || startsWith(name, "active_source_") ||
      grepl("(_input|_input_file|_decisions_input)$", name) ||
      name %in% c("input_rdata_file", "try_export_file")
  ) {
    return("source")
  }

  # File-writing and presentation support targets are publication concerns.
  if (grepl("(_file|_files|_manifest)$", name)) {
    return("publication")
  }

  stop("No workflow stage is registered for target `", name, "`.", call. = FALSE)
}

#' Return the registered target-stage inventory.
#'
#' @return A data frame with one row per target definition.
pipeline_target_stage_table <- function() {
  definitions <- pipeline_target_definitions()
  names <- vapply(definitions, function(target) target$name, character(1))
  if (anyDuplicated(names)) {
    stop("Target definitions must have unique names.", call. = FALSE)
  }
  data.frame(
    target = names,
    stage = vapply(names, pipeline_target_stage, character(1)),
    stringsAsFactors = FALSE
  )
}

#' Select all target definitions registered to one stage.
#'
#' @param stage One stage identifier.
#' @return A list of target-definition objects.
targets_for_stage <- function(stage) {
  definitions <- pipeline_target_definitions()
  inventory <- pipeline_target_stage_table()
  keep <- inventory$stage == stage &
    !inventory$target %in% tangled_scientific_target_names() &
    !inventory$target %in% disconnected_factory_target_names()
  definitions[keep]
}

#' Define source and configuration targets.
#' @return Source-stage target definitions.
source_targets <- function() targets_for_stage("source")

#' Define community import and taxonomy targets.
#' @return Community-stage target definitions.
community_targets <- function() targets_for_stage("community")

#' Define climate, trait, and other covariate targets.
#' @return Covariate-stage target definitions.
covariate_targets <- function() targets_for_stage("covariate")

#' Define AUC, winner/loser, and temporal-response targets.
#' @return Response-stage target definitions.
response_targets <- function() targets_for_stage("response")

#' Define resident climate, species-context, and family targets.
#' @return Resident-model target definitions.
resident_model_targets <- function() targets_for_stage("resident_model")

#' Define trait, functional, and phylogenetic predictor targets.
#' @return Predictor-model target definitions.
predictor_model_targets <- function() targets_for_stage("predictor_model")

#' Define heterogeneity, performance, and diagnostic targets.
#' @return Evaluation-stage target definitions.
evaluation_targets <- function() targets_for_stage("evaluation")

#' Define artifact, figure, and report targets.
#' @return Publication-stage target definitions.
publication_targets <- function() targets_for_stage("publication")

#' Define current release targets with no historical inputs.
#' @return Current release-stage target definitions.
release_targets <- function() targets_for_stage("release")

#' Define temporary historical-certification targets.
#'
#' These targets are the only permitted consumers of frozen historical result
#' fixtures. The entire factory is deleted after explicit Phase 2 acceptance.
#'
#' @return Temporary legacy-certification target definitions.
legacy_certification_targets <- function() targets_for_stage("legacy")
