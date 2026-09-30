# Module: R/targets/definitions/models
# Responsibility: Implement the named scientific workflow stage or shared infrastructure.
# Inputs: Function arguments and explicitly declared target dependencies.
# Outputs: Validated in-memory objects and declared workflow artifacts.
# Side effects: Limited to functions whose contracts explicitly write, copy, or render files.

#' Model Target Definitions.
#'
#' @return An R object described by the function name; its structure is validated before later stages use it.
model_target_definitions <- function() {
  list(
    tar_target(analysis_settings, resident_analysis_spec(config$root)),
    tar_target(predictor_formula_registry, predictor_formulas()),
    tar_target(
      legacy_transition_reference_fixture_file,
      file.path(
        config$root, "tests", "fixtures", "canonical_transition_reference.rds"
      ),
      format = "file"
    ),
    tar_target(
      legacy_transition_reference,
      readRDS(legacy_transition_reference_fixture_file)
    ),
    tar_target(
      legacy_artifact_ledger_file,
      file.path(config$root, "config", "reference_artifacts.csv"),
      format = "file"
    ),
    tar_target(
      legacy_artifact_ledger,
      read_reference_csv(legacy_artifact_ledger_file)
    ),
    tar_target(
      legacy_file_artifacts,
      file.path(
        config$root,
        legacy_artifact_ledger$path[
          legacy_artifact_ledger$checksum_method == "file_sha256"
        ]
      ),
      format = "file"
    ),

    # Compute the downstream boundary from active upstream file targets. The
    # explicit CSV round trip is part of the anchored Step 22 interface: it
    # preserves the exact double representation used by the original runner
    # while keeping every read visible in the target graph.
    tar_target(
      winner_loser_input,
      read_reference_csv(winner_loser_file)
    ),
    tar_target(
      analysis_data_input,
      read_reference_csv(analysis_data_file)
    ),
    tar_target(
      climate_context_input,
      read_reference_csv(climate_context_experiments_file)
    ),
    tar_target(
      resident_auc_unscaled,
      build_resident_auc_data_from_objects(
        winner_loser_input,
        analysis_data_input,
        climate_context_input,
        "max"
      )
    ),
    tar_target(
      resident_scaled,
      scale_resident_auc_data(resident_auc_unscaled)
    ),
    tar_target(
      resident_auc_data,
      resident_scaled$data |>
        dplyr::mutate(TP_z = .data$T_z * .data$P_z)
    ),
    tar_target(resident_model_data, resident_auc_data),
    tar_target(
      resident_climate_scaling,
      resident_scaled$scaling |>
        dplyr::filter(.data$scaled_field %in% c("T_z", "P_z", "W_z")) |>
        dplyr::mutate(
          n_distinct_experiments = dplyr::n_distinct(
            resident_auc_data$experiment_id
          )
        )
    ),
    tar_target(
      resident_sampling_variance_matrix,
      build_sampling_variance_matrix(resident_auc_data, "auc_se")
    ),
    tar_target(
      resident_sampling_variance_check,
      validate_sampling_variance_matrix(
        resident_sampling_variance_matrix, resident_auc_data, "auc_se"
      )
    ),
    tar_target(
      legacy_resident_auc_data,
      legacy_transition_reference$downstream_auc_data
    ),
    tar_target(
      legacy_resident_V_auc,
      legacy_transition_reference$downstream_V_auc
    ),
    tar_target(
      resident_input_equivalence,
      compare_resident_inputs(
        resident_auc_data, resident_sampling_variance_matrix,
        legacy_resident_auc_data, legacy_resident_V_auc
      )
    ),
    tar_target(
      resident_scope,
      validate_resident_scope_artifacts(
        resident_auc_data, resident_sampling_variance_matrix, analysis_settings
      )
    ),
    tar_target(
      resident_sampling_variance_validation,
      resident_sampling_variance_release_check(
        resident_scope$sampling_variance_validation
      )
    ),

    tar_target(legacy_resident_climate_model_signature, legacy_transition_reference$mboth_model_signature),
    tar_target(
      legacy_resident_climate_model_bundle,
      legacy_transition_reference$mboth_bundle
    ),
    tar_target(
      resident_climate_model,
      fit_meta_regression(
        resident_model_data,
        resident_sampling_variance_matrix,
        analysis_settings$climate_formula,
        method = "REML",
        label = "resident climate REML",
        random = analysis_settings$species_random
      )
    ),
    tar_target(
      legacy_resident_climate_model_equivalence,
      compare_model_fit(resident_climate_model, legacy_resident_climate_model_signature)
    ),
    tar_target(
      resident_climate_variance_profile_specification,
      variance_profile_specifications(),
      iteration = "list"
    ),
    tar_target(
      resident_climate_variance_profile_component,
      profile_model_variance(
        resident_climate_model, resident_climate_variance_profile_specification
      ),
      pattern = map(resident_climate_variance_profile_specification),
      iteration = "list"
    ),
    tar_target(
      resident_climate_variance_profiles,
      dplyr::bind_rows(resident_climate_variance_profile_component) |>
        dplyr::mutate(
          warnings = if (all(is.na(.data$warnings) | .data$warnings == "")) {
            rep(NA, dplyr::n())
          } else {
            dplyr::na_if(as.character(.data$warnings), "")
          },
          error = if (all(is.na(.data$error) | .data$error == "")) {
            rep(NA, dplyr::n())
          } else {
            dplyr::na_if(as.character(.data$error), "")
          }
        )
    ),
    tar_target(
      legacy_resident_climate_variance_profile_equivalence,
      assert_reference_table(
        resident_climate_variance_profiles,
        legacy_resident_climate_model_bundle$variance_profiles,
        "M_both variance profiles", "component"
      )
    ),
    tar_target(
      resident_climate_model_results,
      assemble_climate_model_results(
        resident_climate_model, resident_climate_variance_profiles
      )
    ),

    tar_target(
      legacy_species_classification,
      legacy_transition_reference$species_classification
    ),
    tar_target(
      species_context_bundle,
      build_species_context_bundle(
        resident_climate_model,
        resident_model_data,
        analysis_settings$climate_formula
      )
    ),
    tar_target(
      species_classification,
      species_context_bundle$classification
    ),
    tar_target(
      species_classification_equivalence,
      compare_classification_reference(
        species_classification,
        legacy_species_classification,
        "species_id"
      )
    ),
    tar_target(
      species_callout_examples,
      select_species_callout_examples(species_classification)
    ),
    tar_target(
      species_callout_examples_file,
      write_species_callout_examples(
        species_callout_examples,
        file.path(config$output_dir, "tables", "species_context_callouts.csv")
      ),
      format = "file"
    ),
    tar_target(
      species_classification_summary,
      summarize_reference_classification(species_classification, "Species")
    ),
    tar_target(
      legacy_family_classification,
      legacy_transition_reference$family$classification
    ),
    tar_target(
      legacy_family_model_summary,
      legacy_transition_reference$family$model_summary
    ),
    tar_target(
      legacy_family_auc_data,
      legacy_transition_reference$family$auc_data
    ),
    tar_target(
      legacy_family_V_reference_container,
      legacy_transition_reference$family$V
    ),
    tar_target(
      legacy_family_model,
      legacy_transition_reference$family$model_signature
    ),
    tar_target(
      legacy_family_support,
      legacy_transition_reference$family$support
    ),
    tar_target(
      legacy_family_random,
      legacy_transition_reference$family$random
    ),
    tar_target(
      legacy_family_cover_audit,
      legacy_transition_reference$family$cover_audit
    ),

    # The family branch reconstructs family cover and estimates family AUC
    # effects once per experiment. Squared contrast SEs supply the diagonal
    # sampling variances used by the family meta-regression.
    tar_target(family_community_input, readRDS(community_file)),
    tar_target(
      family_experiments_input,
      read_reference_csv(experiments_file)
    ),
    tar_target(
      family_override_table,
      family_overrides(resident_species_overrides)
    ),
    tar_target(
      family_resident_membership,
      build_family_resident_membership(
        winner_loser_input, analysis_data_input,
        family_override_table
      )
    ),
    tar_target(
      family_plot_auc_bundle,
      build_family_plot_auc(
        family_community_input,
        family_experiments_input,
        family_resident_membership
      )
    ),
    tar_target(
      family_auc_bundle,
      build_family_auc_effects(
        family_plot_auc_bundle$plot_auc,
        family_resident_membership,
        resident_model_data
      )
    ),
    tar_target(family_auc_data, family_auc_bundle$data),
    tar_target(
      family_sampling_variance_matrix,
      family_auc_bundle$sampling_variance_matrix
    ),
    tar_target(
      family_auc_data_equivalence,
      assert_reference_table(
        family_auc_data, legacy_family_auc_data,
        "family_auc_data", "effect_id"
      )
    ),
    tar_target(
      family_sampling_variance_equivalence,
      compare_family_sampling_variances(
        family_sampling_variance_matrix, legacy_family_V_reference_container
      )
    ),
    tar_target(
      family_cover_audit_equivalence,
      assert_reference_table(
        family_plot_auc_bundle$cover_audit,
        legacy_family_cover_audit,
        "family_plot_cover_reconstruction_audit", "experiment_id"
      )
    ),
    tar_target(
      family_model_bundle,
      fit_family_aggregated_model(family_auc_bundle)
    ),
    tar_target(
      family_model_equivalence,
      compare_model_fit(
        family_model_bundle, legacy_family_model
      )
    ),
    tar_target(
      family_support,
      family_support_table(family_model_bundle$data) |>
        dplyr::mutate(family_id = as.character(.data$family_id))
    ),
    tar_target(
      family_support_equivalence,
      assert_reference_table(
        family_support, legacy_family_support,
        "family_support", "family_id"
      )
    ),
    tar_target(
      family_random_bundle,
      extract_family_random_coefficients(family_model_bundle, family_support)
    ),
    tar_target(
      family_random_coefficients,
      family_random_bundle$coefficients
    ),
    tar_target(
      family_random_equivalence,
      assert_reference_table(
        family_random_coefficients, legacy_family_random,
        "family_random_coefficients", "family_id"
      )
    ),
    tar_target(
      family_classification,
      classify_family_context(
        family_model_bundle, family_support,
        family_random_bundle
      )
    ),
    tar_target(
      family_classification_equivalence,
      compare_classification_reference(
        family_classification,
        legacy_family_classification,
        "family_id"
      )
    ),
    tar_target(
      family_classification_summary,
      summarize_reference_classification(family_classification, "Family")
    ),
    tar_target(
      family_model_summary,
      summarize_family_model(family_model_bundle)
    ),
    tar_target(
      family_model_summary_equivalence,
      compare_family_model_summary_reference(
        family_model_summary, legacy_family_model_summary
      )
    ),
    tar_target(
      family_bundle_validation,
      validate_family_bundle(
        family_auc_data,
        family_sampling_variance_matrix,
        family_classification,
        family_model_summary
      )
    ),
    tar_target(
      species_family_resolution_comparison,
      compare_classification_resolutions(
        species_classification, family_classification
      )
    ),

    tar_target(
      legacy_trait_scaling,
      legacy_transition_reference$trait_scaling
    ),
    tar_target(
      legacy_trait_species,
      legacy_transition_reference$trait_species
    ),
    tar_target(
      trait_summary_input,
      read_reference_csv(try_summary_file)
    ),
    tar_target(
      trait_species,
      build_trait_species(
        trait_summary_input,
        resident_auc_data$species_id
      )
    ),
    tar_target(
      trait_species_equivalence,
      compare_trait_species_reference(
        trait_species,
        legacy_trait_species
      )
    ),
    tar_target(
      prepared_ldmc,
      prepare_trait_data(
        trait_species, resident_auc_data, resident_sampling_variance_matrix, "LDMC"
      )
    ),
    tar_target(
      prepared_lma,
      prepare_trait_data(
        trait_species, resident_auc_data, resident_sampling_variance_matrix, "LMA"
      )
    ),
    tar_target(
      prepared_height,
      prepare_trait_data(
        trait_species, resident_auc_data, resident_sampling_variance_matrix,
        "vegetative_height"
      )
    ),
    tar_target(
      prepared_trait_models,
      list(
        LDMC = prepared_ldmc,
        LMA = prepared_lma,
        vegetative_height = prepared_height
      )
    ),
    tar_target(
      trait_scaling,
      dplyr::bind_rows(lapply(prepared_trait_models, `[[`, "scaling"))
    ),
    tar_target(
      trait_scaling_equivalence,
      assert_reference_table(
        trait_scaling, legacy_trait_scaling,
        "trait scaling", "trait"
      )
    ),
    tar_target(
      legacy_climate_scaling,
      legacy_transition_reference$climate_scaling
    ),
    tar_target(climate_scaling, resident_climate_scaling),
    tar_target(
      climate_scaling_equivalence,
      assert_reference_table(
        climate_scaling, legacy_climate_scaling,
        "climate scaling", "scaled_field"
      )
    ),
    tar_target(legacy_pev_scaling, legacy_transition_reference$pev_scaling),
    tar_target(
      standardized_pev_species,
      standardize_unique_species_predictors(
        resident_pev_scores_oriented |>
          dplyr::select("species_id", "PEV1", "PEV2", "PEV3"),
        c("PEV1", "PEV2", "PEV3")
      )
    ),
    tar_target(pev_scaling, standardized_pev_species$scaling),
    tar_target(
      pev_scaling_equivalence,
      assert_reference_table(
        pev_scaling, legacy_pev_scaling,
        "PEV scaling", "PEV"
      )
    ),
    tar_target(
      phylogeny_auc_data,
      resident_auc_data |>
        dplyr::left_join(
          standardized_pev_species$data |>
            dplyr::select("species_id", "PEV1_z", "PEV2_z", "PEV3_z"),
          by = "species_id"
        )
    ),
    tar_target(
      legacy_functional_species,
      legacy_transition_reference$functional_species
    ),
    tar_target(
      functional_groups_input,
      read_reference_csv(functional_groups_file)
    ),
    tar_target(
      functional_species,
      build_functional_species(
        functional_groups_input,
        resident_auc_data$species_id,
        resident_species_overrides
      )
    ),
    tar_target(
      functional_source_coverage_audit,
      audit_functional_source_coverage(
        functional_groups_input,
        functional_species,
        resident_species_overrides
      )
    ),
    tar_target(
      functional_species_equivalence,
      compare_functional_species_reference(
        functional_species, legacy_functional_species
      )
    ),
    tar_target(
      functional_auc_data,
      prepare_functional_predictors(resident_auc_data, functional_species)
    ),
    tar_target(
      functional_climate_coverage,
      dplyr::bind_rows(
        summarize_functional_climate_coverage(
          functional_auc_data, "GraminoidStatus", "Graminoid status"
        ),
        summarize_functional_climate_coverage(
          functional_auc_data, "Woodiness", "Woodiness"
        )
      )
    ),
    tar_target(
      ldmc_predictor_specification,
      trait_predictor_specification("LDMC", "LDMC", prepared_ldmc)
    ),
    tar_target(
      lma_predictor_specification,
      trait_predictor_specification("LMA", "LMA", prepared_lma)
    ),
    tar_target(
      height_predictor_specification,
      trait_predictor_specification(
        "Height", "Height", prepared_height
      )
    ),
    tar_target(
      graminoid_predictor_specification,
      functional_predictor_specification(
        "Graminoid", "Graminoid status", functional_auc_data,
        resident_sampling_variance_matrix
      )
    ),
    tar_target(
      woodiness_predictor_specification,
      functional_predictor_specification(
        "Woodiness", "Woodiness", functional_auc_data,
        resident_sampling_variance_matrix
      )
    ),
    tar_target(
      phylogeny_predictor_specification,
      build_phylogeny_predictor_specification(
        phylogeny_auc_data, resident_sampling_variance_matrix
      )
    ),
    tar_target(
      predictor_specifications,
      list(
        LDMC = ldmc_predictor_specification,
        LMA = lma_predictor_specification,
        Height = height_predictor_specification,
        Graminoid = graminoid_predictor_specification,
        Woodiness = woodiness_predictor_specification,
        Phylogeny = phylogeny_predictor_specification
      )
    ),
    tar_target(
      predictor_specification_validation,
      dplyr::bind_rows(lapply(predictor_specifications, function(spec) {
        validate_predictor_specification(
          spec,
          analysis_settings$climate_formula,
          resident_auc_data$effect_id
        )
      }))
    ),
    tar_target(
      legacy_predictor_model_fit,
      legacy_transition_reference$predictor$fit
    ),
    tar_target(
      legacy_predictor_model_comparisons,
      legacy_transition_reference$predictor$comparisons
    ),
    tar_target(legacy_predictor_CR1, legacy_transition_reference$predictor$CR1),
    tar_target(
      legacy_predictor_pair_CR1,
      legacy_predictor_pair_cr1(legacy_predictor_CR1)
    ),
    tar_target(
      resident_climate_model_ML,
      fit_meta_regression(
        resident_auc_data,
        resident_sampling_variance_matrix,
        analysis_settings$climate_formula,
        "ML",
        "resident climate ML"
      )
    ),
    tar_target(
      ldmc_model_pair,
      fit_predictor_model_pair(
        ldmc_predictor_specification,
        analysis_settings$climate_formula
      )
    ),
    tar_target(
      lma_model_pair,
      fit_predictor_model_pair(
        lma_predictor_specification,
        analysis_settings$climate_formula
      )
    ),
    tar_target(
      height_model_pair,
      fit_predictor_model_pair(
        height_predictor_specification,
        analysis_settings$climate_formula
      )
    ),
    tar_target(
      graminoid_model_pair,
      fit_predictor_model_pair(
        graminoid_predictor_specification,
        analysis_settings$climate_formula,
        resident_climate_model_ML,
        resident_climate_model
      )
    ),
    tar_target(
      woodiness_model_pair,
      fit_predictor_model_pair(
        woodiness_predictor_specification,
        analysis_settings$climate_formula,
        resident_climate_model_ML,
        resident_climate_model
      )
    ),
    tar_target(
      phylogeny_model_pair,
      fit_predictor_model_pair(
        phylogeny_predictor_specification,
        analysis_settings$climate_formula,
        resident_climate_model_ML,
        resident_climate_model
      )
    ),
    # Validation-only historical refit. This branch deliberately consumes the
    # frozen transition PEV values, but no active scientific model depends on it.
    tar_target(
      legacy_standardized_pev_species_reference,
      standardize_unique_species_predictors(
        legacy_resident_pev_scores |>
          dplyr::select("species_id", "PEV1", "PEV2", "PEV3"),
        c("PEV1", "PEV2", "PEV3")
      )
    ),
    tar_target(
      legacy_phylogeny_auc_data_reference,
      resident_auc_data |>
        dplyr::left_join(
          legacy_standardized_pev_species_reference$data |>
            dplyr::select("species_id", "PEV1_z", "PEV2_z", "PEV3_z"),
          by = "species_id"
        )
    ),
    tar_target(
      legacy_phylogeny_predictor_specification_reference,
      build_phylogeny_predictor_specification(
        legacy_phylogeny_auc_data_reference,
        resident_sampling_variance_matrix
      )
    ),
    tar_target(
      legacy_phylogeny_model_pair_reference,
      fit_predictor_model_pair(
        legacy_phylogeny_predictor_specification_reference,
        analysis_settings$climate_formula,
        resident_climate_model_ML,
        resident_climate_model
      )
    ),
    tar_target(
      predictor_model_tables,
      assemble_predictor_model_tables(
        ldmc_model_pair,
        lma_model_pair,
        height_model_pair,
        graminoid_model_pair,
        woodiness_model_pair,
        phylogeny_model_pair
      )
    ),
    # Report-facing assemblies stay branch-specific. This prevents a change to
    # one scientific branch from invalidating unrelated reports.
    tar_target(
      trait_model_tables,
      assemble_predictor_model_tables(
        ldmc_model_pair,
        lma_model_pair,
        height_model_pair
      )
    ),
    tar_target(
      functional_model_tables,
      assemble_predictor_model_tables(
        graminoid_model_pair,
        woodiness_model_pair
      )
    ),
    tar_target(
      phylogeny_model_tables,
      assemble_predictor_model_tables(phylogeny_model_pair)
    ),
    tar_target(
      predictor_model_fit,
      predictor_model_tables$fit_statistics
    ),
    tar_target(
      predictor_model_comparisons,
      predictor_model_tables$comparisons
    ),
    tar_target(predictor_CR1, predictor_model_tables$CR1),
    tar_target(
      legacy_predictor_pair_fit,
      legacy_predictor_pair_fit_statistics(legacy_predictor_model_fit)
    ),
    tar_target(
      legacy_predictor_pair_comparisons,
      derive_legacy_predictor_pair_comparisons(legacy_predictor_model_fit)
    ),
    tar_target(
      predictor_model_fit_equivalence,
      validate_reference_table(
        predictor_model_fit, legacy_predictor_pair_fit,
        "predictor ML fit statistics"
      )
    ),
    tar_target(
      predictor_model_comparisons_equivalence,
      validate_reference_table(
        predictor_model_comparisons, legacy_predictor_pair_comparisons,
        "predictor climate-only versus full ML comparisons"
      )
    ),
    tar_target(
      predictor_CR1_equivalence,
      validate_reference_table(
        predictor_CR1, legacy_predictor_pair_CR1,
        "predictor full REML CR1"
      )
    ),
    tar_target(
      trait_model_fit_equivalence,
      assert_reference_table(
        trait_model_tables$fit_statistics,
        dplyr::filter(
          legacy_predictor_pair_fit,
          .data$predictor %in% c("LDMC", "LMA", "Height")
        ),
        "trait-model ML fit statistics"
      )
    ),
    tar_target(
      trait_model_comparisons_equivalence,
      assert_reference_table(
        trait_model_tables$comparisons,
        dplyr::filter(
          legacy_predictor_pair_comparisons,
          .data$predictor %in% c("LDMC", "LMA", "Height")
        ),
        "trait-model ML nested comparisons"
      )
    ),
    tar_target(
      trait_model_CR1_equivalence,
      assert_reference_table(
        trait_model_tables$CR1,
        dplyr::filter(
          legacy_predictor_pair_CR1,
          .data$predictor %in% c("LDMC", "LMA", "Height")
        ),
        "trait-model full REML CR1"
      )
    ),
    tar_target(
      functional_model_fit_equivalence,
      assert_reference_table(
        functional_model_tables$fit_statistics,
        dplyr::filter(
          legacy_predictor_pair_fit,
          .data$predictor %in% c("Graminoid status", "Woodiness")
        ),
        "functional-model ML fit statistics"
      )
    ),
    tar_target(
      functional_model_comparisons_equivalence,
      assert_reference_table(
        functional_model_tables$comparisons,
        dplyr::filter(
          legacy_predictor_pair_comparisons,
          .data$predictor %in% c("Graminoid status", "Woodiness")
        ),
        "functional-model ML nested comparisons"
      )
    ),
    tar_target(
      functional_model_CR1_equivalence,
      assert_reference_table(
        functional_model_tables$CR1,
        dplyr::filter(
          legacy_predictor_pair_CR1,
          .data$predictor %in% c("Graminoid status", "Woodiness")
        ),
        "functional-model full REML CR1"
      )
    ),
    tar_target(
      phylogeny_model_fit_equivalence,
      validate_reference_table(
        phylogeny_model_tables$fit_statistics,
        dplyr::filter(
          legacy_predictor_pair_fit,
          .data$predictor == "Phylogenetic position (PEV1-PEV3)"
        ),
        "phylogeny-model ML fit statistics"
      )
    ),
    tar_target(
      phylogeny_model_comparisons_equivalence,
      validate_reference_table(
        phylogeny_model_tables$comparisons,
        dplyr::filter(
          legacy_predictor_pair_comparisons,
          .data$predictor == "Phylogenetic position (PEV1-PEV3)"
        ),
        "phylogeny-model ML nested comparisons"
      )
    ),
    tar_target(
      phylogeny_model_CR1_equivalence,
      validate_reference_table(
        phylogeny_model_tables$CR1,
        dplyr::filter(
          legacy_predictor_pair_CR1,
          .data$predictor == "Phylogenetic position (PEV1-PEV3)"
        ),
        "phylogeny-model full REML CR1"
      )
    ),

    tar_target(
      legacy_heterogeneity_components,
      legacy_transition_reference$heterogeneity$components
    ),
    tar_target(
      legacy_heterogeneity_synthesis,
      legacy_transition_reference$heterogeneity$synthesis
    ),
    tar_target(
      heterogeneity_bundle,
      assemble_heterogeneity_synthesis(
        ldmc_model_pair,
        lma_model_pair,
        height_model_pair,
        graminoid_model_pair,
        woodiness_model_pair,
        phylogeny_model_pair
      )
    ),
    tar_target(
      heterogeneity_components,
      heterogeneity_bundle$components
    ),
    tar_target(
      heterogeneity_synthesis,
      heterogeneity_bundle$synthesis
    ),
    tar_target(
      legacy_heterogeneity_bundle_reference,
      assemble_heterogeneity_synthesis(
        ldmc_model_pair,
        lma_model_pair,
        height_model_pair,
        graminoid_model_pair,
        woodiness_model_pair,
        legacy_phylogeny_model_pair_reference
      )
    ),
    tar_target(
      phylogeny_migration_equivalence,
      build_phylogeny_migration_equivalence(
        resident_pev_scores_oriented,
        legacy_resident_pev_scores,
        phylogeny_model_pair,
        legacy_phylogeny_model_pair_reference,
        heterogeneity_bundle,
        legacy_heterogeneity_bundle_reference
      )
    ),
    tar_target(
      phylogeny_migration_release_validation,
      phylogeny_migration_release_check(phylogeny_migration_equivalence)
    ),
    tar_target(
      resident_consistency_checks,
      build_resident_consistency_checks(
        resident_auc_data,
        winner_loser_input,
        analysis_data_input,
        analysis_settings,
        family_resident_membership,
        family_plot_auc_bundle$cover_audit,
        resident_climate_scaling,
        prepared_trait_models,
        standardized_pev_species,
        phylogeny_auc_data,
        predictor_specifications,
        predictor_specification_validation,
        functional_species,
        heterogeneity_synthesis,
        species_classification
      )
    ),
    tar_target(
      heterogeneity_components_equivalence,
      validate_reference_table(
        heterogeneity_components, legacy_heterogeneity_components,
        "response-scale heterogeneity components"
      )
    ),
    tar_target(
      heterogeneity_synthesis_equivalence,
      validate_reference_table(
        heterogeneity_synthesis, legacy_heterogeneity_synthesis,
        "heterogeneity synthesis"
      )
    ),

    # Validation-only predictive performance. These folds preserve the
    # canonical specifications and scaling, but none replaces a canonical fit.
    tar_target(
      model_fit_r2,
      assemble_model_fit_r2(
        resident_climate_model,
        ldmc_model_pair,
        lma_model_pair,
        height_model_pair,
        graminoid_model_pair,
        woodiness_model_pair,
        phylogeny_model_pair
      )
    ),
    tar_target(
      model_fixed_effects,
      assemble_model_fixed_effects(
        resident_climate_model_results$CR1,
        predictor_CR1,
        model_fit_r2
      )
    ),
    tar_target(
      model_random_profile_job,
      build_predictor_variance_profile_jobs(
        ldmc_model_pair,
        lma_model_pair,
        height_model_pair,
        graminoid_model_pair,
        woodiness_model_pair,
        phylogeny_model_pair
      ),
      iteration = "list"
    ),
    tar_target(
      model_random_profile_component,
      profile_predictor_variance_component(
        model_random_profile_job
      ),
      pattern = map(model_random_profile_job),
      iteration = "list"
    ),
    tar_target(
      model_random_effect_profiles,
      assemble_model_random_effect_profiles(
        resident_climate_variance_profiles,
        model_random_profile_component
      )
    ),
    tar_target(
      model_random_effect_correlations,
      assemble_model_random_correlations(
        resident_climate_model,
        ldmc_model_pair,
        lma_model_pair,
        height_model_pair,
        graminoid_model_pair,
        woodiness_model_pair,
        phylogeny_model_pair
      )
    ),
    tar_target(
      model_effects_validation,
      validate_model_effects(
        model_fixed_effects,
        model_random_effect_profiles,
        model_random_effect_correlations
      )
    ),
    tar_target(
      model_performance_gradient,
      sort(unique(as.character(resident_auc_data$gradient_id))),
      iteration = "vector"
    ),
    tar_target(
      model_performance_fold,
      run_model_performance_gradient_fold(
        model_performance_gradient,
        predictor_specifications,
        resident_auc_data,
        resident_sampling_variance_matrix,
        analysis_settings$climate_formula
      ),
      pattern = map(model_performance_gradient),
      iteration = "list"
    ),
    tar_target(
      model_performance_cv,
      combine_model_performance_folds(model_performance_fold)
    ),
    tar_target(
      model_cross_validation_predictions,
      model_performance_cv$predictions
    ),
    tar_target(
      model_cross_validation_metrics,
      summarise_model_cross_validation(
        model_performance_cv$predictions,
        model_performance_cv$gaussian_scores,
        bootstrap_replicates = 5000L,
        seed = 20260914L
      )
    ),
    tar_target(
      model_predictive_skill,
      summarise_model_predictive_skill(
        model_performance_cv$predictions,
        model_performance_cv$gaussian_scores,
        bootstrap_replicates = 5000L,
        seed = 20260914L
      )
    ),
    tar_target(
      model_performance_validation,
      validate_model_performance(
        model_performance_fold,
        model_performance_cv,
        model_fit_r2,
        model_cross_validation_metrics,
        model_predictive_skill
      )
    ),
    tar_target(
      model_cross_validation_predictions_file,
      write_csv_atomic(
        model_cross_validation_predictions,
        file.path(
          analysis_settings$output$data,
          "model_cross_validation_predictions.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      model_fit_r2_file,
      write_csv_atomic(
        model_fit_r2,
        file.path(
          analysis_settings$output$tables,
          "model_fit_r2.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      model_cross_validation_metrics_file,
      write_csv_atomic(
        model_cross_validation_metrics,
        file.path(
          analysis_settings$output$tables,
          "model_cross_validation_metrics.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      model_predictive_skill_file,
      write_csv_atomic(
        model_predictive_skill,
        file.path(
          analysis_settings$output$tables,
          "model_predictive_skill.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      model_performance_validation_file,
      write_csv_atomic(
        model_performance_validation,
        file.path(
          analysis_settings$output$audits,
          "model_performance_validation.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      model_fixed_effects_file,
      write_csv_atomic(
        model_fixed_effects,
        file.path(
          analysis_settings$output$tables,
          "model_fixed_effects_robust_inference.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      model_random_effect_profiles_file,
      write_csv_atomic(
        model_random_effect_profiles,
        file.path(
          analysis_settings$output$tables,
          "model_random_effect_profile_intervals.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      model_random_effect_correlations_file,
      write_csv_atomic(
        model_random_effect_correlations,
        file.path(
          analysis_settings$output$tables,
          "model_random_effect_correlations.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      model_effects_validation_file,
      write_csv_atomic(
        model_effects_validation,
        file.path(
          analysis_settings$output$audits,
          "model_effects_validation.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      model_effects_release_validation,
      {
        model_fixed_effects_file
        model_random_effect_profiles_file
        model_random_effect_correlations_file
        model_effects_validation_file
        assert_model_effects_validation(model_effects_validation)
      },
      description = paste(
        "Current fixed-effect CR1 and random-effect profile checks for the",
        "resident climate and six full predictor models"
      )
    ),
    tar_target(
      model_performance_release_validation,
      {
        model_cross_validation_predictions_file
        model_fit_r2_file
        model_cross_validation_metrics_file
        model_predictive_skill_file
        model_performance_validation_file
        assert_model_performance_validation(
          model_performance_validation
        )
      }
    ),

    tar_target(
      legacy_pcoa_eigen_summary,
      legacy_transition_reference$phylogeny$eigen_summary
    ),
    tar_target(
      legacy_resident_phylogeny,
      legacy_transition_reference$phylogeny$tree_rds
    ),
    tar_target(
      phylogeny_placement_decisions_file,
      file.path(config$root, "config", "phylogeny_placement_decisions.csv"),
      format = "file"
    ),
    tar_target(
      phylogeny_placement_decisions,
      read_reference_csv(phylogeny_placement_decisions_file)
    ),
    tar_target(
      resident_phylogeny_rebuild,
      build_phylogeny_from_placement_decisions(
        phylogeny_placement_decisions,
        resident_model_data$species_id
      )
    ),
    tar_target(
      legacy_resident_phylogeny_newick,
      legacy_transition_reference$phylogeny$tree_newick
    ),
    tar_target(
      legacy_phylogeny_species_order,
      legacy_transition_reference$phylogeny$species_order
    ),
    tar_target(
      legacy_resident_phylogeny_newick_reference_restored,
      restore_phylogeny_tip_labels(
        legacy_resident_phylogeny_newick,
        legacy_phylogeny_species_order
      )
    ),
    # The canonical tree is rebuilt from the tracked placement-decision ledger
    # and the pinned V.PhyloMaker2 backbone. The old RDS and Newick remain
    # independent regression authorities until this gate passes.
    tar_target(resident_phylogeny_candidate, resident_phylogeny_rebuild$tree),
    tar_target(
      phylogeny_species_order,
      build_phylogeny_species_order(
        resident_phylogeny_candidate,
        phylogeny_placement_decisions,
        resident_model_data$species_id
      )
    ),
    tar_target(
      legacy_resident_pev_scores,
      legacy_transition_reference$phylogeny$pev_scores
    ),
    tar_target(
      resident_pcoa_candidate,
      recompute_phylogenetic_pcoa(
        resident_phylogeny_candidate,
        metadata = NULL,
        tip_order = phylogeny_species_order$species_id
      )
    ),
    tar_target(
      resident_pev_scores_recomputed,
      resident_pcoa_candidate$scores |>
        dplyr::select("species_id", "PEV1", "PEV2", "PEV3")
    ),
    tar_target(
      resident_pev_scores_oriented,
      resident_pev_scores_recomputed
    ),
    tar_target(
      pcoa_eigen_summary,
      resident_pcoa_candidate$eigen_summary
    ),
    tar_target(
      phylogenetic_positive_variation_PEV1_PEV3,
      validated_phylogenetic_variation(pcoa_eigen_summary)
    ),
    tar_target(
      legacy_phylogenetic_distance,
      legacy_transition_reference$phylogeny$distance
    ),
    tar_target(phylogenetic_distance, resident_pcoa_candidate$distance),
    tar_target(
      zero_distance_audit,
      build_zero_distance_audit(
        phylogenetic_distance,
        phylogeny_placement_decisions
      )
    ),
    tar_target(
      legacy_zero_distance_audit,
      legacy_transition_reference$phylogeny$zero_distance_audit
    ),
    tar_target(
      legacy_phylogeny_reference_artifact_validation,
      {
        legacy_file_artifacts
        phylogeny_package_installation_files
        verify_reference_artifacts(
          config$root, references = legacy_artifact_ledger
        )
      }
    ),
    tar_target(
      phylogeny_package_installation_files,
      {
        package_root <- find.package("V.PhyloMaker2")
        file.path(
          package_root,
          c(
            "DESCRIPTION", "R/V.PhyloMaker2.rdb", "R/V.PhyloMaker2.rdx",
            "data/Rdata.rdb", "data/Rdata.rds", "data/Rdata.rdx"
          )
        )
      },
      format = "file"
    ),
    tar_target(
      phylogeny_package_source_provenance,
      {
        phylogeny_package_installation_files
        phylogeny_package_provenance()
      }
    ),
    tar_target(
      phylogeny_package_validation,
      validate_phylogeny_package_provenance(
        phylogeny_package_source_provenance
      )
    ),
    tar_target(
      phylogeny_reproduction_validation,
      assert_phylogeny_validation(
        resident_phylogeny_rebuild$audit,
        phylogeny_package_validation,
        compare_phylogeny(
          resident_phylogeny_candidate, legacy_resident_phylogeny,
          phylogeny_species_order, legacy_phylogeny_species_order
        ),
        compare_phylogeny(
          legacy_resident_phylogeny_newick_reference_restored,
          legacy_resident_phylogeny,
          legacy_phylogeny_species_order,
          legacy_phylogeny_species_order
        ),
        compare_phylogenetic_distance(
          phylogenetic_distance, legacy_phylogenetic_distance
        ),
        compare_pev_coordinates(
          resident_pev_scores_recomputed,
          legacy_resident_pev_scores
        ),
        compare_pcoa_eigenstructure(
          pcoa_eigen_summary, legacy_pcoa_eigen_summary
        ),
        validate_zero_distance_placements(phylogenetic_distance, zero_distance_audit),
        compare_zero_distance_audit(
          zero_distance_audit,
          legacy_zero_distance_audit
        ),
        tibble::tibble(
          check = "reference artifact hashes",
          pass = all(legacy_phylogeny_reference_artifact_validation$hash_matches),
          detail = paste(nrow(legacy_phylogeny_reference_artifact_validation),
                         "configured artifacts")
        )
      )
    ),
    tar_target(
      current_scientific_release_validation,
      current_release_checks(
        resident_scope,
        predictor_specification_validation,
        predictor_model_tables,
        family_model_summary,
        resident_phylogeny_rebuild$audit,
        zero_distance_audit,
        heterogeneity_components,
        heterogeneity_synthesis
      ),
      description = paste(
        "Current cohort model family phylogeny predictor and heterogeneity",
        "checks with no historical comparison inputs"
      )
    ),
    tar_target(
      phylogeny_package_provenance_file,
      write_csv_atomic(
        phylogeny_package_source_provenance,
        file.path(config$output_dir, "audits", "phylogeny_package_provenance.csv")
      ),
      format = "file"
    ),
    tar_target(
      legacy_phylogeny_reference_artifact_validation_file,
      write_csv_atomic(
        legacy_phylogeny_reference_artifact_validation,
        file.path(
          config$output_dir, "audits",
          "legacy_phylogeny_reference_artifact_validation.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      phylogeny_reproduction_validation_file,
      write_csv_atomic(
        phylogeny_reproduction_validation,
        file.path(config$output_dir, "audits", "phylogeny_reproduction_validation.csv")
      ),
      format = "file"
    ),
    tar_target(
      phylogeny_migration_equivalence_file,
      write_csv_atomic(
        phylogeny_migration_equivalence,
        file.path(config$output_dir, "audits", "phylogeny_migration_equivalence.csv")
      ),
      format = "file"
    ),
    tar_target(
      phylogeny_migration_equivalence_report_file,
      write_text_atomic(
        format_phylogeny_migration_report(phylogeny_migration_equivalence),
        file.path(
          config$root, "reports", "rendered",
          "phylogeny_migration_equivalence.txt"
        )
      ),
      format = "file"
    ),
    # Canonical exports are deliberately scoped to their immediate scientific
    # objects. None depends on the late Step 22 validation assembly.
    tar_target(
      resident_data_file,
      write_csv_atomic(
        resident_auc_data,
        file.path(
          analysis_settings$output$data,
          "resident_effects.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      resident_sampling_variance_file,
      write_rds_atomic(
        resident_sampling_variance_matrix,
        file.path(
          analysis_settings$output$data,
          "resident_sampling_variance_matrix.rds"
        )
      ),
      format = "file"
    ),
    tar_target(
      resident_climate_model_file,
      write_rds_atomic(
        resident_climate_model$fit,
        file.path(
          analysis_settings$output$models,
          "resident_climate_model.rds"
        )
      ),
      format = "file"
    ),
    tar_target(
      resident_climate_fixed_effects_file,
      write_csv_atomic(
        resident_climate_model_results$fixed,
        file.path(
          analysis_settings$output$tables,
          "resident_climate_fixed_effects.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      resident_climate_robust_inference_file,
      write_csv_atomic(
        resident_climate_model_results$CR1,
        file.path(
          analysis_settings$output$tables,
          "resident_climate_robust_inference.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      resident_climate_random_effects_file,
      write_csv_atomic(
        resident_climate_model_results$random,
        file.path(
          analysis_settings$output$tables,
          "resident_climate_random_effects.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      resident_climate_variance_profiles_file,
      write_csv_atomic(
        resident_climate_variance_profiles,
        file.path(
          analysis_settings$output$tables,
          "resident_climate_variance_profiles.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      species_classification_file,
      write_csv_atomic(
        species_classification,
        file.path(
          analysis_settings$output$tables,
          "species_context_classification.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      family_data_file,
      write_csv_atomic(
        family_auc_data,
        file.path(
          analysis_settings$output$data,
          "family_auc_data.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      family_sampling_variance_file,
      write_rds_atomic(
        family_sampling_variance_matrix,
        file.path(
          analysis_settings$output$data,
          "family_sampling_variance_matrix.rds"
        )
      ),
      format = "file"
    ),
    tar_target(
      family_model_file,
      write_rds_atomic(
        family_model_bundle$fit,
        file.path(
          analysis_settings$output$models,
          "family_context_model.rds"
        )
      ),
      format = "file"
    ),
    tar_target(
      family_support_file,
      write_csv_atomic(
        family_support,
        file.path(
          analysis_settings$output$tables,
          "family_support.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      family_random_file,
      write_csv_atomic(
        family_random_coefficients,
        file.path(
          analysis_settings$output$tables,
          "family_random_coefficients.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      family_classification_file,
      write_csv_atomic(
        family_classification,
        file.path(
          analysis_settings$output$tables,
          "family_context_classification.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      family_summary_file,
      write_csv_atomic(
        family_model_summary,
        file.path(
          analysis_settings$output$tables,
          "family_model_diagnostics.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      family_cover_audit_file,
      write_csv_atomic(
        family_plot_auc_bundle$cover_audit,
        file.path(
          analysis_settings$output$audits,
          "family_plot_cover_reconstruction.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      climate_scaling_file,
      write_rds_atomic(
        resident_climate_scaling,
        file.path(
          analysis_settings$output$data,
          "resident_climate_scaling.rds"
        )
      ),
      format = "file"
    ),
    tar_target(
      trait_scaling_file,
      write_csv_atomic(
        trait_scaling,
        file.path(
          analysis_settings$output$tables,
          "trait_scaling.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      functional_species_file,
      write_csv_atomic(
        functional_species,
        file.path(
          analysis_settings$output$data,
          "species_graminoid_woodiness.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      pev_scores_file,
      write_csv_atomic(
        resident_pev_scores_oriented,
        file.path(
          analysis_settings$output$data,
          "phylogeny_predictor_scores.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      pev_scaling_file,
      write_csv_atomic(
        pev_scaling,
        file.path(
          analysis_settings$output$tables,
          "phylogeny_predictor_scaling.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      phylogeny_file,
      write_rds_atomic(
        resident_phylogeny_candidate,
        file.path(
          analysis_settings$output$models,
          "phylogeny_tree.rds"
        )
      ),
      format = "file"
    ),
    tar_target(
      phylogeny_newick_file,
      write_text_atomic(
        ape::write.tree(resident_phylogeny_candidate),
        file.path(
          analysis_settings$output$models,
          "phylogeny_tree.newick"
        )
      ),
      format = "file"
    ),
    tar_target(
      phylogenetic_distance_file,
      write_rds_atomic(
        phylogenetic_distance,
        file.path(
          analysis_settings$output$data,
          "phylogeny_distances.rds"
        )
      ),
      format = "file"
    ),
    tar_target(
      zero_distance_audit_file,
      write_csv_atomic(
        zero_distance_audit,
        file.path(
          analysis_settings$output$audits,
          "zero_distance_pair_audit.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      ldmc_model_file,
      write_rds_atomic(
        ldmc_model_pair$fits$full_REML$fit,
        file.path(
          analysis_settings$output$models,
          "predictor_ldmc_model.rds"
        )
      ),
      format = "file"
    ),
    tar_target(
      lma_model_file,
      write_rds_atomic(
        lma_model_pair$fits$full_REML$fit,
        file.path(
          analysis_settings$output$models,
          "predictor_lma_model.rds"
        )
      ),
      format = "file"
    ),
    tar_target(
      height_model_file,
      write_rds_atomic(
        height_model_pair$fits$full_REML$fit,
        file.path(
          analysis_settings$output$models,
          "predictor_height_model.rds"
        )
      ),
      format = "file"
    ),
    tar_target(
      graminoid_model_file,
      write_rds_atomic(
        graminoid_model_pair$fits$full_REML$fit,
        file.path(
          analysis_settings$output$models,
          "predictor_graminoid_model.rds"
        )
      ),
      format = "file"
    ),
    tar_target(
      woodiness_model_file,
      write_rds_atomic(
        woodiness_model_pair$fits$full_REML$fit,
        file.path(
          analysis_settings$output$models,
          "predictor_woodiness_model.rds"
        )
      ),
      format = "file"
    ),
    tar_target(
      phylogeny_model_file,
      write_rds_atomic(
        phylogeny_model_pair$fits$full_REML$fit,
        file.path(
          analysis_settings$output$models,
          "predictor_phylogeny_model.rds"
        )
      ),
      format = "file"
    ),
    tar_target(
      pcoa_eigen_file,
      write_csv_atomic(
        pcoa_eigen_summary,
        file.path(
          analysis_settings$output$tables,
          "pcoa_eigenvalue_summary.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      predictor_fit_file,
      write_csv_atomic(
        predictor_model_fit,
        file.path(
          analysis_settings$output$tables,
          "predictor_model_fit_statistics.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      predictor_comparison_file,
      write_csv_atomic(
        predictor_model_comparisons,
        file.path(
          analysis_settings$output$tables,
          "predictor_model_comparisons.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      predictor_CR1_file,
      write_csv_atomic(
        predictor_CR1,
        file.path(
          analysis_settings$output$tables,
          "predictor_model_robust_inference.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      predictor_diagnostics_file,
      write_csv_atomic(
        predictor_model_tables$diagnostics,
        file.path(
          analysis_settings$output$audits,
          "predictor_model_random_effects.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      predictor_audit_file,
      write_csv_atomic(
        predictor_model_tables$audit,
        file.path(
          analysis_settings$output$audits,
          "predictor_model_data_audit.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      heterogeneity_components_file,
      write_csv_atomic(
        heterogeneity_components,
        file.path(
          analysis_settings$output$tables,
          "heterogeneity_components.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      heterogeneity_synthesis_file,
      write_csv_atomic(
        heterogeneity_synthesis,
        file.path(
          analysis_settings$output$tables,
          "heterogeneity_synthesis.csv"
        )
      ),
      format = "file"
    ),

    # Temporary Phase 1 assembly: only historical certification depends on it.
    tar_target(
      legacy_comparison_candidate,
      assemble_legacy_comparison_candidate(
        metadata = release_reference_metadata(analysis_settings),
        data_validation = legacy_data_validation_component(
          resident_scope,
          resident_consistency_checks
        ),
        scaling = list(
          climate = resident_climate_scaling,
          traits = trait_scaling,
          PEV = pev_scaling
        ),
        M_both = resident_climate_model_results,
        classification = list(
          species = species_classification_summary,
          family = family_classification_summary
        ),
        family_model = family_model_summary,
        predictor_model_fit = predictor_model_fit,
        predictor_model_comparisons = predictor_model_comparisons,
        predictor_CR1 = predictor_CR1,
        heterogeneity_components = heterogeneity_components,
        heterogeneity_synthesis = heterogeneity_synthesis,
        phylogenetic_positive_variation_percent_PEV1_PEV3 =
          phylogenetic_positive_variation_PEV1_PEV3
      )
    ),
    tar_target(
      frozen_release_reference_file,
      require_existing_file(file.path(
        config$root, "output", "reference",
        "canonical_reference_results.rds"
      )),
      format = "file"
    ),
    tar_target(
      frozen_release_reference_validation,
      validate_frozen_release_reference(
        legacy_comparison_candidate,
        frozen_release_reference_file
      )
    ),
    tar_target(
      legacy_reference_fixture_file,
      file.path(config$root, "tests", "fixtures", "step22_reference_results.rds"),
      format = "file"
    ),
    tar_target(legacy_reference_results, readRDS(legacy_reference_fixture_file)),
    tar_target(
      legacy_reference_validation,
      validate_legacy_reference(legacy_comparison_candidate, legacy_reference_results)
    ),
    tar_target(
      legacy_certification_validation,
      assert_release_validation(assemble_release_validation(list(
        canonical_upstream = community_validation,
        trait_ranges = try_range_validation,
        climate_coordinates = chelsa_coordinate_audit$summary,
        upstream_winner_loser = winner_loser_validation,
        temporal_flows = flow_validation,
        downstream_inputs = resident_input_equivalence,
        downstream_sampling_variance = resident_sampling_variance_validation,
        downstream_consistency = resident_consistency_checks,
        mboth_fit = legacy_resident_climate_model_equivalence,
        mboth_profiles = legacy_resident_climate_variance_profile_equivalence,
        species_classification = species_classification_equivalence,
        family_data = family_auc_data_equivalence,
        family_sampling_variance = family_sampling_variance_equivalence,
        family_cover = family_cover_audit_equivalence,
        family_fit = family_model_equivalence,
        family_support = family_support_equivalence,
        family_random = family_random_equivalence,
        family_classification = family_classification_equivalence,
        family_summary = family_model_summary_equivalence,
        family_bundle = family_bundle_validation,
        trait_species = trait_species_equivalence,
        trait_scaling = trait_scaling_equivalence,
        climate_scaling = climate_scaling_equivalence,
        pev_scaling = pev_scaling_equivalence,
        functional_source = functional_source_coverage_audit,
        functional_species = functional_species_equivalence,
        predictor_specifications = predictor_specification_validation,
        trait_model_fit = trait_model_fit_equivalence,
        trait_model_comparisons = trait_model_comparisons_equivalence,
        trait_model_CR1 = trait_model_CR1_equivalence,
        functional_model_fit = functional_model_fit_equivalence,
        functional_model_comparisons = functional_model_comparisons_equivalence,
        functional_model_CR1 = functional_model_CR1_equivalence,
        phylogeny_reproduction = phylogeny_reproduction_validation,
        phylogeny_migration = phylogeny_migration_release_validation,
        active_downstream_release = current_scientific_release_validation,
        predictive_performance = model_performance_release_validation,
        model_effects = model_effects_release_validation,
        canonical_reference = frozen_release_reference_validation
      )))
    ),
    tar_target(
      legacy_certification_validation_file,
      write_csv_atomic(
        legacy_certification_validation,
        file.path(config$output_dir, "audits", "legacy_certification_validation.csv")
      ),
      format = "file"
    ),
    tar_target(
      legacy_reference_validation_file,
      write_csv_atomic(
        legacy_reference_validation,
        file.path(config$output_dir, "audits", "legacy_reference_validation.csv")
      ),
      format = "file"
    ),
    tar_target(
      run_provenance,
      build_run_provenance(
        config$root, current_scientific_release_validation
      ),
      cue = tar_cue(mode = "always")
    ),
    tar_target(
      run_provenance_file,
      write_csv_atomic(
        run_provenance,
        file.path(config$output_dir, "manifests", "run_provenance.csv")
      ),
      format = "file"
    ),
    tar_target(
      release_artifact_files,
      {
        paths <- c(
        community_file = community_file,
        try_individual_file = try_individual_file,
        auc_by_plot_file = auc_file,
        resident_data_file = resident_data_file,
        resident_sampling_variance_file = resident_sampling_variance_file,
        resident_climate_model_file = resident_climate_model_file,
        resident_climate_fixed_effects_file = resident_climate_fixed_effects_file,
        resident_climate_robust_inference_file = resident_climate_robust_inference_file,
        resident_climate_random_effects_file = resident_climate_random_effects_file,
        resident_climate_variance_profiles_file = resident_climate_variance_profiles_file,
        species_classification_file = species_classification_file,
        family_data_file = family_data_file,
        family_sampling_variance_file = family_sampling_variance_file,
        family_model_file = family_model_file,
        family_support_file = family_support_file,
        family_random_file = family_random_file,
        family_classification_file = family_classification_file,
        family_summary_file = family_summary_file,
        family_cover_audit_file = family_cover_audit_file,
        climate_scaling_file = climate_scaling_file,
        trait_scaling_file = trait_scaling_file,
        functional_species_file = functional_species_file,
        pev_scores_file = pev_scores_file,
        pev_scaling_file = pev_scaling_file,
        phylogeny_file = phylogeny_file,
        phylogeny_newick_file = phylogeny_newick_file,
        phylogenetic_distance_file = phylogenetic_distance_file,
        zero_distance_audit_file = zero_distance_audit_file,
        pcoa_eigen_file = pcoa_eigen_file,
        ldmc_model_file = ldmc_model_file,
        lma_model_file = lma_model_file,
        height_model_file = height_model_file,
        graminoid_model_file = graminoid_model_file,
        woodiness_model_file = woodiness_model_file,
        phylogeny_model_file = phylogeny_model_file,
        predictor_fit_file = predictor_fit_file,
        predictor_comparison_file = predictor_comparison_file,
        predictor_CR1_file = predictor_CR1_file,
        predictor_diagnostics_file = predictor_diagnostics_file,
        predictor_audit_file = predictor_audit_file,
        heterogeneity_components_file = heterogeneity_components_file,
        heterogeneity_synthesis_file = heterogeneity_synthesis_file,
        model_fit_r2_file = model_fit_r2_file,
        model_cross_validation_metrics_file =
          model_cross_validation_metrics_file,
        model_performance_validation_file =
          model_performance_validation_file,
        model_fixed_effects_file = model_fixed_effects_file,
        model_random_effect_profiles_file =
          model_random_effect_profiles_file,
        model_random_effect_correlations_file =
          model_random_effect_correlations_file,
        model_effects_validation_file =
          model_effects_validation_file,
        phylogeny_package_provenance_file = phylogeny_package_provenance_file,
        run_provenance_file = run_provenance_file,
        active_source_inventory_file =
          active_source_inventory_file
        )
        list(
          paths = paths,
          sha256 = stats::setNames(
            vapply(unname(paths), sha256_file, character(1)),
            names(paths)
          )
        )
      }
    ),
    tar_target(
      release_artifact_manifest,
      build_artifact_manifest(
        release_artifact_files$paths,
        config$root,
        paste(
          "current targets; current scientific release checks passed;",
          "gradient-held-out performance validated"
        )
      )
    ),
    tar_target(
      release_artifact_manifest_file,
      write_csv_atomic(
        release_artifact_manifest,
        file.path(config$output_dir, "manifests", "release_artifacts.csv")
      ),
      format = "file"
    )
  )
}
