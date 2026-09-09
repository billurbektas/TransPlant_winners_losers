canonical_analysis_targets <- function() {
  list(
    tar_target(canonical_downstream_settings, canonical_downstream_config(config$root)),
    tar_target(canonical_predictor_formula_registry, canonical_predictor_formulas()),
    tar_target(
      transition_reference_fixture_file,
      file.path(
        config$root, "tests", "fixtures", "canonical_transition_reference.rds"
      ),
      format = "file"
    ),
    tar_target(
      transition_reference,
      readRDS(transition_reference_fixture_file)
    ),
    tar_target(
      reference_artifacts_config_file,
      file.path(config$root, "config", "reference_artifacts.csv"),
      format = "file"
    ),
    tar_target(
      reference_artifact_ledger,
      read_csv_canonical(reference_artifacts_config_file)
    ),
    tar_target(
      reference_file_artifacts,
      file.path(
        config$root,
        reference_artifact_ledger$path[
          reference_artifact_ledger$checksum_method == "file_sha256"
        ]
      ),
      format = "file"
    ),

    # Compute the downstream boundary from active upstream file targets. The
    # explicit CSV round trip is part of the anchored Step 22 interface: it
    # preserves the exact double representation used by the original runner
    # while keeping every read visible in the target graph.
    tar_target(
      canonical_winner_loser_input,
      read_csv_canonical(winner_loser_file)
    ),
    tar_target(
      canonical_analysis_data_input,
      read_csv_canonical(analysis_data_file)
    ),
    tar_target(
      canonical_climate_context_input,
      read_csv_canonical(climate_context_experiments_file)
    ),
    tar_target(
      downstream_auc_unscaled_526,
      build_resident_auc_data_from_objects(
        canonical_winner_loser_input,
        canonical_analysis_data_input,
        canonical_climate_context_input,
        "max"
      )
    ),
    tar_target(
      downstream_scaled_526,
      scale_resident_auc_data(downstream_auc_unscaled_526)
    ),
    tar_target(
      downstream_auc_data,
      downstream_scaled_526$data |>
        dplyr::mutate(TP_z = .data$T_z * .data$P_z)
    ),
    tar_target(downstream_model_data_526, downstream_auc_data),
    tar_target(
      downstream_climate_scaling_526,
      downstream_scaled_526$scaling |>
        dplyr::filter(.data$scaled_field %in% c("T_z", "P_z", "W_z")) |>
        dplyr::mutate(
          n_distinct_experiments = dplyr::n_distinct(
            downstream_auc_data$experiment_id
          )
        )
    ),
    tar_target(
      downstream_sampling_covariance_526,
      reconstruct_auc_covariance(downstream_auc_data, auc_inputs$auc)
    ),
    tar_target(downstream_V_auc, downstream_sampling_covariance_526$V),
    tar_target(
      downstream_auc_data_reference,
      transition_reference$downstream_auc_data
    ),
    tar_target(
      downstream_V_auc_reference,
      transition_reference$downstream_V_auc
    ),
    tar_target(
      downstream_input_equivalence,
      compare_downstream_inputs(
        downstream_auc_data, downstream_V_auc,
        downstream_auc_data_reference, downstream_V_auc_reference
      )
    ),
    tar_target(
      downstream_consistency_checks,
      transition_reference$downstream_consistency_checks
    ),
    tar_target(
      downstream_scope_526,
      validate_downstream_scope_artifacts(
        downstream_auc_data, downstream_V_auc, downstream_consistency_checks,
        canonical_downstream_settings
      )
    ),
    tar_target(
      downstream_covariance_validation,
      downstream_covariance_release_check(downstream_scope_526$V_validation)
    ),

    tar_target(mboth_526_model, transition_reference$mboth_model_signature),
    tar_target(
      mboth_526_bundle,
      transition_reference$mboth_bundle
    ),
    tar_target(
      mboth_526_refit,
      fit_canonical_model(
        downstream_model_data_526,
        downstream_V_auc,
        canonical_downstream_settings$climate_formula,
        method = "REML",
        label = "canonical M_both_526 REML",
        random = canonical_downstream_settings$species_random
      )
    ),
    tar_target(
      mboth_526_refit_equivalence,
      compare_canonical_model_fit(mboth_526_refit, mboth_526_model)
    ),
    tar_target(
      mboth_variance_profile_specification,
      canonical_variance_profile_specifications(),
      iteration = "list"
    ),
    tar_target(
      mboth_variance_profile_component,
      profile_canonical_variance(
        mboth_526_refit, mboth_variance_profile_specification
      ),
      pattern = map(mboth_variance_profile_specification),
      iteration = "list"
    ),
    tar_target(
      mboth_variance_profiles_526,
      dplyr::bind_rows(mboth_variance_profile_component) |>
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
      mboth_variance_profile_equivalence,
      assert_reference_table(
        mboth_variance_profiles_526,
        mboth_526_bundle$variance_profiles,
        "M_both variance profiles", "component"
      )
    ),
    tar_target(
      mboth_526_refit_bundle,
      assemble_mboth_reference_component(
        mboth_526_refit, mboth_variance_profiles_526
      )
    ),

    tar_target(
      species_classification_reference_526,
      transition_reference$species_classification
    ),
    tar_target(
      species_context_bundle_526,
      build_species_context_bundle(
        mboth_526_refit,
        downstream_model_data_526,
        canonical_downstream_settings$climate_formula
      )
    ),
    tar_target(
      species_classification_526,
      species_context_bundle_526$classification
    ),
    tar_target(
      species_classification_equivalence,
      compare_classification_reference(
        species_classification_526,
        species_classification_reference_526,
        "species_id"
      )
    ),
    tar_target(
      species_callout_examples_526,
      select_species_callout_examples(species_classification_526)
    ),
    tar_target(
      species_callout_examples_file,
      write_species_callout_examples(
        species_callout_examples_526,
        file.path(config$output_dir, "tables", "species_minmax_callout_examples_526.csv")
      ),
      format = "file"
    ),
    tar_target(
      species_classification_summary,
      classification_summary_canonical(species_classification_526, "Species")
    ),
    tar_target(
      family_classification_reference_526,
      transition_reference$family$classification
    ),
    tar_target(
      family_model_summary_reference,
      transition_reference$family$model_summary
    ),
    tar_target(
      family_auc_data_reference_526,
      transition_reference$family$auc_data
    ),
    tar_target(
      family_V_reference_container,
      transition_reference$family$V
    ),
    tar_target(
      family_model_reference,
      transition_reference$family$model_signature
    ),
    tar_target(
      family_support_reference_526,
      transition_reference$family$support
    ),
    tar_target(
      family_random_reference_526,
      transition_reference$family$random
    ),
    tar_target(
      family_cover_audit_reference_526,
      transition_reference$family$cover_audit
    ),

    # The family branch reconstructs family cover and sampling covariance from
    # the active upstream file targets. The explicit read boundary mirrors the
    # validated Step 18/22 interface without sourcing either historical runner.
    tar_target(canonical_family_community_input, readRDS(community_file)),
    tar_target(
      canonical_family_experiments_input,
      read_csv_canonical(experiments_file)
    ),
    tar_target(
      canonical_family_override_table,
      canonical_family_overrides(downstream_species_overrides)
    ),
    tar_target(
      family_resident_membership_526,
      build_family_resident_membership(
        canonical_winner_loser_input, canonical_analysis_data_input,
        canonical_family_override_table
      )
    ),
    tar_target(
      family_plot_auc_bundle_526,
      build_family_plot_auc(
        canonical_family_community_input,
        canonical_family_experiments_input,
        family_resident_membership_526
      )
    ),
    tar_target(
      family_auc_bundle_526,
      build_family_auc_covariance(
        family_plot_auc_bundle_526$plot_auc,
        family_resident_membership_526,
        downstream_model_data_526
      )
    ),
    tar_target(family_auc_data_526, family_auc_bundle_526$data),
    tar_target(family_V_auc_526, family_auc_bundle_526$V),
    tar_target(
      family_auc_data_equivalence,
      assert_reference_table(
        family_auc_data_526, family_auc_data_reference_526,
        "family_auc_data", "effect_id"
      )
    ),
    tar_target(
      family_V_equivalence,
      compare_family_covariance(
        family_V_auc_526, family_V_reference_container
      )
    ),
    tar_target(
      family_cover_audit_equivalence,
      assert_reference_table(
        family_plot_auc_bundle_526$cover_audit,
        family_cover_audit_reference_526,
        "family_plot_cover_reconstruction_audit", "experiment_id"
      )
    ),
    tar_target(
      family_model_bundle_526,
      fit_family_aggregated_model(family_auc_bundle_526)
    ),
    tar_target(
      family_model_equivalence,
      compare_canonical_model_fit(
        family_model_bundle_526, family_model_reference
      )
    ),
    tar_target(
      family_support_526,
      family_support_table(family_model_bundle_526$data) |>
        dplyr::mutate(family_id = as.character(.data$family_id))
    ),
    tar_target(
      family_support_equivalence,
      assert_reference_table(
        family_support_526, family_support_reference_526,
        "family_support", "family_id"
      )
    ),
    tar_target(
      family_random_bundle_526,
      family_random_coefficients(family_model_bundle_526, family_support_526)
    ),
    tar_target(
      family_random_coefficients_526,
      family_random_bundle_526$coefficients
    ),
    tar_target(
      family_random_equivalence,
      assert_reference_table(
        family_random_coefficients_526, family_random_reference_526,
        "family_random_coefficients", "family_id"
      )
    ),
    tar_target(
      family_classification_526,
      classify_family_context(
        family_model_bundle_526, family_support_526,
        family_random_bundle_526
      )
    ),
    tar_target(
      family_classification_equivalence,
      compare_classification_reference(
        family_classification_526,
        family_classification_reference_526,
        "family_id"
      )
    ),
    tar_target(
      family_classification_summary,
      classification_summary_canonical(family_classification_526, "Family")
    ),
    tar_target(
      family_model_summary,
      summarize_family_model(family_model_bundle_526)
    ),
    tar_target(
      family_model_summary_equivalence,
      assert_reference_table(
        family_model_summary, family_model_summary_reference,
        "family_model_summary", "model"
      )
    ),
    tar_target(
      family_bundle_validation,
      validate_family_bundle(
        family_auc_data_526,
        family_V_auc_526,
        family_classification_526,
        family_model_summary
      )
    ),
    tar_target(
      species_family_resolution_comparison,
      compare_classification_resolutions(
        species_classification_526, family_classification_526
      )
    ),

    tar_target(
      trait_scaling_reference_526,
      transition_reference$trait_scaling
    ),
    tar_target(
      trait_species_reference_526,
      transition_reference$trait_species
    ),
    tar_target(
      canonical_try_summary_input,
      read_csv_canonical(try_summary_file)
    ),
    tar_target(
      trait_species_526,
      build_canonical_trait_species(
        canonical_try_summary_input,
        downstream_auc_data$species_id
      )
    ),
    tar_target(
      trait_species_equivalence,
      compare_trait_species_reference(
        trait_species_526,
        trait_species_reference_526
      )
    ),
    tar_target(
      prepared_ldmc_526,
      prepare_canonical_trait_data(
        trait_species_526, downstream_auc_data, downstream_V_auc, "LDMC"
      )
    ),
    tar_target(
      prepared_lma_526,
      prepare_canonical_trait_data(
        trait_species_526, downstream_auc_data, downstream_V_auc, "LMA"
      )
    ),
    tar_target(
      prepared_height_526,
      prepare_canonical_trait_data(
        trait_species_526, downstream_auc_data, downstream_V_auc,
        "vegetative_height"
      )
    ),
    tar_target(
      prepared_trait_models_526,
      list(
        LDMC = prepared_ldmc_526,
        LMA = prepared_lma_526,
        vegetative_height = prepared_height_526
      )
    ),
    tar_target(
      trait_scaling_526,
      dplyr::bind_rows(lapply(prepared_trait_models_526, `[[`, "scaling"))
    ),
    tar_target(
      trait_scaling_equivalence,
      assert_reference_table(
        trait_scaling_526, trait_scaling_reference_526,
        "trait scaling", "trait"
      )
    ),
    tar_target(
      climate_scaling_reference_526,
      transition_reference$climate_scaling
    ),
    tar_target(climate_scaling_526, downstream_climate_scaling_526),
    tar_target(
      climate_scaling_equivalence,
      assert_reference_table(
        climate_scaling_526, climate_scaling_reference_526,
        "climate scaling", "scaled_field"
      )
    ),
    tar_target(pev_scaling_reference_526, transition_reference$pev_scaling),
    tar_target(
      standardized_pev_species_526,
      standardize_unique_species_predictors(
        resident_pev_scores_oriented |>
          dplyr::select("species_id", "PEV1", "PEV2", "PEV3"),
        c("PEV1", "PEV2", "PEV3")
      )
    ),
    tar_target(pev_scaling_526, standardized_pev_species_526$scaling),
    tar_target(
      pev_scaling_equivalence,
      assert_reference_table(
        pev_scaling_526, pev_scaling_reference_526,
        "PEV scaling", "PEV"
      )
    ),
    tar_target(
      phylogeny_auc_data_526,
      downstream_auc_data |>
        dplyr::left_join(
          standardized_pev_species_526$data |>
            dplyr::select("species_id", "PEV1_z", "PEV2_z", "PEV3_z"),
          by = "species_id"
        )
    ),
    tar_target(
      functional_species_reference_526,
      transition_reference$functional_species
    ),
    tar_target(
      canonical_functional_groups_input,
      read_csv_canonical(functional_groups_file)
    ),
    tar_target(
      functional_species_526,
      build_canonical_functional_species(
        canonical_functional_groups_input,
        downstream_auc_data$species_id,
        downstream_species_overrides
      )
    ),
    tar_target(
      functional_source_coverage_audit,
      audit_functional_source_coverage(
        canonical_functional_groups_input,
        functional_species_526,
        downstream_species_overrides
      )
    ),
    tar_target(
      functional_species_equivalence,
      compare_functional_species_reference(
        functional_species_526, functional_species_reference_526
      )
    ),
    tar_target(
      functional_auc_data_526,
      prepare_functional_predictors(downstream_auc_data, functional_species_526)
    ),
    tar_target(
      functional_climate_coverage_526,
      dplyr::bind_rows(
        functional_climate_coverage(
          functional_auc_data_526, "GraminoidStatus", "Graminoid status"
        ),
        functional_climate_coverage(
          functional_auc_data_526, "Woodiness", "Woodiness"
        )
      )
    ),
    tar_target(
      ldmc_predictor_specification_526,
      canonical_trait_predictor_specification("LDMC", "LDMC", prepared_ldmc_526)
    ),
    tar_target(
      lma_predictor_specification_526,
      canonical_trait_predictor_specification("LMA", "LMA", prepared_lma_526)
    ),
    tar_target(
      height_predictor_specification_526,
      canonical_trait_predictor_specification(
        "Height", "Height", prepared_height_526
      )
    ),
    tar_target(
      graminoid_predictor_specification_526,
      canonical_functional_predictor_specification(
        "Graminoid", "Graminoid status", functional_auc_data_526,
        downstream_V_auc
      )
    ),
    tar_target(
      woodiness_predictor_specification_526,
      canonical_functional_predictor_specification(
        "Woodiness", "Woodiness", functional_auc_data_526,
        downstream_V_auc
      )
    ),
    tar_target(
      phylogeny_predictor_specification_526,
      canonical_phylogeny_predictor_specification(
        phylogeny_auc_data_526, downstream_V_auc
      )
    ),
    tar_target(
      predictor_specifications_526,
      list(
        LDMC = ldmc_predictor_specification_526,
        LMA = lma_predictor_specification_526,
        Height = height_predictor_specification_526,
        Graminoid = graminoid_predictor_specification_526,
        Woodiness = woodiness_predictor_specification_526,
        Phylogeny = phylogeny_predictor_specification_526
      )
    ),
    tar_target(
      predictor_specification_validation,
      dplyr::bind_rows(lapply(predictor_specifications_526, function(spec) {
        validate_predictor_specification(
          spec,
          canonical_downstream_settings$climate_formula,
          downstream_auc_data$effect_id
        )
      }))
    ),
    tar_target(
      predictor_model_fit_reference,
      transition_reference$predictor$fit
    ),
    tar_target(
      predictor_model_comparisons_reference,
      transition_reference$predictor$comparisons
    ),
    tar_target(predictor_CR1_reference, transition_reference$predictor$CR1),
    tar_target(
      full_526_M0_ML,
      fit_canonical_model(
        downstream_auc_data,
        downstream_V_auc,
        canonical_downstream_settings$climate_formula,
        "ML",
        "Full 526-species M0 ML"
      )
    ),
    tar_target(
      ldmc_model_sequence_526,
      fit_predictor_sequence(
        ldmc_predictor_specification_526,
        canonical_downstream_settings$climate_formula
      )
    ),
    tar_target(
      lma_model_sequence_526,
      fit_predictor_sequence(
        lma_predictor_specification_526,
        canonical_downstream_settings$climate_formula
      )
    ),
    tar_target(
      height_model_sequence_526,
      fit_predictor_sequence(
        height_predictor_specification_526,
        canonical_downstream_settings$climate_formula
      )
    ),
    tar_target(
      graminoid_model_sequence_526,
      fit_predictor_sequence(
        graminoid_predictor_specification_526,
        canonical_downstream_settings$climate_formula,
        full_526_M0_ML,
        mboth_526_refit
      )
    ),
    tar_target(
      woodiness_model_sequence_526,
      fit_predictor_sequence(
        woodiness_predictor_specification_526,
        canonical_downstream_settings$climate_formula,
        full_526_M0_ML,
        mboth_526_refit
      )
    ),
    tar_target(
      phylogeny_model_sequence_526,
      fit_predictor_sequence(
        phylogeny_predictor_specification_526,
        canonical_downstream_settings$climate_formula,
        full_526_M0_ML,
        mboth_526_refit
      )
    ),
    tar_target(
      predictor_sequence_tables_526,
      assemble_predictor_sequence_tables(
        ldmc_model_sequence_526,
        lma_model_sequence_526,
        height_model_sequence_526,
        graminoid_model_sequence_526,
        woodiness_model_sequence_526,
        phylogeny_model_sequence_526
      )
    ),
    # Report-facing assemblies stay branch-specific. This prevents a change to
    # one scientific branch from invalidating unrelated reports.
    tar_target(
      trait_model_tables_526,
      assemble_predictor_sequence_tables(
        ldmc_model_sequence_526,
        lma_model_sequence_526,
        height_model_sequence_526
      )
    ),
    tar_target(
      functional_model_tables_526,
      assemble_predictor_sequence_tables(
        graminoid_model_sequence_526,
        woodiness_model_sequence_526
      )
    ),
    tar_target(
      phylogeny_model_tables_526,
      assemble_predictor_sequence_tables(phylogeny_model_sequence_526)
    ),
    tar_target(
      predictor_model_fit,
      predictor_sequence_tables_526$fit_statistics
    ),
    tar_target(
      predictor_model_comparisons,
      predictor_sequence_tables_526$comparisons
    ),
    tar_target(predictor_CR1, predictor_sequence_tables_526$CR1),
    tar_target(
      predictor_model_fit_equivalence,
      assert_reference_table(
        predictor_model_fit, predictor_model_fit_reference,
        "predictor ML fit statistics"
      )
    ),
    tar_target(
      predictor_model_comparisons_equivalence,
      assert_reference_table(
        predictor_model_comparisons,
        predictor_model_comparisons_reference,
        "predictor ML nested comparisons"
      )
    ),
    tar_target(
      predictor_CR1_equivalence,
      assert_reference_table(
        predictor_CR1, predictor_CR1_reference,
        "predictor full REML CR1"
      )
    ),
    tar_target(
      trait_model_fit_equivalence,
      assert_reference_table(
        trait_model_tables_526$fit_statistics,
        dplyr::filter(
          predictor_model_fit_reference,
          .data$predictor %in% c("LDMC", "LMA", "Height")
        ),
        "trait-model ML fit statistics"
      )
    ),
    tar_target(
      trait_model_comparisons_equivalence,
      assert_reference_table(
        trait_model_tables_526$comparisons,
        dplyr::filter(
          predictor_model_comparisons_reference,
          .data$predictor %in% c("LDMC", "LMA", "Height")
        ),
        "trait-model ML nested comparisons"
      )
    ),
    tar_target(
      trait_model_CR1_equivalence,
      assert_reference_table(
        trait_model_tables_526$CR1,
        dplyr::filter(
          predictor_CR1_reference,
          .data$predictor %in% c("LDMC", "LMA", "Height")
        ),
        "trait-model full REML CR1"
      )
    ),
    tar_target(
      functional_model_fit_equivalence,
      assert_reference_table(
        functional_model_tables_526$fit_statistics,
        dplyr::filter(
          predictor_model_fit_reference,
          .data$predictor %in% c("Graminoid status", "Woodiness")
        ),
        "functional-model ML fit statistics"
      )
    ),
    tar_target(
      functional_model_comparisons_equivalence,
      assert_reference_table(
        functional_model_tables_526$comparisons,
        dplyr::filter(
          predictor_model_comparisons_reference,
          .data$predictor %in% c("Graminoid status", "Woodiness")
        ),
        "functional-model ML nested comparisons"
      )
    ),
    tar_target(
      functional_model_CR1_equivalence,
      assert_reference_table(
        functional_model_tables_526$CR1,
        dplyr::filter(
          predictor_CR1_reference,
          .data$predictor %in% c("Graminoid status", "Woodiness")
        ),
        "functional-model full REML CR1"
      )
    ),
    tar_target(
      phylogeny_model_fit_equivalence,
      assert_reference_table(
        phylogeny_model_tables_526$fit_statistics,
        dplyr::filter(
          predictor_model_fit_reference,
          .data$predictor == "Phylogenetic position (PEV1-PEV3)"
        ),
        "phylogeny-model ML fit statistics"
      )
    ),
    tar_target(
      phylogeny_model_comparisons_equivalence,
      assert_reference_table(
        phylogeny_model_tables_526$comparisons,
        dplyr::filter(
          predictor_model_comparisons_reference,
          .data$predictor == "Phylogenetic position (PEV1-PEV3)"
        ),
        "phylogeny-model ML nested comparisons"
      )
    ),
    tar_target(
      phylogeny_model_CR1_equivalence,
      assert_reference_table(
        phylogeny_model_tables_526$CR1,
        dplyr::filter(
          predictor_CR1_reference,
          .data$predictor == "Phylogenetic position (PEV1-PEV3)"
        ),
        "phylogeny-model full REML CR1"
      )
    ),

    tar_target(
      heterogeneity_components_reference,
      transition_reference$heterogeneity$components
    ),
    tar_target(
      heterogeneity_synthesis_reference,
      transition_reference$heterogeneity$synthesis
    ),
    tar_target(
      heterogeneity_bundle_526,
      assemble_heterogeneity_synthesis(
        ldmc_model_sequence_526,
        lma_model_sequence_526,
        height_model_sequence_526,
        graminoid_model_sequence_526,
        woodiness_model_sequence_526,
        phylogeny_model_sequence_526
      )
    ),
    tar_target(
      heterogeneity_components,
      heterogeneity_bundle_526$components
    ),
    tar_target(
      heterogeneity_synthesis,
      heterogeneity_bundle_526$synthesis
    ),
    tar_target(
      heterogeneity_components_equivalence,
      assert_reference_table(
        heterogeneity_components, heterogeneity_components_reference,
        "response-scale heterogeneity components"
      )
    ),
    tar_target(
      heterogeneity_synthesis_equivalence,
      assert_reference_table(
        heterogeneity_synthesis, heterogeneity_synthesis_reference,
        "heterogeneity synthesis"
      )
    ),

    tar_target(
      pcoa_eigen_summary_reference,
      transition_reference$phylogeny$eigen_summary
    ),
    tar_target(
      resident_phylogeny_reference,
      transition_reference$phylogeny$tree_rds
    ),
    tar_target(
      phylogeny_placement_decisions_file,
      file.path(config$root, "config", "phylogeny_placement_decisions.csv"),
      format = "file"
    ),
    tar_target(
      phylogeny_placement_decisions,
      read_csv_canonical(phylogeny_placement_decisions_file)
    ),
    tar_target(
      resident_phylogeny_rebuild,
      build_phylogeny_from_placement_decisions(
        phylogeny_placement_decisions,
        downstream_model_data_526$species_id
      )
    ),
    tar_target(
      resident_phylogeny_newick,
      transition_reference$phylogeny$tree_newick
    ),
    tar_target(
      phylogeny_species_order,
      transition_reference$phylogeny$species_order
    ),
    tar_target(
      resident_phylogeny_newick_candidate,
      restore_phylogeny_tip_labels(resident_phylogeny_newick, phylogeny_species_order)
    ),
    # The canonical tree is rebuilt from the tracked placement-decision ledger
    # and the pinned V.PhyloMaker2 backbone. The old RDS and Newick remain
    # independent regression authorities until this gate passes.
    tar_target(resident_phylogeny_candidate, resident_phylogeny_rebuild$tree),
    tar_target(
      resident_pev_scores_reference,
      transition_reference$phylogeny$pev_scores
    ),
    tar_target(
      resident_pcoa_candidate,
      recompute_phylogenetic_pcoa(
        resident_phylogeny_candidate,
        resident_pev_scores_reference,
        phylogeny_species_order$species_id
      )
    ),
    tar_target(
      resident_pev_scores_recomputed,
      resident_pcoa_candidate$scores
    ),
    tar_target(
      pev_modeling_interface_526,
      build_anchored_pev_modeling_interface(
        resident_pev_scores_recomputed,
        resident_pev_scores_reference
      )
    ),
    tar_target(
      resident_pev_scores_oriented,
      pev_modeling_interface_526$data
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
      phylogenetic_distance_reference,
      transition_reference$phylogeny$distance
    ),
    tar_target(phylogenetic_distance, resident_pcoa_candidate$distance),
    tar_target(
      zero_distance_audit,
      transition_reference$phylogeny$zero_distance_audit
    ),
    tar_target(
      phylogeny_reference_artifact_validation,
      {
        reference_file_artifacts
        phylogeny_package_installation_files
        verify_reference_artifacts(
          config$root, references = reference_artifact_ledger
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
          resident_phylogeny_candidate, resident_phylogeny_reference,
          phylogeny_species_order, phylogeny_species_order
        ),
        compare_phylogeny(
          resident_phylogeny_newick_candidate, resident_phylogeny_reference,
          phylogeny_species_order, phylogeny_species_order
        ),
        compare_phylogenetic_distance(
          phylogenetic_distance, phylogenetic_distance_reference
        ),
        compare_pev_coordinates(
          resident_pev_scores_recomputed,
          resident_pev_scores_reference
        ),
        compare_pcoa_eigenstructure(
          pcoa_eigen_summary, pcoa_eigen_summary_reference
        ),
        validate_zero_distance_placements(phylogenetic_distance, zero_distance_audit),
        tibble::tibble(
          check = "reference artifact hashes",
          pass = all(phylogeny_reference_artifact_validation$hash_matches),
          detail = paste(nrow(phylogeny_reference_artifact_validation),
                         "configured artifacts")
        )
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
      phylogeny_reference_artifact_validation_file,
      write_csv_atomic(
        phylogeny_reference_artifact_validation,
        file.path(
          config$output_dir, "audits",
          "phylogeny_reference_artifact_validation.csv"
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
      pev_modeling_interface_audit_file,
      write_csv_atomic(
        pev_modeling_interface_526$audit,
        file.path(config$output_dir, "audits", "pev_modeling_interface_audit.csv")
      ),
      format = "file"
    ),

    # Canonical exports are deliberately scoped to their immediate scientific
    # objects. None depends on the late Step 22 validation assembly.
    tar_target(
      canonical_downstream_data_file,
      write_csv_atomic(
        downstream_auc_data,
        file.path(
          canonical_downstream_settings$canonical_output$data,
          "downstream_auc_data_526.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_downstream_V_file,
      write_rds_atomic(
        downstream_V_auc,
        file.path(
          canonical_downstream_settings$canonical_output$data,
          "downstream_V_auc_526.rds"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_mboth_model_file,
      write_rds_atomic(
        mboth_526_refit$fit,
        file.path(
          canonical_downstream_settings$canonical_output$models,
          "M_both_526.rds"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_mboth_fixed_file,
      write_csv_atomic(
        mboth_526_refit_bundle$fixed,
        file.path(
          canonical_downstream_settings$canonical_output$tables,
          "M_both_526_fixed_coefficients.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_mboth_CR1_file,
      write_csv_atomic(
        mboth_526_refit_bundle$CR1,
        file.path(
          canonical_downstream_settings$canonical_output$tables,
          "M_both_526_CR1.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_mboth_random_file,
      write_csv_atomic(
        mboth_526_refit_bundle$random,
        file.path(
          canonical_downstream_settings$canonical_output$tables,
          "M_both_526_random_diagnostics.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_mboth_profiles_file,
      write_csv_atomic(
        mboth_variance_profiles_526,
        file.path(
          canonical_downstream_settings$canonical_output$tables,
          "M_both_526_variance_profile_CIs.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_species_classification_file,
      write_csv_atomic(
        species_classification_526,
        file.path(
          canonical_downstream_settings$canonical_output$tables,
          "species_minmax_classification_526.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_family_data_file,
      write_csv_atomic(
        family_auc_data_526,
        file.path(
          canonical_downstream_settings$canonical_output$data,
          "family_auc_data_526.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_family_V_file,
      write_rds_atomic(
        family_V_auc_526,
        file.path(
          canonical_downstream_settings$canonical_output$data,
          "family_V_auc_526.rds"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_family_model_file,
      write_rds_atomic(
        family_model_bundle_526$fit,
        file.path(
          canonical_downstream_settings$canonical_output$models,
          "M_family_aggregated_526.rds"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_family_support_file,
      write_csv_atomic(
        family_support_526,
        file.path(
          canonical_downstream_settings$canonical_output$tables,
          "family_support_526.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_family_random_file,
      write_csv_atomic(
        family_random_coefficients_526,
        file.path(
          canonical_downstream_settings$canonical_output$tables,
          "family_random_coefficients_526.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_family_classification_file,
      write_csv_atomic(
        family_classification_526,
        file.path(
          canonical_downstream_settings$canonical_output$tables,
          "family_minmax_classification_526.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_family_summary_file,
      write_csv_atomic(
        family_model_summary,
        file.path(
          canonical_downstream_settings$canonical_output$tables,
          "family_model_diagnostics_526.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_family_cover_audit_file,
      write_csv_atomic(
        family_plot_auc_bundle_526$cover_audit,
        file.path(
          canonical_downstream_settings$canonical_output$audits,
          "family_plot_cover_reconstruction_526.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_trait_scaling_file,
      write_csv_atomic(
        trait_scaling_526,
        file.path(
          canonical_downstream_settings$canonical_output$tables,
          "trait_scaling_526.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_functional_species_file,
      write_csv_atomic(
        functional_species_526,
        file.path(
          canonical_downstream_settings$canonical_output$data,
          "species_graminoid_woodiness_526.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_pev_scores_file,
      write_csv_atomic(
        resident_pev_scores_oriented,
        file.path(
          canonical_downstream_settings$canonical_output$data,
          "species_pev_scores_526.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_pev_scaling_file,
      write_csv_atomic(
        pev_scaling_526,
        file.path(
          canonical_downstream_settings$canonical_output$tables,
          "PEV1_PEV3_scaling_526.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_phylogeny_file,
      write_rds_atomic(
        resident_phylogeny_candidate,
        file.path(
          canonical_downstream_settings$canonical_output$models,
          "resident_phylogeny_526.rds"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_phylogeny_newick_file,
      write_text_atomic(
        ape::write.tree(resident_phylogeny_candidate),
        file.path(
          canonical_downstream_settings$canonical_output$models,
          "resident_phylogeny_526.newick"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_phylogenetic_distance_file,
      write_rds_atomic(
        phylogenetic_distance,
        file.path(
          canonical_downstream_settings$canonical_output$data,
          "phylogenetic_distance_matrix_526.rds"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_zero_distance_audit_file,
      write_csv_atomic(
        zero_distance_audit,
        file.path(
          canonical_downstream_settings$canonical_output$audits,
          "zero_distance_pair_audit_526.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_ldmc_model_file,
      write_rds_atomic(
        ldmc_model_sequence_526$fits$full_REML$fit,
        file.path(
          canonical_downstream_settings$canonical_output$models,
          "LDMC_full_REML_526.rds"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_lma_model_file,
      write_rds_atomic(
        lma_model_sequence_526$fits$full_REML$fit,
        file.path(
          canonical_downstream_settings$canonical_output$models,
          "LMA_full_REML_526.rds"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_height_model_file,
      write_rds_atomic(
        height_model_sequence_526$fits$full_REML$fit,
        file.path(
          canonical_downstream_settings$canonical_output$models,
          "Height_full_REML_526.rds"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_graminoid_model_file,
      write_rds_atomic(
        graminoid_model_sequence_526$fits$full_REML$fit,
        file.path(
          canonical_downstream_settings$canonical_output$models,
          "Graminoid_full_REML_526.rds"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_woodiness_model_file,
      write_rds_atomic(
        woodiness_model_sequence_526$fits$full_REML$fit,
        file.path(
          canonical_downstream_settings$canonical_output$models,
          "Woodiness_full_REML_526.rds"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_phylogeny_model_file,
      write_rds_atomic(
        phylogeny_model_sequence_526$fits$full_REML$fit,
        file.path(
          canonical_downstream_settings$canonical_output$models,
          "Phylogeny_full_REML_526.rds"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_pcoa_eigen_file,
      write_csv_atomic(
        pcoa_eigen_summary,
        file.path(
          canonical_downstream_settings$canonical_output$tables,
          "pcoa_eigenvalue_summary_526.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_predictor_fit_file,
      write_csv_atomic(
        predictor_model_fit,
        file.path(
          canonical_downstream_settings$canonical_output$tables,
          "predictor_ML_fit_statistics_526.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_predictor_comparison_file,
      write_csv_atomic(
        predictor_model_comparisons,
        file.path(
          canonical_downstream_settings$canonical_output$tables,
          "predictor_ML_nested_comparisons_526.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_predictor_CR1_file,
      write_csv_atomic(
        predictor_CR1,
        file.path(
          canonical_downstream_settings$canonical_output$tables,
          "predictor_full_REML_CR1_526.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_predictor_diagnostics_file,
      write_csv_atomic(
        predictor_sequence_tables_526$diagnostics,
        file.path(
          canonical_downstream_settings$canonical_output$audits,
          "predictor_random_diagnostics_526.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_predictor_audit_file,
      write_csv_atomic(
        predictor_sequence_tables_526$audit,
        file.path(
          canonical_downstream_settings$canonical_output$audits,
          "predictor_matched_data_audit_526.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_heterogeneity_components_file,
      write_csv_atomic(
        heterogeneity_components,
        file.path(
          canonical_downstream_settings$canonical_output$tables,
          "response_scale_heterogeneity_components_526.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      canonical_heterogeneity_synthesis_file,
      write_csv_atomic(
        heterogeneity_synthesis,
        file.path(
          canonical_downstream_settings$canonical_output$tables,
          "heterogeneity_resolution_synthesis_526.csv"
        )
      ),
      format = "file"
    ),

    # Late assembly only: no ordinary analysis/report/figure depends on this.
    tar_target(
      candidate_step22_results,
      assemble_step22_candidate(
        metadata = validated_step22_metadata(canonical_downstream_settings),
        data_validation = step22_data_validation_component(downstream_scope_526),
        scaling = list(
          climate = downstream_climate_scaling_526,
          traits = trait_scaling_526,
          PEV = pev_scaling_526
        ),
        M_both = mboth_526_refit_bundle,
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
      step22_reference_fixture_file,
      canonical_downstream_settings$fixture_path,
      format = "file"
    ),
    tar_target(step22_reference_results, readRDS(step22_reference_fixture_file)),
    tar_target(
      step22_reference_validation,
      assert_step22_reference(
        validate_step22_reference(candidate_step22_results, step22_reference_results)
      )
    ),
    tar_target(
      numerical_release_validation,
      assert_release_validation(assemble_release_validation(list(
        canonical_upstream = canonical_validation,
        trait_ranges = try_range_validation,
        climate_coordinates = chelsa_coordinate_audit$summary,
        upstream_winner_loser = winner_loser_validation,
        temporal_flows = flow_validation,
        downstream_inputs = downstream_input_equivalence,
        downstream_covariance = downstream_covariance_validation,
        downstream_consistency = downstream_scope_526$consistency_checks,
        mboth_fit = mboth_526_refit_equivalence,
        mboth_profiles = mboth_variance_profile_equivalence,
        species_classification = species_classification_equivalence,
        family_data = family_auc_data_equivalence,
        family_covariance = family_V_equivalence,
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
        phylogeny_model_fit = phylogeny_model_fit_equivalence,
        phylogeny_model_comparisons = phylogeny_model_comparisons_equivalence,
        phylogeny_model_CR1 = phylogeny_model_CR1_equivalence,
        all_predictor_fit = predictor_model_fit_equivalence,
        all_predictor_comparisons = predictor_model_comparisons_equivalence,
        all_predictor_CR1 = predictor_CR1_equivalence,
        heterogeneity_components = heterogeneity_components_equivalence,
        heterogeneity_synthesis = heterogeneity_synthesis_equivalence,
        phylogeny_reproduction = phylogeny_reproduction_validation,
        step22 = step22_reference_validation
      )))
    ),
    tar_target(
      numerical_release_validation_file,
      write_csv_atomic(
        numerical_release_validation,
        file.path(config$output_dir, "audits", "numerical_release_validation.csv")
      ),
      format = "file"
    ),
    tar_target(
      step22_reference_validation_file,
      write_csv_atomic(
        step22_reference_validation,
        file.path(config$output_dir, "audits", "step22_reference_validation.csv")
      ),
      format = "file"
    ),
    tar_target(
      canonical_run_provenance,
      build_canonical_run_provenance(
        config$root, numerical_release_validation
      ),
      cue = tar_cue(mode = "always")
    ),
    tar_target(
      canonical_run_provenance_file,
      write_csv_atomic(
        canonical_run_provenance,
        file.path(config$output_dir, "manifests", "canonical_run_provenance.csv")
      ),
      format = "file"
    ),
    tar_target(
      canonical_numerical_export_files,
      c(
        canonical_downstream_data_file = canonical_downstream_data_file,
        canonical_downstream_V_file = canonical_downstream_V_file,
        canonical_mboth_model_file = canonical_mboth_model_file,
        canonical_mboth_fixed_file = canonical_mboth_fixed_file,
        canonical_mboth_CR1_file = canonical_mboth_CR1_file,
        canonical_mboth_random_file = canonical_mboth_random_file,
        canonical_mboth_profiles_file = canonical_mboth_profiles_file,
        canonical_species_classification_file = canonical_species_classification_file,
        canonical_family_data_file = canonical_family_data_file,
        canonical_family_V_file = canonical_family_V_file,
        canonical_family_model_file = canonical_family_model_file,
        canonical_family_support_file = canonical_family_support_file,
        canonical_family_random_file = canonical_family_random_file,
        canonical_family_classification_file = canonical_family_classification_file,
        canonical_family_summary_file = canonical_family_summary_file,
        canonical_family_cover_audit_file = canonical_family_cover_audit_file,
        canonical_trait_scaling_file = canonical_trait_scaling_file,
        canonical_functional_species_file = canonical_functional_species_file,
        canonical_pev_scores_file = canonical_pev_scores_file,
        canonical_pev_scaling_file = canonical_pev_scaling_file,
        canonical_phylogeny_file = canonical_phylogeny_file,
        canonical_phylogeny_newick_file = canonical_phylogeny_newick_file,
        canonical_phylogenetic_distance_file = canonical_phylogenetic_distance_file,
        canonical_zero_distance_audit_file = canonical_zero_distance_audit_file,
        canonical_pcoa_eigen_file = canonical_pcoa_eigen_file,
        canonical_ldmc_model_file = canonical_ldmc_model_file,
        canonical_lma_model_file = canonical_lma_model_file,
        canonical_height_model_file = canonical_height_model_file,
        canonical_graminoid_model_file = canonical_graminoid_model_file,
        canonical_woodiness_model_file = canonical_woodiness_model_file,
        canonical_phylogeny_model_file = canonical_phylogeny_model_file,
        canonical_predictor_fit_file = canonical_predictor_fit_file,
        canonical_predictor_comparison_file = canonical_predictor_comparison_file,
        canonical_predictor_CR1_file = canonical_predictor_CR1_file,
        canonical_predictor_diagnostics_file = canonical_predictor_diagnostics_file,
        canonical_predictor_audit_file = canonical_predictor_audit_file,
        canonical_heterogeneity_components_file = canonical_heterogeneity_components_file,
        canonical_heterogeneity_synthesis_file = canonical_heterogeneity_synthesis_file,
        phylogeny_package_provenance_file = phylogeny_package_provenance_file,
        phylogeny_reference_artifact_validation_file = phylogeny_reference_artifact_validation_file,
        phylogeny_reproduction_validation_file = phylogeny_reproduction_validation_file,
        pev_modeling_interface_audit_file = pev_modeling_interface_audit_file,
        numerical_release_validation_file = numerical_release_validation_file,
        step22_reference_validation_file = step22_reference_validation_file,
        canonical_run_provenance_file = canonical_run_provenance_file,
        canonical_scientific_source_inventory_file =
          canonical_scientific_source_inventory_file
      )
    ),
    tar_target(
      canonical_output_manifest,
      build_artifact_manifest(
        canonical_numerical_export_files,
        config$root,
        "canonical targets; Step 22 numerically validated"
      )
    ),
    tar_target(
      canonical_output_manifest_file,
      write_csv_atomic(
        canonical_output_manifest,
        file.path(config$output_dir, "manifests", "canonical_outputs.csv")
      ),
      format = "file"
    )
  )
}
