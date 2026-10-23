# Module: R/targets/definitions/models
# Responsibility: Declare supporting fitted-model and evaluation targets.
# Inputs: QMD model specifications and prepared scientific datasets.
# Outputs: targets definitions for model summaries and diagnostics.
# Side effects: None when sourced; target commands fit models during builds.

#' Model Target Definitions.
#'
#' @return An R object described by the function name; its structure is validated before later stages use it.
model_target_definitions <- function() {
  list(

    # Profile intervals for the fitted resident climate random effects.
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
      species_callout_examples,
      select_species_callout_examples(species_context_bundle$classification)
    ),
    tar_target(
      species_callout_examples_file,
      write_species_callout_examples(
        species_callout_examples,
        file.path(config$output_dir, "tables", "species_context_callouts.csv")
      ),
      format = "file"
    ),

    # The family branch reconstructs family cover and estimates family AUC
    # effects once per experiment. Squared contrast SEs supply the diagonal
    # sampling variances used by the family meta-regression.
    tar_target(
      family_support,
      family_support_table(family_model_bundle$data) |>
        dplyr::mutate(family_id = as.character(.data$family_id))
    ),
    tar_target(
      family_random_bundle,
      extract_family_random_coefficients(family_model_bundle, family_support)
    ),
    tar_target(
      family_classification,
      classify_family_context(
        family_model_bundle, family_support,
        family_random_bundle
      )
    ),
    tar_target(
      family_classification_summary,
      summarize_context_classification(family_classification, "Family")
    ),
    tar_target(
      family_bundle_validation,
      validate_family_bundle(
        family_auc_bundle$data,
        family_auc_bundle$sampling_variance_matrix,
        family_classification,
        family_model_summary
      )
    ),
    tar_target(
      species_family_resolution_comparison,
      compare_classification_resolutions(
        species_context_bundle$classification, family_classification
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
      functional_source_coverage_audit,
      audit_functional_source_coverage(
        functional_groups_input,
        functional_species,
        resident_species_overrides
      )
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


    # Fixed and random effect summaries for model diagnostics.
    tar_target(
      model_fixed_effects,
      assemble_model_fixed_effects(
        resident_climate_model_results$CR1,
        predictor_model_tables$CR1,
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

    # The tree is rebuilt from reviewed placements and the pinned backbone.
    tar_target(
      resident_pev_scores_recomputed,
      resident_pcoa_candidate$scores |>
        dplyr::select("species_id", "PEV1", "PEV2", "PEV3")
    ),
    tar_target(
      phylogenetic_positive_variation_PEV1_PEV3,
      validated_phylogenetic_variation(resident_pcoa_candidate$eigen_summary)
    ),
    tar_target(
      zero_distance_audit,
      build_zero_distance_audit(
        resident_pcoa_candidate$distance,
        phylogeny_placement_decisions
      )
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
      current_scientific_release_validation,
      current_release_checks(
        resident_scope,
        predictor_specification_validation,
        predictor_model_tables,
        family_model_summary,
        resident_phylogeny_rebuild$audit,
        zero_distance_audit,
        heterogeneity_bundle$components,
        heterogeneity_bundle$synthesis
      ),
      description = paste(
        "Current cohort model family phylogeny predictor and heterogeneity",
        "checks"
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
    # Export each declared result from its scientific target.
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
        species_context_bundle$classification,
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
        family_auc_bundle$data,
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
        family_auc_bundle$sampling_variance_matrix,
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
        family_random_bundle$coefficients,
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
        standardized_pev_species$scaling,
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
        resident_phylogeny_rebuild$tree,
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
        ape::write.tree(resident_phylogeny_rebuild$tree),
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
        resident_pcoa_candidate$distance,
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
        resident_pcoa_candidate$eigen_summary,
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
        predictor_model_tables$fit_statistics,
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
        predictor_model_tables$comparisons,
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
        predictor_model_tables$CR1,
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
        heterogeneity_bundle$components,
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
        heterogeneity_bundle$synthesis,
        file.path(
          analysis_settings$output$tables,
          "heterogeneity_synthesis.csv"
        )
      ),
      format = "file"
    ),

    # Record provenance and release artifact paths for the current workflow.
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
