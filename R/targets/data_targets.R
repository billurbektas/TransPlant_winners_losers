# Canonical upstream data, temporal, and descriptive targets. Target names and
# expressions are preserved from the pre-refactor graph.
canonical_data_targets <- function() {
list(
  tar_target(config, transplant_config()),
  tar_target(
    canonical_scientific_source_files_input,
    canonical_scientific_source_files(config$root),
    format = "file"
  ),
  tar_target(
    canonical_scientific_source_inventory,
    build_scientific_source_inventory(
      canonical_scientific_source_files_input, config$root
    )
  ),
  tar_target(
    canonical_scientific_source_inventory_file,
    write_csv_atomic(
      canonical_scientific_source_inventory,
      file.path(
        config$root, "output", "manifests",
        "active_scientific_source_inventory.csv"
      )
    ),
    format = "file"
  ),
  tar_target(
    downstream_species_overrides_input,
    file.path(config$root, "config", "downstream_species_overrides.csv"),
    format = "file"
  ),
  tar_target(
    downstream_species_overrides,
    readr::read_csv(
      downstream_species_overrides_input,
      show_col_types = FALSE,
      col_types = readr::cols(.default = readr::col_character())
    )
  ),
  tar_target(input_rdata_file, config$input_rdata, format = "file"),
  tar_target(current_input, read_current_transplant_input(input_rdata_file)),
  tar_target(input_inventory, current_input_inventory(current_input)),
  tar_target(
    input_inventory_file,
    write_csv_atomic(input_inventory, file.path(config$output_dir, "tables", "input_inventory.csv")),
    format = "file"
  ),
  tar_target(site_metadata, site_metadata_from_regions(current_input$regional)),
  tar_target(species_mapping_decisions_input, file.path(config$root, "config", "species_mapping_decisions.csv"), format = "file"),
  tar_target(
    species_mapping_decisions,
    readr::read_csv(species_mapping_decisions_input, show_col_types = FALSE,
                    col_types = readr::cols(.default = readr::col_character()))
  ),
  tar_target(
    species_mapping_review,
    local_species_mapping_review(current_input$cleanedspecies, species_mapping_decisions)
  ),
  tar_target(
    species_mapping_review_file,
    write_csv_atomic(species_mapping_review, file.path(config$output_dir, "tables", "species_mapping_review.csv")),
    format = "file"
  ),
  tar_target(species_mapping, resolve_local_species_mappings(current_input$cleanedspecies, species_mapping_decisions)),
  tar_target(community, canonical_community(current_input$dat, species_mapping, site_metadata)),
  tar_target(sites, build_sites(community, current_input$regional)),
  tar_target(experiments, build_experiments(community, sites)),
  tar_target(canonical_validation, validate_canonical_snapshot(sites, community, experiments)),
  tar_target(
    sites_file,
    write_csv_atomic(sites, file.path(config$output_dir, "tables", "sites.csv")),
    format = "file"
  ),
  tar_target(
    community_file,
    write_rds_atomic(community, file.path(config$output_dir, "data", "community.rds")),
    format = "file"
  ),
  tar_target(
    experiments_file,
    write_csv_atomic(experiments, file.path(config$output_dir, "tables", "experiments.csv")),
    format = "file"
  ),
  tar_target(
    canonical_validation_file,
    write_csv_atomic(canonical_validation, file.path(config$output_dir, "tables", "canonical_validation.csv")),
    format = "file"
  ),

  tar_target(taxonomy_overrides_input, file.path(config$root, "config", "taxonomy_overrides.csv"), format = "file"),
  tar_target(
    taxonomy_overrides,
    readr::read_csv(taxonomy_overrides_input, show_col_types = FALSE,
                    col_types = readr::cols(.default = readr::col_character()))
  ),
  tar_target(vegetation_species, sort(unique(community$species))),
  tar_target(community_taxonomy_response, taxonomy_response(vegetation_species)),
  tar_target(
    community_taxonomy_response_file,
    write_csv_atomic(community_taxonomy_response, file.path(config$output_dir, "tables", "tnrs_vegetation_full.csv")),
    format = "file"
  ),
  tar_target(
    community_taxonomy_review,
    apply_taxonomy_review(community_taxonomy_response, taxonomy_overrides, config$taxonomy_min_score)
  ),
  tar_target(
    community_taxonomy_review_file,
    write_csv_atomic(community_taxonomy_review, file.path(config$output_dir, "tables", "tnrs_vegetation_review.csv")),
    format = "file"
  ),
  tar_target(community_taxonomy, validate_taxonomy_review(community_taxonomy_review, "vegetation taxonomy")),
  tar_target(vegetation_taxonomic_resolution, taxonomic_resolution_review(community_taxonomy)),
  tar_target(
    vegetation_taxonomic_resolution_file,
    write_csv_atomic(
      vegetation_taxonomic_resolution,
      file.path(config$output_dir, "tables", "vegetation_taxonomic_resolution.csv")
    ),
    format = "file"
  ),
  tar_target(
    omitted_non_species_taxa_file,
    write_csv_atomic(
      dplyr::filter(vegetation_taxonomic_resolution, !.data$include_species_level),
      file.path(config$output_dir, "tables", "omitted_non_species_taxa.csv")
    ),
    format = "file"
  ),
  tar_target(community_taxonomy_species, filter_species_level_taxonomy(vegetation_taxonomic_resolution)),
  tar_target(
    included_vegetation_species_file,
    write_csv_atomic(
      community_taxonomy_species |>
        dplyr::select(dplyr::any_of(c(
          "name_original", "accepted_name_reviewed", "score", "source",
          "identification_level", "resolution_decision"
        ))),
      file.path(config$output_dir, "tables", "included_vegetation_species.csv")
    ),
    format = "file"
  ),
  tar_target(analysis_community, filter_community_to_reviewed_taxa(community, community_taxonomy_species)),

  tar_target(try_export_file, config$try_file, format = "file"),
  tar_target(try_dictionary_input, file.path(config$root, "config", "try_trait_map.csv"), format = "file"),
  tar_target(try_dictionary, readr::read_csv(try_dictionary_input, show_col_types = FALSE)),
  tar_target(
    try_filter_names,
    unique(c(community_taxonomy_species$name_original, community_taxonomy_species$accepted_name_reviewed))
  ),
  tar_target(
    try_cleaned,
    clean_try_export(try_export_file, try_filter_names, try_dictionary, config$error_risk_max)
  ),
  tar_target(try_taxonomy_response, taxonomy_response(try_side_species(try_cleaned))),
  tar_target(
    try_taxonomy_response_file,
    write_csv_atomic(try_taxonomy_response, file.path(config$output_dir, "tables", "tnrs_try_full.csv")),
    format = "file"
  ),
  tar_target(
    try_taxonomy_review,
    apply_taxonomy_review(try_taxonomy_response, taxonomy_overrides, config$taxonomy_min_score)
  ),
  tar_target(
    try_taxonomy_review_file,
    write_csv_atomic(try_taxonomy_review, file.path(config$output_dir, "tables", "tnrs_try_review.csv")),
    format = "file"
  ),
  tar_target(try_taxonomy, validate_taxonomy_review(try_taxonomy_review, "TRY taxonomy")),
  tar_target(
    try_climate_filtered,
    filter_try_climate(try_cleaned, try_taxonomy, config$excluded_try_climates)
  ),
  tar_target(try_range_validation, validate_trait_ranges(try_climate_filtered$observations)),
  tar_target(try_trait_outputs, summarise_try_traits(try_climate_filtered$observations)),
  tar_target(
    try_audit_file,
    write_csv_atomic(try_climate_filtered$audit, file.path(config$output_dir, "tables", "try_filter_audit.csv")),
    format = "file"
  ),
  tar_target(
    try_individual_file,
    write_rds_atomic(try_trait_outputs$long, file.path(config$output_dir, "data", "try_traits_individual.rds")),
    format = "file"
  ),
  tar_target(
    try_summary_file,
    write_csv_atomic(try_trait_outputs$summary, file.path(config$output_dir, "tables", "try_species_traits.csv")),
    format = "file"
  ),
  tar_target(
    try_range_validation_file,
    write_csv_atomic(try_range_validation, file.path(config$output_dir, "tables", "try_range_validation.csv")),
    format = "file"
  ),

  tar_target(chelsa_coordinate_audit, audit_chelsa_coordinates(sites, experiments)),
  tar_target(
    chelsa_coordinate_validation_file,
    write_csv_atomic(
      chelsa_coordinate_audit$summary,
      file.path(config$output_dir, "tables", "chelsa_coordinate_validation.csv")
    ),
    format = "file"
  ),
  tar_target(
    chelsa_duplicate_coordinates_file,
    write_csv_atomic(
      chelsa_coordinate_audit$duplicate_coordinates,
      file.path(config$output_dir, "tables", "chelsa_duplicate_coordinates.csv")
    ),
    format = "file"
  ),
  tar_target(
    chelsa_coordinate_reuse_file,
    write_csv_atomic(
      chelsa_coordinate_audit$reused_coordinates,
      file.path(config$output_dir, "tables", "chelsa_coordinate_reuse_across_experiments.csv")
    ),
    format = "file"
  ),
  tar_target(
    chelsa_site_coordinates,
    {
      chelsa_coordinate_audit
      validate_site_coordinates(sites)
    }
  ),
  tar_target(chelsa_job_list, chelsa_jobs(config$climate_download_years, config$climate_variables), iteration = "list"),
  tar_target(chelsa_job, chelsa_job_list, pattern = map(chelsa_job_list), iteration = "list"),
  tar_target(
    chelsa_monthly_file,
    run_chelsa_job(
      chelsa_job, chelsa_site_coordinates, config$chelsa_cache,
      config$chelsa_retries, config$chelsa_backoff_seconds
    ),
    pattern = map(chelsa_job),
    format = "file"
  ),
  tar_target(chelsa_consolidated, consolidate_chelsa(chelsa_monthly_file, sites, config$climate_download_years)),
  tar_target(site_year_climate, summarise_site_climate(chelsa_consolidated$monthly)),
  tar_target(experiment_climate, summarise_experiment_climate(site_year_climate, experiments)),
  tar_target(
    experiment_climate_coverage,
    experiment_climate |>
      dplyr::select(
        experiment_id, first_year, last_year, years_expected,
        origin_years_available, destination_years_available,
        origin_missing_years, destination_missing_years,
        climate_complete, climate_summary_status, climate_coverage_note
      )
  ),
  tar_target(
    site_year_climate_file,
    write_csv_atomic(site_year_climate, file.path(config$output_dir, "tables", "site_year_climate.csv")),
    format = "file"
  ),
  tar_target(
    experiment_climate_file,
    write_csv_atomic(experiment_climate, file.path(config$output_dir, "tables", "experiment_climate.csv")),
    format = "file"
  ),
  tar_target(
    experiment_climate_coverage_file,
    write_csv_atomic(
      experiment_climate_coverage,
      file.path(config$output_dir, "tables", "experiment_climate_coverage.csv")
    ),
    format = "file"
  ),

  tar_target(auc_inputs, calculate_auc_inputs(analysis_community, experiments, config$horizons)),
  tar_target(auc_contrasts, fit_auc_contrasts(auc_inputs$auc)),
  tar_target(winner_loser, classify_winners_losers(auc_contrasts, auc_inputs$pools)),
  tar_target(winner_loser_validation, validate_winner_loser_outputs(winner_loser, auc_inputs$pools)),
  tar_target(
    legacy_auc_validation,
    compare_legacy_auc(
      winner_loser,
      file.path(config$root, "tests", "fixtures", "legacy_wl_determined.rds"),
      tolerance = 1e-8
    )
  ),
  tar_target(
    winner_loser_file,
    write_csv_atomic(winner_loser, file.path(config$output_dir, "tables", "winner_loser_results.csv")),
    format = "file"
  ),
  tar_target(
    auc_file,
    write_rds_atomic(auc_inputs$auc, file.path(config$output_dir, "data", "auc_by_plot.rds")),
    format = "file"
  ),

  tar_target(temporal_transition_rows, temporal_transitions(winner_loser)),
  tar_target(temporal_linked_rows, temporal_linked_transitions(winner_loser)),
  tar_target(temporal_linked_summary, linked_transition_summary(temporal_linked_rows)),
  tar_target(temporal_category_persistence, summarise_temporal_category_persistence(temporal_linked_rows)),
  tar_target(temporal_transition_matrix, build_temporal_transition_matrix(temporal_linked_rows)),
  tar_target(temporal_category_composition_rows, temporal_category_composition(winner_loser)),
  tar_target(temporal_summary, transition_summary(temporal_transition_rows)),
  tar_target(temporal_kappa, kappa_cluster_bootstrap(temporal_transition_rows, seed = config$seed)),
  tar_target(temporal_strength_models, fit_temporal_strength_models(winner_loser)),
  tar_target(flow_validation, validate_flow_conservation(temporal_transition_rows, temporal_summary$counts)),
  tar_target(
    temporal_transitions_file,
    write_csv_atomic(temporal_summary$counts, file.path(config$output_dir, "tables", "temporal_transitions.csv")),
    format = "file"
  ),
  tar_target(
    temporal_persistence_file,
    write_csv_atomic(temporal_summary$persistence, file.path(config$output_dir, "tables", "temporal_persistence.csv")),
    format = "file"
  ),
  tar_target(
    temporal_state_counts_file,
    write_csv_atomic(temporal_summary$state_counts, file.path(config$output_dir, "tables", "temporal_state_counts.csv")),
    format = "file"
  ),
  tar_target(
    temporal_linked_flows_file,
    write_csv_atomic(temporal_linked_summary, file.path(config$output_dir, "tables", "temporal_linked_flows.csv")),
    format = "file"
  ),
  tar_target(
    temporal_category_persistence_file,
    write_csv_atomic(
      temporal_category_persistence$by_category,
      file.path(config$output_dir, "tables", "temporal_category_persistence.csv")
    ),
    format = "file"
  ),
  tar_target(
    temporal_overall_persistence_file,
    write_csv_atomic(
      temporal_category_persistence$overall,
      file.path(config$output_dir, "tables", "temporal_overall_persistence.csv")
    ),
    format = "file"
  ),
  tar_target(
    temporal_transition_matrix_file,
    write_csv_atomic(
      temporal_transition_matrix,
      file.path(config$output_dir, "tables", "temporal_transition_matrix.csv")
    ),
    format = "file"
  ),
  tar_target(
    temporal_category_composition_file,
    write_csv_atomic(
      temporal_category_composition_rows,
      file.path(config$output_dir, "tables", "temporal_category_composition.csv")
    ),
    format = "file"
  ),
  tar_target(
    temporal_tests_file,
    write_csv_atomic(temporal_kappa, file.path(config$output_dir, "tables", "temporal_kappa.csv")),
    format = "file"
  ),

  tar_target(taxonomy_family, tnrs_family_lookup(community_taxonomy_species$accepted_name_reviewed)),
  tar_target(functional_group_overrides_input, file.path(config$root, "config", "functional_group_overrides.csv"), format = "file"),
  tar_target(
    functional_group_overrides,
    readr::read_csv(functional_group_overrides_input, show_col_types = FALSE,
                    col_types = readr::cols(.default = readr::col_character()))
  ),
  tar_target(
    functional_groups,
    derive_functional_groups(
      taxonomy_family, functional_group_overrides,
      downstream_species_overrides
    )
  ),
  tar_target(
    functional_groups_file,
    write_csv_atomic(functional_groups, file.path(config$output_dir, "tables", "functional_groups.csv")),
    format = "file"
  ),
  tar_target(
    descriptive_analysis_data,
    analysis_dataset(
      winner_loser, community_taxonomy_species, try_trait_outputs$summary,
      functional_groups, experiments
    )
  ),
  tar_target(
    analysis_data,
    join_experiment_climate(descriptive_analysis_data, experiment_climate) |>
      transform_model_predictors()
  ),
  tar_target(climate_context_result, reviewed_five_climate_contexts(analysis_data)),
  tar_target(climate_context_experiments, climate_context_result$assigned),
  tar_target(
    climate_context_gate,
    validate_five_climate_contexts(climate_context_result),
    error = "continue"
  ),
  tar_target(
    climate_context_experiments_file,
    write_csv_atomic(
      climate_context_experiments,
      file.path(config$output_dir, "tables", "climate_context_experiments.csv")
    ),
    format = "file"
  ),
  tar_target(
    climate_context_review_file,
    write_csv_atomic(
      climate_context_result$review,
      file.path(config$output_dir, "tables", "climate_context_review.csv")
    ),
    format = "file"
  ),
  tar_target(
    analysis_data_file,
    write_csv_atomic(analysis_data, file.path(config$output_dir, "tables", "analysis_dataset.csv")),
    format = "file"
  ),
  tar_target(
    descriptive_analysis_data_file,
    write_csv_atomic(
      descriptive_analysis_data,
      file.path(config$output_dir, "tables", "winner_loser_annotated.csv")
    ),
    format = "file"
  ),
  tar_target(
    figure_destination_climate_warming_gradient,
    plot_destination_climate_warming_gradient(
      analysis_data,
      climate_context_gate,
      file.path(
        config$root, "figures", "manuscript",
        "destination_climate_warming_gradient.pdf"
      )
    ),
    format = "file"
  ),
  tar_target(
    figure_sensitivity,
    plot_sensitivity(
      winner_loser,
      file.path(
        config$root, "figures", "manuscript",
        "compare_rel.abund_vs.sqrt.rel.abund.pdf"
      )
    ),
    format = "file"
  ),
  tar_target(figures_manifest, figure_manifest(config$root), cue = tar_cue(mode = "always")),
  tar_target(
    figures_manifest_file,
    write_csv_atomic(figures_manifest, file.path(config$output_dir, "manifests", "figures.csv")),
    format = "file"
  )
)
}
