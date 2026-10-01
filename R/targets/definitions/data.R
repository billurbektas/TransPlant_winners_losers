# Module: R/targets/definitions/data
# Responsibility: Declare supporting file and data targets not defined in the QMD.
# Inputs: Configuration, current inputs, and QMD scientific targets.
# Outputs: targets definitions for scientific data and audits.
# Side effects: None when sourced; target commands write declared artifacts during builds.

# File exports and support targets for community, climate, and response data.
#' Data Target Definitions.
#'
#' @return An R object described by the function name; its structure is validated before later stages use it.
data_target_definitions <- function() {
list(
  tar_target(
    active_source_files_input,
    active_scientific_source_files(config$root),
    format = "file"
  ),
  tar_target(
    active_source_inventory,
    build_scientific_source_inventory(
      active_source_files_input, config$root
    )
  ),
  tar_target(
    active_source_inventory_file,
    write_csv_atomic(
      active_source_inventory,
      file.path(
        config$root, "output", "manifests",
        "active_scientific_source_inventory.csv"
      )
    ),
    format = "file"
  ),
  tar_target(input_inventory, current_input_inventory(current_input)),
  tar_target(
    input_inventory_file,
    write_csv_atomic(input_inventory, file.path(config$output_dir, "tables", "input_inventory.csv")),
    format = "file"
  ),
  tar_target(
    species_mapping_review_file,
    write_csv_atomic(species_mapping_review, file.path(config$output_dir, "tables", "species_mapping_review.csv")),
    format = "file"
  ),
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
    community_validation_file,
    write_csv_atomic(community_validation, file.path(config$output_dir, "tables", "community_validation.csv")),
    format = "file"
  ),

  tar_target(
    community_taxonomy_response_file,
    write_csv_atomic(community_taxonomy_response, file.path(config$output_dir, "tables", "tnrs_vegetation_full.csv")),
    format = "file"
  ),
  tar_target(
    community_taxonomy_review_file,
    write_csv_atomic(community_taxonomy_review, file.path(config$output_dir, "tables", "tnrs_vegetation_review.csv")),
    format = "file"
  ),
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

  tar_target(
    try_taxonomy_response_file,
    write_csv_atomic(try_taxonomy_response, file.path(config$output_dir, "tables", "tnrs_try_full.csv")),
    format = "file"
  ),
  tar_target(
    try_taxonomy_review_file,
    write_csv_atomic(try_taxonomy_review, file.path(config$output_dir, "tables", "tnrs_try_review.csv")),
    format = "file"
  ),
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

  tar_target(
    functional_groups_file,
    write_csv_atomic(functional_groups, file.path(config$output_dir, "tables", "functional_groups.csv")),
    format = "file"
  ),
  tar_target(
    climate_context_gate,
    validate_five_climate_contexts(climate_context_result),
    error = "continue"
  ),
  tar_target(
    climate_context_experiments_file,
    write_csv_atomic(
      climate_context_result$assigned,
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
        config$root, "figures", "candidates",
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
