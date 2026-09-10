canonical_publication_targets <- function() {
  list(
    tar_target(
      figure01_festuca_data,
      prepare_figure01_festuca_data(
        winner_loser,
        species_classification_526,
        auc_inputs$auc,
        analysis_community,
        experiments
      )
    ),
    tar_target(
      figure01_winners_losers,
      build_figure01_winners_losers(figure01_festuca_data)
    ),
    tar_target(
      figure01_winners_losers_pdf,
      write_vector_figure(
        figure01_winners_losers,
        file.path(
          config$root, "figures", "manuscript",
          "fig01_winners_losers.pdf"
        ),
        width = 21,
        height = 14.25
      ),
      format = "file"
    ),
    tar_target(
      figure01_winners_losers_png,
      save_plot_atomic(
        figure01_winners_losers,
        file.path(
          config$root, "figures", "manuscript",
          "fig01_winners_losers.png"
        ),
        width = 21,
        height = 14.25,
        dpi = 300,
        bg = "white"
      ),
      format = "file"
    ),
    tar_target(
      figure02_temporal_strength_correlation_inset,
      build_temporal_strength_correlation_inset(winner_loser)
    ),
    tar_target(
      figure02_temporal_category_flows_file,
      plot_temporal_flows(
        temporal_linked_rows,
        file.path(
          config$root, "figures", "manuscript",
          "fig02_temporal_category_flows.pdf"
        ),
        include_colonizers = FALSE,
        base_size = 22,
        flow_label_size = 5.4,
        stratum_label_size = 5.0,
        axis_text_size = 22,
        legend_text_size = 18,
        legend_position = "none",
        direct_right_labels = TRUE,
        direct_label_size = 5.8,
        width = 18,
        height = 11,
        bg = "white"
      ),
      format = "file"
    ),
    tar_target(
      figure02_temporal_category_flows_png,
      plot_temporal_flows(
        temporal_linked_rows,
        file.path(
          config$root, "figures", "manuscript",
          "fig02_temporal_category_flows.png"
        ),
        include_colonizers = FALSE,
        base_size = 22,
        flow_label_size = 5.4,
        stratum_label_size = 5.0,
        axis_text_size = 22,
        legend_text_size = 18,
        legend_position = "none",
        direct_right_labels = TRUE,
        direct_label_size = 5.8,
        width = 18,
        height = 11,
        bg = "white"
      ),
      format = "file"
    ),
    tar_target(
      appendix_temporal_winning_losing_strength_correlations_file,
      save_plot_atomic(
        figure02_temporal_strength_correlation_inset,
        file.path(
          config$root, "figures", "manuscript",
          "appendix_temporal_winning_losing_strength_correlations.pdf"
        ),
        width = 12,
        height = 5.5,
        bg = "white"
      ),
      format = "file"
    ),
    tar_target(
      appendix_temporal_winning_losing_strength_correlations_png,
      save_plot_atomic(
        figure02_temporal_strength_correlation_inset,
        file.path(
          config$root, "figures", "manuscript",
          "appendix_temporal_winning_losing_strength_correlations.png"
        ),
        width = 12,
        height = 5.5,
        dpi = 300,
        bg = "white"
      ),
      format = "file"
    ),
    tar_target(
      appendix_destination_climate_warming_gradient_file,
      plot_destination_climate_warming_gradient_continuous(
        analysis_data,
        file.path(
          config$root, "figures", "manuscript",
          "appendix_destination_climate_warming_gradient.pdf"
        ),
        base_size = 25,
        experiment_label_size = 5.0,
        axis_title_size = 25,
        axis_text_size = 22,
        legend_title_size = 23,
        legend_text_size = 21,
        width = 18,
        height = 12
      ),
      format = "file"
    ),
    tar_target(
      appendix_destination_climate_warming_gradient_png,
      plot_destination_climate_warming_gradient_continuous(
        analysis_data,
        file.path(
          config$root, "figures", "manuscript",
          "appendix_destination_climate_warming_gradient.png"
        ),
        base_size = 25,
        experiment_label_size = 5.0,
        axis_title_size = 25,
        axis_text_size = 22,
        legend_title_size = 23,
        legend_text_size = 21,
        width = 18,
        height = 12
      ),
      format = "file"
    ),
    tar_target(
      appendix_temporal_transition_matrices_file,
      plot_temporal_transition_matrices(
        temporal_transition_matrix,
        temporal_category_composition_rows,
        file.path(
          config$root, "figures", "manuscript",
          "appendix_temporal_transition_matrices.pdf"
        ),
        include_colonizers = FALSE,
        resident_classified = winner_loser
      ),
      format = "file"
    ),
    tar_target(
      appendix_temporal_transition_matrices_png,
      plot_temporal_transition_matrices(
        temporal_transition_matrix,
        temporal_category_composition_rows,
        file.path(
          config$root, "figures", "manuscript",
          "appendix_temporal_transition_matrices.png"
        ),
        include_colonizers = FALSE,
        resident_classified = winner_loser
      ),
      format = "file"
    ),
    tar_target(
      species_context_figure,
      build_species_context_callout_figure(
        mboth_526_refit$fit,
        downstream_model_data_526,
        species_classification_526,
        species_context_bundle_526$random_coefficients,
        downstream_climate_scaling_526,
        species_callout_examples_526,
        classification_base_size = 14,
        classification_label_size = 4.8,
        inset_base_size = 14
      )
    ),
    tar_target(
      figure04_species_minmax_classification_pdf,
      {
        path <- write_vector_figure(
          species_context_figure,
          file.path(
            config$root, "figures", "manuscript",
            "fig04_species_minmax_classification_526.pdf"
          ),
          width = 16,
          height = 12.5
        )
        writeLines(
          format(Sys.time(), tz = "UTC", usetz = TRUE),
          paste0(path, ".generated")
        )
        path
      },
      format = "file"
    ),
    tar_target(
      figure04_species_minmax_classification_png,
      save_plot_atomic(
        species_context_figure,
        file.path(
          config$root, "figures", "manuscript",
          "fig04_species_minmax_classification_526.png"
        ),
        width = 16,
        height = 12.5,
        dpi = 300,
        bg = "white"
      ),
      format = "file"
    ),
    tar_target(
      retained_figure_files,
      copy_retained_figures_verified(config$root),
      format = "file"
    ),
    tar_target(
      figure03_repeated_species_consistency_pdf,
      plot_signed_consistency(
        descriptive_analysis_data,
        file.path(
          config$root, "figures", "manuscript",
          "fig03_repeated_species_consistency.pdf"
        ),
        include_colonizers = FALSE
      ),
      format = "file"
    ),
    tar_target(
      figure03_repeated_species_consistency_png,
      plot_signed_consistency(
        descriptive_analysis_data,
        file.path(
          config$root, "figures", "manuscript",
          "fig03_repeated_species_consistency.png"
        ),
        include_colonizers = FALSE
      ),
      format = "file"
    ),
    tar_target(
      retained_figure_copy_audit,
      {
        retained_figure_files
        audit_retained_figure_copies(config$root)
      }
    ),
    tar_target(
      retained_figure_validation,
      retained_figure_release_check(retained_figure_copy_audit)
    ),
    tar_target(
      retained_figure_copy_audit_file,
      write_csv_atomic(
        retained_figure_copy_audit,
        file.path(
          config$root, "output", "manifests",
          "retained_figure_copy_audit.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      figure_provenance,
      {
        figure01_winners_losers_pdf
        figure01_winners_losers_png
        figure02_temporal_category_flows_file
        figure02_temporal_category_flows_png
        appendix_temporal_winning_losing_strength_correlations_file
        appendix_temporal_winning_losing_strength_correlations_png
        figure03_repeated_species_consistency_pdf
        figure03_repeated_species_consistency_png
        appendix_destination_climate_warming_gradient_file
        appendix_destination_climate_warming_gradient_png
        appendix_temporal_transition_matrices_file
        appendix_temporal_transition_matrices_png
        figure04_species_minmax_classification_pdf
        figure04_species_minmax_classification_png
        woodiness_figure_file
        retained_figure_files
        build_figure_provenance(config$root)
      },
      cue = tar_cue(mode = "always")
    ),
    tar_target(
      figure_provenance_file,
      write_csv_atomic(
        figure_provenance,
        file.path(config$root, "output", "manifests", "figure_provenance.csv")
      ),
      format = "file"
    ),
    tar_target(
      woodiness_support_mask_specification,
      list(
        status = "PLOTTING_REFINEMENT_IMPLEMENTED",
        needs_group_specific_observed_climate_support_mask = TRUE,
        groups = c("non_woody", "woody"),
        model_refit_permitted = FALSE,
        invalidates = c(
          "woodiness_support_mask", "woodiness_figure",
          "report_13_functional_information", "figure_provenance"
        )
      )
    ),
    tar_target(
      woodiness_support_mask,
      build_group_specific_climate_support(
        functional_auc_data_526, "Woodiness", points = 65L
      )
    ),
    tar_target(
      woodiness_support_mask_audit,
      woodiness_support_mask |>
        dplyr::summarise(
          grid_cells = dplyr::n(),
          supported_cells = sum(.data$supported),
          support_fraction = mean(.data$supported),
          support_method = paste(unique(.data$support_method), collapse = "; "),
          .by = "group"
        )
    ),
    tar_target(
      woodiness_figure,
      build_woodiness_response_surface_figure(
        woodiness_model_sequence_526,
        woodiness_support_mask,
        functional_auc_data_526,
        downstream_climate_scaling_526
      )
    ),
    tar_target(
      woodiness_figure_file,
      write_vector_figure(
        woodiness_figure,
        file.path(
          config$root, "figures", "manuscript",
          "woodiness_response_surfaces_526.pdf"
        ),
        width = 9.5,
        height = 7.2
      ),
      format = "file"
    ),
    tar_target(
      current_downstream_model_report,
      build_current_downstream_model_report(
        canonical_downstream_settings,
        downstream_scope_526,
        downstream_climate_scaling_526,
        mboth_526_refit_bundle,
        species_classification_summary,
        family_classification_summary,
        trait_scaling_526,
        predictor_model_comparisons,
        heterogeneity_synthesis,
        downstream_consistency_checks,
        active_downstream_release_validation
      )
    ),
    tar_target(
      current_downstream_model_results_file,
      write_text_atomic(
        current_downstream_model_report,
        file.path(
          config$root, "reports", "rendered",
          "current_downstream_model_results.txt"
        )
      ),
      format = "file"
    ),
    tar_target(
      release_validation,
      assert_release_validation(assemble_release_validation(list(
        numerical = numerical_release_validation,
        retained_figures = retained_figure_validation
      )))
    ),
    tar_target(
      release_validation_file,
      write_csv_atomic(
        release_validation,
        file.path(config$output_dir, "audits", "release_validation.csv")
      ),
      format = "file"
    ),
    tarchetypes::tar_quarto(
      report_01_inputs_and_scope, "reports/01_inputs_and_scope.qmd",
      output_file = "01_inputs_and_scope.html", working_directory = "."
    ),
    tarchetypes::tar_quarto(
      report_02_community_cleaning, "reports/02_community_cleaning.qmd",
      output_file = "02_community_cleaning.html", working_directory = "."
    ),
    tarchetypes::tar_quarto(
      report_03_taxonomy, "reports/03_taxonomy.qmd",
      output_file = "03_taxonomy.html", working_directory = "."
    ),
    tarchetypes::tar_quarto(
      report_04_climate, "reports/04_climate.qmd",
      output_file = "04_climate.html", working_directory = "."
    ),
    tarchetypes::tar_quarto(
      report_05_trait_preparation, "reports/05_trait_preparation.qmd",
      output_file = "05_trait_preparation.html", working_directory = "."
    ),
    tarchetypes::tar_quarto(
      report_06_auc_validation, "reports/06_auc_validation.qmd",
      output_file = "06_auc_validation.html", working_directory = "."
    ),
    tarchetypes::tar_quarto(
      report_07_temporal_winner_loser, "reports/07_temporal_winner_loser.qmd",
      output_file = "07_temporal_winner_loser.html", working_directory = "."
    ),
    tarchetypes::tar_quarto(
      report_08_downstream_scope_526, "reports/08_downstream_scope_526.qmd",
      output_file = "08_downstream_scope_526.html", working_directory = "."
    ),
    tarchetypes::tar_quarto(
      report_09_mboth_526, "reports/09_mboth_526.qmd",
      output_file = "09_mboth_526.html", working_directory = "."
    ),
    tarchetypes::tar_quarto(
      report_10_species_context, "reports/10_species_context.qmd",
      output_file = "10_species_context.html", working_directory = "."
    ),
    tarchetypes::tar_quarto(
      report_11_family_resolution, "reports/11_family_resolution.qmd",
      output_file = "11_family_resolution.html", working_directory = "."
    ),
    tarchetypes::tar_quarto(
      report_12_trait_information, "reports/12_trait_information.qmd",
      output_file = "12_trait_information.html", working_directory = "."
    ),
    tarchetypes::tar_quarto(
      report_13_functional_information, "reports/13_functional_information.qmd",
      output_file = "13_functional_information.html", working_directory = "."
    ),
    tarchetypes::tar_quarto(
      report_14_phylogeny, "reports/14_phylogeny.qmd",
      output_file = "14_phylogeny.html", working_directory = "."
    ),
    tarchetypes::tar_quarto(
      report_15_heterogeneity_synthesis, "reports/15_heterogeneity_synthesis.qmd",
      output_file = "15_heterogeneity_synthesis.html", working_directory = "."
    ),
    tarchetypes::tar_quarto(
      report_16_reference_and_release, "reports/16_reference_and_release.qmd",
      output_file = "16_reference_and_release.html", working_directory = "."
    ),
    tar_target(
      published_report_01_inputs_and_scope,
      copy_file_atomic_verified(
        report_01_inputs_and_scope,
        file.path(config$root, "reports", "rendered", "01_inputs_and_scope.html")
      ), format = "file"
    ),
    tar_target(
      published_report_02_community_cleaning,
      copy_file_atomic_verified(
        report_02_community_cleaning,
        file.path(config$root, "reports", "rendered", "02_community_cleaning.html")
      ), format = "file"
    ),
    tar_target(
      published_report_03_taxonomy,
      copy_file_atomic_verified(
        report_03_taxonomy,
        file.path(config$root, "reports", "rendered", "03_taxonomy.html")
      ), format = "file"
    ),
    tar_target(
      published_report_04_climate,
      copy_file_atomic_verified(
        report_04_climate,
        file.path(config$root, "reports", "rendered", "04_climate.html")
      ), format = "file"
    ),
    tar_target(
      published_report_05_trait_preparation,
      copy_file_atomic_verified(
        report_05_trait_preparation,
        file.path(config$root, "reports", "rendered", "05_trait_preparation.html")
      ), format = "file"
    ),
    tar_target(
      published_report_06_auc_validation,
      copy_file_atomic_verified(
        report_06_auc_validation,
        file.path(config$root, "reports", "rendered", "06_auc_validation.html")
      ), format = "file"
    ),
    tar_target(
      published_report_07_temporal_winner_loser,
      copy_file_atomic_verified(
        report_07_temporal_winner_loser,
        file.path(config$root, "reports", "rendered", "07_temporal_winner_loser.html")
      ), format = "file"
    ),
    tar_target(
      published_report_08_downstream_scope_526,
      copy_file_atomic_verified(
        report_08_downstream_scope_526,
        file.path(config$root, "reports", "rendered", "08_downstream_scope_526.html")
      ), format = "file"
    ),
    tar_target(
      published_report_09_mboth_526,
      copy_file_atomic_verified(
        report_09_mboth_526,
        file.path(config$root, "reports", "rendered", "09_mboth_526.html")
      ), format = "file"
    ),
    tar_target(
      published_report_10_species_context,
      copy_file_atomic_verified(
        report_10_species_context,
        file.path(config$root, "reports", "rendered", "10_species_context.html")
      ), format = "file"
    ),
    tar_target(
      published_report_11_family_resolution,
      copy_file_atomic_verified(
        report_11_family_resolution,
        file.path(config$root, "reports", "rendered", "11_family_resolution.html")
      ), format = "file"
    ),
    tar_target(
      published_report_12_trait_information,
      copy_file_atomic_verified(
        report_12_trait_information,
        file.path(config$root, "reports", "rendered", "12_trait_information.html")
      ), format = "file"
    ),
    tar_target(
      published_report_13_functional_information,
      copy_file_atomic_verified(
        report_13_functional_information,
        file.path(config$root, "reports", "rendered", "13_functional_information.html")
      ), format = "file"
    ),
    tar_target(
      published_report_14_phylogeny,
      copy_file_atomic_verified(
        report_14_phylogeny,
        file.path(config$root, "reports", "rendered", "14_phylogeny.html")
      ), format = "file"
    ),
    tar_target(
      published_report_15_heterogeneity_synthesis,
      copy_file_atomic_verified(
        report_15_heterogeneity_synthesis,
        file.path(config$root, "reports", "rendered", "15_heterogeneity_synthesis.html")
      ), format = "file"
    ),
    tar_target(
      published_report_16_reference_and_release,
      copy_file_atomic_verified(
        report_16_reference_and_release,
        file.path(config$root, "reports", "rendered", "16_reference_and_release.html")
      ), format = "file"
    )
  )
}
