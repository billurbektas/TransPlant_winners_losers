# Module: R/targets/definitions/publication
# Responsibility: Implement the named scientific workflow stage or shared infrastructure.
# Inputs: Function arguments and explicitly declared target dependencies.
# Outputs: Validated in-memory objects and declared workflow artifacts.
# Side effects: Limited to functions whose contracts explicitly write, copy, or render files.

#' Publication Target Definitions.
#'
#' @return An R object described by the function name; its structure is validated before later stages use it.
publication_target_definitions <- function() {
  list(
    tar_target(
      figure01_festuca_data,
      prepare_figure01_festuca_data(
        winner_loser,
        species_classification,
        auc_inputs$auc,
        analysis_community,
        experiments
      )
    ),
    tar_target(
      figure01_winners_losers,
      build_figure01_winners_losers(figure01_festuca_data),
      description = paste(
        "Main winner-loser figure showing the focal Festuca example and",
        "network-wide response classification"
      )
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
      figure03_temporal_strength_correlation_inset,
      build_temporal_strength_correlation_inset(winner_loser)
    ),
    tar_target(
      figure03_temporal_category_flows_file,
      plot_temporal_flows(
        temporal_linked_rows,
        file.path(
          config$root, "figures", "manuscript",
          "fig03_temporal_category_flows.pdf"
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
      figure03_temporal_category_flows_png,
      plot_temporal_flows(
        temporal_linked_rows,
        file.path(
          config$root, "figures", "manuscript",
          "fig03_temporal_category_flows.png"
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
        figure03_temporal_strength_correlation_inset,
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
        figure03_temporal_strength_correlation_inset,
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
        resident_climate_model$fit,
        resident_model_data,
        species_classification,
        species_context_bundle$random_coefficients,
        resident_climate_scaling,
        species_callout_examples,
        classification_base_size = 14,
        classification_label_size = 4.8,
        inset_base_size = 14
      ),
      description = paste(
        "Species winner loser and context-dependent classifications with",
        "selected fitted climate-response examples"
      )
    ),
    tar_target(
      figure04_species_minmax_classification_pdf,
      {
        path <- write_vector_figure(
          species_context_figure,
          file.path(
            config$root, "figures", "manuscript",
            "fig04_species_minmax_classification.pdf"
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
          "fig04_species_minmax_classification.png"
        ),
        width = 16,
        height = 12.5,
        dpi = 300,
        bg = "white"
      ),
      format = "file"
    ),
    tar_target(
      figure05_climate_only_fixed_effect_renderer,
      file.path(
        config$root, "scripts",
        "render_figure.R"
      ),
      format = "file"
    ),
    tar_target(
      figure05_climate_only_fixed_effect_files,
      {
        invisible(list(
          resident_climate_model_file,
          resident_climate_fixed_effects_file,
          climate_scaling_file,
          experiments_file,
          sites_file
        ))
        output_stem <- file.path(
          config$root, "figures", "manuscript",
          "fig05_climate_only_fixed_effect_prediction_surfaces"
        )
        status <- system2(
          file.path(R.home("bin"), "Rscript"),
          args = c(
            shQuote(figure05_climate_only_fixed_effect_renderer),
            "climate-response",
            shQuote(output_stem),
            shQuote(config$output_dir)
          )
        )
        if (!identical(status, 0L)) {
          stop("Figure 5 renderer failed with status ", status, ".")
        }
        output_name <- basename(output_stem)
        paths <- c(
          paste0(output_stem, c(".pdf", ".png")),
          file.path(
            config$output_dir, "data",
            paste0(output_name, "_data.csv")
          ),
          file.path(
            config$output_dir, "tables",
            paste0(output_name, "_quantiles.csv")
          ),
          file.path(
            config$output_dir, "audits",
            paste0(output_name, "_validation.csv")
          )
        )
        if (any(!file.exists(paths))) {
          stop("Figure 5 renderer did not create every declared output.")
        }
        paths
      },
      format = "file"
    ),
    tar_target(
      figure05_climate_only_fixed_effect_pdf,
      figure05_climate_only_fixed_effect_files[
        grepl("[.]pdf$", figure05_climate_only_fixed_effect_files)
      ],
      format = "file"
    ),
    tar_target(
      figure05_climate_only_fixed_effect_png,
      figure05_climate_only_fixed_effect_files[
        grepl("[.]png$", figure05_climate_only_fixed_effect_files)
      ],
      format = "file"
    ),
    tar_target(
      figure06_predictor_synthesis_data,
      prepare_predictor_synthesis(
        predictor_fit_file,
        predictor_comparison_file,
        heterogeneity_synthesis_file
      )
    ),
    tar_target(
      figure06_predictor_synthesis_data_file,
      write_csv_atomic(
        figure06_predictor_synthesis_data,
        file.path(
          analysis_settings$output$tables,
          "predictor_synthesis_support_and_heterogeneity_data.csv"
        )
      ),
      format = "file"
    ),
    tar_target(
      figure06_predictor_synthesis,
      {
        figure06_predictor_synthesis_data_file
        build_predictor_synthesis_figure(
          figure06_predictor_synthesis_data
        )
      },
      description = paste(
        "Direct climate-only versus full predictor comparisons and changes",
        "in fitted response heterogeneity"
      )
    ),
    tar_target(
      figure06_predictor_synthesis_pdf,
      {
        path <- write_vector_figure(
          figure06_predictor_synthesis,
          file.path(
            config$root, "figures", "manuscript",
            "fig06_predictor_synthesis_support_and_heterogeneity.pdf"
          ),
          width = 16,
          height = 8.5
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
      figure06_predictor_synthesis_png,
      save_plot_atomic(
        figure06_predictor_synthesis,
        file.path(
          config$root, "figures", "manuscript",
          "fig06_predictor_synthesis_support_and_heterogeneity.png"
        ),
        width = 16,
        height = 8.5,
        dpi = 300,
        bg = "white"
      ),
      format = "file"
    ),
    tar_target(
      appendix_model_predictive_performance,
      {
        model_performance_release_validation
        build_model_performance_figure(
          model_fit_r2,
          model_cross_validation_metrics
        )
      },
      description = paste(
        "Appendix Figure S7 with equal-gradient marginal and conditional R2",
        "and leave-gradient-out predictive R2 and RMSE"
      )
    ),
    tar_target(
      appendix_model_predictive_performance_pdf,
      {
        path <- write_vector_figure(
          appendix_model_predictive_performance,
          file.path(
            config$root, "figures", "manuscript",
            "appendix_model_predictive_performance.pdf"
          ),
          width = 18,
          height = 12
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
      appendix_model_predictive_performance_png,
      save_plot_atomic(
        appendix_model_predictive_performance,
        file.path(
          config$root, "figures", "manuscript",
          "appendix_model_predictive_performance.png"
        ),
        width = 18,
        height = 12,
        dpi = 300,
        bg = "white"
      ),
      format = "file"
    ),
    tar_target(
      appendix_model_fixed_and_random_effects,
      {
        model_effects_release_validation
        build_model_effects_figure(
          model_fixed_effects,
          model_random_effect_profiles
        )
      }
    ),
    tar_target(
      appendix_model_fixed_and_random_effects_pdf,
      {
        path <- write_vector_figure(
          appendix_model_fixed_and_random_effects,
          file.path(
            config$root, "figures", "manuscript",
            "appendix_model_fixed_and_random_effects.pdf"
          ),
          width = 20,
          height = 28
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
      appendix_model_fixed_and_random_effects_png,
      save_plot_atomic(
        appendix_model_fixed_and_random_effects,
        file.path(
          config$root, "figures", "manuscript",
          "appendix_model_fixed_and_random_effects.png"
        ),
        width = 20,
        height = 28,
        dpi = 300,
        bg = "white"
      ),
      format = "file"
    ),
    tar_target(
      figure07_predictor_response_atlas_renderer,
      file.path(
        config$root, "scripts", "render_figure.R"
      ),
      format = "file"
    ),
    tar_target(
      figure07_predictor_response_atlas_files,
      {
        invisible(list(
          ldmc_model_file,
          lma_model_file,
          height_model_file,
          graminoid_model_file,
          woodiness_model_file,
          phylogeny_model_file,
          trait_scaling_file,
          pev_scaling_file,
          pev_scores_file,
          climate_scaling_file
        ))
        output_stem <- file.path(
          config$root, "figures", "manuscript",
          "fig07_predictor_response_atlas"
        )
        status <- system2(
          file.path(R.home("bin"), "Rscript"),
          args = c(
            shQuote(figure07_predictor_response_atlas_renderer),
            "predictor-atlas",
            shQuote(output_stem),
            shQuote(climate_scaling_file),
            shQuote(config$output_dir)
          )
        )
        if (!identical(status, 0L)) {
          stop("Figure 7 renderer failed with status ", status, ".")
        }
        output_name <- basename(output_stem)
        paths <- c(
          paste0(output_stem, c(".pdf", ".png", "_caption.txt")),
          file.path(
            config$output_dir, "data",
            paste0(output_name, "_data.csv")
          ),
          file.path(
            config$output_dir, "tables",
            paste0(output_name, "_quantiles.csv")
          ),
          file.path(
            config$output_dir, "audits",
            paste0(output_name, "_validation.csv")
          )
        )
        if (any(!file.exists(paths))) {
          stop("Figure 7 renderer did not create every declared output.")
        }
        paths
      },
      format = "file"
    ),
    tar_target(
      figure07_predictor_response_atlas_pdf,
      figure07_predictor_response_atlas_files[
        grepl("[.]pdf$", figure07_predictor_response_atlas_files)
      ],
      format = "file"
    ),
    tar_target(
      figure07_predictor_response_atlas_png,
      figure07_predictor_response_atlas_files[
        grepl("[.]png$", figure07_predictor_response_atlas_files)
      ],
      format = "file"
    ),
    tar_target(
      appendix_family_minmax_renderer,
      file.path(
        config$root, "scripts",
        "render_figure.R"
      ),
      format = "file"
    ),
    tar_target(
      appendix_family_minmax_files,
      {
        invisible(list(
          family_model_file,
          family_data_file,
          family_classification_file,
          family_random_file,
          climate_scaling_file
        ))
        output_stem <- file.path(
          config$root, "figures", "manuscript",
          "appendix_family_minmax_classification"
        )
        status <- system2(
          file.path(R.home("bin"), "Rscript"),
          args = c(
            shQuote(appendix_family_minmax_renderer),
            "family-context",
            shQuote(output_stem),
            shQuote(climate_scaling_file),
            shQuote(config$output_dir)
          )
        )
        if (!identical(status, 0L)) {
          stop("Family appendix renderer failed with status ", status, ".")
        }
        paths <- c(
          paste0(output_stem, c(".pdf", ".png")),
          file.path(
            config$output_dir, "tables",
            "appendix_family_minmax_classification_selected_families.csv"
          ),
          file.path(
            config$output_dir, "audits",
            "appendix_family_minmax_classification_validation.csv"
          )
        )
        if (any(!file.exists(paths))) {
          stop("Family appendix renderer did not create every declared output.")
        }
        paths
      },
      format = "file"
    ),
    tar_target(
      appendix_family_minmax_pdf,
      appendix_family_minmax_files[
        grepl("[.]pdf$", appendix_family_minmax_files)
      ],
      format = "file"
    ),
    tar_target(
      appendix_family_minmax_png,
      appendix_family_minmax_files[
        grepl("[.]png$", appendix_family_minmax_files)
      ],
      format = "file"
    ),
    tar_target(
      appendix_pev2_pev3_atlas_renderer,
      file.path(
        config$root, "scripts",
        "render_figure.R"
      ),
      format = "file"
    ),
    tar_target(
      appendix_pev2_pev3_atlas_files,
      {
        invisible(list(
          phylogeny_model_file,
          pev_scaling_file,
          pev_scores_file,
          climate_scaling_file,
          figure07_predictor_response_atlas_files
        ))
        output_stem <- file.path(
          config$root, "figures", "manuscript",
          "appendix_pev2_pev3_prediction_atlas"
        )
        figure07_data_file <- figure07_predictor_response_atlas_files[
          grepl("_data[.]csv$", figure07_predictor_response_atlas_files)
        ]
        status <- system2(
          file.path(R.home("bin"), "Rscript"),
          args = c(
            shQuote(appendix_pev2_pev3_atlas_renderer),
            "phylogeny-axes",
            shQuote(output_stem),
            shQuote(climate_scaling_file),
            shQuote(figure07_data_file),
            shQuote(config$output_dir)
          )
        )
        if (!identical(status, 0L)) {
          stop("PEV2/PEV3 appendix renderer failed with status ", status, ".")
        }
        paths <- c(
          paste0(output_stem, c(".pdf", ".png", "_caption.txt")),
          file.path(
            config$output_dir, "data",
            "appendix_pev2_pev3_prediction_atlas_data.csv"
          ),
          file.path(
            config$output_dir, "tables",
            "appendix_pev2_pev3_prediction_atlas_quantiles.csv"
          ),
          file.path(
            config$output_dir, "audits",
            "appendix_pev2_pev3_prediction_atlas_validation.csv"
          )
        )
        if (any(!file.exists(paths))) {
          stop("PEV2/PEV3 appendix renderer did not create every declared output.")
        }
        paths
      },
      format = "file"
    ),
    tar_target(
      appendix_pev2_pev3_atlas_pdf,
      appendix_pev2_pev3_atlas_files[
        grepl("[.]pdf$", appendix_pev2_pev3_atlas_files)
      ],
      format = "file"
    ),
    tar_target(
      appendix_pev2_pev3_atlas_png,
      appendix_pev2_pev3_atlas_files[
        grepl("[.]png$", appendix_pev2_pev3_atlas_files)
      ],
      format = "file"
    ),
    tar_target(
      appendix_phylogeny_pcoa_summary_renderer,
      file.path(
        config$root, "scripts",
        "render_figure.R"
      ),
      format = "file"
    ),
    tar_target(
      appendix_phylogeny_pcoa_summary_files,
      {
        invisible(list(
          pev_scores_file,
          pcoa_eigen_file
        ))
        output_stem <- file.path(
          config$root, "figures", "manuscript",
          "appendix_phylogeny_pcoa_summary"
        )
        status <- system2(
          file.path(R.home("bin"), "Rscript"),
          args = c(
            shQuote(appendix_phylogeny_pcoa_summary_renderer),
            "phylogeny-summary",
            shQuote(output_stem),
            shQuote(config$output_dir)
          )
        )
        if (!identical(status, 0L)) {
          stop("Phylogeny PCoA appendix renderer failed with status ", status, ".")
        }
        paths <- c(
          paste0(output_stem, c(".pdf", ".png")),
          file.path(
            config$output_dir, "data",
            "appendix_phylogeny_pcoa_summary_data.csv"
          ),
          file.path(
            config$output_dir, "tables",
            "appendix_phylogeny_pcoa_first10_axes.csv"
          ),
          file.path(
            config$output_dir, "audits",
            "appendix_phylogeny_pcoa_summary_validation.csv"
          )
        )
        if (any(!file.exists(paths))) {
          stop("Phylogeny PCoA appendix renderer did not create every declared output.")
        }
        paths
      },
      format = "file"
    ),
    tar_target(
      appendix_phylogeny_pcoa_summary_pdf,
      appendix_phylogeny_pcoa_summary_files[
        grepl("[.]pdf$", appendix_phylogeny_pcoa_summary_files)
      ],
      format = "file"
    ),
    tar_target(
      appendix_phylogeny_pcoa_summary_png,
      appendix_phylogeny_pcoa_summary_files[
        grepl("[.]png$", appendix_phylogeny_pcoa_summary_files)
      ],
      format = "file"
    ),
    tar_target(
      retained_figure_files,
      copy_retained_figures_verified(config$root),
      format = "file"
    ),
    tar_target(
      figure02_repeated_species_consistency_pdf,
      plot_signed_consistency(
        descriptive_analysis_data,
        file.path(
          config$root, "figures", "manuscript",
          "fig02_repeated_species_consistency.pdf"
        ),
        include_colonizers = FALSE
      ),
      format = "file"
    ),
    tar_target(
      figure02_repeated_species_consistency_png,
      plot_signed_consistency(
        descriptive_analysis_data,
        file.path(
          config$root, "figures", "manuscript",
          "fig02_repeated_species_consistency.png"
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
        figure02_repeated_species_consistency_pdf
        figure02_repeated_species_consistency_png
        appendix_temporal_winning_losing_strength_correlations_file
        appendix_temporal_winning_losing_strength_correlations_png
        figure03_temporal_category_flows_file
        figure03_temporal_category_flows_png
        appendix_destination_climate_warming_gradient_file
        appendix_destination_climate_warming_gradient_png
        appendix_temporal_transition_matrices_file
        appendix_temporal_transition_matrices_png
        figure04_species_minmax_classification_pdf
        figure04_species_minmax_classification_png
        figure05_climate_only_fixed_effect_pdf
        figure05_climate_only_fixed_effect_png
        figure06_predictor_synthesis_pdf
        figure06_predictor_synthesis_png
        appendix_model_predictive_performance_pdf
        appendix_model_predictive_performance_png
        appendix_model_fixed_and_random_effects_pdf
        appendix_model_fixed_and_random_effects_png
        figure07_predictor_response_atlas_pdf
        figure07_predictor_response_atlas_png
        appendix_family_minmax_pdf
        appendix_family_minmax_png
        appendix_pev2_pev3_atlas_pdf
        appendix_pev2_pev3_atlas_png
        appendix_phylogeny_pcoa_summary_pdf
        appendix_phylogeny_pcoa_summary_png
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
          "woodiness_support_mask", "report_predictor_models"
        )
      )
    ),
    tar_target(
      woodiness_support_mask,
      build_group_specific_climate_support(
        functional_auc_data, "Woodiness", points = 65L
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
      current_model_report,
      build_current_model_report(
        analysis_settings,
        resident_scope,
        resident_climate_scaling,
        resident_climate_model_results,
        species_classification_summary,
        family_classification_summary,
        trait_scaling,
        predictor_model_comparisons,
        heterogeneity_synthesis,
        current_scientific_release_validation
      )
    ),
    tar_target(
      current_model_results_file,
      write_text_atomic(
        current_model_report,
        file.path(
          config$root, "reports", "rendered",
          "current_model_results.txt"
        )
      ),
      format = "file"
    ),
    tar_target(
      release_validation,
      assert_release_validation(assemble_release_validation(list(
        scientific = current_scientific_release_validation,
        predictive_performance = model_performance_release_validation,
        model_effects = model_effects_release_validation,
        retained_figures = retained_figure_validation
      ))),
      description = paste(
        "Final current-workflow checks for scientific results Figure S7",
        "model diagnostics and registered retained figures"
      )
    ),
    tar_target(
      release_validation_file,
      write_csv_atomic(
        release_validation,
        file.path(config$output_dir, "audits", "release_validation.csv")
      ),
      format = "file"
    ),
    tar_target(
      refactor_certification,
      assert_release_validation(assemble_release_validation(list(
        current_release = release_validation,
        historical_equivalence = legacy_certification_validation,
        historical_fixture_equivalence = legacy_reference_validation
      )))
    ),
    tarchetypes::tar_quarto(
      report_analysis_workflow,
      "analysis_workflow.qmd",
      output_file = "analysis_workflow.html",
      working_directory = ".",
      description = paste(
        "Linear reviewer-facing analysis with cached data checks, explicit",
        "model specifications, figures, and targets inspection"
      )
    ),
    tarchetypes::tar_quarto(
      report_data_foundation, "reports/data_foundation.qmd",
      output_file = "data_foundation.html", working_directory = "."
    ),
    tarchetypes::tar_quarto(
      report_response_estimation, "reports/response_estimation.qmd",
      output_file = "response_estimation.html", working_directory = "."
    ),
    tarchetypes::tar_quarto(
      report_resident_climate_model, "reports/resident_climate_model.qmd",
      output_file = "resident_climate_model.html", working_directory = "."
    ),
    tarchetypes::tar_quarto(
      report_family_context_model, "reports/family_context_model.qmd",
      output_file = "family_context_model.html", working_directory = "."
    ),
    tarchetypes::tar_quarto(
      report_predictor_models, "reports/predictor_models.qmd",
      output_file = "predictor_models.html", working_directory = "."
    ),
    tarchetypes::tar_quarto(
      report_model_synthesis, "reports/model_synthesis.qmd",
      output_file = "model_synthesis.html", working_directory = ".",
      quarto_args = c("--execute-dir", ".")
    ),
    tarchetypes::tar_quarto(
      report_release_validation, "reports/release_validation.qmd",
      output_file = "release_validation.html", working_directory = "."
    ),
    tar_target(
      published_report_data_foundation,
      copy_file_atomic_verified(
        report_data_foundation,
        file.path(config$root, "reports", "rendered", "data_foundation.html")
      ), format = "file"
    ),
    tar_target(
      published_report_response_estimation,
      copy_file_atomic_verified(
        report_response_estimation,
        file.path(config$root, "reports", "rendered", "response_estimation.html")
      ), format = "file"
    ),
    tar_target(
      published_report_resident_climate_model,
      copy_file_atomic_verified(
        report_resident_climate_model,
        file.path(config$root, "reports", "rendered", "resident_climate_model.html")
      ), format = "file"
    ),
    tar_target(
      published_report_family_context_model,
      copy_file_atomic_verified(
        report_family_context_model,
        file.path(config$root, "reports", "rendered", "family_context_model.html")
      ), format = "file"
    ),
    tar_target(
      published_report_predictor_models,
      copy_file_atomic_verified(
        report_predictor_models,
        file.path(config$root, "reports", "rendered", "predictor_models.html")
      ), format = "file"
    ),
    tar_target(
      published_report_model_synthesis,
      copy_file_atomic_verified(
        report_model_synthesis,
        file.path(config$root, "reports", "rendered", "model_synthesis.html")
      ), format = "file"
    ),
    tar_target(
      published_report_release_validation,
      copy_file_atomic_verified(
        report_release_validation,
        file.path(config$root, "reports", "rendered", "release_validation.html")
      ), format = "file"
    )
  )
}
