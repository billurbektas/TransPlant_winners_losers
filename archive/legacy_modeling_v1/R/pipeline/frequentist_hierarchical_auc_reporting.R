# Visual reporting for the resident continuous-AUC analysis -----------------

auc_reporting_required_packages <- function() {
  c("dplyr", "ggplot2", "knitr", "patchwork", "readr", "rmarkdown", "scales", "tidyr")
}

assert_auc_reporting_packages <- function() {
  missing <- auc_reporting_required_packages()[
    !vapply(
      auc_reporting_required_packages(), requireNamespace, logical(1), quietly = TRUE
    )
  ]
  if (length(missing)) {
    stop(
      "Install the resident AUC reporting dependencies before rendering: ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  if (!rmarkdown::pandoc_available()) {
    stop(
      "Pandoc is required to render the resident AUC HTML reports. ",
      "Run the standalone renderer after installing Pandoc.", call. = FALSE
    )
  }
  invisible(TRUE)
}

auc_report_palette <- function() {
  list(
    negative = "#0072B2", zero = "#7A7A7A", positive = "#D55E00",
    low = "#0072B2", mid = "#F7F7F7", high = "#D55E00",
    accent = "#009E73", purple = "#CC79A7", gold = "#E69F00"
  )
}

auc_report_theme <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title.position = "plot", plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )
}

save_auc_report_plot <- function(plot, stem, width = 11, height = 7) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  png_path <- paste0(stem, ".png")
  pdf_path <- paste0(stem, ".pdf")
  ggplot2::ggsave(
    png_path, plot = plot, width = width, height = height,
    units = "in", dpi = 180, bg = "white"
  )
  ggplot2::ggsave(
    pdf_path, plot = plot, width = width, height = height,
    units = "in", device = grDevices::cairo_pdf
  )
  c(png = png_path, pdf = pdf_path)
}

plot_auc_audit_effects <- function(data, species_replication, trait_coverage) {
  palette <- auc_report_palette()
  effect_data <- data |>
    dplyr::mutate(horizon = factor(
      as.character(.data$horizon), levels = c("max", "3"),
      labels = c("Maximum duration", "First 3 years")
    ))
  p_distribution <- ggplot2::ggplot(
    effect_data, ggplot2::aes(x = .data$auc_effect, fill = .data$horizon)
  ) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.4) +
    ggplot2::geom_histogram(bins = 55, alpha = 0.72, position = "identity") +
    ggplot2::facet_wrap(~ horizon, ncol = 1, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c(palette$purple, palette$accent), guide = "none") +
    auc_report_theme() +
    ggplot2::labs(
      title = "Continuous warmed-minus-origin AUC",
      subtitle = "One observation per resident species × experiment × horizon",
      x = "Raw AUC contrast", y = "Number of effects"
    )

  sign_data <- effect_data |>
    dplyr::mutate(
      response_sign = factor(
        dplyr::case_when(
          .data$auc_effect < 0 ~ "Negative", .data$auc_effect > 0 ~ "Positive",
          TRUE ~ "Zero"
        ),
        levels = c("Negative", "Zero", "Positive")
      )
    ) |>
    dplyr::count(.data$horizon, .data$response_sign, name = "n") |>
    tidyr::complete(.data$horizon, .data$response_sign, fill = list(n = 0L))
  p_sign <- ggplot2::ggplot(
    sign_data, ggplot2::aes(x = .data$response_sign, y = .data$n, fill = .data$response_sign)
  ) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::geom_text(ggplot2::aes(label = .data$n), vjust = -0.25, size = 3.3) +
    ggplot2::facet_wrap(~ horizon) +
    ggplot2::scale_fill_manual(values = c(
      Negative = palette$negative, Zero = palette$zero, Positive = palette$positive
    ), guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    auc_report_theme() +
    ggplot2::labs(
      title = "Direction, not significance category", x = NULL, y = "Number of effects"
    )

  replication_data <- species_replication |>
    dplyr::mutate(horizon = factor(
      as.character(.data$horizon), levels = c("max", "3"),
      labels = c("Maximum duration", "First 3 years")
    ))
  p_replication <- ggplot2::ggplot(
    replication_data, ggplot2::aes(x = .data$n_experiments, fill = .data$horizon)
  ) +
    ggplot2::geom_histogram(binwidth = 1, boundary = 0.5, colour = "white") +
    ggplot2::geom_vline(xintercept = 3, linetype = 2, colour = "grey30") +
    ggplot2::facet_wrap(~ horizon, ncol = 1, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c(palette$purple, palette$accent), guide = "none") +
    auc_report_theme() +
    ggplot2::labs(
      title = "Repeated-species frequency",
      subtitle = "Dashed line: three-experiment display threshold",
      x = "Experiments represented per species", y = "Number of species"
    )

  coverage_data <- trait_coverage |>
    dplyr::filter(as.character(.data$horizon) == "max") |>
    dplyr::select(
      "trait", effect_coverage = "effect_coverage",
      species_coverage = "species_coverage"
    ) |>
    tidyr::pivot_longer(
      -dplyr::all_of("trait"), names_to = "level", values_to = "coverage"
    ) |>
    dplyr::mutate(
      level = factor(
        .data$level, levels = c("effect_coverage", "species_coverage"),
        labels = c("Effects", "Species")
      )
    )
  p_traits <- ggplot2::ggplot(
    coverage_data, ggplot2::aes(x = .data$trait, y = .data$coverage, fill = .data$level)
  ) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.78), width = 0.7) +
    ggplot2::geom_text(
      ggplot2::aes(label = scales::percent(.data$coverage, accuracy = 1)),
      position = ggplot2::position_dodge(width = 0.78), vjust = -0.25, size = 3.2
    ) +
    ggplot2::scale_fill_manual(values = c(Effects = palette$gold, Species = palette$accent)) +
    ggplot2::scale_y_continuous(
      labels = scales::percent, limits = c(0, 1.08),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    auc_report_theme() +
    ggplot2::labs(
      title = "Trait coverage at the primary horizon", x = NULL, y = "Available", fill = NULL
    )

  patchwork::wrap_plots(
    p_distribution, p_sign, p_replication, p_traits, ncol = 2
  ) +
    patchwork::plot_annotation(tag_levels = "A")
}

plot_auc_audit_experiments <- function(experiment_audit, primary_horizon = "max") {
  palette <- auc_report_palette()
  primary <- experiment_audit |>
    dplyr::filter(as.character(.data$horizon) == as.character(.env$primary_horizon))
  horizon_data <- experiment_audit |>
    dplyr::mutate(horizon = factor(
      as.character(.data$horizon), levels = c("max", "3"),
      labels = c("Maximum duration", "First 3 years")
    ))
  p_climate <- ggplot2::ggplot(
    primary,
    ggplot2::aes(
      x = .data$T, y = .data$P_raw, colour = .data$W, size = .data$n_species
    )
  ) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::scale_colour_viridis_c(option = "C", end = 0.9) +
    auc_report_theme() +
    ggplot2::labs(
      title = "Observed experiment climate support",
      subtitle = "One point per experiment; size is resident-species count",
      x = "Destination JJA temperature (°C)",
      y = "Destination JJA precipitation", colour = "Warming (°C)", size = "Species"
    )

  p_duration <- ggplot2::ggplot(
    horizon_data, ggplot2::aes(x = .data$duration, fill = .data$horizon)
  ) +
    ggplot2::geom_histogram(binwidth = 1, boundary = 0.5, colour = "white") +
    ggplot2::facet_wrap(~ horizon, ncol = 1, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c(palette$purple, palette$accent), guide = "none") +
    auc_report_theme() +
    ggplot2::labs(
      title = "Experimental duration", x = "Inclusive analyzed years", y = "Experiments"
    )

  p_area <- ggplot2::ggplot(primary, ggplot2::aes(x = .data$plot_area)) +
    ggplot2::geom_histogram(bins = 18, fill = palette$gold, colour = "white") +
    ggplot2::scale_x_log10() + auc_report_theme() +
    ggplot2::labs(
      title = "Plot area", subtitle = "Primary horizon; logarithmic x-axis",
      x = "Plot area", y = "Experiments"
    )

  p_plots <- ggplot2::ggplot(
    primary, ggplot2::aes(x = .data$origin_plots, y = .data$warmed_plots)
  ) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey45") +
    ggplot2::geom_point(colour = palette$accent, size = 2.3, alpha = 0.85) +
    ggplot2::coord_equal() + auc_report_theme() +
    ggplot2::labs(
      title = "Plot replication by experiment",
      subtitle = "Dashed line indicates equal replication",
      x = "Distinct origin-control plots", y = "Distinct warmed plots"
    )

  patchwork::wrap_plots(p_climate, p_duration, p_area, p_plots, ncol = 2) +
    patchwork::plot_annotation(tag_levels = "A")
}

plot_auc_audit_missingness <- function(missingness) {
  traits <- c("LDMC", "LMA", "vegetative_height")
  display <- missingness |>
    dplyr::mutate(
      variable_group = dplyr::if_else(.data$column %in% traits, "Trait", "Core/design"),
      horizon = factor(
        as.character(.data$horizon), levels = c("max", "3"),
        labels = c("Maximum duration", "First 3 years")
      ),
      column = stats::reorder(.data$column, .data$proportion_missing)
    )
  ggplot2::ggplot(
    display,
    ggplot2::aes(x = .data$proportion_missing, y = .data$column, colour = .data$variable_group)
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = .data$proportion_missing, yend = .data$column),
      linewidth = 0.7
    ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~ horizon) +
    ggplot2::scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
    ggplot2::scale_colour_manual(values = c("Core/design" = "#0072B2", Trait = "#D55E00")) +
    auc_report_theme() +
    ggplot2::labs(
      title = "Missingness before model fitting",
      subtitle = "All inferential fields; one denominator per horizon",
      x = "Rows missing", y = NULL, colour = NULL
    )
}

render_auc_rmarkdown <- function(template, output_path, params) {
  assert_auc_reporting_packages()
  if (!file.exists(template)) stop("Missing report template: ", template, call. = FALSE)
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  rendered <- rmarkdown::render(
    input = template, output_file = basename(output_path),
    output_dir = dirname(output_path), params = params,
    envir = new.env(parent = globalenv()), quiet = TRUE,
    intermediates_dir = dirname(output_path), clean = TRUE
  )
  normalizePath(rendered, mustWork = TRUE)
}

auc_validate_saved_table <- function(output_dir, relative_path, required_columns) {
  path <- file.path(output_dir, relative_path)
  if (!file.exists(path)) {
    stop("Missing report input: ", relative_path, call. = FALSE)
  }
  columns <- names(readr::read_csv(path, n_max = 0, show_col_types = FALSE))
  missing <- setdiff(required_columns, columns)
  if (length(missing)) {
    stop(
      "Report input ", relative_path, " is missing required columns: ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  invisible(path)
}

render_frequentist_hierarchical_auc_data_audit <- function(
    output_dir = file.path("output", "frequentist_hierarchical_auc"),
    root = here::here()) {
  output_dir <- if (grepl("^/", output_dir)) output_dir else file.path(root, output_dir)
  auc_validate_saved_table(
    output_dir, "tables/data_audit.csv", c("horizon", "metric", "value")
  )
  auc_validate_saved_table(
    output_dir, "tables/column_map.csv",
    c("analysis_field", "source_field", "definition")
  )
  auc_validate_saved_table(
    output_dir, "tables/missingness.csv",
    c("horizon", "column", "n_missing", "n_rows", "proportion_missing")
  )
  auc_validate_saved_table(
    output_dir, "audit/experiment_audit.csv",
    c(
      "horizon", "experiment_id", "gradient_id", "T", "P_raw", "W",
      "duration", "integration_span", "plot_area", "n_species",
      "origin_plots", "warmed_plots"
    )
  )
  auc_validate_saved_table(
    output_dir, "audit/trait_coverage.csv",
    c(
      "horizon", "trait", "available_effects", "total_effects",
      "effect_coverage", "available_species", "total_species", "species_coverage"
    )
  )
  template <- file.path(
    root, "R", "reporting", "frequentist_hierarchical_auc_data_audit.Rmd"
  )
  render_auc_rmarkdown(
    template, file.path(output_dir, "audit", "data_audit.html"),
    list(output_dir = normalizePath(output_dir, mustWork = TRUE), root = root)
  )
}

write_auc_pre_model_audit <- function(
    data, horizons, auc_by_plot, scaled, output_dir, root = here::here(),
    render_html = TRUE) {
  message("Validating and writing the pre-model resident AUC audit ...")
  validation <- validate_auc_pre_model_data(data, horizons)
  audits <- purrr::map_dfr(
    as.character(horizons), ~ auc_data_audit(data, .x, auc_by_plot)
  )
  missingness <- auc_missingness_audit(data, horizons)
  species_replication <- auc_species_replication(data) |>
    dplyr::filter(as.character(.data$horizon) %in% as.character(.env$horizons))
  trait_coverage <- auc_trait_coverage(data, horizons)
  experiment_audit <- purrr::map_dfr(
    as.character(horizons), ~ auc_experiment_audit(data, .x, auc_by_plot)
  )
  trait_scaling <- purrr::map_dfr(
    c("LDMC", "LMA", "vegetative_height"),
    ~ trait_complete_data(data, .x)$scaling
  )

  write_auc_table(validation, output_dir, "audit", "pre_model_validation.csv")
  write_auc_table(auc_column_map(), output_dir, "tables", "column_map.csv")
  write_auc_table(audits, output_dir, "tables", "data_audit.csv")
  write_auc_table(missingness, output_dir, "tables", "missingness.csv")
  write_auc_table(scaled$scaling, output_dir, "tables", "climate_design_scaling.csv")
  write_auc_table(trait_scaling, output_dir, "tables", "trait_scaling.csv")
  write_auc_table(species_replication, output_dir, "tables", "species_replication.csv")
  write_auc_table(trait_coverage, output_dir, "audit", "trait_coverage.csv")
  write_auc_table(experiment_audit, output_dir, "audit", "experiment_audit.csv")
  write_auc_table(
    auc_context_support(data), output_dir, "tables", "climate_context_support.csv"
  )

  figure_dir <- file.path(output_dir, "audit", "figures")
  figures <- list(
    effects = save_auc_report_plot(
      plot_auc_audit_effects(data, species_replication, trait_coverage),
      file.path(figure_dir, "effect_and_biological_coverage"), 12, 9
    ),
    experiments = save_auc_report_plot(
      plot_auc_audit_experiments(experiment_audit, as.character(horizons[[1]])),
      file.path(figure_dir, "experiment_and_design_coverage"), 12, 9
    ),
    missingness = save_auc_report_plot(
      plot_auc_audit_missingness(missingness),
      file.path(figure_dir, "inferential_field_missingness"), 11, 8
    )
  )
  audit_report <- if (isTRUE(render_html)) {
    render_frequentist_hierarchical_auc_data_audit(output_dir, root)
  } else {
    NA_character_
  }
  invisible(list(
    validation = validation, data_audit = audits, missingness = missingness,
    species_replication = species_replication, trait_coverage = trait_coverage,
    experiment_audit = experiment_audit, trait_scaling = trait_scaling,
    figures = figures, report = audit_report
  ))
}

auc_require_report_files <- function(output_dir, relative_paths) {
  paths <- file.path(output_dir, relative_paths)
  missing <- relative_paths[!file.exists(paths)]
  if (length(missing)) {
    stop(
      "Cannot render the resident AUC visual report; missing outputs: ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  invisible(paths)
}

render_frequentist_hierarchical_auc_visual_report <- function(
    output_dir = file.path("output", "frequentist_hierarchical_auc"),
    report_file = "visual_report.html", root = here::here()) {
  output_dir <- if (grepl("^/", output_dir)) output_dir else file.path(root, output_dir)
  required <- c(
    "audit/data_audit.html", "audit/figures/effect_and_biological_coverage.png",
    "audit/figures/experiment_and_design_coverage.png",
    "audit/figures/inferential_field_missingness.png",
    "tables/data_audit.csv", "tables/model_robustness_ledger.csv",
    "tables/climate_cluster_robust_block_tests.csv",
    "tables/climate_coefficients.csv", "tables/auc_uncertainty_audit.csv",
    "tables/trait_ML_comparisons.csv", "tables/trait_cluster_robust_block_tests.csv",
    "tables/trait_species_slope_variance_reduction.csv",
    "tables/species_random_effects.csv", "tables/species_support.csv",
    "tables/horizon_comparison.csv", "tables/horizon_paired_effects.csv",
    "predictions/climate_surface_horizon_max.csv",
    "predictions/warming_slope_surface_horizon_max.csv",
    "diagnostics/influential_experiments.csv",
    "diagnostics/variance_profile_diagnostics.csv"
  )
  auc_require_report_files(output_dir, required)
  schemas <- list(
    "tables/data_audit.csv" = c("horizon", "metric", "value"),
    "tables/model_robustness_ledger.csv" = c(
      "horizon", "analysis", "model_id", "formula", "random_structure",
      "stable", "selected", "robust_method"
    ),
    "tables/climate_cluster_robust_block_tests.csv" = c(
      "horizon", "block", "statistic", "numerator_df", "denominator_df", "p_value"
    ),
    "tables/trait_ML_comparisons.csv" = c("trait", "comparison", "p_value"),
    "tables/trait_cluster_robust_block_tests.csv" = c("trait", "block", "p_value"),
    "tables/species_random_effects.csv" = c(
      "species_id", "component", "estimate", "conf_low", "conf_high"
    ),
    "tables/species_support.csv" = c("species_id", "display_eligible"),
    "predictions/climate_surface_horizon_max.csv" = c(
      "T_z", "P_z", "warming_level", "supported", "estimate", "robust_se"
    ),
    "predictions/warming_slope_surface_horizon_max.csv" = c(
      "T_z", "P_z", "warming_level", "supported", "estimate"
    )
  )
  purrr::iwalk(
    schemas,
    ~ auc_validate_saved_table(output_dir, .y, .x)
  )
  template <- file.path(
    root, "R", "reporting", "frequentist_hierarchical_auc_visual_report.Rmd"
  )
  output_path <- if (grepl("^/", report_file)) report_file else file.path(output_dir, report_file)
  render_auc_rmarkdown(
    template, output_path,
    list(output_dir = normalizePath(output_dir, mustWork = TRUE), root = root)
  )
}
