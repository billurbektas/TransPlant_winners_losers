# Module: R/figures/predictor_synthesis
# Responsibility: Implement the named scientific workflow stage or shared infrastructure.
# Inputs: Function arguments and explicitly declared target dependencies.
# Outputs: Validated in-memory objects and declared workflow artifacts.
# Side effects: Limited to functions whose contracts explicitly write, copy, or render files.

# Figure 6: predictor model support and residual heterogeneity synthesis.
#
# These helpers read exported summary tables only. They do not load
# model objects, refit models, or recompute earlier scientific results.

#' Predictor Synthesis Key.
#'
#' @return An R object described by the function name; its structure is validated before later stages use it.
predictor_synthesis_key <- function() {
  tibble::tribble(
    ~source_predictor,                   ~predictor,
    "LDMC",                              "LDMC",
    "LMA",                               "LMA",
    "Height",                            "Vegetative height",
    "Graminoid status",                  "Graminoid status",
    "Woodiness",                         "Woodiness",
    "Phylogenetic position (PEV1-PEV3)", "Phylogeny"
  )
}

#' Assert Predictor Synthesis Columns.
#'
#' @param data A data frame whose required columns are checked before use.
#' @param required A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param label A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
assert_predictor_synthesis_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(
      label, " lacks required column(s): ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

#' Prepare Predictor Synthesis.
#'
#' @param fit_file A project-relative or absolute filesystem path.
#' @param comparison_file A project-relative or absolute filesystem path.
#' @param heterogeneity_file A project-relative or absolute filesystem path.
#' @return An R object described by the function name; its structure is validated before later stages use it.
prepare_predictor_synthesis <- function(
    fit_file, comparison_file, heterogeneity_file) {
  input_files <- c(fit_file, comparison_file, heterogeneity_file)
  missing_inputs <- input_files[!file.exists(input_files)]
  if (length(missing_inputs) > 0L) {
    stop(
      "Required exported table(s) missing: ",
      paste(missing_inputs, collapse = ", "),
      call. = FALSE
    )
  }

  fit_statistics <- readr::read_csv(fit_file, show_col_types = FALSE)
  comparisons <- readr::read_csv(comparison_file, show_col_types = FALSE)
  heterogeneity <- readr::read_csv(
    heterogeneity_file, show_col_types = FALSE
  )
  predictor_key <- predictor_synthesis_key()

  required_fit_columns <- c("predictor", "model", "method", "AIC")
  required_comparison_columns <- c(
    "predictor", "comparison", "likelihood_ratio", "df", "p_value",
    "valid_nested_comparison", "same_rows_and_sampling_variances", "method"
  )
  required_heterogeneity_columns <- c(
    "Predictor",
    "Overall heterogeneity resolved (%)",
    "Warming heterogeneity resolved (%)",
    "Context heterogeneity resolved (%)",
    "Total heterogeneity resolved (%)"
  )
  assert_predictor_synthesis_columns(
    fit_statistics, required_fit_columns, basename(fit_file)
  )
  assert_predictor_synthesis_columns(
    comparisons, required_comparison_columns, basename(comparison_file)
  )
  assert_predictor_synthesis_columns(
    heterogeneity, required_heterogeneity_columns,
    basename(heterogeneity_file)
  )

  requested_heterogeneity <- heterogeneity |>
    dplyr::filter(.data$Predictor %in% predictor_key$source_predictor)
  if (nrow(requested_heterogeneity) != nrow(predictor_key) ||
      anyDuplicated(requested_heterogeneity$Predictor)) {
    stop(
      "Heterogeneity synthesis must contain exactly one row per predictor.",
      call. = FALSE
    )
  }

  # Support compares the full predictor-by-climate model directly with the
  # matched climate-only model. Negative delta AIC favours the full model.
  support_data <- predictor_key |>
    dplyr::left_join(
      fit_statistics |>
        dplyr::filter(
          .data$method == "ML", .data$model == "climate only"
        ) |>
        dplyr::select(
          source_predictor = "predictor",
          AIC_climate = "AIC"
        ),
      by = "source_predictor"
    ) |>
    dplyr::left_join(
      fit_statistics |>
        dplyr::filter(
          .data$method == "ML", .data$model == "full predictor-climate"
        ) |>
        dplyr::select(
          source_predictor = "predictor",
          AIC_full = "AIC"
        ),
      by = "source_predictor"
    ) |>
    dplyr::mutate(
      comparison = "climate only vs full predictor-climate",
      delta_AIC_full_minus_climate = .data$AIC_full - .data$AIC_climate
    ) |>
    dplyr::left_join(
      comparisons |>
        dplyr::filter(.data$method == "ML") |>
        dplyr::select(
          source_predictor = "predictor",
          comparison = "comparison",
          likelihood_ratio = "likelihood_ratio",
          lrt_df = "df",
          lrt_p_value = "p_value",
          valid_nested_comparison = "valid_nested_comparison",
          same_rows_and_sampling_variances =
            "same_rows_and_sampling_variances"
        ),
      by = c("source_predictor", "comparison")
    )

  if (nrow(support_data) != nrow(predictor_key) ||
      anyDuplicated(support_data$source_predictor)) {
    stop(
      "Support-table joins did not return exactly one row per predictor.",
      call. = FALSE
    )
  }
  if (any(!is.finite(support_data$AIC_climate)) ||
      any(!is.finite(support_data$AIC_full)) ||
      any(!is.finite(support_data$lrt_p_value))) {
    stop(
      "One or more requested support statistics are missing or non-finite.",
      call. = FALSE
    )
  }
  if (any(!support_data$valid_nested_comparison) ||
      any(!support_data$same_rows_and_sampling_variances)) {
    stop(
      "A requested support comparison is not a valid matched nested comparison.",
      call. = FALSE
    )
  }

  plot_data <- support_data |>
    dplyr::left_join(
      requested_heterogeneity |>
        dplyr::transmute(
          source_predictor = .data$Predictor,
          heterogeneity_overall_percent =
            .data$`Overall heterogeneity resolved (%)`,
          heterogeneity_warming_percent =
            .data$`Warming heterogeneity resolved (%)`,
          heterogeneity_background_climate_percent =
            .data$`Context heterogeneity resolved (%)`,
          heterogeneity_total_percent =
            .data$`Total heterogeneity resolved (%)`
        ),
      by = "source_predictor"
    ) |>
    dplyr::mutate(
      lrt_significance = dplyr::case_when(
        .data$lrt_p_value < 0.001 ~ "***",
        .data$lrt_p_value < 0.01 ~ "**",
        .data$lrt_p_value < 0.05 ~ "*",
        TRUE ~ ""
      )
    ) |>
    dplyr::select(
      "predictor", "source_predictor", "comparison",
      "AIC_climate", "AIC_full", "delta_AIC_full_minus_climate",
      "likelihood_ratio", "lrt_df", "lrt_p_value", "lrt_significance",
      "heterogeneity_overall_percent", "heterogeneity_warming_percent",
      "heterogeneity_background_climate_percent",
      "heterogeneity_total_percent"
    )

  heterogeneity_columns <- c(
    "heterogeneity_overall_percent",
    "heterogeneity_warming_percent",
    "heterogeneity_background_climate_percent",
    "heterogeneity_total_percent"
  )
  if (any(!is.finite(as.matrix(plot_data[heterogeneity_columns])))) {
    stop(
      "One or more requested heterogeneity percentages are missing or non-finite.",
      call. = FALSE
    )
  }
  plot_data
}

#' Build Predictor Synthesis Figure.
#'
#' @param plot_data A data frame whose required columns are checked before use.
#' @return A plot object, or the declared figure path when rendering to disk.
build_predictor_synthesis_figure <- function(plot_data) {
  predictor_order <- predictor_synthesis_key()$predictor
  heterogeneity_columns <- c(
    "heterogeneity_overall_percent",
    "heterogeneity_warming_percent",
    "heterogeneity_background_climate_percent",
    "heterogeneity_total_percent"
  )
  assert_predictor_synthesis_columns(
    plot_data,
    c(
      "predictor", "delta_AIC_full_minus_climate", "lrt_p_value",
      heterogeneity_columns
    ),
    "Figure 6 plotting data"
  )

  support_plot_data <- plot_data |>
    dplyr::mutate(
      predictor = factor(.data$predictor, levels = rev(predictor_order)),
      lrt_significance = dplyr::case_when(
        .data$lrt_p_value < 0.001 ~ "***",
        .data$lrt_p_value < 0.01 ~ "**",
        .data$lrt_p_value < 0.05 ~ "*",
        TRUE ~ ""
      ),
      lrt_p_label = paste0(
        "p = ", formatC(.data$lrt_p_value, format = "f", digits = 3),
        ifelse(
          .data$lrt_significance == "", "",
          paste0(" ", .data$lrt_significance)
        )
      )
    )
  support_span <- diff(range(c(
    0, support_plot_data$delta_AIC_full_minus_climate
  )))
  label_offset <- 0.035 * support_span
  support_plot_data <- support_plot_data |>
    dplyr::mutate(
      label_x = .data$delta_AIC_full_minus_climate +
        dplyr::if_else(
          .data$delta_AIC_full_minus_climate >= 0,
          label_offset,
          -label_offset
        ),
      label_hjust = dplyr::if_else(
        .data$delta_AIC_full_minus_climate >= 0, 0, 1
      )
    )

  panel_a <- ggplot2::ggplot(
    support_plot_data,
    ggplot2::aes(
      x = .data$delta_AIC_full_minus_climate,
      y = .data$predictor
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0, colour = "grey25", linewidth = 0.75
    ) +
    ggplot2::geom_col(
      width = 0.62, fill = "#3B4A6B", colour = "white", linewidth = 0.25
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = .data$label_x,
        label = .data$lrt_p_label,
        hjust = .data$label_hjust
      ),
      size = 5.3,
      fontface = "plain"
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.52, 0.45))
    ) +
    ggplot2::labs(title = "a", x = "ΔAIC (full − climate only)", y = NULL) +
    ggplot2::theme_bw(base_size = 20) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", size = 22, hjust = 0,
        margin = ggplot2::margin(b = 8)
      ),
      axis.title.x = ggplot2::element_text(size = 20),
      axis.text.y = ggplot2::element_text(size = 18, colour = "black"),
      axis.text.x = ggplot2::element_text(size = 16, colour = "black"),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(
        colour = "grey88", linewidth = 0.5
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        colour = "grey25", linewidth = 0.8
      ),
      plot.margin = ggplot2::margin(10, 8, 10, 10)
    )

  heterogeneity_long <- plot_data |>
    dplyr::mutate(
      predictor = factor(.data$predictor, levels = rev(predictor_order))
    ) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(heterogeneity_columns),
      names_to = "component",
      values_to = "percent_resolved"
    ) |>
    dplyr::mutate(
      component = dplyr::recode(
        .data$component,
        heterogeneity_overall_percent =
          "Response at mean climate\n(random species intercept)",
        heterogeneity_warming_percent =
          "Warming sensitivity\n(random warming slope)",
        heterogeneity_background_climate_percent =
          "Climate-context sensitivity\n(random T/P/T×P slope)",
        heterogeneity_total_percent = "Total response heterogeneity"
      ),
      component = factor(
        .data$component,
        levels = c(
          "Response at mean climate\n(random species intercept)",
          "Warming sensitivity\n(random warming slope)",
          "Climate-context sensitivity\n(random T/P/T×P slope)",
          "Total response heterogeneity"
        )
      )
    )
  component_colours <- c(
    "Response at mean climate\n(random species intercept)" = "#C06C9D",
    "Warming sensitivity\n(random warming slope)" = "#5A78A6",
    "Climate-context sensitivity\n(random T/P/T×P slope)" = "#8B6FB1",
    "Total response heterogeneity" = "#2F4B7C"
  )

  panel_b <- ggplot2::ggplot(
    heterogeneity_long,
    ggplot2::aes(
      x = .data$percent_resolved,
      y = .data$predictor,
      fill = .data$component
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0, colour = "grey25", linewidth = 0.75
    ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge2(
        width = 0.82, preserve = "single", padding = 0.08,
        reverse = TRUE
      ),
      width = 0.78,
      colour = "white",
      linewidth = 0.2
    ) +
    ggplot2::scale_fill_manual(values = component_colours, drop = FALSE) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.08, 0.08))
    ) +
    ggplot2::labs(
      title = "b", x = "% heterogeneity resolved", y = NULL, fill = NULL
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(nrow = 1, byrow = TRUE)
    ) +
    ggplot2::theme_bw(base_size = 20) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", size = 22, hjust = 0,
        margin = ggplot2::margin(b = 8)
      ),
      axis.title.x = ggplot2::element_text(size = 20),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 16, colour = "black"),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(
        colour = "grey88", linewidth = 0.5
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        colour = "grey25", linewidth = 0.8
      ),
      legend.position = "bottom",
      legend.justification = "center",
      legend.text = ggplot2::element_text(size = 16, lineheight = 0.9),
      legend.key.height = grid::unit(6, "mm"),
      legend.key.width = grid::unit(10, "mm"),
      legend.spacing.x = grid::unit(3, "mm"),
      plot.margin = ggplot2::margin(10, 10, 10, 8)
    )

  legend <- cowplot::get_legend(
    panel_b + ggplot2::theme(
      legend.position = "bottom", legend.justification = "center"
    )
  )
  panels <- cowplot::plot_grid(
    panel_a,
    panel_b + ggplot2::theme(legend.position = "none"),
    nrow = 1L,
    rel_widths = c(1.05, 1.35),
    align = "h",
    axis = "tb"
  )
  cowplot::plot_grid(
    panels, legend,
    ncol = 1L,
    rel_heights = c(1, 0.2)
  )
}
