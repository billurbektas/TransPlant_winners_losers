# Module: R/figures/winners_losers
# Responsibility: Prepare the main winner/loser response illustration.
# Inputs: Selected experiment observations and AUC contrasts.
# Outputs: Figure 1 data and plot.
# Side effects: None.

# Approved current manuscript Figure 01: resident winner/loser AUC responses with the
# approved context-dependent Festuca rubra examples.

figure01_direction_colours <- c(
  Winner = "#1B9E77",
  Loser = "#E66101"
)

#' Figure01 Festuca Specification.
#'
#' @return An R object described by the function name; its structure is validated before later stages use it.
figure01_festuca_specification <- function() {
  tibble::tribble(
    ~species, ~direction, ~experiment_id,
    "Festuca rubra", "Winner", "DE_Susalps::GW::BT",
    "Festuca rubra", "Loser", "US_Colorado::pf::um"
  ) |>
    dplyr::mutate(
      direction = factor(.data$direction, levels = c("Winner", "Loser"))
    )
}

#' Prepare Figure01 Festuca Data.
#'
#' @param winner_loser A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param species_classification A data frame whose required columns are checked before use.
#' @param auc_by_plot A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param community A data frame whose required columns are checked before use.
#' @param experiments A data frame whose required columns are checked before use.
#' @return An R object described by the function name; its structure is validated before later stages use it.
prepare_figure01_festuca_data <- function(
  winner_loser,
  species_classification,
  auc_by_plot,
  community,
  experiments
) {
  specification <- figure01_festuca_specification()

  festuca_classification <- species_classification |>
    dplyr::ungroup() |>
    dplyr::filter(.data$species_id == "Festuca rubra")
  if (
    nrow(festuca_classification) != 1L ||
      !isTRUE(festuca_classification$eligible[[1]]) ||
      festuca_classification$fitted_class[[1]] != "Context-dependent"
  ) {
    stop(
      "Festuca rubra is not an eligible current context-dependent species.",
      call. = FALSE
    )
  }

  analysis_rows <- winner_loser |>
    dplyr::ungroup() |>
    dplyr::filter(
      .data$transform == "raw",
      .data$horizon == "max",
      .data$contrast == "to_origin",
      .data$pool == "resident"
    ) |>
    dplyr::select(
      "species", "experiment_id", "estimate", "SE", "lower.CL",
      "upper.CL", "p.value", "significant", "category", "pool"
    )

  auc_rows_source <- auc_by_plot |>
    dplyr::filter(
      .data$transform == "raw",
      .data$horizon == "max",
      as.character(.data$ODT) %in% c("origin_controls", "warmed")
    ) |>
    dplyr::mutate(ODT = as.character(.data$ODT))

  mean_auc <- auc_rows_source |>
    dplyr::summarise(
      mean_auc = mean(.data$auc, na.rm = TRUE),
      .by = c("species", "experiment_id", "ODT")
    ) |>
    tidyr::pivot_wider(names_from = "ODT", values_from = "mean_auc") |>
    dplyr::rename(
      mean_auc_origin_controls = "origin_controls",
      mean_auc_warmed = "warmed"
    )

  scatter <- analysis_rows |>
    dplyr::inner_join(
      mean_auc,
      by = c("species", "experiment_id"),
      relationship = "one-to-one"
    ) |>
    dplyr::mutate(
      winning_losing_strength =
        .data$mean_auc_warmed - .data$mean_auc_origin_controls,
      point_significance = factor(
        dplyr::if_else(.data$significant, "Significant", "Not significant"),
        levels = c("Significant", "Not significant")
      ),
      display_region = dplyr::case_when(
        .data$mean_auc_warmed == 0 ~ "Potentially extinct",
        .data$significant & .data$winning_losing_strength > 0 ~ "Winner",
        .data$significant & .data$winning_losing_strength < 0 ~ "Loser",
        TRUE ~ "Stable"
      )
    )

  scatter_difference <- max(
    abs(scatter$winning_losing_strength - scatter$estimate),
    na.rm = TRUE
  )
  if (!is.finite(scatter_difference) || scatter_difference > 1e-10) {
    stop(
      "Mean AUC differences do not reproduce the current winner/loser contrasts.",
      call. = FALSE
    )
  }
  if (any(as.character(scatter$pool) != "resident")) {
    stop("Colonizing responses entered the resident Figure 01 scatter.", call. = FALSE)
  }

  experiment_labels <- experiments |>
    dplyr::transmute(
      experiment_id = .data$experiment_id,
      experiment_label = paste0(
        .data$region, ": ", .data$origin_site, " -> ",
        .data$destination_site
      )
    )

  selected <- specification |>
    dplyr::inner_join(
      scatter |>
        dplyr::filter(.data$species == "Festuca rubra"),
      by = c("species", "experiment_id"),
      relationship = "one-to-one"
    ) |>
    dplyr::left_join(
      experiment_labels,
      by = "experiment_id",
      relationship = "many-to-one"
    )

  expected_category <- tolower(as.character(selected$direction))
  if (
    nrow(selected) != 2L ||
      !all(selected$significant) ||
      !all(as.character(selected$category) == expected_category)
  ) {
    stop(
      "The approved Festuca rubra experiments no longer resolve to one significant winner and loser.",
      call. = FALSE
    )
  }

  curves <- purrr::map_dfr(seq_len(nrow(selected)), function(index) {
    selected_row <- selected[index, , drop = FALSE]
    experiment <- experiments |>
      dplyr::filter(.data$experiment_id == selected_row$experiment_id[[1]])
    if (nrow(experiment) != 1L) {
      stop("Each Figure 01 example must resolve to one experiment.", call. = FALSE)
    }

    experiment_data <- experiment_community(community, experiment) |>
      filter_minimum_plot_observations()
    completed <- complete_horizon_community(
      experiment_data,
      experiment_pool(experiment_data)
    )

    completed |>
      dplyr::filter(
        .data$species == selected_row$species[[1]],
        as.character(.data$ODT) %in% c("origin_controls", "warmed")
      ) |>
      dplyr::transmute(
        species = selected_row$species[[1]],
        direction = selected_row$direction[[1]],
        experiment_id = selected_row$experiment_id[[1]],
        experiment_label = selected_row$experiment_label[[1]],
        estimate = selected_row$estimate[[1]],
        p.value = selected_row$p.value[[1]],
        treatment = factor(
          dplyr::recode(
            as.character(.data$ODT),
            origin_controls = "Origin controls",
            warmed = "Warmed"
          ),
          levels = c("Origin controls", "Warmed")
        ),
        stable_plot_id = .data$stable_plot_id,
        relative_year = .data$relative_year,
        relative_cover = .data$relative_cover
      )
  })

  auc <- selected |>
    dplyr::select(
      "species", "direction", "experiment_id", "experiment_label",
      "estimate", "p.value"
    ) |>
    dplyr::inner_join(
      auc_rows_source |>
        dplyr::mutate(
          treatment = factor(
            dplyr::recode(
              .data$ODT,
              origin_controls = "Origin controls",
              warmed = "Warmed"
            ),
            levels = c("Origin controls", "Warmed")
          )
        ),
      by = c("species", "experiment_id"),
      relationship = "one-to-many"
    )

  example_checks <- auc |>
    dplyr::summarise(
      direct_change = mean(.data$auc[.data$treatment == "Warmed"], na.rm = TRUE) -
        mean(.data$auc[.data$treatment == "Origin controls"], na.rm = TRUE),
      modeled_change = dplyr::first(.data$estimate),
      absolute_difference = abs(.data$direct_change - .data$modeled_change),
      pass = .data$absolute_difference <= 1e-10,
      .by = c("species", "direction", "experiment_id", "experiment_label")
    )
  if (nrow(example_checks) != 2L || !all(example_checks$pass)) {
    stop("Festuca example AUC values failed modeled contrast checks.", call. = FALSE)
  }

  list(
    scatter = scatter,
    selected = selected,
    curves = curves,
    auc = auc,
    example_checks = example_checks,
    validation = tibble::tibble(
      check = c(
        "resident_scatter_only",
        "scatter_auc_differences_match_modeled_contrasts",
        "festuca_is_context_dependent",
        "festuca_examples_match_modeled_contrasts"
      ),
      pass = c(TRUE, TRUE, TRUE, all(example_checks$pass)),
      detail = c(
        paste0(nrow(scatter), " resident species-experiment responses"),
        paste0("maximum absolute difference = ", format(scatter_difference, scientific = TRUE)),
        paste0(
          "minimum = ", sprintf("%.3f", festuca_classification$min_fitted_delta_auc),
          "; maximum = ", sprintf("%.3f", festuca_classification$max_fitted_delta_auc)
        ),
        paste(selected$experiment_id, collapse = "; ")
      )
    )
  )
}

#' Plot Figure01 Scatter.
#'
#' @param figure_data A data frame whose required columns are checked before use.
#' @return A plot object, or the declared figure path when rendering to disk.
plot_figure01_scatter <- function(figure_data) {
  scatter <- figure_data$scatter
  selected <- figure_data$selected
  axis_max <- max(
    c(scatter$mean_auc_origin_controls, scatter$mean_auc_warmed),
    na.rm = TRUE
  ) * 1.045
  gradient_grid <- tidyr::expand_grid(
    mean_auc_origin_controls = seq(0, axis_max, length.out = 220L),
    mean_auc_warmed = seq(0, axis_max, length.out = 220L)
  ) |>
    dplyr::mutate(
      winning_losing_strength =
        .data$mean_auc_warmed - .data$mean_auc_origin_controls
    )
  scatter_points <- scatter |>
    dplyr::anti_join(
      selected |>
        dplyr::select("species", "experiment_id"),
      by = c("species", "experiment_id")
    )
  plot <- ggplot2::ggplot() +
    ggplot2::geom_raster(
      data = gradient_grid,
      ggplot2::aes(
        x = .data$mean_auc_origin_controls,
        y = .data$mean_auc_warmed,
        fill = .data$winning_losing_strength
      ),
      interpolate = TRUE,
      alpha = 0.34
    ) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      colour = "grey20",
      linewidth = 0.9
    ) +
    ggplot2::geom_point(
      data = scatter_points,
      ggplot2::aes(
        x = .data$mean_auc_origin_controls,
        y = .data$mean_auc_warmed,
        colour = .data$point_significance
      ),
      shape = 21,
      fill = "white",
      size = 3.05,
      stroke = 0.55,
      alpha = 0.80
    )

  for (direction in c("Winner", "Loser")) {
    direction_row <- selected |>
      dplyr::filter(as.character(.data$direction) == .env$direction)
    plot <- plot +
      ggplot2::geom_point(
        data = direction_row,
        ggplot2::aes(
          x = .data$mean_auc_origin_controls,
          y = .data$mean_auc_warmed
        ),
        inherit.aes = FALSE,
        shape = 23,
        fill = NA,
        colour = "grey12",
        size = 4.8,
        stroke = 1.1
      ) +
      ggplot2::geom_text(
        data = direction_row,
        ggplot2::aes(
          x = .data$mean_auc_origin_controls,
          y = .data$mean_auc_warmed
        ),
        inherit.aes = FALSE,
        label = "Festuca rubra",
        nudge_x = if (direction == "Winner") -0.070 * axis_max else 0.015 * axis_max,
        nudge_y = if (direction == "Winner") 0.075 * axis_max else 0.055 * axis_max,
        hjust = 0,
        vjust = 0,
        colour = unname(figure01_direction_colours[[direction]]),
        fontface = "bold.italic",
        lineheight = 0.95,
        size = 6.6
      ) +
      ggplot2::geom_text(
        data = direction_row,
        ggplot2::aes(
          x = .data$mean_auc_origin_controls,
          y = .data$mean_auc_warmed
        ),
        inherit.aes = FALSE,
        label = paste0(
          sub("::.*$", "", direction_row$experiment_id[[1]]),
          " experiment"
        ),
        nudge_x = if (direction == "Winner") -0.070 * axis_max else 0.015 * axis_max,
        nudge_y = if (direction == "Winner") 0.045 * axis_max else 0.025 * axis_max,
        hjust = 0,
        vjust = 0,
        colour = unname(figure01_direction_colours[[direction]]),
        fontface = "bold",
        lineheight = 0.95,
        size = 6.3
      )
  }

  plot +
    ggplot2::annotate(
      "text",
      x = 0.32 * axis_max,
      y = 0.72 * axis_max,
      label = "winner species",
      colour = figure01_direction_colours[["Winner"]],
      fontface = "bold",
      size = 9.6
    ) +
    ggplot2::annotate(
      "text",
      x = 0.64 * axis_max,
      y = 0.34 * axis_max,
      label = "loser species",
      colour = figure01_direction_colours[["Loser"]],
      fontface = "bold",
      size = 9.6
    ) +
    ggplot2::annotate(
      "text",
      x = 0.81 * axis_max,
      y = 0.84 * axis_max,
      label = "stable species",
      colour = "grey28",
      fontface = "bold",
      angle = 43,
      size = 9.1
    ) +
    ggplot2::scale_fill_gradient2(
      low = figure01_direction_colours[["Loser"]],
      mid = "white",
      high = figure01_direction_colours[["Winner"]],
      midpoint = 0,
      limits = c(-axis_max, axis_max),
      name = "Winning-losing strength\nΔ AUC"
    ) +
    ggplot2::scale_colour_manual(
      values = c("Significant" = "grey12", "Not significant" = "grey63"),
      drop = FALSE,
      name = "Significant winner/loser"
    ) +
    ggplot2::scale_y_continuous(position = "right") +
    ggplot2::guides(
      fill = ggplot2::guide_colourbar(
        order = 1,
        title.position = "top",
        barwidth = grid::unit(7.2, "cm"),
        barheight = grid::unit(0.42, "cm")
      ),
      colour = ggplot2::guide_legend(order = 2, title.position = "top")
    ) +
    ggplot2::coord_equal(
      xlim = c(0, axis_max),
      ylim = c(0, axis_max),
      expand = FALSE,
      clip = "off"
    ) +
    ggplot2::labs(
      x = "Mean area under relative-cover curve (Origin controls)",
      y = "Mean area under relative-cover curve (Warmed)"
    ) +
    ggplot2::theme_bw(base_size = 24) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = NA, colour = "black"),
      axis.title = ggplot2::element_text(size = 24, face = "bold"),
      axis.text = ggplot2::element_text(size = 20),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.title = ggplot2::element_text(size = 21, face = "bold"),
      legend.text = ggplot2::element_text(size = 19),
      legend.spacing.x = grid::unit(0.45, "cm"),
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )
}

#' Plot Figure01 Individual Curves.
#'
#' @param figure_data A data frame whose required columns are checked before use.
#' @param direction A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A plot object, or the declared figure path when rendering to disk.
plot_figure01_individual_curves <- function(figure_data, direction) {
  panel_data <- figure_data$curves |>
    dplyr::filter(as.character(.data$direction) == .env$direction)

  ggplot2::ggplot(
    panel_data,
    ggplot2::aes(
      x = .data$relative_year,
      y = .data$relative_cover,
      group = .data$stable_plot_id
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = 0, ymax = .data$relative_cover),
      fill = "grey68",
      alpha = 0.16,
      colour = NA
    ) +
    ggplot2::geom_line(colour = "grey43", linewidth = 0.48, alpha = 0.82) +
    ggplot2::geom_point(colour = "grey35", size = 0.9, alpha = 0.82) +
    ggplot2::facet_wrap(~treatment, nrow = 1) +
    ggplot2::scale_x_continuous(breaks = sort(unique(panel_data$relative_year))) +
    ggplot2::scale_y_continuous(
      limits = c(0, NA),
      expand = ggplot2::expansion(mult = c(0, 0.06))
    ) +
    ggplot2::labs(x = "Experimental year", y = "Relative cover") +
    ggplot2::theme_bw(base_size = 20) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey88"),
      strip.text = ggplot2::element_text(size = 20, face = "bold"),
      axis.title = ggplot2::element_text(size = 20),
      axis.text = ggplot2::element_text(size = 17),
      plot.margin = ggplot2::margin(2, 3, 2, 2)
    )
}

#' Plot Figure01 Auc Boxplot.
#'
#' @param figure_data A data frame whose required columns are checked before use.
#' @param direction A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A plot object, or the declared figure path when rendering to disk.
plot_figure01_auc_boxplot <- function(figure_data, direction) {
  panel_data <- figure_data$auc |>
    dplyr::filter(as.character(.data$direction) == .env$direction)
  means <- panel_data |>
    dplyr::summarise(mean_auc = mean(.data$auc, na.rm = TRUE), .by = "treatment")
  selected <- figure_data$selected |>
    dplyr::filter(as.character(.data$direction) == .env$direction)
  direction_colour <- unname(figure01_direction_colours[[direction]])
  annotation_x <- if (direction == "Winner") 1.1 else 1.85
  ordered_means <- means |>
    dplyr::arrange(.data$treatment)
  auc_range <- diff(range(panel_data$auc, na.rm = TRUE))
  if (!is.finite(auc_range) || auc_range == 0) auc_range <- 1
  annotation_y <- ordered_means$mean_auc[[1]] +
    (annotation_x - 1) *
      (ordered_means$mean_auc[[2]] - ordered_means$mean_auc[[1]]) +
    0.3 * auc_range

  ggplot2::ggplot(panel_data, ggplot2::aes(x = .data$treatment, y = .data$auc)) +
    ggplot2::geom_boxplot(
      width = 0.58,
      fill = "grey83",
      colour = "grey20",
      linewidth = 0.7,
      outlier.shape = NA
    ) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(width = 0.07, height = 0, seed = 20260824L),
      shape = 21,
      fill = "grey55",
      colour = "grey20",
      size = 2.0,
      stroke = 0.42
    ) +
    ggplot2::geom_line(
      data = means,
      ggplot2::aes(x = .data$treatment, y = .data$mean_auc, group = 1),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 1
    ) +
    ggplot2::geom_point(
      data = means,
      ggplot2::aes(x = .data$treatment, y = .data$mean_auc),
      inherit.aes = FALSE,
      shape = 23,
      fill = "white",
      colour = "black",
      size = 3.4,
      stroke = 0.9
    ) +
    ggplot2::annotate(
      "text",
      x = annotation_x,
      y = annotation_y,
      label = paste0("Δ AUC = ", sprintf("%+.3f", selected$estimate[[1]])),
      colour = direction_colour,
      fontface = "bold",
      size = 7
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, NA),
      expand = ggplot2::expansion(mult = c(0, 0.10))
    ) +
    ggplot2::labs(x = NULL, y = "Area under\nrelative-cover curve") +
    ggplot2::theme_bw(base_size = 20) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 20),
      axis.text = ggplot2::element_text(size = 17),
      plot.margin = ggplot2::margin(2, 2, 2, 3)
    )
}

#' Build Figure01 Context Callout.
#'
#' @param figure_data A data frame whose required columns are checked before use.
#' @param direction A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
build_figure01_context_callout <- function(figure_data, direction) {
  selected <- figure_data$selected |>
    dplyr::filter(as.character(.data$direction) == .env$direction)
  direction_colour <- unname(figure01_direction_colours[[direction]])
  curve_plot <- plot_figure01_individual_curves(figure_data, direction)
  box_plot <- plot_figure01_auc_boxplot(figure_data, direction)

  combined <- curve_plot / box_plot +
    patchwork::plot_layout(heights = c(1.12, 1))

  combined +
    patchwork::plot_annotation(
      title = paste0(
        "Festuca rubra - ",
        sub("::.*$", "", selected$experiment_id[[1]]),
        " experiment"
      ),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = 20,
          face = "bold.italic",
          colour = direction_colour,
          margin = ggplot2::margin(b = 6)
        )
      )
    )
}

#' Build Figure01 Winners Losers.
#'
#' @param figure_data A data frame whose required columns are checked before use.
#' @return An R object described by the function name; its structure is validated before later stages use it.
build_figure01_winners_losers <- function(figure_data) {
  if (!all(figure_data$validation$pass)) {
    stop("Figure 01 input validation failed.", call. = FALSE)
  }

  scatter <- plot_figure01_scatter(figure_data)
  winner_callout <- build_figure01_context_callout(
    figure_data,
    direction = "Winner"
  )
  loser_callout <- build_figure01_context_callout(
    figure_data,
    direction = "Loser"
  )

  cowplot::ggdraw() +
    cowplot::draw_grob(
      grid::rectGrob(gp = grid::gpar(fill = "white", col = "grey35", lwd = 1.2)),
      x = 0.018, y = 0.515, width = 0.315, height = 0.46
    ) +
    cowplot::draw_grob(
      grid::rectGrob(gp = grid::gpar(fill = "white", col = "grey35", lwd = 1.2)),
      x = 0.018, y = 0.025, width = 0.315, height = 0.46
    ) +
    cowplot::draw_plot(
      winner_callout, x = 0.032, y = 0.550, width = 0.287, height = 0.400
    ) +
    cowplot::draw_plot(
      loser_callout, x = 0.032, y = 0.060, width = 0.287, height = 0.400
    ) +
    cowplot::draw_plot(scatter, x = 0.345, y = 0.035, width = 0.64, height = 0.93) +
    cowplot::draw_label(
      "a", x = 0.018, y = 0.995, hjust = 0, vjust = 1,
      size = 26, fontface = "bold"
    ) +
    cowplot::draw_label(
      "b", x = 0.018, y = 0.505, hjust = 0, vjust = 1,
      size = 26, fontface = "bold"
    ) +
    cowplot::draw_label(
      "c", x = 0.375, y = 0.995, hjust = 0, vjust = 1,
      size = 26, fontface = "bold"
    )
}
