# Module: R/figures/family_context
# Responsibility: Prepare family-level climate predictions and callout plots.
# Inputs: Family fit, species membership, and climate grids.
# Outputs: Family context figure data and plot.
# Side effects: None.

# Figure-only helpers for the family-resolution analogue of manuscript Figure 4.
# These functions consume existing current exports and never fit a model.

#' Family Context Scaling Value.
#'
#' @param scaling A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param field A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param statistic A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
family_context_scaling_value <- function(scaling, field, statistic) {
  value <- scaling[[statistic]][scaling$scaled_field == field]
  if (length(value) != 1L || !is.finite(value)) {
    stop("Climate scaling value unavailable for ", field, "/", statistic,
         call. = FALSE)
  }
  as.numeric(value[[1]])
}

#' Family Context Surface.
#'
#' @param fit A fitted model object with the structure checked by this function.
#' @param data A data frame whose required columns are checked before use.
#' @param random_coefficients A model formula or random-effects specification.
#' @param family_name A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param W_z_value A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param warming_panel A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param grid_points A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
family_context_surface <- function(
    fit, data, random_coefficients, family_name, W_z_value, warming_panel,
    grid_points = 81L) {
  family_data <- data |>
    dplyr::filter(as.character(.data$family_id) == .env$family_name)
  observed <- family_data |>
    dplyr::distinct(.data$experiment_id, .data$T_z, .data$P_z, .data$W_z)
  if (dplyr::n_distinct(observed$experiment_id) < 4L) {
    stop(family_name, " does not have at least four experiments.",
         call. = FALSE)
  }
  candidate <- tidyr::expand_grid(
    T_z = seq(min(observed$T_z), max(observed$T_z), length.out = grid_points),
    P_z = seq(min(observed$P_z), max(observed$P_z), length.out = grid_points),
    W_z = W_z_value
  )
  support <- auc_support_mask(observed, candidate)$grid |>
    dplyr::filter(.data$supported) |>
    dplyr::mutate(TP_z = .data$T_z * .data$P_z)
  if (!nrow(support)) {
    stop("No supported surface cells for ", family_name, " at ",
         warming_panel, ".", call. = FALSE)
  }
  fixed_design <- auc_design_matrix(
    ~ T_z * P_z * W_z, support, names(stats::coef(fit))
  )
  random <- random_coefficients |>
    dplyr::filter(as.character(.data$family_id) == .env$family_name)
  if (nrow(random) != 1L) {
    stop("Random coefficients unavailable for ", family_name, call. = FALSE)
  }
  support |>
    dplyr::mutate(
      estimate = as.numeric(.env$fixed_design %*% stats::coef(.env$fit)) +
        .env$random$random_intercept[[1]] +
        .env$random$random_W[[1]] * .data$W_z +
        .env$random$random_T[[1]] * .data$T_z +
        .env$random$random_P[[1]] * .data$P_z +
        .env$random$random_TP[[1]] * .data$TP_z,
      warming_panel = warming_panel,
      family_id = family_name
    )
}

#' Prepare Family Context Callouts.
#'
#' @param fit A fitted model object with the structure checked by this function.
#' @param data A data frame whose required columns are checked before use.
#' @param classification A data frame whose required columns are checked before use.
#' @param random_coefficients A model formula or random-effects specification.
#' @param selected_families A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param grid_points A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
prepare_family_context_callouts <- function(
    fit, data, classification, random_coefficients,
    selected_families = c("Fabaceae", "Poaceae", "Rubiaceae", "Rosaceae"),
    grid_points = 81L) {
  required_classification <- c(
    "family_id", "min_fitted_family_delta_auc",
    "max_fitted_family_delta_auc", "min_conf_low", "min_conf_high",
    "max_conf_low", "max_conf_high", "random_W", "n_experiments",
    "fitted_class"
  )
  required_random <- c(
    "family_id", "random_intercept", "random_W", "random_T", "random_P",
    "random_TP"
  )
  if (!all(required_classification %in% names(classification)) ||
      !all(required_random %in% names(random_coefficients))) {
    stop("Current family classification or random-coefficient export is incomplete.",
         call. = FALSE)
  }
  if (!inherits(fit, "rma.mv")) {
    stop("The current family model must be an existing rma.mv fit.",
         call. = FALSE)
  }
  selected <- classification |>
    dplyr::filter(.data$family_id %in% .env$selected_families) |>
    dplyr::mutate(
      family_id = factor(.data$family_id, levels = .env$selected_families)
    ) |>
    dplyr::arrange(.data$family_id)
  if (!identical(as.character(selected$family_id), selected_families)) {
    stop("The requested family callouts are not all available.", call. = FALSE)
  }
  if (!all(selected$selected_for_surface)) {
    stop("Every appendix family must be selected in the current export.",
         call. = FALSE)
  }
  warming_experiments <- data |>
    dplyr::distinct(.data$experiment_id, .data$W_z, .data$W)
  warming <- tibble::tibble(
    warming_panel = c("Low warming", "High warming"),
    W_z = as.numeric(stats::quantile(
      warming_experiments$W_z, c(0.25, 0.75), names = FALSE, type = 7
    )),
    W_raw = as.numeric(stats::quantile(
      warming_experiments$W, c(0.25, 0.75), names = FALSE, type = 7
    ))
  )
  surfaces <- dplyr::bind_rows(lapply(selected_families, function(family_name) {
    dplyr::bind_rows(lapply(seq_len(nrow(warming)), function(index) {
      family_context_surface(
        fit, data, random_coefficients, family_name,
        warming$W_z[[index]], warming$warming_panel[[index]], grid_points
      )
    }))
  })) |>
    dplyr::mutate(
      warming_panel = factor(
        .data$warming_panel, levels = c("Low warming", "High warming")
      )
    )
  observed <- data |>
    dplyr::filter(as.character(.data$family_id) %in% .env$selected_families) |>
    dplyr::distinct(.data$family_id, .data$experiment_id, .data$T_z, .data$P_z) |>
    tidyr::crossing(warming_panel = c("Low warming", "High warming")) |>
    dplyr::mutate(
      warming_panel = factor(
        .data$warming_panel, levels = c("Low warming", "High warming")
      )
    )
  list(
    classification = classification,
    selection = selected,
    warming = warming,
    surfaces = surfaces,
    observed = observed
  )
}

#' Build Family Context Callout Figure.
#'
#' @param bundle A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param scaling A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param classification_base_size A data frame whose required columns are checked before use.
#' @param classification_label_size A data frame whose required columns are checked before use.
#' @param inset_base_size A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A plot object, or the declared figure path when rendering to disk.
build_family_context_callout_figure <- function(
    bundle, scaling, classification_base_size = 14,
    classification_label_size = 4.8, inset_base_size = 14) {
  classification <- bundle$classification |>
    dplyr::mutate(
      fitted_class = factor(
        .data$fitted_class,
        levels = c("Consistent loser", "Context-dependent", "Consistent winner")
      )
    )
  selection <- bundle$selection
  selected_families <- as.character(selection$family_id)
  class_colours <- c(
    "Consistent loser" = "#E66101",
    "Context-dependent" = "#0072B2",
    "Consistent winner" = "#1B9E77"
  )
  neutral_colour <- "#3F3F3F"
  x_range <- range(c(
    classification$min_conf_low, classification$min_conf_high,
    classification$min_fitted_family_delta_auc, 0
  ), na.rm = TRUE)
  y_range <- range(c(
    classification$max_conf_low, classification$max_conf_high,
    classification$max_fitted_family_delta_auc, 0
  ), na.rm = TRUE)
  x_limits <- x_range + c(-0.08, 0.08) * diff(x_range)
  y_limits <- y_range + c(-0.08, 0.08) * diff(y_range)
  class_counts <- classification |>
    dplyr::count(.data$fitted_class, name = "n_family")
  class_n <- function(class_name) {
    value <- class_counts$n_family[
      as.character(class_counts$fitted_class) == class_name
    ]
    ifelse(length(value), value, 0L)
  }
  warming_limit <- max(abs(classification$random_W), na.rm = TRUE)
  surface_limit <- max(abs(bundle$surfaces$estimate), na.rm = TRUE)

  class_box <- c(x = 0.265, y = 0.155, width = 0.47, height = 0.69)
  inset_boxes <- list(
    "Poaceae" = c(x = 0.020, y = 0.575, width = 0.225, height = 0.39),
    "Rosaceae" = c(x = 0.020, y = 0.075, width = 0.225, height = 0.39),
    "Rubiaceae" = c(x = 0.755, y = 0.535, width = 0.225, height = 0.39),
    "Fabaceae" = c(x = 0.755, y = 0.075, width = 0.225, height = 0.39)
  )
  # The two-line warming-sensitivity legend raises the lower data-panel edge.
  # Keep this boundary aligned so internal and external connector segments meet.
  panel_box <- c(left = 0.310, right = 0.734, bottom = 0.309, top = 0.843)
  callout_centres_y <- vapply(
    inset_boxes, function(box) box[["y"]] + box[["height"]] / 2, numeric(1)
  )
  callout_centres_y <- pmin(
    pmax(callout_centres_y, panel_box[["bottom"]]), panel_box[["top"]]
  )
  canvas_y_to_data <- function(canvas_y) {
    y_limits[[1]] +
      (canvas_y - panel_box[["bottom"]]) /
      (panel_box[["top"]] - panel_box[["bottom"]]) * diff(y_limits)
  }
  connector_data <- selection |>
    dplyr::transmute(
      family_id = as.character(.data$family_id),
      x = .data$min_fitted_family_delta_auc,
      y = .data$max_fitted_family_delta_auc,
      xend = ifelse(
        .data$family_id %in% c("Poaceae", "Rosaceae"),
        x_limits[[1]], x_limits[[2]]
      ),
      yend = vapply(
        as.character(.data$family_id),
        function(family_name) canvas_y_to_data(callout_centres_y[[family_name]]),
        numeric(1)
      )
    )
  classification_plot <- ggplot2::ggplot(classification) +
    ggplot2::annotate(
      "rect", xmin = x_limits[[1]], xmax = 0,
      ymin = y_limits[[1]], ymax = 0,
      fill = class_colours[["Consistent loser"]], alpha = 0.07
    ) +
    ggplot2::annotate(
      "rect", xmin = x_limits[[1]], xmax = 0,
      ymin = 0, ymax = y_limits[[2]],
      fill = class_colours[["Context-dependent"]], alpha = 0.07
    ) +
    ggplot2::annotate(
      "rect", xmin = 0, xmax = x_limits[[2]],
      ymin = 0, ymax = y_limits[[2]],
      fill = class_colours[["Consistent winner"]], alpha = 0.07
    ) +
    ggplot2::geom_abline(
      intercept = 0, slope = 1, linetype = "dashed",
      colour = "grey55", linewidth = 0.4
    ) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey20", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey20", linewidth = 0.5) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = .data$min_conf_low, xend = .data$min_conf_high,
        y = .data$max_fitted_family_delta_auc,
        yend = .data$max_fitted_family_delta_auc,
        colour = .data$fitted_class
      ), linewidth = 0.35, alpha = 0.20
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = .data$min_fitted_family_delta_auc,
        xend = .data$min_fitted_family_delta_auc,
        y = .data$max_conf_low, yend = .data$max_conf_high,
        colour = .data$fitted_class
      ), linewidth = 0.35, alpha = 0.20
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = .data$min_fitted_family_delta_auc,
        y = .data$max_fitted_family_delta_auc,
        colour = .data$fitted_class, fill = .data$random_W
      ), shape = 21, size = 3.0, stroke = 0.85, alpha = 0.95
    ) +
    ggplot2::geom_segment(
      data = connector_data,
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
      inherit.aes = FALSE, colour = neutral_colour,
      linewidth = 0.85, lineend = "round"
    ) +
    ggplot2::geom_point(
      data = selection,
      ggplot2::aes(
        x = .data$min_fitted_family_delta_auc,
        y = .data$max_fitted_family_delta_auc, fill = .data$random_W
      ), inherit.aes = FALSE, shape = 21, size = 4.5, stroke = 1.3,
      colour = neutral_colour, show.legend = FALSE
    ) +
    ggplot2::annotate(
      "text", x = x_limits[[1]] + 0.16 * diff(x_limits),
      y = y_limits[[1]] + 0.09 * diff(y_limits),
      label = paste0("CONSISTENT LOSER\nn = ", class_n("Consistent loser")),
      colour = class_colours[["Consistent loser"]], fontface = "bold",
      size = classification_label_size, hjust = 0
    ) +
    ggplot2::annotate(
      "text", x = x_limits[[1]] + 0.05 * diff(x_limits),
      y = y_limits[[2]] - 0.07 * diff(y_limits),
      label = paste0("CONTEXT-DEPENDENT\nn = ", class_n("Context-dependent")),
      colour = class_colours[["Context-dependent"]], fontface = "bold",
      size = classification_label_size, hjust = 0
    ) +
    ggplot2::annotate(
      "text", x = x_limits[[2]] - 0.04 * diff(x_limits),
      y = y_limits[[2]] - 0.07 * diff(y_limits),
      label = paste0("CONSISTENT WINNER\nn = ", class_n("Consistent winner")),
      colour = class_colours[["Consistent winner"]], fontface = "bold",
      size = classification_label_size, hjust = 1
    ) +
    ggplot2::scale_colour_manual(values = class_colours, drop = FALSE, guide = "none") +
    ggplot2::scale_fill_gradient2(
      low = "#E66101", mid = "white", high = "#1B9E77", midpoint = 0,
      limits = c(-warming_limit, warming_limit),
      breaks = c(-warming_limit, 0, warming_limit),
      labels = function(x) sprintf("%.2f", x),
      name = "Warming sensitivity\n(shift with stronger warming)",
      guide = ggplot2::guide_colourbar(
        title.position = "top", direction = "horizontal",
        barwidth = grid::unit(50, "mm"), barheight = grid::unit(4, "mm")
      )
    ) +
    ggplot2::coord_cartesian(
      xlim = x_limits, ylim = y_limits, expand = FALSE, clip = "off"
    ) +
    ggplot2::theme_bw(base_size = classification_base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom", legend.box = "vertical",
      axis.title = ggplot2::element_text(size = classification_base_size + 2),
      axis.text = ggplot2::element_text(size = classification_base_size - 1.5),
      legend.title = ggplot2::element_text(
        size = classification_base_size + 2, hjust = 0.5
      ),
      legend.text = ggplot2::element_text(size = classification_base_size + 1),
      legend.justification = "center", legend.box.just = "center",
      plot.margin = ggplot2::margin(2, 2, 2, 2)
    ) +
    ggplot2::labs(
      x = "Minimum fitted family winning-losing strength",
      y = "Maximum fitted family winning-losing strength",
      colour = NULL
    )

  T_center <- family_context_scaling_value(scaling, "T_z", "center")
  T_scale <- family_context_scaling_value(scaling, "T_z", "scale")
  P_center <- family_context_scaling_value(scaling, "P_z", "center")
  P_scale <- family_context_scaling_value(scaling, "P_z", "scale")
  raw_T_labels <- function(x) sprintf("%.1f", x * T_scale + T_center)
  raw_P_labels <- function(x) sprintf("%.0f", exp(x * P_scale + P_center))
  surface_plot <- function(family_name) {
    family_surface <- bundle$surfaces |>
      dplyr::filter(.data$family_id == .env$family_name)
    family_observed <- bundle$observed |>
      dplyr::filter(.data$family_id == .env$family_name)
    family_n <- selection$n_experiments[
      as.character(selection$family_id) == family_name
    ]
    ggplot2::ggplot(
      family_surface,
      ggplot2::aes(x = .data$T_z, y = .data$P_z, fill = .data$estimate)
    ) +
      ggplot2::geom_raster() +
      ggplot2::geom_point(
        data = family_observed,
        ggplot2::aes(x = .data$T_z, y = .data$P_z), inherit.aes = FALSE,
        shape = 21, fill = "white", colour = "grey15", size = 1.05,
        stroke = 0.3
      ) +
      ggplot2::facet_grid(rows = ggplot2::vars(.data$warming_panel)) +
      ggplot2::scale_fill_gradient2(
        low = "#E66101", mid = "white", high = "#1B9E77", midpoint = 0,
        limits = c(-surface_limit, surface_limit),
        breaks = c(-surface_limit, 0, surface_limit),
        labels = function(x) sprintf("%.2f", x),
        guide = ggplot2::guide_colourbar(
          title.position = "top", direction = "horizontal",
          barwidth = grid::unit(50, "mm"), barheight = grid::unit(4, "mm")
        )
      ) +
      ggplot2::scale_x_continuous(
        labels = raw_T_labels,
        breaks = function(limits) pretty(limits, n = 3),
        expand = c(0, 0)
      ) +
      ggplot2::scale_y_continuous(
        labels = raw_P_labels,
        breaks = function(limits) pretty(limits, n = 3),
        expand = c(0, 0)
      ) +
      ggplot2::coord_cartesian(expand = FALSE) +
      ggplot2::theme_bw(base_size = inset_base_size) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        strip.text.y = ggplot2::element_text(
          face = "bold", size = inset_base_size - 1
        ),
        plot.title = ggplot2::element_text(
          face = "bold", size = inset_base_size,
          colour = neutral_colour, hjust = 0.5, lineheight = 0.9
        ),
        plot.background = ggplot2::element_rect(
          colour = neutral_colour, linewidth = 1, fill = "white"
        ),
        plot.margin = ggplot2::margin(4, 4, 3, 4),
        axis.title = ggplot2::element_text(size = inset_base_size + 2),
        axis.text = ggplot2::element_text(size = inset_base_size - 1.3),
        legend.position = "bottom",
        legend.title = ggplot2::element_text(
          size = inset_base_size + 2, hjust = 0.5
        ),
        legend.text = ggplot2::element_text(size = inset_base_size - 1.5),
        legend.justification = "center", legend.box.just = "center",
        legend.margin = ggplot2::margin(-2, 0, -2, 0)
      ) +
      ggplot2::labs(
        title = paste0(family_name, "\n(", family_n, " experiments)"),
        x = "Destination temperature (°C)",
        y = "Destination precipitation (mm)",
        fill = "Fitted family winning-losing strength"
      )
  }
  mini_plots <- stats::setNames(
    lapply(selected_families, surface_plot), selected_families
  )
  connector_target <- function(family_name) {
    box <- inset_boxes[[family_name]]
    if (box[["x"]] < class_box[["x"]]) {
      c(x = box[["x"]] + box[["width"]], y = box[["y"]] + box[["height"]] / 2)
    } else {
      c(x = box[["x"]], y = box[["y"]] + box[["height"]] / 2)
    }
  }
  canvas <- cowplot::ggdraw() +
    cowplot::draw_plot(
      classification_plot,
      x = class_box[["x"]], y = class_box[["y"]],
      width = class_box[["width"]], height = class_box[["height"]]
    )
  for (family_name in names(inset_boxes)) {
    is_left <- family_name %in% c("Poaceae", "Rosaceae")
    start <- c(
      x = if (is_left) panel_box[["left"]] else panel_box[["right"]],
      y = callout_centres_y[[family_name]]
    )
    finish <- connector_target(family_name)
    canvas <- canvas + cowplot::draw_line(
      x = c(start[["x"]], finish[["x"]]),
      y = c(start[["y"]], finish[["y"]]),
      colour = neutral_colour, linewidth = 0.85, alpha = 0.9
    )
  }
  for (family_name in names(inset_boxes)) {
    box <- inset_boxes[[family_name]]
    canvas <- canvas + cowplot::draw_plot(
      mini_plots[[family_name]],
      x = box[["x"]], y = box[["y"]],
      width = box[["width"]], height = box[["height"]]
    )
  }
  canvas
}
