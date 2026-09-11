# Plotting data helpers. Model fitting is deliberately absent from this file.

functional_climate_coverage <- function(data, variable, label) {
  data |>
    dplyr::distinct(
      .data$species_id, .data$experiment_id, .data$T_z, .data$P_z, .data$W_z,
      group = .data[[variable]]
    ) |>
    dplyr::summarise(
      n_species = dplyr::n_distinct(.data$species_id),
      n_experiments = dplyr::n_distinct(.data$experiment_id),
      T_z_min = min(.data$T_z), T_z_max = max(.data$T_z),
      P_z_min = min(.data$P_z), P_z_max = max(.data$P_z),
      W_z_min = min(.data$W_z), W_z_max = max(.data$W_z),
      .by = "group"
    ) |>
    dplyr::mutate(predictor = label, .before = 1)
}

# Support is calculated separately within each functional group. This is a
# response-independent plotting mask and never changes the fitted model.
build_group_specific_climate_support <- function(data, group_column,
                                                 points = 65L) {
  groups <- unique(as.character(data[[group_column]]))
  result <- lapply(groups, function(group) {
    group_data <- data[as.character(data[[group_column]]) == group, , drop = FALSE]
    support <- build_auc_climate_grid(group_data, points = points)$grid
    support[[group_column]] <- group
    support
  })
  dplyr::bind_rows(result) |>
    dplyr::mutate(
      group = as.character(.data[[group_column]]),
      support_scope = "group-specific observed three-dimensional climate hull"
    )
}

build_woodiness_response_surface_figure <- function(
    model_sequence, support_mask, analysis_data, scaling) {
  fit <- model_sequence$fits$full_REML$fit
  group_levels <- levels(model_sequence$data$Woodiness)
  expected_levels <- c("non_woody", "woody")
  if (!identical(group_levels, expected_levels)) {
    stop("Woodiness model levels are not the canonical binary levels.",
         call. = FALSE)
  }
  required_mask <- c(
    "T_z", "P_z", "W_z", "warming_level", "supported", "group",
    "support_scope"
  )
  if (!all(required_mask %in% names(support_mask))) {
    stop("Woodiness support mask is incomplete.", call. = FALSE)
  }
  if (!all(grepl("group-specific", support_mask$support_scope, fixed = TRUE))) {
    stop("Woodiness figure requires group-specific observed climate support.",
         call. = FALSE)
  }

  warming_levels <- c("low (25th percentile)", "high (75th percentile)")
  grid <- support_mask |>
    dplyr::filter(
      .data$supported,
      as.character(.data$warming_level) %in% .env$warming_levels
    ) |>
    dplyr::mutate(
      Woodiness = factor(.data$group, levels = .env$group_levels),
      TP_z = .data$T_z * .data$P_z,
      warming_panel = factor(
        dplyr::recode(
          as.character(.data$warming_level),
          `low (25th percentile)` = "Low warming\n(group 25th percentile)",
          `high (75th percentile)` = "High warming\n(group 75th percentile)"
        ),
        levels = c(
          "Low warming\n(group 25th percentile)",
          "High warming\n(group 75th percentile)"
        )
      ),
      woodiness_panel = factor(
        dplyr::recode(.data$group, non_woody = "Non-woody", woody = "Woody"),
        levels = c("Non-woody", "Woody")
      )
    ) |>
    dplyr::group_by(.data$group) |>
    dplyr::mutate(
      tile_width = min(diff(sort(unique(.data$T_z)))) * 1.025,
      tile_height = min(diff(sort(unique(.data$P_z)))) * 1.025
    ) |>
    dplyr::ungroup()
  design <- auc_design_matrix(
    canonical_predictor_formulas()$Woodiness_full,
    grid,
    names(stats::coef(fit))
  )
  grid$predicted_delta_auc <- as.numeric(design %*% stats::coef(fit))

  observed <- analysis_data |>
    dplyr::transmute(
      group = as.character(.data$Woodiness),
      experiment_id = as.character(.data$experiment_id),
      T_z = .data$T_z,
      P_z = .data$P_z
    ) |>
    dplyr::distinct() |>
    tidyr::crossing(warming_panel = levels(grid$warming_panel)) |>
    dplyr::mutate(
      warming_panel = factor(
        .data$warming_panel, levels = levels(grid$warming_panel)
      ),
      woodiness_panel = factor(
        dplyr::recode(.data$group, non_woody = "Non-woody", woody = "Woody"),
        levels = levels(grid$woodiness_panel)
      )
    )

  scaling_value <- function(field, statistic) {
    value <- scaling[[statistic]][scaling$scaled_field == field]
    if (length(value) != 1L || !is.finite(value)) {
      stop("Climate scaling value unavailable for ", field, "/", statistic,
           call. = FALSE)
    }
    value[[1]]
  }
  T_center <- scaling_value("T_z", "center")
  T_scale <- scaling_value("T_z", "scale")
  P_center <- scaling_value("P_z", "center")
  P_scale <- scaling_value("P_z", "scale")
  limit <- max(abs(grid$predicted_delta_auc), na.rm = TRUE)

  ggplot2::ggplot(
    grid,
    ggplot2::aes(
      x = .data$T_z, y = .data$P_z, fill = .data$predicted_delta_auc
    )
  ) +
    ggplot2::geom_tile(
      ggplot2::aes(width = .data$tile_width, height = .data$tile_height)
    ) +
    ggplot2::geom_point(
      data = observed,
      ggplot2::aes(x = .data$T_z, y = .data$P_z),
      inherit.aes = FALSE, shape = 21, fill = "white", colour = "grey25",
      size = 1, stroke = 0.3
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(.data$warming_panel),
      cols = ggplot2::vars(.data$woodiness_panel)
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#E66101", mid = "white", high = "#1B9E77", midpoint = 0,
      limits = c(-limit, limit), name = "Predicted\nDelta AUC"
    ) +
    ggplot2::scale_x_continuous(
      labels = function(x) sprintf("%.1f", x * T_scale + T_center)
    ) +
    ggplot2::scale_y_continuous(
      labels = function(x) sprintf("%.0f", exp(x * P_scale + P_center))
    ) +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::labs(
      title = "Woodiness: anchored population-level response surfaces",
      subtitle = paste(
        "Blank cells lie outside each group's own observed three-dimensional",
        "temperature-precipitation-warming support."
      ),
      x = "Destination JJA temperature (°C)",
      y = "Destination JJA precipitation",
      caption = "Plotting-mask refinement only; the Step 22 Woodiness fit is unchanged."
    )
}

write_species_callout_examples <- function(selection, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(selection, path, na = "")
  path
}

species_callout_surface <- function(fit, data, random_coefficients,
                                    species_name, W_z_value,
                                    warming_panel, grid_points = 81L) {
  observed <- data |>
    dplyr::filter(as.character(.data$species_id) == .env$species_name) |>
    dplyr::distinct(.data$experiment_id, .data$T_z, .data$P_z)
  design <- cbind(1, observed$T_z, observed$P_z, observed$T_z * observed$P_z)
  if (nrow(observed) < 4L || qr(design)$rank != 4L) {
    stop(species_name, " does not have full-rank observed T-P support.",
         call. = FALSE)
  }
  grid <- tidyr::expand_grid(
    T_z = seq(min(observed$T_z), max(observed$T_z), length.out = grid_points),
    P_z = seq(min(observed$P_z), max(observed$P_z), length.out = grid_points)
  )
  hull <- geometry::convhulln(
    as.matrix(observed[, c("T_z", "P_z")]), options = "Qt"
  )
  grid <- grid |>
    dplyr::filter(geometry::inhulln(
      hull, as.matrix(dplyr::select(grid, "T_z", "P_z"))
    )) |>
    dplyr::mutate(W_z = W_z_value, TP_z = .data$T_z * .data$P_z)
  random <- random_coefficients |>
    dplyr::filter(as.character(.data$species_id) == .env$species_name)
  if (nrow(random) != 1L) {
    stop("Random coefficients unavailable for ", species_name, call. = FALSE)
  }
  fixed_design <- auc_design_matrix(
    ~ T_z * P_z * W_z, grid, names(stats::coef(fit))
  )
  predicted <- grid |>
    dplyr::mutate(
      estimate = as.numeric(fixed_design %*% stats::coef(fit)) +
        random$random_intercept[[1]] + random$random_W[[1]] * .data$W_z +
        random$random_T[[1]] * .data$T_z + random$random_P[[1]] * .data$P_z +
        random$random_TP[[1]] * .data$TP_z,
      tile_width = diff(range(observed$T_z)) / (grid_points - 1L) * 1.025,
      tile_height = diff(range(observed$P_z)) / (grid_points - 1L) * 1.025,
      warming_panel = warming_panel
    )
  list(
    grid = predicted,
    observed = dplyr::mutate(observed, warming_panel = warming_panel)
  )
}

build_species_context_callout_figure <- function(
    fit, data, classification, random_coefficients, scaling, selection,
    grid_points = 81L, classification_base_size = 14,
    classification_label_size = 4.8, inset_base_size = 14) {
  requested_species <- canonical_callout_species()
  if (!identical(as.character(selection$species_id), requested_species)) {
    stop("Callout examples are not in the canonical display order.", call. = FALSE)
  }
  classification <- classification |>
    dplyr::mutate(
      species_id = as.character(.data$species_id),
      fitted_class = factor(
        as.character(.data$fitted_class),
        levels = c("Consistent loser", "Context-dependent", "Consistent winner")
      )
    )
  warming_experiments <- data |>
    dplyr::distinct(.data$experiment_id, .data$W_z, .data$W)
  warming_lookup <- tibble::tibble(
    warming_name = c("Low warming", "High warming"),
    W_z = as.numeric(stats::quantile(warming_experiments$W_z, c(0.25, 0.75))),
    W_C = as.numeric(stats::quantile(warming_experiments$W, c(0.25, 0.75)))
  ) |>
    dplyr::mutate(
      warming_panel = .data$warming_name
    )
  panel_levels <- warming_lookup$warming_panel
  surface_sets <- stats::setNames(
    lapply(requested_species, function(species_name) {
      rows <- lapply(seq_len(nrow(warming_lookup)), function(index) {
        result <- species_callout_surface(
          fit, data, random_coefficients, species_name,
          warming_lookup$W_z[[index]], panel_levels[[index]], grid_points
        )
        list(
          grid = dplyr::mutate(
            result$grid,
            warming_panel = factor(.data$warming_panel, levels = panel_levels)
          ),
          observed = dplyr::mutate(
            result$observed,
            warming_panel = factor(.data$warming_panel, levels = panel_levels)
          )
        )
      })
      list(
        grid = dplyr::bind_rows(lapply(rows, `[[`, "grid")),
        observed = dplyr::bind_rows(lapply(rows, `[[`, "observed"))
      )
    }),
    requested_species
  )
  surface_limit <- max(abs(unlist(lapply(
    surface_sets, function(x) x$grid$estimate
  ))), na.rm = TRUE)
  warming_limit <- max(abs(classification$random_W), na.rm = TRUE)
  class_colours <- c(
    "Consistent loser" = "#E66101",
    "Context-dependent" = "#0072B2",
    "Consistent winner" = "#1B9E77"
  )
  neutral_colour <- "#3F3F3F"
  x_range <- range(c(
    classification$min_CI_low, classification$min_CI_high,
    classification$min_fitted_delta_auc, 0
  ), na.rm = TRUE)
  y_range <- range(c(
    classification$max_CI_low, classification$max_CI_high,
    classification$max_fitted_delta_auc, 0
  ), na.rm = TRUE)
  x_limits <- x_range + c(-0.08, 0.08) * diff(x_range)
  y_limits <- y_range + c(-0.08, 0.08) * diff(y_range)
  class_counts <- classification |>
    dplyr::count(.data$fitted_class, name = "n_species")
  class_n <- function(class_name) {
    class_counts$n_species[as.character(class_counts$fitted_class) == class_name]
  }

  class_box <- c(x = 0.265, y = 0.155, width = 0.47, height = 0.69)
  inset_boxes <- list(
    "Elymus repens" = c(x = 0.012, y = 0.575, width = 0.225, height = 0.39),
    "Silene acaulis" = c(x = 0.012, y = 0.075, width = 0.225, height = 0.39),
    "Potentilla erecta" = c(x = 0.763, y = 0.535, width = 0.225, height = 0.39),
    "Plantago lanceolata" = c(x = 0.763, y = 0.075, width = 0.225, height = 0.39)
  )
  # Canvas coordinates of the classification plot's data panel after the
  # enlarged manuscript axes and legend are laid out.
  panel_box <- c(left = 0.317, right = 0.734, bottom = 0.290, top = 0.843)
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
      species_id = as.character(.data$species_id),
      x = .data$min_fitted_delta_auc,
      y = .data$max_fitted_delta_auc,
      xend = ifelse(
        .data$species_id %in% c("Elymus repens", "Silene acaulis"),
        x_limits[[1]], x_limits[[2]]
      ),
      yend = vapply(
        .data$species_id,
        function(species_name) canvas_y_to_data(callout_centres_y[[species_name]]),
        numeric(1)
      )
    )
  classification_plot <- ggplot2::ggplot(classification) +
    ggplot2::annotate(
      "rect", xmin = x_limits[[1]], xmax = 0, ymin = y_limits[[1]], ymax = 0,
      fill = class_colours[["Consistent loser"]], alpha = 0.07
    ) +
    ggplot2::annotate(
      "rect", xmin = x_limits[[1]], xmax = 0, ymin = 0, ymax = y_limits[[2]],
      fill = class_colours[["Context-dependent"]], alpha = 0.07
    ) +
    ggplot2::annotate(
      "rect", xmin = 0, xmax = x_limits[[2]], ymin = 0, ymax = y_limits[[2]],
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
        x = .data$min_CI_low, xend = .data$min_CI_high,
        y = .data$max_fitted_delta_auc, yend = .data$max_fitted_delta_auc,
        colour = .data$fitted_class
      ), linewidth = 0.35, alpha = 0.17
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = .data$min_fitted_delta_auc, xend = .data$min_fitted_delta_auc,
        y = .data$max_CI_low, yend = .data$max_CI_high,
        colour = .data$fitted_class
      ), linewidth = 0.35, alpha = 0.17
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = .data$min_fitted_delta_auc, y = .data$max_fitted_delta_auc,
        colour = .data$fitted_class, fill = .data$random_W
      ), shape = 21, size = 2.5, stroke = 0.75, alpha = 0.94
    ) +
    ggplot2::geom_segment(
      data = connector_data,
      ggplot2::aes(
        x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend
      ), inherit.aes = FALSE, colour = neutral_colour,
      linewidth = 0.85, lineend = "round"
    ) +
    ggplot2::geom_point(
      data = selection,
      ggplot2::aes(
        x = .data$min_fitted_delta_auc, y = .data$max_fitted_delta_auc,
        fill = .data$random_W
      ), inherit.aes = FALSE, shape = 21, size = 4.2, stroke = 1.25,
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
      legend.justification = "center",
      legend.box.just = "center",
      plot.margin = ggplot2::margin(2, 2, 2, 2)
    ) +
    ggplot2::labs(
      x = "Minimum fitted winning-losing strength",
      y = "Maximum fitted winning-losing strength",
      colour = NULL
    )

  scaling_value <- function(field, statistic) {
    value <- scaling[[statistic]][scaling$scaled_field == field]
    if (length(value) != 1L || !is.finite(value)) {
      stop("Climate scaling value unavailable for ", field, "/", statistic,
           call. = FALSE)
    }
    value[[1]]
  }
  T_center <- scaling_value("T_z", "center")
  T_scale <- scaling_value("T_z", "scale")
  P_center <- scaling_value("P_z", "center")
  P_scale <- scaling_value("P_z", "scale")
  raw_T_labels <- function(x) sprintf("%.1f", x * T_scale + T_center)
  raw_P_labels <- function(x) sprintf("%.0f", exp(x * P_scale + P_center))
  surface_plot <- function(species_name) {
    ggplot2::ggplot(
      surface_sets[[species_name]]$grid,
      ggplot2::aes(x = .data$T_z, y = .data$P_z, fill = .data$estimate)
    ) +
      ggplot2::geom_tile(
        ggplot2::aes(width = .data$tile_width, height = .data$tile_height)
      ) +
      ggplot2::geom_point(
        data = surface_sets[[species_name]]$observed,
        ggplot2::aes(x = .data$T_z, y = .data$P_z), inherit.aes = FALSE,
        shape = 21, fill = "white", colour = "grey15", size = 1.05, stroke = 0.3
      ) +
      ggplot2::facet_grid(rows = ggplot2::vars(.data$warming_panel)) +
      ggplot2::scale_fill_gradient2(
        low = "#E66101", mid = "white", high = "#1B9E77", midpoint = 0,
        limits = c(-surface_limit, surface_limit),
        guide = ggplot2::guide_colourbar(
          title.position = "top", direction = "horizontal",
          barwidth = grid::unit(50, "mm"), barheight = grid::unit(4, "mm")
        )
      ) +
      ggplot2::scale_x_continuous(labels = raw_T_labels) +
      ggplot2::scale_y_continuous(labels = raw_P_labels) +
      ggplot2::theme_bw(base_size = inset_base_size) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        strip.text.y = ggplot2::element_text(
          face = "bold", size = inset_base_size - 1
        ),
        plot.title = ggplot2::element_text(
          face = "bold.italic", size = inset_base_size,
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
        legend.justification = "center",
        legend.box.just = "center",
        legend.margin = ggplot2::margin(-2, 0, -2, 0)
      ) +
      ggplot2::labs(
        title = paste0(
          species_name, "\n(",
          dplyr::filter(selection, .data$species_id == .env$species_name)$n_experiments,
          " experiments)"
        ),
        x = "Destination temperature (°C)",
        y = "Destination precipitation (mm)",
        fill = "Fitted winning/losing strength"
      )
  }
  mini_plots <- stats::setNames(
    lapply(requested_species, surface_plot), requested_species
  )
  connector_target <- function(species_name) {
    box <- inset_boxes[[species_name]]
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
  for (species_name in names(inset_boxes)) {
    is_left <- species_name %in% c("Elymus repens", "Silene acaulis")
    start <- c(
      x = if (is_left) panel_box[["left"]] else panel_box[["right"]],
      y = callout_centres_y[[species_name]]
    )
    finish <- connector_target(species_name)
    canvas <- canvas + cowplot::draw_line(
      x = c(start[["x"]], finish[["x"]]),
      y = c(start[["y"]], finish[["y"]]),
      colour = neutral_colour, linewidth = 0.85, alpha = 0.9
    )
  }
  for (species_name in names(inset_boxes)) {
    box <- inset_boxes[[species_name]]
    canvas <- canvas + cowplot::draw_plot(
      mini_plots[[species_name]],
      x = box[["x"]], y = box[["y"]],
      width = box[["width"]], height = box[["height"]]
    )
  }
  canvas
}

write_vector_figure <- function(figure, path, width, height) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = path, plot = figure, width = width, height = height,
    device = grDevices::cairo_pdf
  )
  path
}
