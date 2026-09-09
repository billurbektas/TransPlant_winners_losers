# Canonical temporal figure functions.
# Each function receives explicit scientific objects and an explicit output path.

plot_transient_strength <- function(classified, path) {
  data <- classified |>
    dplyr::ungroup() |>
    dplyr::filter(.data$transform == "raw") |>
    dplyr::select(experiment_id, species, horizon, estimate) |>
    tidyr::pivot_wider(names_from = horizon, values_from = estimate, names_prefix = "h_") |>
    tidyr::pivot_longer(c(h_4, h_6), names_to = "later", values_to = "later_estimate") |>
    dplyr::mutate(later = factor(.data$later, levels = c("h_4", "h_6"), labels = c("First 4 years", "First 6 years")))
  correlations <- data |>
    dplyr::filter(is.finite(.data$h_2), is.finite(.data$later_estimate)) |>
    dplyr::summarise(
      correlation = stats::cor(.data$h_2, .data$later_estimate, method = "pearson"),
      p_value = stats::cor.test(
        .data$h_2, .data$later_estimate, method = "pearson"
      )$p.value,
      .by = "later"
    ) |>
    dplyr::mutate(
      p_label = dplyr::if_else(
        .data$p_value < 0.001,
        "p < 0.001",
        sprintf("p = %.3f", .data$p_value)
      ),
      label = sprintf("Pearson r = %.2f\n%s", .data$correlation, .data$p_label)
    )
  plot <- ggplot2::ggplot(data, ggplot2::aes(.data$h_2, .data$later_estimate)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey75") +
    ggplot2::geom_vline(xintercept = 0, colour = "grey75") +
    ggplot2::geom_abline(slope = 1, colour = "grey60") +
    ggplot2::geom_point(alpha = 0.35, size = 0.7, na.rm = TRUE) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, colour = "#2C3E50", na.rm = TRUE) +
    ggplot2::geom_text(
      data = correlations,
      ggplot2::aes(x = -Inf, y = Inf, label = .data$label),
      inherit.aes = FALSE, hjust = -0.08, vjust = 1.35,
      fontface = "bold", size = 4
    ) +
    ggplot2::facet_wrap(~ later, scales = "free") +
    transplant_figure_theme() +
    ggplot2::labs(x = "Strength over the first 2 years", y = "Strength over the later horizon")
  save_plot_atomic(plot, path, 9, 4.5)
}

temporal_pairwise_flow_summary <- function(linked_transitions) {
  flow_weight <- if ("n" %in% names(linked_transitions)) linked_transitions$n else 1
  weighted <- linked_transitions |>
    dplyr::mutate(.flow_weight = .env$flow_weight)
  summarise_pair <- function(from_column, to_column, source_x, destination_x, pair_id,
                             comparison) {
    weighted |>
      dplyr::transmute(
        source_state = .data[[from_column]],
        destination_state = .data[[to_column]],
        n = .data$.flow_weight
      ) |>
      dplyr::group_by(.data$source_state, .data$destination_state) |>
      dplyr::summarise(n = sum(.data$n), .groups = "drop") |>
      dplyr::mutate(
        source_x = .env$source_x, destination_x = .env$destination_x,
        pair_id = .env$pair_id, comparison = .env$comparison
      )
  }
  dplyr::bind_rows(
    summarise_pair(
      "state_2", "state_4", 1, 2, "2_to_4", "First 2 → First 4"
    ),
    summarise_pair(
      "state_4", "state_6", 2, 3, "4_to_6", "First 4 → First 6"
    )
  ) |>
    dplyr::mutate(
      visible = .data$source_state != "follow_up_unavailable" &
        .data$destination_state != "follow_up_unavailable"
    ) |>
    dplyr::group_by(.data$pair_id, .data$source_state) |>
    dplyr::mutate(
      observed_source_total = sum(.data$n[.data$visible]),
      row_proportion = dplyr::if_else(
        .data$visible, .data$n / .data$observed_source_total, NA_real_
      )
    ) |>
    dplyr::ungroup()
}

plot_temporal_flows <- function(linked_transitions, path) {
  state_levels <- c(
    "loser", "non_significant_loser", "non_significant_winner",
    "winner", "colonizer", "non_significant_colonizer",
    "follow_up_unavailable"
  )
  visible_levels <- setdiff(state_levels, "follow_up_unavailable")
  flow_colours <- c(category_palette, follow_up_unavailable = "white")
  horizon_columns <- c("state_2", "state_4", "state_6")
  flow_weight <- if ("n" %in% names(linked_transitions)) linked_transitions$n else 1
  weighted_transitions <- linked_transitions |>
    dplyr::mutate(.flow_weight = .env$flow_weight)
  lodes <- weighted_transitions |>
    tidyr::pivot_longer(
      dplyr::all_of(horizon_columns), names_to = "horizon", values_to = "state"
    ) |>
    dplyr::mutate(
      x = match(.data$horizon, .env$horizon_columns),
      state = factor(.data$state, levels = .env$state_levels),
      stratum_alpha = dplyr::if_else(
        as.character(.data$state) == "follow_up_unavailable", 0, 1
      )
    )
  pair_flows <- temporal_pairwise_flow_summary(linked_transitions) |>
    dplyr::mutate(
      flow_id = paste(.data$pair_id, .data$source_state, .data$destination_state, sep = "::"),
      flow_alpha = dplyr::if_else(
        .data$source_state != "follow_up_unavailable" &
          .data$destination_state != "follow_up_unavailable",
        0.52, 0
      ),
      source_state = factor(.data$source_state, levels = .env$state_levels),
      source_fill = .data$source_state
    )
  pair_lodes <- pair_flows |>
    tidyr::pivot_longer(
      c("source_state", "destination_state"),
      names_to = "endpoint", values_to = "state"
    ) |>
    dplyr::mutate(
      x = dplyr::if_else(
        .data$endpoint == "source_state", .data$source_x, .data$destination_x
      ),
      state = factor(.data$state, levels = .env$state_levels)
    )

  flow_layer <- function(pair_id) {
    ggalluvial::geom_flow(
      data = dplyr::filter(pair_lodes, .data$pair_id == .env$pair_id),
      mapping = ggplot2::aes(
        x = .data$x, stratum = .data$state, alluvium = .data$flow_id,
        y = .data$n, fill = .data$source_fill, alpha = .data$flow_alpha
      ),
      inherit.aes = FALSE, width = 0.30, curve_type = "quintic",
      colour = scales::alpha("white", 0), linewidth = 0.1
    )
  }

  base <- ggplot2::ggplot() +
    flow_layer("2_to_4") +
    flow_layer("4_to_6") +
    ggalluvial::geom_stratum(
      data = lodes,
      mapping = ggplot2::aes(
        x = .data$x, stratum = .data$state, y = .data$.flow_weight,
        fill = .data$state, alpha = .data$stratum_alpha
      ),
      inherit.aes = FALSE,
      width = 0.30, colour = "white", linewidth = 0.7
    ) +
    ggplot2::scale_alpha_identity(guide = "none") +
    ggplot2::scale_fill_manual(
      values = flow_colours,
      breaks = visible_levels,
      labels = clean_figure_label(unname(category_display_labels[visible_levels])),
      drop = FALSE,
      name = NULL
    )

  built <- ggplot2::ggplot_build(base)$data
  flow_positions <- purrr::map2_dfr(
    built[1:2], c("2_to_4", "4_to_6"),
    function(layer, pair_id) {
      tibble::as_tibble(layer) |>
        dplyr::transmute(
          pair_id = .env$pair_id, PANEL = .data$PANEL, group = .data$group,
          x = .data$x, flow = .data$flow, stratum = .data$stratum,
          n = .data$ymax - .data$ymin, y = .data$y
        )
    }
  )
  flow_labels <- flow_positions |>
    dplyr::group_by(.data$pair_id, .data$PANEL, .data$group) |>
    dplyr::filter(dplyr::n() == 2L) |>
    dplyr::summarise(
      source_x = .data$x[.data$flow == "from"],
      source_y = .data$y[.data$flow == "from"],
      destination_y = .data$y[.data$flow == "to"],
      source_state = as.character(.data$stratum[.data$flow == "from"]),
      destination_state = as.character(.data$stratum[.data$flow == "to"]),
      n = dplyr::first(.data$n),
      .groups = "drop"
    ) |>
    dplyr::filter(
      .data$source_state != "follow_up_unavailable",
      .data$destination_state != "follow_up_unavailable"
    ) |>
    dplyr::mutate(
      source_response = dplyr::case_when(
        .data$source_state %in% c("loser", "non_significant_loser") ~ "loser",
        .data$source_state %in% c("winner", "non_significant_winner") ~ "winner",
        TRUE ~ "other"
      ),
      destination_response = dplyr::case_when(
        .data$destination_state %in% c("loser", "non_significant_loser") ~ "loser",
        .data$destination_state %in% c("winner", "non_significant_winner") ~ "winner",
        TRUE ~ "other"
      )
    ) |>
    dplyr::group_by(.data$source_x, .data$source_state) |>
    dplyr::mutate(
      transition_proportion = .data$n / sum(.data$n),
      percentage_label = dplyr::if_else(
        100 * .data$transition_proportion < 1,
        sprintf("%.1f%%", 100 * .data$transition_proportion),
        sprintf("%.0f%%", 100 * .data$transition_proportion)
      ),
      label = sprintf(
        "%d (%s)", round(.data$n), .data$percentage_label
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(
      .data$source_response %in% c("loser", "winner"),
      .data$destination_response %in% c("loser", "winner"),
      .data$source_response != .data$destination_response
    ) |>
    dplyr::mutate(
      label_x = .data$source_x + 0.78,
      label_y = .data$destination_y
    )
  stratum_labels <- tibble::as_tibble(built[[3]]) |>
    dplyr::transmute(
      x = .data$x, y = .data$y, state = as.character(.data$stratum),
      n = .data$ymax - .data$ymin, stratum_height = .data$ymax - .data$ymin
    ) |>
    dplyr::filter(.data$state != "follow_up_unavailable", .data$n > 0) |>
    dplyr::mutate(
      label = format(round(.data$n), scientific = FALSE, trim = TRUE)
    )

  plot <- base +
    ggrepel::geom_label_repel(
      data = flow_labels,
      ggplot2::aes(.data$label_x, .data$label_y, label = .data$label),
      inherit.aes = FALSE, size = 3.0, fontface = "bold", linewidth = 0,
      label.padding = grid::unit(0.06, "lines"),
      fill = scales::alpha("white", 0.82), colour = "grey15",
      direction = "y", box.padding = 0.12, point.padding = 0,
      min.segment.length = Inf, max.overlaps = Inf, seed = 47163
    ) +
    ggplot2::geom_label(
      data = stratum_labels,
      ggplot2::aes(.data$x, .data$y, label = .data$label),
      inherit.aes = FALSE, size = 3.6, fontface = "bold", lineheight = 0.88,
      label.padding = grid::unit(0.12, "lines"), linewidth = 0,
      fill = scales::alpha("white", 0.90), colour = "grey10"
    ) +
    ggplot2::scale_x_continuous(
      breaks = 1:3, labels = c("First 2 years", "First 4 years", "First 6 years"),
      limits = c(0.78, 3.22), expand = c(0, 0), position = "top"
    ) +
    ggplot2::scale_y_continuous(breaks = NULL, expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::theme_void(base_size = 16) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        colour = "#1A1A1A", size = 16, face = "bold", margin = ggplot2::margin(b = 12)
      ),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 13),
      legend.key.width = grid::unit(1.5, "lines"),
      plot.margin = ggplot2::margin(12, 32, 8, 32)
    )
  save_plot_atomic(plot, path, 18, 11)
}

plot_temporal_transition_matrices <- function(transitions, composition, path) {
  category_levels <- c(
    "loser", "non_significant_loser", "non_significant_winner", "winner",
    "colonizer", "non_significant_colonizer"
  )
  category_labels <- c(
    loser = "Loser",
    non_significant_loser = "Loser\n(non-significant)",
    non_significant_winner = "Winner\n(non-significant)",
    winner = "Winner",
    colonizer = "Colonizer",
    non_significant_colonizer = "Colonizer\n(non-significant)"
  )
  temporal_colours <- c(
    loser = "#B2182B", non_significant_loser = "#EF8A62",
    non_significant_winner = "#7FBF7B", winner = "#006837",
    colonizer = "#542788", non_significant_colonizer = "#B2ABD2"
  )
  comparison_levels <- c("First 2 → First 4", "First 4 → First 6", "First 2 → First 6")
  comparison_labels <- c("First 2 to First 4", "First 4 to First 6", "First 2 to First 6")
  matrix_data <- transitions |>
    dplyr::filter(
      .data$from_category != "no_change",
      .data$to_category != "no_change"
    ) |>
    dplyr::mutate(
      comparison = factor(.data$comparison, levels = comparison_levels, labels = comparison_labels),
      from_category_raw = .data$from_category,
      to_category_raw = .data$to_category,
      from_index = match(.data$from_category_raw, rev(category_levels)),
      to_index = match(.data$to_category_raw, category_levels),
      from_pool_display = dplyr::if_else(
        .data$from_category_raw %in% c("colonizer", "non_significant_colonizer"),
        "colonizer", "resident"
      ),
      to_pool_display = dplyr::if_else(
        .data$to_category_raw %in% c("colonizer", "non_significant_colonizer"),
        "colonizer", "resident"
      ),
      valid_transition = .data$from_pool_display == .data$to_pool_display,
      tile_category = dplyr::if_else(
        .data$valid_transition, .data$to_category_raw, "structurally_not_applicable"
      ),
      label = dplyr::case_when(
        !.data$valid_transition ~ "",
        .data$n == 0L ~ "0%",
        TRUE ~ sprintf("%d\n%.1f%%", .data$n, 100 * .data$row_proportion)
      ),
      from_response = dplyr::case_when(
        .data$from_category_raw %in% c("loser", "non_significant_loser") ~ "loser",
        .data$from_category_raw %in% c("winner", "non_significant_winner") ~ "winner",
        TRUE ~ "colonizer"
      ),
      to_response = dplyr::case_when(
        .data$to_category_raw %in% c("loser", "non_significant_loser") ~ "loser",
        .data$to_category_raw %in% c("winner", "non_significant_winner") ~ "winner",
        TRUE ~ "colonizer"
      ),
      reversal = .data$valid_transition &
        .data$from_response %in% c("loser", "winner") &
        .data$to_response %in% c("loser", "winner") &
        .data$from_response != .data$to_response,
      dark_tile = .data$valid_transition &
        .data$to_category_raw %in% c("loser", "winner", "colonizer")
    )
  x_axis_labels <- tidyr::crossing(
    comparison = factor(comparison_labels, levels = comparison_labels),
    category = category_levels
  ) |>
    dplyr::mutate(
      x = match(.data$category, category_levels), y = 0.32,
      label = unname(category_labels[.data$category])
    )
  y_axis_labels <- tibble::tibble(
    comparison = factor(comparison_labels[[1]], levels = comparison_labels),
    category = category_levels
  ) |>
    dplyr::mutate(
      x = 0.24, y = match(.data$category, rev(category_levels)),
      label = unname(category_labels[.data$category])
    )
  matrix_plot <- ggplot2::ggplot(
    matrix_data,
    ggplot2::aes(.data$to_index, .data$from_index)
  ) +
    ggplot2::geom_tile(
      ggplot2::aes(fill = .data$tile_category),
      colour = "grey72", linewidth = 0.4
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(matrix_data, !.data$reversal, !.data$dark_tile),
      ggplot2::aes(label = .data$label),
      colour = "black", size = 3.8, lineheight = 0.92
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(matrix_data, !.data$reversal, .data$dark_tile),
      ggplot2::aes(label = .data$label),
      colour = "white", size = 3.8, lineheight = 0.92
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(matrix_data, .data$reversal, !.data$dark_tile),
      ggplot2::aes(label = .data$label),
      colour = "black", size = 4.6, lineheight = 0.92, fontface = "bold"
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(matrix_data, .data$reversal, .data$dark_tile),
      ggplot2::aes(label = .data$label),
      colour = "white", size = 4.6, lineheight = 0.92, fontface = "bold"
    ) +
    ggplot2::geom_text(
      data = x_axis_labels,
      ggplot2::aes(.data$x, .data$y, label = .data$label, colour = .data$category),
      inherit.aes = FALSE, angle = 45, hjust = 1, vjust = 1, size = 4.0,
      lineheight = 0.9
    ) +
    ggplot2::geom_text(
      data = y_axis_labels,
      ggplot2::aes(.data$x, .data$y, label = .data$label, colour = .data$category),
      inherit.aes = FALSE, hjust = 1, size = 4.0, lineheight = 0.9
    ) +
    ggplot2::facet_wrap(~ comparison, nrow = 1, drop = FALSE) +
    ggplot2::scale_fill_manual(
      values = c(temporal_colours, structurally_not_applicable = "grey80"),
      drop = FALSE
    ) +
    ggplot2::scale_colour_manual(values = temporal_colours, guide = "none") +
    ggplot2::scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = NULL, expand = c(0, 0)) +
    ggplot2::coord_fixed(xlim = c(0.5, 6.5), ylim = c(0.5, 6.5), clip = "off") +
    transplant_figure_theme() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 15),
      strip.text = ggplot2::element_text(face = "bold", size = 15),
      panel.grid = ggplot2::element_blank(),
      legend.position = "none",
      panel.spacing.x = grid::unit(1.6, "lines"),
      plot.margin = ggplot2::margin(8, 8, 76, 104)
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL
    )

  required_composition_columns <- c(
    "year", "category", "proportion", "total_species_experiment",
    "n_species", "n_experiments"
  )
  missing_composition_columns <- setdiff(required_composition_columns, names(composition))
  if (length(missing_composition_columns)) {
    stop(
      "Temporal composition is missing columns: ",
      paste(missing_composition_columns, collapse = ", "),
      call. = FALSE
    )
  }
  composition_data <- composition |>
    dplyr::filter(.data$category != "no_change") |>
    dplyr::mutate(
      category = factor(.data$category, levels = category_levels),
      percentage = 100 * .data$proportion,
      significance = factor(
        dplyr::if_else(
          as.character(.data$category) %in% c("loser", "winner", "colonizer"),
          "Significant", "Non-significant"
        ),
        levels = c("Significant", "Non-significant")
      )
    )
  composition_years <- sort(unique(composition_data$year))
  year_metadata <- composition |>
    dplyr::distinct(
      .data$year, .data$total_species_experiment,
      .data$n_species, .data$n_experiments
    ) |>
    dplyr::arrange(.data$year)
  if (anyDuplicated(year_metadata$year)) {
    stop("Temporal sample-size metadata is inconsistent within at least one year.", call. = FALSE)
  }
  year_labels <- stats::setNames(
    sprintf(
      "%d\nC=%s  S=%s  E=%s",
      year_metadata$year,
      year_metadata$total_species_experiment,
      year_metadata$n_species,
      year_metadata$n_experiments
    ),
    year_metadata$year
  )
  composition_plot <- ggplot2::ggplot(
    composition_data,
    ggplot2::aes(.data$year, .data$percentage, colour = .data$category, group = .data$category)
  ) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_point(size = 3.5) +
    ggplot2::facet_grid(ggplot2::vars(.data$significance), scales = "free_y") +
    ggplot2::scale_colour_manual(values = temporal_colours, labels = category_labels, drop = FALSE) +
    ggplot2::scale_x_continuous(
      breaks = composition_years,
      labels = year_labels[as.character(composition_years)],
      limits = range(composition_years)
    ) +
    ggplot2::scale_y_continuous(
      labels = function(x) paste0(round(x), "%"),
      expand = ggplot2::expansion(mult = c(0.06, 0.08))
    ) +
    transplant_figure_theme() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 13),
      axis.text.x = ggplot2::element_text(lineheight = 0.95),
      axis.title = ggplot2::element_text(size = 15),
      strip.text.y = ggplot2::element_text(face = "bold", size = 14),
      legend.position = "none"
    ) +
    ggplot2::labs(
      x = "Year\n(C = combinations; S = distinct species; E = distinct experiments)",
      y = "Percentage of species × experiment combinations",
      colour = NULL
    )

  matrix_with_y_title <- patchwork::wrap_plots(
    patchwork::wrap_elements(
      full = grid::textGrob("Starting category", rot = 90, gp = grid::gpar(fontsize = 16))
    ),
    matrix_plot,
    nrow = 1,
    widths = c(0.035, 1)
  )
  matrix_labeled <- patchwork::wrap_plots(
    matrix_with_y_title,
    patchwork::wrap_elements(
      full = grid::textGrob("Later category", gp = grid::gpar(fontsize = 16))
    ),
    ncol = 1,
    heights = c(1, 0.045)
  )
  combined <- patchwork::wrap_plots(
    composition_plot, matrix_labeled,
    ncol = 1,
    heights = c(1.0, 1.25)
  )
  save_plot_atomic(combined, path, 17, 13)
}
