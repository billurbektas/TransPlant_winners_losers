category_palette <- c(
  loser = "#E66101",
  non_significant_loser = "#FDBB84",
  no_change = "#D9D9D9",
  non_significant_winner = "#99D8C9",
  winner = "#1B9E77",
  colonizer = "#5E3C99",
  non_significant_colonizer = "#B7A6D8"
)

category_display_labels <- c(
  loser = "Loser",
  winner = "Winner",
  colonizer = "Colonizer",
  non_significant_loser = "Loser\n(non-significant)",
  non_significant_winner = "Winner\n(non-significant)",
  non_significant_colonizer = "Colonizer\n(non-significant)",
  no_change = "No significant\nchange"
)

clean_figure_label <- function(x) {
  stringr::str_replace_all(as.character(x), "_", " ")
}

transplant_figure_theme <- function(base_size = 15) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 2),
      strip.text = ggplot2::element_text(size = base_size - 1, face = "bold"),
      legend.text = ggplot2::element_text(size = base_size - 2),
      legend.title = ggplot2::element_text(size = base_size - 1)
    )
}

save_plot_atomic <- function(plot, path, width, height, dpi = 300) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  extension <- tools::file_ext(path)
  temporary <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path), fileext = paste0(".", extension))
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  ggplot2::ggsave(temporary, plot = plot, width = width, height = height, dpi = dpi)
  if (!file.rename(temporary, path)) stop("Could not save figure atomically: ", path, call. = FALSE)
  receipt <- paste0(path, ".generated")
  receipt_temporary <- tempfile(pattern = paste0(basename(receipt), "."), tmpdir = dirname(path))
  writeLines(format(Sys.time(), tz = "UTC", usetz = TRUE), receipt_temporary)
  if (!file.rename(receipt_temporary, receipt)) {
    stop("Could not save figure-generation receipt: ", receipt, call. = FALSE)
  }
  path
}

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

experiment_category_proportions <- function(data, categories = unique(data$category)) {
  data |>
    dplyr::filter(.data$category %in% .env$categories) |>
    dplyr::count(.data$experiment_id, .data$category, name = "n") |>
    dplyr::group_by(.data$experiment_id) |>
    dplyr::mutate(proportion = .data$n / sum(.data$n)) |>
    dplyr::ungroup()
}

category_density_data_xy <- function(data, x_column, y_column, grid_size = 120L) {
  categories <- c("loser", "winner", "colonizer")
  coordinates <- data |>
    dplyr::transmute(
      experiment_id = .data$experiment_id,
      climate_x = .data[[x_column]],
      climate_y = .data[[y_column]]
    ) |>
    dplyr::distinct()
  coordinate_counts <- coordinates |>
    dplyr::count(.data$experiment_id, name = "n_coordinates")
  if (any(coordinate_counts$n_coordinates != 1L)) {
    stop(
      "Each experiment must have exactly one pair of climate coordinates for ",
      x_column, " and ", y_column, ".",
      call. = FALSE
    )
  }
  coordinates <- coordinates |>
    dplyr::filter(
      is.finite(.data$climate_x),
      is.finite(.data$climate_y)
    )
  counts <- data |>
    dplyr::filter(.data$category %in% .env$categories) |>
    dplyr::count(.data$experiment_id, .data$category, name = "n") |>
    tidyr::complete(
      experiment_id = coordinates$experiment_id,
      category = categories,
      fill = list(n = 0L)
    ) |>
    dplyr::left_join(coordinates, by = "experiment_id")
  count_check <- counts |>
    dplyr::summarise(total = sum(.data$n), .by = "experiment_id")
  if (any(!is.finite(count_check$total)) || any(count_check$total <= 0)) {
    stop("Every experiment must have at least one significant response.", call. = FALSE)
  }
  x_range <- range(coordinates$climate_x)
  y_range <- range(coordinates$climate_y)
  x_grid <- seq(x_range[[1]], x_range[[2]], length.out = grid_size)
  y_grid <- seq(y_range[[1]], y_range[[2]], length.out = grid_size)
  bandwidth_x <- stats::bw.nrd0(coordinates$climate_x)
  bandwidth_y <- stats::bw.nrd0(coordinates$climate_y)
  if (!is.finite(bandwidth_x) || bandwidth_x <= 0 ||
      !is.finite(bandwidth_y) || bandwidth_y <= 0) {
    stop("Destination climate bandwidths must be positive and finite.", call. = FALSE)
  }
  density_grid <- purrr::map_dfr(categories, function(focal_category) {
    category_data <- counts |>
      dplyr::filter(.data$category == .env$focal_category, .data$n > 0)
    grid <- tidyr::expand_grid(
      climate_x = x_grid,
      climate_y = y_grid
    )
    density <- numeric(nrow(grid))
    for (index in seq_len(nrow(category_data))) {
      density <- density + category_data$n[[index]] *
        stats::dnorm(
          grid$climate_x,
          mean = category_data$climate_x[[index]], sd = bandwidth_x
        ) *
        stats::dnorm(
          grid$climate_y,
          mean = category_data$climate_y[[index]], sd = bandwidth_y
        )
    }
    dplyr::mutate(grid, category = focal_category, density = density)
  })
  maximum_density <- max(density_grid$density)
  density_grid <- density_grid |>
    dplyr::mutate(
      relative_density = .data$density / maximum_density,
      relative_density = dplyr::if_else(.data$relative_density < 0.01, 0, .data$relative_density)
    )
  list(
    density = density_grid,
    counts = counts,
    proportions = counts,
    bandwidth = c(temperature = bandwidth_x, precipitation = bandwidth_y)
  )
}

destination_category_density_data <- function(data, grid_size = 120L) {
  density_data <- category_density_data_xy(
    data, "destination_T_JJA", "destination_P_JJA", grid_size
  )
  density_data$density <- density_data$density |>
    dplyr::rename(
      destination_T_JJA = "climate_x",
      destination_P_JJA = "climate_y"
    )
  density_data$counts <- density_data$counts |>
    dplyr::rename(
      destination_T_JJA = "climate_x",
      destination_P_JJA = "climate_y"
    )
  density_data$proportions <- density_data$counts
  density_data
}




climate_response_surface_plot <- function(
    predictions, estimands, estimand_labels, value_column = "estimate",
    fill_title = "Predicted probability", fill_limits = NULL,
    show_experiments = TRUE, show_contours = FALSE) {
  warming_labels <- as.character(predictions$warming_levels$warming_label)
  surface_data <- predictions$grid |>
    dplyr::filter(.data$estimand %in% .env$estimands) |>
    dplyr::mutate(
      warming_label = factor(
        as.character(.data$warming_label), levels = .env$warming_labels
      ),
      estimand_label = factor(
        unname(.env$estimand_labels[.data$estimand]),
        levels = unname(.env$estimand_labels[.env$estimands])
      ),
      display_value = dplyr::if_else(
        .data$supported, .data[[value_column]], NA_real_
      )
    )
  if (is.null(fill_limits)) {
    fill_limits <- range(surface_data$display_value, na.rm = TRUE)
  }
  if (!all(is.finite(fill_limits)) || diff(fill_limits) <= 0) {
    stop("Climate-response surface limits must span a finite range.", call. = FALSE)
  }
  cutoffs <- surface_data |>
    dplyr::distinct(
      .data$temperature_cutoff_low, .data$temperature_cutoff_high,
      .data$precipitation_cutoff_dry, .data$precipitation_cutoff_wet
    )
  if (nrow(cutoffs) != 1L) stop("Prediction cutoffs must be unique.", call. = FALSE)
  experiments <- predictions$experiment_support |>
    dplyr::mutate(
      warming_label = factor(
        as.character(.data$warming_label), levels = .env$warming_labels
      )
    )
  plot <- ggplot2::ggplot(surface_data) +
    ggplot2::geom_raster(
      ggplot2::aes(
        .data$destination_T_JJA, .data$destination_P_JJA,
        fill = .data$display_value
      ),
      interpolate = TRUE, na.rm = TRUE
    )
  if (show_contours) {
    plot <- plot + ggplot2::geom_contour(
      ggplot2::aes(
        .data$destination_T_JJA, .data$destination_P_JJA,
        z = .data$display_value
      ),
      bins = 5, colour = scales::alpha("white", 0.75), linewidth = 0.32,
      na.rm = TRUE, show.legend = FALSE
    )
  }
  plot <- plot +
    ggplot2::geom_vline(
      xintercept = c(
        cutoffs$temperature_cutoff_low, cutoffs$temperature_cutoff_high
      ),
      colour = "grey30", linetype = "dashed", linewidth = 0.35
    ) +
    ggplot2::geom_hline(
      yintercept = c(
        cutoffs$precipitation_cutoff_dry, cutoffs$precipitation_cutoff_wet
      ),
      colour = "grey30", linetype = "dashed", linewidth = 0.35
    )
  if (show_experiments) {
    plot <- plot + ggplot2::geom_point(
      data = experiments,
      ggplot2::aes(.data$destination_T_JJA, .data$destination_P_JJA),
      inherit.aes = FALSE, shape = 21, size = 1.25, stroke = 0.35,
      fill = scales::alpha("white", 0.75), colour = "grey20"
    )
  }
  plot +
    ggplot2::facet_grid(
      rows = ggplot2::vars(.data$estimand_label),
      cols = ggplot2::vars(.data$warming_label), drop = FALSE
    ) +
    ggplot2::scale_fill_viridis_c(
      option = "C", limits = fill_limits, oob = scales::squish,
      na.value = "white", name = fill_title
    ) +
    transplant_figure_theme(base_size = 13) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 12),
      legend.position = "bottom",
      axis.text = ggplot2::element_text(size = 10)
    ) +
    ggplot2::labs(
      x = "Destination JJA mean temperature (°C)",
      y = "Destination JJA precipitation (mm)"
    )
}




# =============================================================================
# Pool-specific climate-response figures
# =============================================================================

pool_climate_term_labels <- c(
  z_destination_T = "Destination temperature",
  z_log_destination_P = "Destination precipitation",
  z_warming = "Warming",
  `z_destination_T:z_log_destination_P` = "Temperature × precipitation",
  `z_destination_T:z_warming` = "Temperature × warming",
  `z_log_destination_P:z_warming` = "Precipitation × warming"
)

pool_climate_category_labels <- c(
  loser = "Loser | resident pool",
  winner = "Winner | resident pool",
  colonizer = "Colonizer | nonresident pool"
)

pool_review_odds_ratio_transform <- function() {
  odds_anchors <- c(0.25, 0.5, 1, 1.5, 3, 4)
  display_anchors <- seq_along(odds_anchors) - 2
  piecewise_linear <- function(values, source, target) {
    result <- rep(NA_real_, length(values))
    valid <- !is.na(values)
    interval <- findInterval(values[valid], source, all.inside = TRUE)
    interval <- pmin(interval, length(source) - 1L)
    slope <- (target[interval + 1L] - target[interval]) /
      (source[interval + 1L] - source[interval])
    result[valid] <- target[interval] +
      slope * (values[valid] - source[interval])
    result
  }
  forward <- function(values) {
    piecewise_linear(values, odds_anchors, display_anchors)
  }
  inverse <- function(values) {
    piecewise_linear(values, display_anchors, odds_anchors)
  }
  scales::new_transform(
    name = "review_odds_ratio",
    transform = forward,
    inverse = inverse,
    domain = c(0.01, Inf)
  )
}

pool_climate_odds_ratio_plot <- function(coefficients, facet_ncol = 1L,
                                         sensitivity = FALSE,
                                         review_axis = FALSE) {
  plot_data <- if ("response" %in% names(coefficients)) {
    coefficients |> dplyr::mutate(category = .data$response)
  } else {
    coefficients |>
      dplyr::mutate(category = .data$category_vs_non_significant)
  }
  plot_data <- plot_data |>
    dplyr::filter(.data$term != "(Intercept)") |>
    dplyr::mutate(
      term_label = factor(
        unname(pool_climate_term_labels[.data$term]),
        levels = rev(unname(pool_climate_term_labels))
      ),
      category_label = factor(
        unname(pool_climate_category_labels[.data$category]),
        levels = unname(pool_climate_category_labels)
      )
    )
  if (sensitivity) {
    model_lookup <- plot_data |>
      dplyr::distinct(.data$category, .data$selected_model) |>
      dplyr::mutate(
        model_label = dplyr::if_else(
          .data$selected_model == "species_random_binomial",
          "species random intercept", "fixed fallback"
        ),
        category_label = factor(
          paste0(
            unname(pool_climate_category_labels[.data$category]),
            "\n", .data$model_label
          ),
          levels = paste0(
            unname(pool_climate_category_labels), "\n",
            dplyr::if_else(
              unname(pool_climate_category_labels) ==
                pool_climate_category_labels[["colonizer"]],
              "fixed fallback", "species random intercept"
            )
          )
        )
      )
    plot_data <- plot_data |>
      dplyr::select(-"category_label") |>
      dplyr::left_join(
        dplyr::select(model_lookup, "category", "category_label"),
        by = "category"
      )
  }
  colours <- c(
    loser = category_palette[["loser"]],
    winner = category_palette[["winner"]],
    colonizer = category_palette[["colonizer"]]
  )
  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(.data$odds_ratio, .data$term_label, colour = .data$category)
  ) +
    ggplot2::geom_vline(xintercept = 1, colour = "grey55", linewidth = 0.65) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$odds_low, xmax = .data$odds_high),
      orientation = "y", width = 0, linewidth = 0.8
    ) +
    ggplot2::geom_point(
      ggplot2::aes(shape = .data$raw_significant),
      size = 3.5, stroke = 1.1
    ) +
    ggplot2::facet_wrap(~ category_label, ncol = facet_ncol, scales = "free_y") +
    ggplot2::scale_colour_manual(values = colours, guide = "none") +
    ggplot2::scale_shape_manual(
      values = c(`FALSE` = 1, `TRUE` = 16),
      labels = c(`FALSE` = "p >= 0.05", `TRUE` = "p < 0.05"), name = NULL
    ) +
    transplant_figure_theme(base_size = 15) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 14),
      legend.position = "bottom",
      panel.spacing = grid::unit(0.45, "cm"),
      plot.margin = ggplot2::margin(12, 8, 8, 8)
    ) +
    ggplot2::labs(
      x = "Odds ratio versus pool-specific non-significant response", y = NULL
    )
  if (review_axis) {
    plot + ggplot2::scale_x_continuous(
      transform = pool_review_odds_ratio_transform(),
      breaks = c(0.5, 1, 1.5, 3),
      labels = c("0.5", "1", "1.5", "3"),
      minor_breaks = NULL
    )
  } else {
    plot + ggplot2::scale_x_log10()
  }
}


plot_pool_binary_sensitivity_odds_ratios <- function(coefficients, path) {
  save_plot_atomic(
    pool_climate_odds_ratio_plot(coefficients, sensitivity = TRUE),
    path, 10, 10
  )
}

pool_climate_surface_row <- function(
    predictions, estimand, value_column, fill_limits, fill_title,
    low_colour, high_colour, show_x = TRUE) {
  warming_labels <- as.character(predictions$warming_levels$warming_label)
  label <- unname(pool_climate_category_labels[[estimand]])
  plot_data <- predictions$grid |>
    dplyr::filter(.data$estimand == .env$estimand) |>
    dplyr::mutate(
      warming_label = factor(
        as.character(.data$warming_label), levels = .env$warming_labels
      ),
      display_value = dplyr::if_else(
        .data$supported, .data[[value_column]], NA_real_
      ),
      category_label = factor(.env$label, levels = .env$label)
    )
  cutoffs <- plot_data |>
    dplyr::distinct(
      .data$temperature_cutoff_low, .data$temperature_cutoff_high,
      .data$precipitation_cutoff_dry, .data$precipitation_cutoff_wet
    )
  experiments <- predictions$experiment_support |>
    dplyr::mutate(
      warming_label = factor(
        as.character(.data$warming_label), levels = .env$warming_labels
      )
    )
  ggplot2::ggplot(plot_data) +
    ggplot2::geom_raster(
      ggplot2::aes(
        .data$destination_T_JJA, .data$destination_P_JJA,
        fill = .data$display_value
      ),
      interpolate = TRUE, na.rm = TRUE
    ) +
    ggplot2::geom_vline(
      xintercept = c(
        cutoffs$temperature_cutoff_low, cutoffs$temperature_cutoff_high
      ),
      colour = "grey35", linetype = "dashed", linewidth = 0.35
    ) +
    ggplot2::geom_hline(
      yintercept = c(
        cutoffs$precipitation_cutoff_dry, cutoffs$precipitation_cutoff_wet
      ),
      colour = "grey35", linetype = "dashed", linewidth = 0.35
    ) +
    ggplot2::geom_point(
      data = experiments,
      ggplot2::aes(.data$destination_T_JJA, .data$destination_P_JJA),
      inherit.aes = FALSE, shape = 21, size = 1.25, stroke = 0.35,
      fill = scales::alpha("white", 0.75), colour = "grey20"
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(.data$category_label),
      cols = ggplot2::vars(.data$warming_label), drop = FALSE
    ) +
    ggplot2::scale_fill_gradient(
      low = low_colour, high = high_colour, limits = fill_limits,
      oob = scales::squish, na.value = "white", name = fill_title
    ) +
    transplant_figure_theme(base_size = 13) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 12),
      legend.position = "bottom",
      axis.text = ggplot2::element_text(size = 10),
      axis.title.x = ggplot2::element_text(
        colour = if (show_x) "black" else NA
      ),
      axis.text.x = ggplot2::element_text(
        colour = if (show_x) "black" else "transparent"
      ),
      axis.ticks.x = ggplot2::element_line(
        colour = if (show_x) "black" else "transparent"
      )
    ) +
    ggplot2::labs(
      x = if (show_x) "Destination JJA mean temperature (°C)" else NULL,
      y = "Destination JJA precipitation (mm)"
    )
}

pool_climate_probability_plot <- function(
    predictions, value_column = "estimate", fill_title = "Probability") {
  estimands <- c("loser", "winner", "colonizer")
  supported <- predictions$grid |>
    dplyr::filter(.data$supported, .data$estimand %in% .env$estimands)
  fill_limits <- range(supported[[value_column]], na.rm = TRUE)
  colours <- list(
    loser = c("#FEE8D0", category_palette[["loser"]]),
    winner = c("#D9F0E8", category_palette[["winner"]]),
    colonizer = c("#E5DDF2", category_palette[["colonizer"]])
  )
  rows <- purrr::imap(estimands, function(estimand, index) {
    pool_climate_surface_row(
      predictions, estimand, value_column, fill_limits, fill_title,
      low_colour = colours[[estimand]][[1]],
      high_colour = colours[[estimand]][[2]],
      show_x = index == length(estimands)
    )
  })
  patchwork::wrap_plots(rows, ncol = 1)
}



pool_climate_review_density_plot <- function(data, climate_contexts, predictions) {
  density_data <- category_density_data_xy(
    data, "destination_T_JJA", "destination_P_JJA"
  )
  categories <- c("loser", "winner", "colonizer")
  density_grid <- density_data$density |>
    dplyr::mutate(category = factor(.data$category, levels = categories)) |>
    dplyr::group_by(.data$category) |>
    dplyr::mutate(within_category_density = .data$density / max(.data$density)) |>
    dplyr::ungroup()
  counts <- density_data$counts |>
    dplyr::filter(.data$n > 0) |>
    dplyr::left_join(
      predictions$experiment_support |>
        dplyr::transmute(
          experiment_id = as.character(.data$experiment_id),
          warming = .data$warming
        ),
      by = "experiment_id"
    ) |>
    dplyr::mutate(category = factor(.data$category, levels = categories))
  experiments <- predictions$experiment_support |>
    dplyr::transmute(
      experiment_id = as.character(.data$experiment_id),
      climate_x = .data$destination_T_JJA,
      climate_y = .data$destination_P_JJA,
      warming = .data$warming
    )
  temperature_edges <- sort(unique(c(
    climate_contexts$temperature_cutoff_low,
    climate_contexts$temperature_cutoff_high
  )))
  precipitation_edges <- sort(unique(c(
    climate_contexts$precipitation_cutoff_dry,
    climate_contexts$precipitation_cutoff_wet
  )))
  density_palette <- category_palette[c("loser", "winner", "colonizer")]
  ggplot2::ggplot() +
    purrr::map(categories, function(category) {
      ggplot2::geom_raster(
        data = dplyr::filter(
          density_grid, as.character(.data$category) == .env$category
        ),
        ggplot2::aes(
          .data$climate_x, .data$climate_y, alpha = .data$relative_density
        ),
        fill = unname(density_palette[[category]]), interpolate = TRUE
      )
    }) +
    ggplot2::geom_contour(
      data = density_grid,
      ggplot2::aes(
        .data$climate_x, .data$climate_y,
        z = .data$within_category_density, colour = .data$category
      ),
      breaks = c(0.25, 0.5, 0.75), linewidth = 0.35,
      alpha = 0.7, show.legend = FALSE
    ) +
    ggplot2::geom_vline(
      xintercept = temperature_edges,
      colour = "grey45", linetype = "dashed", linewidth = 0.4
    ) +
    ggplot2::geom_hline(
      yintercept = precipitation_edges,
      colour = "grey45", linetype = "dashed", linewidth = 0.4
    ) +
    ggplot2::geom_point(
      data = experiments,
      ggplot2::aes(.data$climate_x, .data$climate_y),
      shape = 21, size = 1.5, stroke = 0.35,
      fill = scales::alpha("white", 0.5), colour = "grey55"
    ) +
    ggplot2::geom_point(
      data = counts,
      ggplot2::aes(
        .data$climate_x, .data$climate_y,
        size = .data$n, fill = .data$warming
      ),
      shape = 21, colour = "white", stroke = 0.55, alpha = 0.95
    ) +
    ggplot2::facet_wrap(
      ~ category, ncol = 1,
      labeller = ggplot2::as_labeller(c(
        loser = "Observed significant losers",
        winner = "Observed significant winners",
        colonizer = "Observed significant colonizers"
      )), drop = FALSE
    ) +
    ggplot2::scale_fill_gradientn(
      colours = c("#FEE5D9", "#FCAE91", "#FB6A4A", "#CB181D", "#67000D"),
      name = "Experimental warming (°C)"
    ) +
    ggplot2::scale_colour_manual(values = density_palette, guide = "none") +
    ggplot2::scale_alpha_continuous(
      limits = c(0, 1), range = c(0, 0.9), guide = "none"
    ) +
    ggplot2::scale_size_area(max_size = 6, name = "Species × experiment count") +
    ggplot2::guides(
      fill = ggplot2::guide_colourbar(
        title.position = "top", barwidth = grid::unit(4.5, "cm")
      ),
      size = ggplot2::guide_legend(title.position = "top", nrow = 1)
    ) +
    transplant_figure_theme(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom", legend.box = "vertical",
      legend.title = ggplot2::element_text(hjust = 0.5),
      strip.text = ggplot2::element_text(size = 11)
    ) +
    ggplot2::labs(
      x = "Destination JJA mean temperature (°C)",
      y = "Destination JJA precipitation (mm)"
    )
}

pool_climate_review_probability_plot <- function(
    predictions, value_column = "estimate",
    fill_title = "Probability within pool") {
  estimands <- c("loser", "winner", "colonizer")
  retained_warming <- predictions$warming_levels |>
    dplyr::filter(.data$warming_quantile %in% c(0.25, 0.75)) |>
    dplyr::arrange(.data$warming_quantile)
  warming_labels <- as.character(retained_warming$warming_level)
  category_labels <- unname(pool_climate_category_labels[estimands])
  plot_data <- predictions$grid |>
    dplyr::filter(
      .data$estimand %in% .env$estimands,
      .data$warming_quantile %in% c(0.25, 0.75)
    ) |>
    dplyr::mutate(
      warming_label = factor(
        as.character(.data$warming_level), levels = .env$warming_labels
      ),
      category_label = factor(
        unname(pool_climate_category_labels[.data$estimand]),
        levels = .env$category_labels
      ),
      display_value = .data[[value_column]]
    )
  fill_limits <- c(0, max(plot_data$display_value, na.rm = TRUE))
  category_gradients <- list(
    loser = c("#FEE8D0", category_palette[["loser"]]),
    winner = c("#D9F0E8", category_palette[["winner"]]),
    colonizer = c("#E5DDF2", category_palette[["colonizer"]])
  )
  colour_ramps <- purrr::map(category_gradients, scales::colour_ramp)
  plot_data <- plot_data |>
    dplyr::mutate(
      colour_position = pmax(.data$display_value, 0) / fill_limits[[2]],
      display_colour = purrr::map2_chr(
        as.character(.data$estimand), .data$colour_position,
        function(category, position) colour_ramps[[category]](position)
      )
    )
  cutoffs <- plot_data |>
    dplyr::distinct(
      .data$temperature_cutoff_low, .data$temperature_cutoff_high,
      .data$precipitation_cutoff_dry, .data$precipitation_cutoff_wet
    )
  experiments <- tidyr::crossing(
    predictions$experiment_support |>
      dplyr::select(
        "experiment_id", "destination_T_JJA", "destination_P_JJA"
      ) |>
      dplyr::distinct(),
    warming_label = factor(warming_labels, levels = warming_labels),
    category_label = factor(category_labels, levels = category_labels)
  )
  surface_plot <- ggplot2::ggplot(plot_data) +
    ggplot2::geom_raster(
      ggplot2::aes(
        .data$destination_T_JJA, .data$destination_P_JJA,
        fill = .data$display_colour
      ),
      interpolate = TRUE
    ) +
    ggplot2::geom_blank(
      ggplot2::aes(colour = .data$display_value),
      show.legend = TRUE
    ) +
    ggplot2::geom_vline(
      xintercept = c(
        cutoffs$temperature_cutoff_low, cutoffs$temperature_cutoff_high
      ),
      colour = "grey45", linetype = "dashed", linewidth = 0.35
    ) +
    ggplot2::geom_hline(
      yintercept = c(
        cutoffs$precipitation_cutoff_dry, cutoffs$precipitation_cutoff_wet
      ),
      colour = "grey45", linetype = "dashed", linewidth = 0.35
    ) +
    ggplot2::geom_point(
      data = experiments,
      ggplot2::aes(.data$destination_T_JJA, .data$destination_P_JJA),
      inherit.aes = FALSE, shape = 21, size = 2.2, stroke = 0.5,
      fill = scales::alpha("white", 0.82), colour = "grey10"
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(.data$category_label),
      cols = ggplot2::vars(.data$warming_label), drop = FALSE
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_colour_gradient(
      low = "white", high = "grey15", limits = fill_limits,
      name = fill_title
    ) +
    transplant_figure_theme(base_size = 13) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 12),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(hjust = 0.5),
      axis.text = ggplot2::element_text(size = 10)
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_colourbar(
        title.position = "top", barwidth = grid::unit(7, "cm")
      )
    ) +
    ggplot2::labs(
      x = "Destination JJA mean temperature (°C)",
      y = "Destination JJA precipitation (mm)"
    )
  surface_plot
}

plot_pool_climate_response_review_prototype <- function(
    data, climate_contexts, coefficients, predictions, path) {
  density_plot <- pool_climate_review_density_plot(
    data, climate_contexts, predictions
  )
  odds_plot <- pool_climate_odds_ratio_plot(
    coefficients, review_axis = TRUE
  ) +
    ggplot2::labs(x = "Odds ratio")
  probability_plot <- pool_climate_review_probability_plot(
    predictions, value_column = "estimate",
    fill_title = "Probability within pool"
  )
  combined <- patchwork::wrap_plots(
    density_plot, odds_plot, probability_plot,
    nrow = 1, widths = c(1.2, 0.95, 2.0)
  )
  save_plot_atomic(combined, path, 26, 14)
}



destination_climate_experiment_map_data <- function(data) {
  required <- c(
    "experiment_id", "destination_T_JJA", "destination_P_JJA",
    "diff_T_JJA", "diff_P_JJA"
  )
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(
      "Experiment climate map is missing columns: ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  if (!"climate_summary_status" %in% names(data)) {
    data$climate_summary_status <- "complete"
  }
  plot_data <- data |>
    dplyr::select(dplyr::all_of(c(required, "climate_summary_status"))) |>
    dplyr::distinct()
  duplicated_experiments <- plot_data |>
    dplyr::count(.data$experiment_id) |>
    dplyr::filter(.data$n != 1L)
  if (nrow(duplicated_experiments)) {
    stop("Every experiment must have one unique climate-map record.", call. = FALSE)
  }
  climate_columns <- c(
    "destination_T_JJA", "destination_P_JJA", "diff_T_JJA", "diff_P_JJA"
  )
  if (any(!is.finite(as.matrix(plot_data[, climate_columns])))) {
    stop("Every experiment requires finite destination and change climate values.", call. = FALSE)
  }
  plot_data |>
    dplyr::mutate(
      experiment_label = sub(
        "^([^:]+)::([^:]+)::(.+)$", "\\1\n\\2 -> \\3", .data$experiment_id
      ),
      experiment_label = dplyr::if_else(
        .data$climate_summary_status == "partial",
        paste0(.data$experiment_label, "*"),
        .data$experiment_label
      ),
      experiment_label = clean_figure_label(.data$experiment_label)
    )
}




destination_climate_display_context <- function(experiment_id, grid_context) {
  if (length(experiment_id) != length(grid_context)) {
    stop("Experiment IDs and equal-width grid contexts must have equal length.", call. = FALSE)
  }
  dplyr::case_when(
    experiment_id == "NO_Gudmedalen::Rambera::Arhelleren" ~
      "cold_intermediate",
    experiment_id == "NO_Skjellingahaugen::Veskre::Ovstedal" ~
      "cold_wet",
    stringr::str_starts(
      experiment_id, stringr::regex("DE_(Grainau|Susalps)::")
    ) ~ "warm_dry",
    grid_context %in% c(
      "cold_dry", "intermediate_dry", "warm_dry",
      "cold_intermediate", "cold_wet"
    ) ~ grid_context,
    TRUE ~ NA_character_
  )
}

plot_destination_climate_warming_gradient <- function(data, climate_contexts, path) {
  plot_data <- data |>
    dplyr::distinct(
      .data$experiment_id, .data$origin_T_JJA, .data$destination_T_JJA,
      .data$destination_P_JJA, .data$diff_T_JJA, .data$climate_summary_status
    ) |>
    dplyr::filter(
      is.finite(.data$origin_T_JJA),
      is.finite(.data$destination_T_JJA),
      is.finite(.data$destination_P_JJA),
      is.finite(.data$diff_T_JJA)
    ) |>
    dplyr::mutate(
      experiment_label = dplyr::if_else(
        .data$climate_summary_status == "partial",
        paste0(.data$experiment_id, "*"),
        .data$experiment_id
      ),
      experiment_label = sub(
        "^([^:]+)::([^:]+)::(.+)$", "\\1\n\\2 -> \\3", .data$experiment_label
      ),
      experiment_label = clean_figure_label(.data$experiment_label)
    )
  if (!isTRUE(all.equal(
    plot_data$diff_T_JJA,
    plot_data$destination_T_JJA - plot_data$origin_T_JJA,
    tolerance = 1e-10
  ))) {
    stop("diff_T_JJA must equal destination_T_JJA - origin_T_JJA.", call. = FALSE)
  }
  x_range <- range(plot_data$destination_T_JJA)
  y_range <- range(plot_data$destination_P_JJA)
  temperature_edges <- seq(x_range[[1]], x_range[[2]], length.out = 4L)
  precipitation_edges <- seq(y_range[[1]], y_range[[2]], length.out = 4L)
  temperature_names <- c("Cold", "Intermediate", "Warm")
  precipitation_names <- c("Dry", "Intermediate", "Wet")
  plot_data <- plot_data |>
    dplyr::mutate(
      temperature_context = cut(
        .data$destination_T_JJA,
        breaks = .env$temperature_edges,
        labels = .env$temperature_names,
        include.lowest = TRUE
      ),
      precipitation_context = cut(
        .data$destination_P_JJA,
        breaks = .env$precipitation_edges,
        labels = .env$precipitation_names,
        include.lowest = TRUE
      ),
      grid_context = paste0(
        tolower(as.character(.data$temperature_context)), "_",
        tolower(as.character(.data$precipitation_context))
      ),
      display_context = destination_climate_display_context(
        .data$experiment_id, .data$grid_context
      )
    )
  if (anyNA(plot_data$temperature_context) || anyNA(plot_data$precipitation_context)) {
    stop("Every experiment must fall inside the equal-width climate grid.", call. = FALSE)
  }
  if (anyNA(plot_data$display_context)) {
    stop(
      "An experiment in a sparse equal-width climate cell needs a reviewed display assignment.",
      call. = FALSE
    )
  }
  x_padding <- diff(x_range) * 0.03
  y_padding <- diff(y_range) * 0.03
  x_min <- x_range[[1]] - x_padding
  x_max <- x_range[[2]] + x_padding
  y_min <- y_range[[1]] - y_padding
  y_max <- y_range[[2]] + y_padding
  x_bounds <- c(x_min, temperature_edges[2:3], x_max)
  y_bounds <- c(y_min, precipitation_edges[2:3], y_max)
  context_counts <- plot_data |>
    dplyr::count(
      .data$temperature_context, .data$precipitation_context,
      name = "n_experiments"
    )
  context_rectangles <- tidyr::expand_grid(
    temperature_index = seq_along(temperature_names),
    precipitation_index = seq_along(precipitation_names)
  ) |>
    dplyr::mutate(
      temperature_context = factor(
        .env$temperature_names[.data$temperature_index],
        levels = .env$temperature_names
      ),
      precipitation_context = factor(
        .env$precipitation_names[.data$precipitation_index],
        levels = .env$precipitation_names
      ),
      xmin = .env$x_bounds[.data$temperature_index],
      xmax = .env$x_bounds[.data$temperature_index + 1L],
      ymin = .env$y_bounds[.data$precipitation_index],
      ymax = .env$y_bounds[.data$precipitation_index + 1L]
    ) |>
    dplyr::left_join(
      context_counts,
      by = c("temperature_context", "precipitation_context")
    ) |>
    dplyr::mutate(
      n_experiments = tidyr::replace_na(.data$n_experiments, 0L),
      climate_context = paste0(
        tolower(as.character(.data$temperature_context)), "_",
        tolower(as.character(.data$precipitation_context))
      ),
      display_context = dplyr::case_when(
        .data$climate_context == "cold_dry" ~ "cold_dry",
        .data$climate_context == "cold_intermediate" ~ "cold_intermediate",
        .data$climate_context == "intermediate_dry" ~ "intermediate_dry",
        .data$climate_context %in% c("warm_dry", "warm_intermediate") ~
          "warm_dry",
        .data$climate_context == "intermediate_intermediate" ~
          "cold_intermediate",
        .data$climate_context == "intermediate_wet" ~ "cold_wet",
        .data$climate_context == "cold_wet" ~ "cold_wet",
        .data$climate_context == "warm_wet" ~ "unoccupied",
        TRUE ~ NA_character_
      )
    )
  if (nrow(context_rectangles) != 9L ||
      sum(context_rectangles$n_experiments) != nrow(plot_data) ||
      any(abs(diff(temperature_edges) - diff(temperature_edges)[[1]]) > 1e-10) ||
      any(abs(diff(precipitation_edges) - diff(precipitation_edges)[[1]]) > 1e-10)) {
    stop("The equal-width 3 x 3 climate grid failed validation.", call. = FALSE)
  }
  display_counts <- plot_data |>
    dplyr::summarise(
      n_experiments = dplyr::n(),
      .by = "display_context"
    )
  context_annotations <- tibble::tribble(
    ~display_context, ~x_mid, ~y_mid, ~display_label,
    "cold_dry",
      (x_bounds[[1]] + x_bounds[[2]]) / 2,
      (y_bounds[[1]] + y_bounds[[2]]) / 2,
      "Cold/Dry",
    "intermediate_dry",
      (x_bounds[[2]] + x_bounds[[3]]) / 2,
      (y_bounds[[1]] + y_bounds[[2]]) / 2,
      "Intermediate/Dry",
    "warm_dry",
      (x_bounds[[3]] + x_bounds[[4]]) / 2,
      (y_bounds[[1]] + y_bounds[[2]]) / 2,
      "Warm/Dry",
    "cold_intermediate",
      (x_bounds[[1]] + x_bounds[[2]]) / 2,
      (y_bounds[[2]] + y_bounds[[3]]) / 2,
      "Cold/Intermediate",
    "cold_wet",
      (x_bounds[[1]] + x_bounds[[2]]) / 2,
      (y_bounds[[3]] + y_bounds[[4]]) / 2,
      "Cold/Wet",
    "unoccupied",
      (x_bounds[[3]] + x_bounds[[4]]) / 2,
      (y_bounds[[3]] + y_bounds[[4]]) / 2,
      "Unoccupied"
  ) |>
    dplyr::left_join(display_counts, by = "display_context") |>
    dplyr::mutate(
      n_experiments = tidyr::replace_na(.data$n_experiments, 0L),
      context_label = paste0(.data$display_label, "\n", .data$n_experiments)
    )
  expected_display_counts <- c(
    cold_dry = 14L, intermediate_dry = 9L, warm_dry = 11L,
    cold_intermediate = 3L, cold_wet = 5L, unoccupied = 0L
  )
  observed_display_counts <- stats::setNames(
    context_annotations$n_experiments, context_annotations$display_context
  )
  if (!identical(
    as.integer(observed_display_counts[names(expected_display_counts)]),
    as.integer(expected_display_counts)
  )) {
    stop(
      "The reviewed merged climate assignments no longer match the 42-experiment snapshot.",
      call. = FALSE
    )
  }
  context_colours <- c(
    cold_dry = "#9CCFE0", intermediate_dry = "#F7E6B5", warm_dry = "#F7C5B8",
    cold_intermediate = "#A7DDB2", cold_wet = "#5AA6C8",
    unoccupied = "#F2F2F2"
  )
  warming_colours <- c("#8C6D00", "#D8A600", "#F16913", "#BD0026", "#67000D")
  plot <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = context_rectangles,
      ggplot2::aes(
        xmin = .data$xmin, xmax = .data$xmax,
        ymin = .data$ymin, ymax = .data$ymax
      ),
      fill = unname(context_colours[context_rectangles$display_context]),
      alpha = 0.18, colour = "white", linewidth = 0.7, show.legend = FALSE
    ) +
    ggplot2::geom_vline(
      xintercept = temperature_edges[2:3],
      colour = "grey55", linetype = "dashed", linewidth = 0.45
    ) +
    ggplot2::geom_hline(
      yintercept = precipitation_edges[2:3],
      colour = "grey55", linetype = "dashed", linewidth = 0.45
    ) +
    ggplot2::geom_text(
      data = context_annotations,
      ggplot2::aes(.data$x_mid, .data$y_mid, label = .data$context_label),
      colour = "grey25", size = 4.2, fontface = "bold", lineheight = 0.9
    ) +
    ggplot2::geom_point(
      data = plot_data,
      ggplot2::aes(
        .data$destination_T_JJA, .data$destination_P_JJA,
        fill = .data$diff_T_JJA
      ),
      shape = 21, size = 5.2, stroke = 0.6, colour = "white"
    ) +
    ggrepel::geom_text_repel(
      data = plot_data,
      ggplot2::aes(
        .data$destination_T_JJA, .data$destination_P_JJA,
        label = .data$experiment_label, colour = .data$diff_T_JJA
      ),
      size = 3, max.overlaps = Inf, min.segment.length = 0,
      box.padding = 0.5, point.padding = 0.3, seed = 20260827,
      force = 4, force_pull = 0.15, max.time = 10, max.iter = 200000,
      segment.colour = "grey55", segment.size = 0.35
    ) +
    ggplot2::scale_fill_gradientn(
      colours = warming_colours,
      name = "JJA warming\n(destination - origin, °C)"
    ) +
    ggplot2::scale_colour_gradientn(colours = warming_colours, guide = "none") +
    ggplot2::scale_x_continuous(
      breaks = sort(unique(c(
        pretty(x_range), temperature_edges[2:3]
      ))),
      labels = scales::label_number(accuracy = 0.01, big.mark = ""),
      expand = ggplot2::expansion(mult = 0.025)
    ) +
    ggplot2::scale_y_continuous(
      breaks = sort(unique(c(
        pretty(y_range), precipitation_edges[2:3]
      ))),
      labels = scales::label_number(accuracy = 0.01, big.mark = ""),
      expand = ggplot2::expansion(mult = 0.025)
    ) +
    transplant_figure_theme() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right"
    ) +
    ggplot2::labs(
      x = "Destination JJA mean temperature (°C)",
      y = "Destination JJA precipitation (mm)"
    )
  save_plot_atomic(plot, path, 18, 12)
}




signed_species_consistency <- function(data) {
  counts <- data |>
    dplyr::transmute(
      species = dplyr::coalesce(.data$accepted_name, .data$species),
      .data$experiment_id,
      category = as.character(.data$category)
    ) |>
    dplyr::distinct(.data$species, .data$experiment_id, .data$category) |>
    dplyr::group_by(.data$species) |>
    dplyr::filter(dplyr::n_distinct(.data$experiment_id) > 3L) |>
    dplyr::count(.data$species, .data$category, name = "n") |>
    tidyr::pivot_wider(names_from = "category", values_from = "n", values_fill = 0) |>
    dplyr::ungroup()
  for (column in c(
    "loser", "non_significant_loser", "no_change", "non_significant_winner",
    "winner", "colonizer", "non_significant_colonizer"
  )) {
    if (!column %in% names(counts)) counts[[column]] <- 0L
  }
  counts |>
    dplyr::mutate(
      significant_total = .data$loser + .data$winner + .data$colonizer,
      neutral_total = .data$non_significant_loser + .data$no_change +
        .data$non_significant_winner + .data$non_significant_colonizer,
      neutral_low = -.data$neutral_total / 2,
      neutral_high = .data$neutral_total / 2,
      significant_categories = (.data$loser > 0L) + (.data$winner > 0L) +
        (.data$colonizer > 0L),
      consistency_group = dplyr::case_when(
        .data$loser > 0L & .data$winner == 0L & .data$colonizer == 0L ~ "Loser only",
        .data$significant_categories > 1L ~ "Experiment-dependent responses",
        .data$winner > 0L & .data$loser == 0L & .data$colonizer == 0L ~ "Winner only",
        .data$colonizer > 0L & .data$loser == 0L & .data$winner == 0L ~ "Colonizer only",
        TRUE ~ NA_character_
      ),
      loser_low = .data$neutral_low - .data$loser,
      winner_high = .data$neutral_high + .data$winner,
      colonizer_high = .data$winner_high + .data$colonizer,
      net = .data$winner + .data$colonizer - .data$loser,
      within_group_order = dplyr::case_when(
        .data$consistency_group == "Loser only" ~ -.data$loser,
        .data$consistency_group == "Experiment-dependent responses" ~ .data$net,
        .data$consistency_group == "Winner only" ~ .data$winner,
        .data$consistency_group == "Colonizer only" ~ .data$colonizer,
        TRUE ~ 0
      ),
      consistency_group = factor(
        .data$consistency_group,
        levels = c(
          "Loser only", "Experiment-dependent responses", "Winner only", "Colonizer only"
        )
      )
    ) |>
    dplyr::filter(.data$significant_total > 0L, !is.na(.data$consistency_group)) |>
    dplyr::arrange(
      .data$consistency_group, .data$within_group_order,
      .data$significant_total, .data$species
    ) |>
    dplyr::mutate(species_position = dplyr::row_number())
}

plot_signed_consistency <- function(data, path) {
  signed <- signed_species_consistency(data)
  rectangles <- dplyr::bind_rows(
    signed |>
      dplyr::transmute(
        .data$species_position, .data$species, .data$consistency_group,
        response = "Loser", ymin = .data$loser_low, ymax = .data$neutral_low
      ),
    signed |>
      dplyr::transmute(
        .data$species_position, .data$species, .data$consistency_group,
        response = "No significant change", ymin = .data$neutral_low,
        ymax = .data$neutral_high
      ),
    signed |>
      dplyr::transmute(
        .data$species_position, .data$species, .data$consistency_group,
        response = "Winner", ymin = .data$neutral_high, ymax = .data$winner_high
      ),
    signed |>
      dplyr::transmute(
        .data$species_position, .data$species, .data$consistency_group,
        response = "Colonizer", ymin = .data$winner_high, ymax = .data$colonizer_high
      )
  ) |>
    dplyr::filter(.data$ymax > .data$ymin) |>
    dplyr::mutate(
      response = factor(
        .data$response,
        levels = c("Loser", "No significant change", "Winner", "Colonizer")
      )
    )
  response_palette <- c(
    Loser = category_palette[["loser"]],
    `No significant change` = "#CFCFCF",
    Winner = category_palette[["winner"]],
    Colonizer = category_palette[["colonizer"]]
  )
  significant_labels <- rectangles |>
    dplyr::filter(.data$response %in% c("Loser", "Winner", "Colonizer")) |>
    dplyr::mutate(
      experiment_count = .data$ymax - .data$ymin,
      label_y = (.data$ymin + .data$ymax) / 2,
      label = sprintf("%d", round(.data$experiment_count))
    )
  if (any(significant_labels$experiment_count <= 0) ||
      any(abs(
        significant_labels$experiment_count -
          round(significant_labels$experiment_count)
      ) > 1e-7) ||
      abs(sum(significant_labels$experiment_count) -
            sum(signed$significant_total)) > 1e-7 ||
      any(significant_labels$response == "No significant change")) {
    stop("Significant species-segment labels are invalid.", call. = FALSE)
  }
  plot <- ggplot2::ggplot(rectangles) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = .data$species_position - 0.42, xmax = .data$species_position + 0.42,
        ymin = .data$ymin, ymax = .data$ymax, fill = .data$response
      ),
      colour = "white", linewidth = 0.2
    ) +
    ggplot2::geom_text(
      data = significant_labels,
      mapping = ggplot2::aes(
        x = .data$species_position, y = .data$label_y,
        label = .data$label
      ),
      inherit.aes = FALSE, colour = "white", fontface = "bold", size = 2.55
    ) +
    ggplot2::geom_hline(yintercept = 0, colour = "white", linewidth = 0.5) +
    ggplot2::facet_grid(
      . ~ consistency_group, scales = "free_x", space = "free_x", drop = TRUE
    ) +
    ggplot2::scale_x_continuous(
      breaks = signed$species_position, labels = signed$species,
      expand = ggplot2::expansion(add = 0.6)
    ) +
    ggplot2::scale_fill_manual(values = response_palette, drop = FALSE, name = NULL) +
    ggplot2::scale_y_continuous(
      breaks = function(limits) {
        breaks <- unique(round(pretty(limits, n = 7)))
        breaks[breaks >= ceiling(limits[[1]]) & breaks <= floor(limits[[2]])]
      },
      labels = scales::label_number(accuracy = 1)
    ) +
    transplant_figure_theme() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.spacing.x = grid::unit(0.25, "cm"),
      strip.text = ggplot2::element_text(face = "bold", size = 13),
      axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, size = 8),
      legend.position = "bottom"
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1)) +
    ggplot2::labs(
      x = NULL,
      y = "Number of experiments\n(losers below zero)"
    )
  save_plot_atomic(plot, path, 14, 8.5)
}

signed_family_consistency <- function(data) {
  if (!"family" %in% names(data)) {
    stop("Family consistency requires a reviewed family column.", call. = FALSE)
  }
  observations <- data |>
    dplyr::transmute(
      family = trimws(as.character(.data$family)),
      species = dplyr::coalesce(.data$accepted_name, .data$species),
      .data$experiment_id,
      category = as.character(.data$category)
    ) |>
    dplyr::mutate(
      family = dplyr::if_else(
        is.na(.data$family) | .data$family == "", "Unresolved", .data$family
      )
    ) |>
    dplyr::distinct(.data$family, .data$species, .data$experiment_id, .data$category)
  unexpected_categories <- setdiff(unique(observations$category), names(category_palette))
  if (length(unexpected_categories)) {
    stop(
      "Family consistency encountered unexpected categories: ",
      paste(unexpected_categories, collapse = ", "), call. = FALSE
    )
  }
  sampling <- observations |>
    dplyr::summarise(
      total_combinations = dplyr::n(),
      n_species = dplyr::n_distinct(.data$species),
      n_experiments = dplyr::n_distinct(.data$experiment_id),
      .by = "family"
    )
  counts <- observations |>
    dplyr::count(.data$family, .data$category, name = "n") |>
    tidyr::pivot_wider(names_from = "category", values_from = "n", values_fill = 0) |>
    dplyr::left_join(sampling, by = "family")
  for (column in names(category_palette)) {
    if (!column %in% names(counts)) counts[[column]] <- 0L
  }
  counts |>
    dplyr::mutate(
      significant_total = .data$loser + .data$winner + .data$colonizer,
      neutral_total = .data$non_significant_loser + .data$no_change +
        .data$non_significant_winner + .data$non_significant_colonizer,
      loser_percent = 100 * .data$loser / .data$total_combinations,
      neutral_percent = 100 * .data$neutral_total / .data$total_combinations,
      winner_percent = 100 * .data$winner / .data$total_combinations,
      colonizer_percent = 100 * .data$colonizer / .data$total_combinations,
      neutral_low = -.data$neutral_percent / 2,
      neutral_high = .data$neutral_percent / 2,
      neutral_count_low = -.data$neutral_total / 2,
      neutral_count_high = .data$neutral_total / 2,
      significant_categories = (.data$loser > 0L) + (.data$winner > 0L) +
        (.data$colonizer > 0L),
      consistency_group = dplyr::case_when(
        .data$loser > 0L & .data$winner == 0L & .data$colonizer == 0L ~ "Loser only",
        .data$significant_categories > 1L ~ "Mixed significant responses",
        .data$winner > 0L & .data$loser == 0L & .data$colonizer == 0L ~ "Winner only",
        .data$colonizer > 0L & .data$loser == 0L & .data$winner == 0L ~ "Colonizer only",
        TRUE ~ NA_character_
      ),
      loser_low = .data$neutral_low - .data$loser_percent,
      winner_high = .data$neutral_high + .data$winner_percent,
      colonizer_high = .data$winner_high + .data$colonizer_percent,
      loser_count_low = .data$neutral_count_low - .data$loser,
      winner_count_high = .data$neutral_count_high + .data$winner,
      colonizer_count_high = .data$winner_count_high + .data$colonizer,
      net_percent = .data$winner_percent + .data$colonizer_percent -
        .data$loser_percent,
      within_group_order = dplyr::case_when(
        .data$consistency_group == "Loser only" ~ -.data$loser_percent,
        .data$consistency_group == "Mixed significant responses" ~ .data$net_percent,
        .data$consistency_group == "Winner only" ~ .data$winner_percent,
        .data$consistency_group == "Colonizer only" ~ .data$colonizer_percent,
        TRUE ~ 0
      ),
      consistency_group = factor(
        .data$consistency_group,
        levels = c(
          "Loser only", "Mixed significant responses", "Winner only", "Colonizer only"
        )
      )
    ) |>
    dplyr::filter(.data$significant_total > 0L, !is.na(.data$consistency_group)) |>
    dplyr::arrange(
      .data$consistency_group, .data$within_group_order,
      .data$significant_total, .data$family
    ) |>
    dplyr::mutate(
      family_position = dplyr::row_number(),
      family_label = .data$family
    )
}

family_consistency_plot_data <- function(signed, metric = c("percentage", "count")) {
  metric <- match.arg(metric)
  if (metric == "percentage") {
    values <- signed |>
      dplyr::transmute(
        .data$family, .data$family_label, .data$consistency_group,
        .data$total_combinations,
        loser = -.data$loser_percent,
        neutral_negative = -.data$neutral_percent / 2,
        neutral_positive = .data$neutral_percent / 2,
        winner = .data$winner_percent,
        colonizer = .data$colonizer_percent
      )
  } else {
    values <- signed |>
      dplyr::transmute(
        .data$family, .data$family_label, .data$consistency_group,
        .data$total_combinations,
        loser = -.data$loser,
        neutral_negative = -.data$neutral_total / 2,
        neutral_positive = .data$neutral_total / 2,
        winner = .data$winner,
        colonizer = .data$colonizer
      )
  }
  stack_levels <- c(
    "colonizer", "winner", "neutral_positive", "loser", "neutral_negative"
  )
  values |>
    tidyr::pivot_longer(
      dplyr::all_of(stack_levels), names_to = "stack_category", values_to = "value"
    ) |>
    dplyr::mutate(
      family_axis = factor(.data$family_label, levels = signed$family_label),
      stack_category = factor(.data$stack_category, levels = stack_levels)
    ) |>
    dplyr::filter(.data$value != 0)
}




species_climate_fingerprint_data <- function(data, climate_contexts, minimum_experiments = 10L) {
  if (length(minimum_experiments) != 1L || is.na(minimum_experiments) ||
      minimum_experiments < 1L) {
    stop("minimum_experiments must be one positive integer.", call. = FALSE)
  }
  required_context <- c(
    "experiment_id", "destination_T_JJA", "destination_P_JJA",
    "temperature_bin_lower", "temperature_bin_upper",
    "precipitation_bin_lower", "precipitation_bin_upper"
  )
  missing_context <- setdiff(required_context, names(climate_contexts))
  if (length(missing_context)) {
    stop(
      "Species climate fingerprints are missing climate columns: ",
      paste(missing_context, collapse = ", "), call. = FALSE
    )
  }
  experiment_climate <- climate_contexts |>
    dplyr::select(dplyr::all_of(required_context)) |>
    dplyr::distinct()
  if (anyDuplicated(experiment_climate$experiment_id)) {
    stop("Each experiment must have one unique destination climate for fingerprints.", call. = FALSE)
  }
  if (any(!is.finite(experiment_climate$destination_T_JJA)) ||
      any(!is.finite(experiment_climate$destination_P_JJA))) {
    stop("Every fingerprint experiment requires finite destination JJA climate.", call. = FALSE)
  }
  observations <- data |>
    dplyr::transmute(
      species = dplyr::coalesce(.data$accepted_name, .data$species),
      .data$experiment_id,
      category = as.character(.data$category)
    ) |>
    dplyr::distinct(.data$species, .data$experiment_id, .data$category)
  duplicated_cells <- observations |>
    dplyr::count(.data$species, .data$experiment_id) |>
    dplyr::filter(.data$n != 1L)
  if (nrow(duplicated_cells)) {
    stop("Accepted species do not uniquely identify fingerprint experiment cells.", call. = FALSE)
  }
  signed <- signed_species_consistency(data)
  eligible <- observations |>
    dplyr::count(.data$species, name = "n_experiments") |>
    dplyr::filter(.data$n_experiments >= as.integer(minimum_experiments)) |>
    dplyr::semi_join(signed, by = "species")
  species_summary <- signed |>
    dplyr::inner_join(eligible, by = "species") |>
    dplyr::arrange(
      .data$consistency_group, .data$within_group_order,
      .data$significant_total, .data$species
    ) |>
    dplyr::mutate(
      species_order = dplyr::row_number(),
      facet_label = sprintf(
        "%s\nN=%d; L=%d, W=%d, C=%d",
        .data$species, .data$n_experiments,
        .data$loser, .data$winner, .data$colonizer
      )
    )
  if (!nrow(species_summary)) {
    stop("No significant-response species meet the fingerprint replication threshold.", call. = FALSE)
  }
  observations <- observations |>
    dplyr::semi_join(species_summary, by = "species") |>
    dplyr::left_join(
      experiment_climate |>
        dplyr::select("experiment_id", "destination_T_JJA", "destination_P_JJA"),
      by = "experiment_id"
    ) |>
    dplyr::mutate(
      response = dplyr::case_when(
        .data$category == "loser" ~ "Loser",
        .data$category == "winner" ~ "Winner",
        .data$category == "colonizer" ~ "Colonizer",
        TRUE ~ "No significant change"
      ),
      response = factor(
        .data$response,
        levels = c("Loser", "No significant change", "Winner", "Colonizer")
      ),
      species = factor(.data$species, levels = species_summary$species)
    )
  if (any(!is.finite(observations$destination_T_JJA)) ||
      any(!is.finite(observations$destination_P_JJA))) {
    stop("Every selected species-experiment cell requires destination JJA climate.", call. = FALSE)
  }
  temperature_boundaries <- sort(unique(c(
    experiment_climate$temperature_bin_lower,
    experiment_climate$temperature_bin_upper
  )))
  precipitation_boundaries <- sort(unique(c(
    experiment_climate$precipitation_bin_lower,
    experiment_climate$precipitation_bin_upper
  )))
  list(
    observations = observations,
    species_summary = species_summary,
    experiment_climate = experiment_climate,
    temperature_cutoffs = temperature_boundaries[-c(1L, length(temperature_boundaries))],
    precipitation_cutoffs = precipitation_boundaries[-c(1L, length(precipitation_boundaries))]
  )
}





plot_sensitivity <- function(classified, path) {
  data <- classified |>
    dplyr::ungroup() |>
    dplyr::filter(.data$horizon == "max") |>
    dplyr::select(experiment_id, species, transform, estimate) |>
    tidyr::pivot_wider(names_from = transform, values_from = estimate)
  plot <- ggplot2::ggplot(data, ggplot2::aes(.data$raw, .data$sqrt)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey75") +
    ggplot2::geom_vline(xintercept = 0, colour = "grey75") +
    ggplot2::geom_abline(slope = 1, colour = "grey60") +
    ggplot2::geom_point(alpha = 0.35, size = 0.7) +
    ggplot2::geom_smooth(method = "lm", colour = "#D55E00", se = FALSE) +
    transplant_figure_theme() +
    ggplot2::labs(x = "Raw relative-cover strength", y = "Square-root relative-cover strength")
  save_plot_atomic(plot, path, 8, 6)
}

figure_manifest <- function(root) {
  curated <- tibble::tribble(
    ~filename, ~question, ~source_data, ~target, ~width, ~height, ~status,
    "figures/manuscript/transient_winning_losing_strength.pdf", "Temporal consistency of continuous strength with Pearson correlations and p-values", "winner_loser_results.csv", "figure_transient_strength", 9, 4.5, "CURRENT_CANONICAL_ANALYSIS_FIGURE",
    "figures/manuscript/temporal_category_flows.pdf", "Observed temporal category transitions with source-conditional percentages", "temporal_linked_flows.csv", "figure_temporal_flows", 18, 11, "CURRENT_CANONICAL_ANALYSIS_FIGURE",
    "figures/manuscript/temporal_transition_matrices.pdf", "Category persistence, transitions, and composition through time", "temporal_transition_matrix.csv; temporal_category_composition.csv", "figure_temporal_transition_matrices", 17, 13, "CURRENT_CANONICAL_ANALYSIS_FIGURE",
    "figures/manuscript/destination_climate_warming_gradient.pdf", "Experimental warming across equal-width climate cells with reviewed experiment-level display assignments", "analysis_dataset.csv; climate_context_experiments.csv", "figure_destination_climate_warming_gradient", 18, 12, "CURRENT_CANONICAL_ANALYSIS_FIGURE",
    "figures/manuscript/repeated_species_consistency.pdf", "Significant-response consistency among species in more than three experiments", "winner_loser_annotated.csv", "figure_species_consistency", 14, 8.5, "CURRENT_CANONICAL_ANALYSIS_FIGURE",
    "figures/manuscript/compare_rel.abund_vs.sqrt.rel.abund.pdf", "Raw-versus-square-root sensitivity", "winner_loser_results.csv", "figure_sensitivity", 8, 6, "CURRENT_CANONICAL_ANALYSIS_FIGURE",
    "plot/climate_response_data_model_probability_prototype.pdf", "Observed densities, pool-specific odds ratios, and probability surfaces", "analysis_dataset.csv; climate_context_experiments.csv; reconstructed model predictions", "R/04_climate_response_model.R", 26, 14, "provisional modeling checkpoint; source outputs reconstructable",
    "plot/response_to_climate_alluvial.pdf", "Significant responses linked to destination-climate groups with validated counts and percentages", "analysis_dataset.csv", "R/10_response_grouping_flows_gradients_review.R", 17, 11, "current; standalone alluvial",
    "plot/response_by_climate_percentages.pdf", "Response-conditional loser, winner, and colonizer percentages across destination climates", "analysis_dataset.csv", "R/10_response_grouping_flows_gradients_review.R", 12, 7, "current; standalone descriptive figure",
    "plot/response_by_climate_percentages_response_axis.pdf", "Response-conditional destination-climate composition shown within loser, winner, and colonizer categories", "analysis_dataset.csv", "R/10_response_grouping_flows_gradients_review.R", 12, 7, "current; standalone descriptive figure",
    "plot/response_to_family_alluvial.pdf", "Significant response pools linked to responses, functional groups, and families", "analysis_dataset.csv", "R/10_response_grouping_flows_gradients_review.R", 28, 20, "current; standalone alluvial",
    "plot/trait_categories_with_overall.pdf", "Significant loser, winner, and colonizer trait distributions across destination-climate backgrounds", "analysis_dataset.csv; climate_context_experiments.csv", "R/07_climate_faceted_descriptive_review.R", 23, 8.5, "current; standalone descriptive figure",
    "plot/eg.sp.1.pdf", "Legacy worked species example", "Legacy winner-loser analysis objects", NA_character_, 7, 9, "legacy artifact; retained unchanged",
    "plot/eg.sp.2.pdf", "Legacy worked species example", "Legacy winner-loser analysis objects", NA_character_, 5, 4, "legacy artifact; retained unchanged",
    "plot/winners_losers.pdf", "Legacy winner-loser summary", "Legacy winner-loser analysis objects", NA_character_, 9, 7, "legacy artifact; retained unchanged",
    "plot/figures_presentation/focal4_trait_phenology_warming.png", "Legacy four-trait phenology panel", "source code unavailable", NA_character_, NA_real_, NA_real_, "legacy artifact; retained unchanged"
  )
  presentation_files <- file.path(
    "plot", "figures_presentation",
    list.files(
      file.path(root, "plot", "figures_presentation"),
      pattern = "[.](pdf|png)$", recursive = FALSE,
      full.names = FALSE, ignore.case = TRUE
    )
  )
  uncatalogued <- tibble::tibble(
    filename = setdiff(presentation_files, curated$filename)
  ) |>
    dplyr::mutate(
      question = "Legacy exploratory or presentation artifact",
      source_data = "See legacy R scripts; exact source may be unavailable",
      target = NA_character_,
      width = NA_real_,
      height = NA_real_,
      status = "legacy artifact; retained unchanged"
    )
  dplyr::bind_rows(curated, uncatalogued) |>
    dplyr::arrange(.data$filename) |>
    dplyr::mutate(
      file_exists = file.exists(file.path(root, .data$filename)),
      generation_receipt = file.exists(paste0(file.path(root, .data$filename), ".generated")),
      generated = dplyr::if_else(
        is.na(.data$target), .data$file_exists,
        .data$file_exists & .data$generation_receipt
      ),
      status = dplyr::case_when(
        .data$generated ~ .data$status,
        !is.na(.data$target) ~ "implemented target; not generated",
        TRUE ~ "documented artifact; file missing"
      ),
      filename = file.path(root, .data$filename)
    )
}
