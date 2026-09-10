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

prepare_temporal_strength_pairs <- function(classified) {
  horizon_rows <- classified |>
    dplyr::ungroup() |>
    dplyr::filter(
      .data$transform == "raw",
      as.character(.data$contrast) == "to_origin",
      .data$pool == "resident",
      .data$horizon %in% c("2", "4", "6")
    ) |>
    dplyr::select(
      "experiment_id", "species", "horizon", "estimate", "category"
    )
  duplicated_rows <- horizon_rows |>
    dplyr::count(.data$experiment_id, .data$species, .data$horizon) |>
    dplyr::filter(.data$n != 1L)
  if (nrow(duplicated_rows)) {
    stop("Temporal strength input is not unique by experiment, species, and horizon.",
         call. = FALSE)
  }

  wide <- horizon_rows |>
    tidyr::pivot_wider(
      names_from = "horizon",
      values_from = c("estimate", "category"),
      names_sep = "_"
    )
  pairs <- dplyr::bind_rows(
    wide |>
      dplyr::transmute(
        comparison = "First 2 vs First 4 years",
        earlier_strength = .data$estimate_2,
        later_strength = .data$estimate_4,
        source_category = as.character(.data$category_2)
      ),
    wide |>
      dplyr::transmute(
        comparison = "First 4 vs First 6 years",
        earlier_strength = .data$estimate_4,
        later_strength = .data$estimate_6,
        source_category = as.character(.data$category_4)
      )
  ) |>
    dplyr::filter(
      is.finite(.data$earlier_strength),
      is.finite(.data$later_strength)
    ) |>
    dplyr::mutate(
      direction = dplyr::case_when(
        .data$source_category %in% c("loser", "non_significant_loser") ~
          "Loser",
        .data$source_category %in% c("winner", "non_significant_winner") ~
          "Winner",
        TRUE ~ NA_character_
      ),
      comparison = factor(
        .data$comparison,
        levels = c(
          "First 2 vs First 4 years",
          "First 4 vs First 6 years"
        )
      ),
      direction = factor(.data$direction, levels = c("Loser", "Winner"))
    )
  if (!nrow(pairs) || any(is.na(pairs$direction))) {
    stop("Temporal strength pairs contain no usable resident winner/loser directions.",
         call. = FALSE)
  }
  pairs
}

summarise_temporal_strength_correlations <- function(pairs) {
  grouped <- dplyr::bind_rows(
    pairs |> dplyr::mutate(correlation_group = "All"),
    pairs |>
      dplyr::mutate(correlation_group = as.character(.data$direction))
  )
  grouped |>
    dplyr::summarise(
      n = dplyr::n(),
      correlation = stats::cor(.data$earlier_strength, .data$later_strength),
      p_value = stats::cor.test(
        .data$earlier_strength, .data$later_strength,
        method = "pearson"
      )$p.value,
      .by = c("comparison", "correlation_group")
    ) |>
    dplyr::mutate(
      correlation_group = factor(
        .data$correlation_group,
        levels = c("All", "Loser", "Winner")
      ),
      label = paste0(
        as.character(.data$correlation_group),
        ": r = ",
        sprintf("%.2f", .data$correlation),
        dplyr::if_else(.data$p_value < 0.05, "*", "")
      ),
      label_order = match(
        as.character(.data$correlation_group), c("All", "Loser", "Winner")
      ),
      label_vjust = 1.15 + 1.15 * (.data$label_order - 1)
    )
}

build_temporal_strength_correlation_inset <- function(classified) {
  pairs <- prepare_temporal_strength_pairs(classified)
  correlations <- summarise_temporal_strength_correlations(pairs)
  correlation_colours <- c(
    All = "grey30",
    Loser = unname(category_palette[["loser"]]),
    Winner = unname(category_palette[["winner"]])
  )
  build_panel <- function(comparison, x_label, y_label) {
    panel_pairs <- pairs |>
      dplyr::filter(as.character(.data$comparison) == .env$comparison)
    panel_correlations <- correlations |>
      dplyr::filter(as.character(.data$comparison) == .env$comparison)

    ggplot2::ggplot(
      panel_pairs,
      ggplot2::aes(.data$earlier_strength, .data$later_strength)
    ) +
      ggplot2::geom_hline(
        yintercept = 0, colour = "grey78", linewidth = 0.35
      ) +
      ggplot2::geom_vline(
        xintercept = 0, colour = "grey78", linewidth = 0.35
      ) +
      ggplot2::geom_abline(
        slope = 1, intercept = 0, colour = "grey70", linewidth = 0.45
      ) +
      ggplot2::geom_point(
        ggplot2::aes(colour = .data$direction),
        alpha = 0.28, size = 0.8
      ) +
      ggplot2::geom_smooth(
        ggplot2::aes(group = 1),
        method = "lm", formula = y ~ x, se = FALSE,
        colour = correlation_colours[["All"]], linewidth = 1
      ) +
      ggplot2::geom_smooth(
        ggplot2::aes(colour = .data$direction),
        method = "lm", formula = y ~ x, se = FALSE,
        linewidth = 1
      ) +
      ggplot2::geom_text(
        data = panel_correlations,
        ggplot2::aes(
          x = -Inf, y = Inf, label = .data$label,
          colour = .data$correlation_group,
          vjust = .data$label_vjust
        ),
        inherit.aes = FALSE,
        hjust = -0.06,
        size = 5.3,
        fontface = "bold"
      ) +
      ggplot2::scale_colour_manual(
        values = correlation_colours, guide = "none"
      ) +
      ggplot2::labs(title = comparison, x = x_label, y = y_label) +
      ggplot2::theme_bw(base_size = 16) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_line(
          colour = "grey90", linewidth = 0.3
        ),
        plot.title = ggplot2::element_text(
          size = 16, face = "bold", hjust = 0.5,
          margin = ggplot2::margin(b = 5)
        ),
        axis.title = ggplot2::element_text(size = 15),
        axis.text = ggplot2::element_text(size = 13),
        plot.margin = ggplot2::margin(6, 6, 6, 6)
      )
  }

  first_panel <- build_panel(
    "First 2 vs First 4 years",
    expression("First 2 years " * Delta * " AUC"),
    expression("First 4 years " * Delta * " AUC")
  )
  second_panel <- build_panel(
    "First 4 vs First 6 years",
    expression("First 4 years " * Delta * " AUC"),
    expression("First 6 years " * Delta * " AUC")
  )
  cowplot::plot_grid(
    first_panel, second_panel,
    ncol = 2, align = "h", axis = "tb"
  )
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

plot_temporal_flows <- function(
    linked_transitions,
    path,
    include_colonizers = TRUE,
    base_size = 16,
    flow_label_size = 3.0,
    stratum_label_size = 3.6,
    axis_text_size = 16,
    legend_text_size = 13,
    legend_position = "bottom",
    direct_right_labels = FALSE,
    direct_label_size = 5.5,
    strength_inset = NULL,
    width = 18,
    height = 11,
    bg = NULL) {
  horizon_columns <- c("state_2", "state_4", "state_6")
  colonizer_levels <- c("colonizer", "non_significant_colonizer")
  if (!include_colonizers) {
    linked_transitions <- linked_transitions |>
      dplyr::filter(
        !dplyr::if_any(
          dplyr::all_of(.env$horizon_columns),
          ~ .x %in% .env$colonizer_levels
        )
      )
  }
  state_levels <- c(
    "loser", "non_significant_loser", "non_significant_winner",
    "winner", "colonizer", "non_significant_colonizer",
    "follow_up_unavailable"
  )
  if (!include_colonizers) {
    state_levels <- setdiff(state_levels, colonizer_levels)
  }
  visible_levels <- setdiff(state_levels, "follow_up_unavailable")
  flow_colours <- c(category_palette, follow_up_unavailable = "white")
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
    ) +
    ggplot2::scale_colour_manual(values = flow_colours, guide = "none")

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
  right_labels <- stratum_labels |>
    dplyr::filter(.data$x == max(.data$x)) |>
    dplyr::mutate(
      label = dplyr::recode(
        .data$state,
        loser = "Loser",
        non_significant_loser = "Loser (ns)",
        non_significant_winner = "Winner (ns)",
        winner = "Winner"
      )
    )

  plot <- base +
    ggrepel::geom_label_repel(
      data = flow_labels,
      ggplot2::aes(
        .data$label_x,
        .data$label_y,
        label = .data$label,
        colour = .data$source_state
      ),
      inherit.aes = FALSE, size = flow_label_size, fontface = "bold", linewidth = 0,
      label.padding = grid::unit(0.06, "lines"),
      fill = scales::alpha("white", 0.82),
      direction = "y", box.padding = 0.12, point.padding = 0,
      min.segment.length = Inf, max.overlaps = Inf, seed = 47163
    ) +
    ggplot2::geom_label(
      data = stratum_labels,
      ggplot2::aes(.data$x, .data$y, label = .data$label),
      inherit.aes = FALSE, size = stratum_label_size, fontface = "bold", lineheight = 0.88,
      label.padding = grid::unit(0.12, "lines"), linewidth = 0,
      fill = scales::alpha("white", 0.90), colour = "grey10"
    ) +
    ggplot2::scale_x_continuous(
      breaks = 1:3, labels = c("First 2 years", "First 4 years", "First 6 years"),
      limits = c(0.78, 3.22), expand = c(0, 0), position = "top"
    ) +
    ggplot2::scale_y_continuous(breaks = NULL, expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        colour = "#1A1A1A", size = axis_text_size, face = "bold",
        margin = ggplot2::margin(b = 12)
      ),
      legend.position = legend_position,
      legend.text = ggplot2::element_text(size = legend_text_size),
      legend.key.width = grid::unit(1.5, "lines"),
      plot.margin = ggplot2::margin(12, 32, 8, 32)
    )
  if (direct_right_labels) {
    plot <- plot +
      ggplot2::geom_text(
        data = right_labels,
        ggplot2::aes(
          x = 3.165,
          y = .data$y,
          label = .data$label,
          colour = .data$state
        ),
        inherit.aes = FALSE,
        hjust = 0,
        size = direct_label_size,
        fontface = "bold"
      ) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::theme(plot.margin = ggplot2::margin(12, 145, 8, 32))
  }
  if (!is.null(strength_inset)) {
    plot <- cowplot::ggdraw() +
      cowplot::draw_plot(plot, x = 0, y = 0, width = 1, height = 1) +
      cowplot::draw_plot(
        strength_inset,
        x = 0.695, y = 0.025, width = 0.290, height = 0.605
      )
  }
  save_plot_atomic(plot, path, width, height, bg = bg)
}

plot_temporal_transition_matrices <- function(
    transitions,
    composition,
    path,
    include_colonizers = TRUE,
    resident_classified = NULL) {
  colonizer_levels <- c("colonizer", "non_significant_colonizer")
  category_levels <- c(
    "loser", "non_significant_loser", "non_significant_winner", "winner",
    "colonizer", "non_significant_colonizer"
  )
  if (!include_colonizers) {
    category_levels <- setdiff(category_levels, colonizer_levels)
  }
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
      .data$to_category != "no_change",
      .data$from_category %in% .env$category_levels,
      .data$to_category %in% .env$category_levels
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
    ggplot2::coord_fixed(
      xlim = c(0.5, length(category_levels) + 0.5),
      ylim = c(0.5, length(category_levels) + 0.5),
      clip = "off"
    ) +
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
    dplyr::filter(
      .data$category != "no_change",
      .data$category %in% .env$category_levels
    )
  if (!include_colonizers) {
    composition_data <- composition_data |>
      dplyr::mutate(
        proportion = .data$n / sum(.data$n),
        total_species_experiment = sum(.data$n),
        .by = "year"
      )
  }
  composition_data <- composition_data |>
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
  if (!include_colonizers && !is.null(resident_classified)) {
    year_metadata <- resident_classified |>
      dplyr::ungroup() |>
      dplyr::filter(
        .data$transform == "raw",
        as.character(.data$contrast) == "to_origin",
        .data$pool == "resident",
        .data$horizon %in% as.character(.env$composition_years)
      ) |>
      dplyr::summarise(
        total_species_experiment = dplyr::n(),
        n_species = dplyr::n_distinct(.data$species),
        n_experiments = dplyr::n_distinct(.data$experiment_id),
        .by = "horizon"
      ) |>
      dplyr::transmute(
        year = as.integer(.data$horizon),
        .data$total_species_experiment,
        .data$n_species,
        .data$n_experiments
      ) |>
      dplyr::arrange(.data$year)
  } else {
    year_metadata <- composition |>
      dplyr::distinct(
        .data$year, .data$total_species_experiment,
        .data$n_species, .data$n_experiments
      ) |>
      dplyr::arrange(.data$year)
  }
  if (anyDuplicated(year_metadata$year)) {
    stop("Temporal sample-size metadata is inconsistent within at least one year.", call. = FALSE)
  }
  if (!identical(year_metadata$year, composition_years)) {
    stop("Temporal sample-size metadata does not cover every plotted year.", call. = FALSE)
  }
  year_labels <- stats::setNames(
    sprintf(
      "%d\nC=%s\nS=%s  E=%s",
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
