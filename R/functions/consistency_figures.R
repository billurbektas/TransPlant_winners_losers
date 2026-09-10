# Canonical consistency, sensitivity, and figure-ledger functions.

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
    "figures/retained/response_to_family_alluvial.pdf", "Significant response pools linked to responses, functional groups, and families", "analysis_dataset.csv", NA_character_, 28, 20, "KEEP_FOR_MANUSCRIPT",
    "figures/retained/trait_categories_with_overall.pdf", "Significant loser, winner, and colonizer trait distributions across destination-climate backgrounds", "analysis_dataset.csv; climate_context_experiments.csv", NA_character_, 23, 8.5, "KEEP_FOR_MANUSCRIPT",
    "figures/retained/eg.sp.1.pdf", "Legacy worked species example", "Legacy winner-loser analysis objects", NA_character_, 7, 9, "KEEP_FOR_MANUSCRIPT",
    "figures/retained/eg.sp.2.pdf", "Legacy worked species example", "Legacy winner-loser analysis objects", NA_character_, 5, 4, "KEEP_FOR_MANUSCRIPT",
    "figures/retained/winners_losers.pdf", "Legacy winner-loser summary", "Legacy winner-loser analysis objects", NA_character_, 9, 7, "KEEP_FOR_MANUSCRIPT"
  )
  curated |>
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

# Family-level consistency helpers retained for audited descriptive outputs.

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
