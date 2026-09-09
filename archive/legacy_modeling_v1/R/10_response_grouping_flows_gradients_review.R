# ==============================================================================
# Script: 10_response_grouping_flows_gradients_review.R
#
# Purpose:
#   Generate the retained significant-response alluvials linking responses to
#   destination climate and to functional-group/family classifications. The
#   optional continuous summaries are written under output/review only when
#   this script is run and are absent from the clean checkpoint.
#
# Notes:
#   - Every alluvial ribbon represents a significant species x experiment
#     response.
#   - In the climate alluvial, response bars show counts only, climate-bar
#     percentages use all significant rows, and every path percentage uses its
#     response total so the five paths sum to 100% within each response.
#   - The companion climate bar charts use the same response-conditional
#     percentages and print the underlying count above each bar. One places
#     climate on the x-axis; the alternative places response on the x-axis.
#   - Families use the outcome-independent >= 50-record rule plus Other.
#   - Family-path percentages use the complete response totals. No additional
#     response counts or percentages are printed to the right of the families.
#   - Continuous summaries use the exact experiment-level destination JJA
#     temperature, destination JJA precipitation, and JJA warming values.
#   - Loser/winner percentages use resident denominators; colonizer
#     percentages use nonresident denominators.
# ============================================================================== 

source("config/analysis.R")
source("R/pipeline/utils.R")
source("R/pipeline/figures.R")

config <- transplant_config()
analysis_path <- file.path(config$output_dir, "tables", "analysis_dataset.csv")
figure_dir <- config$plot_dir
review_output_dir <- file.path(config$output_dir, "review")

if (!file.exists(analysis_path)) {
  stop("The saved analysis dataset is missing.", call. = FALSE)
}
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(review_output_dir, recursive = TRUE, showWarnings = FALSE)

response_levels <- c(
  "Loser", "Winner", "Colonizer",
  "Loser (non-significant)",
  "Winner (non-significant)",
  "Colonizer (non-significant)"
)
response_palette <- c(
  `Loser` = category_palette[["loser"]],
  `Winner` = category_palette[["winner"]],
  `Colonizer` = category_palette[["colonizer"]],
  `Loser (non-significant)` = category_palette[["non_significant_loser"]],
  `Winner (non-significant)` = category_palette[["non_significant_winner"]],
  `Colonizer (non-significant)` =
    category_palette[["non_significant_colonizer"]]
)
significant_levels <- c("Loser", "Winner", "Colonizer")
significant_palette <- response_palette[significant_levels]

climate_levels <- c(
  "Cold/Dry", "Intermediate/Dry", "Warm/Dry",
  "Cold/Intermediate", "Cold/Wet"
)
climate_palette <- c(
  `Cold/Dry` = "#9CCFE0",
  `Intermediate/Dry` = "#F7E6B5",
  `Warm/Dry` = "#F7C5B8",
  `Cold/Intermediate` = "#A7DDB2",
  `Cold/Wet` = "#5AA6C8"
)

analysis_data <- readr::read_csv(analysis_path, show_col_types = FALSE) |>
  dplyr::mutate(
    category = as.character(.data$analysis_category),
    species = dplyr::coalesce(.data$accepted_name, .data$species),
    experiment_id = as.character(.data$experiment_id),
    family = dplyr::if_else(
      is.na(.data$family) | trimws(.data$family) == "",
      "Unresolved", trimws(.data$family)
    ),
    response_state = factor(
      dplyr::recode(
        .data$category,
        loser = "Loser",
        winner = "Winner",
        colonizer = "Colonizer",
        non_significant_loser = "Loser (non-significant)",
        non_significant_winner = "Winner (non-significant)",
        non_significant_colonizer = "Colonizer (non-significant)"
      ),
      levels = response_levels
    ),
    functional_group_display = factor(
      dplyr::recode(
        .data$functional_group,
        forb = "Forb", graminoid = "Graminoid",
        other_unresolved = "Other/unresolved"
      ),
      levels = c("Forb", "Graminoid", "Other/unresolved")
    )
  )

if (nrow(analysis_data) != 2239L ||
    anyDuplicated(dplyr::select(analysis_data, "species", "experiment_id")) ||
    anyNA(analysis_data$response_state)) {
  stop("The species x experiment response snapshot changed.", call. = FALSE)
}

# Shared climate data and reviewed categorical assignments ---------------------
experiment_climate <- analysis_data |>
  dplyr::distinct(
    .data$experiment_id,
    .data$destination_T_JJA, .data$destination_P_JJA, .data$diff_T_JJA
  ) |>
  dplyr::filter(
    is.finite(.data$destination_T_JJA),
    is.finite(.data$destination_P_JJA),
    is.finite(.data$diff_T_JJA)
  )
if (nrow(experiment_climate) != 42L) {
  stop("Expected 42 experiment climate profiles.", call. = FALSE)
}

temperature_edges3 <- seq(
  min(experiment_climate$destination_T_JJA),
  max(experiment_climate$destination_T_JJA), length.out = 4L
)
precipitation_edges3 <- seq(
  min(experiment_climate$destination_P_JJA),
  max(experiment_climate$destination_P_JJA), length.out = 4L
)
experiment_climate <- experiment_climate |>
  dplyr::mutate(
    temperature_context = cut(
      .data$destination_T_JJA, .env$temperature_edges3,
      labels = c("cold", "intermediate", "warm"), include.lowest = TRUE
    ),
    precipitation_context = cut(
      .data$destination_P_JJA, .env$precipitation_edges3,
      labels = c("dry", "intermediate", "wet"), include.lowest = TRUE
    ),
    grid_context = paste0(
      as.character(.data$temperature_context), "_",
      as.character(.data$precipitation_context)
    ),
    display_context = destination_climate_display_context(
      .data$experiment_id, .data$grid_context
    ),
    climate_context = factor(
      dplyr::recode(
        .data$display_context,
        cold_dry = "Cold/Dry",
        intermediate_dry = "Intermediate/Dry",
        warm_dry = "Warm/Dry",
        cold_intermediate = "Cold/Intermediate",
        cold_wet = "Cold/Wet"
      ),
      levels = climate_levels
    )
  )

expected_climate_counts <- c(14L, 9L, 11L, 3L, 5L)
if (anyNA(experiment_climate$climate_context) ||
    !identical(
      as.integer(table(experiment_climate$climate_context)),
      expected_climate_counts
    )) {
  stop("The reviewed climate display assignments changed.", call. = FALSE)
}

analysis_data <- analysis_data |>
  dplyr::left_join(
    dplyr::select(experiment_climate, "experiment_id", "climate_context"),
    by = "experiment_id"
  )

# Outcome-independent family display groups -----------------------------------
common_families <- analysis_data |>
  dplyr::count(.data$family, name = "records") |>
  dplyr::filter(.data$records >= 50L) |>
  dplyr::arrange(dplyr::desc(.data$records), .data$family)
if (nrow(common_families) != 13L) {
  stop("Expected 13 families with at least 50 records.", call. = FALSE)
}
family_levels <- c(common_families$family, "Other")
analysis_data <- analysis_data |>
  dplyr::mutate(
    family_display = factor(
      dplyr::if_else(
        .data$family %in% common_families$family,
        .data$family, "Other"
      ),
      levels = family_levels
    )
  )

# Alluvial figures -------------------------------------------------------------
make_alluvial <- function(data, destination, destination_label,
                          destination_levels, destination_palette,
                          path, width, height, label_size = 4,
                          caption = NULL) {
  destination_symbol <- rlang::sym(destination)
  active_response_levels <- response_levels[
    response_levels %in% as.character(unique(data$response_state))
  ]
  flow_data <- data |>
    dplyr::count(
      .data$response_state,
      destination_group = !!destination_symbol,
      name = "species_experiment_count"
    ) |>
    dplyr::filter(!is.na(.data$destination_group)) |>
    dplyr::mutate(
      response_state = factor(
        .data$response_state, levels = active_response_levels
      ),
      destination_group = factor(
        as.character(.data$destination_group), levels = destination_levels
      ),
      flow_id = dplyr::row_number()
    )
  if (sum(flow_data$species_experiment_count) != nrow(data)) {
    stop(destination_label, " flow does not conserve records.", call. = FALSE)
  }
  combined_palette <- c(response_palette, destination_palette)
  response_totals <- data |>
    dplyr::count(.data$response_state, name = "n")
  response_legend_labels <- stats::setNames(
    active_response_levels,
    active_response_levels
  )
  plot <- ggplot2::ggplot(
    flow_data,
    ggplot2::aes(
      y = .data$species_experiment_count,
      axis1 = .data$response_state,
      axis2 = .data$destination_group
    )
  ) +
    ggalluvial::geom_alluvium(
      ggplot2::aes(fill = .data$response_state),
      width = 0.16, alpha = 0.48, knot.pos = 0.45,
      colour = "white", linewidth = 0.12
    ) +
    ggalluvial::geom_stratum(
      ggplot2::aes(fill = ggplot2::after_stat(stratum)),
      width = 0.18, colour = "white", linewidth = 0.6,
      show.legend = FALSE
    ) +
    ggplot2::scale_x_discrete(
      limits = c("Response", destination_label),
      expand = ggplot2::expansion(add = c(0.12, 0.12))
    ) +
    ggplot2::scale_fill_manual(
      values = combined_palette,
      breaks = active_response_levels,
      labels = response_legend_labels,
      name = NULL
    ) +
    transplant_figure_theme(base_size = 16) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      legend.position = "bottom",
      plot.caption = ggplot2::element_text(
        size = 11, hjust = 0, colour = "grey30"
      )
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        nrow = if (length(active_response_levels) <= 3L) 1 else 2,
        byrow = TRUE, override.aes = list(alpha = 1)
      )
    ) +
    ggplot2::labs(
      x = NULL, y = "Significant species × experiment count",
      caption = caption
    )
  built <- ggplot2::ggplot_build(plot)$data
  stratum_positions <- tibble::as_tibble(built[[2]]) |>
    dplyr::transmute(
      x = .data$x,
      state = as.character(.data$stratum),
      y = .data$y,
      n = .data$ymax - .data$ymin
    )
  response_bar_labels <- stratum_positions |>
    dplyr::filter(.data$x == 1) |>
    dplyr::mutate(label = sprintf("%d", round(.data$n)))
  destination_bar_labels <- stratum_positions |>
    dplyr::filter(.data$x == 2) |>
    dplyr::mutate(
      percentage = .data$n / sum(.data$n),
      label = sprintf(
        "%s\n%d (%.0f%%)", .data$state, round(.data$n),
        100 * .data$percentage
      )
    )
  flow_end_positions <- tibble::as_tibble(built[[1]]) |>
    dplyr::filter(.data$x == 2) |>
    dplyr::transmute(
      flow_id = as.integer(.data$alluvium),
      endpoint_y = .data$y
    )
  path_labels <- flow_data |>
    dplyr::mutate(
      destination_group = as.character(.data$destination_group),
      response_state = factor(
        as.character(.data$response_state), levels = active_response_levels
      )
    ) |>
    dplyr::left_join(flow_end_positions, by = "flow_id") |>
    dplyr::group_by(.data$response_state) |>
    dplyr::mutate(
      response_total = sum(.data$species_experiment_count),
      percentage = .data$species_experiment_count / .data$response_total,
      percentage_label = dplyr::case_when(
        .data$species_experiment_count == 0L ~ "0%",
        100 * .data$percentage < 1 ~ sprintf(
          "%.1f%%", 100 * .data$percentage
        ),
        TRUE ~ sprintf("%.0f%%", 100 * .data$percentage)
      ),
      label = sprintf(
        "%d (%s)", .data$species_experiment_count,
        .data$percentage_label
      ),
      label_x = 1.78
    ) |>
    dplyr::ungroup()
  if (nrow(response_bar_labels) != 3L ||
      abs(sum(response_bar_labels$n) - nrow(data)) > 1e-7 ||
      nrow(destination_bar_labels) != length(destination_levels) ||
      abs(sum(destination_bar_labels$n) - nrow(data)) > 1e-7 ||
      abs(sum(destination_bar_labels$percentage) - 1) > 1e-7 ||
      nrow(path_labels) !=
        length(destination_levels) * length(active_response_levels) ||
      anyNA(path_labels$endpoint_y) ||
      sum(path_labels$species_experiment_count) != nrow(data) ||
      any(abs(
        path_labels |>
          dplyr::summarise(total = sum(.data$percentage),
                           .by = "response_state") |>
          dplyr::pull("total") - 1
      ) > 1e-7)) {
    stop(destination_label, " annotation denominators changed.", call. = FALSE)
  }
  plot <- plot +
    ggplot2::geom_label(
      data = response_bar_labels,
      mapping = ggplot2::aes(
        x = .data$x, y = .data$y, label = .data$label
      ),
      inherit.aes = FALSE, size = label_size, fontface = "bold",
      lineheight = 0.9, label.padding = grid::unit(0.11, "lines"),
      linewidth = 0, fill = scales::alpha("white", 0.90),
      colour = "grey10"
    ) +
    ggplot2::geom_label(
      data = destination_bar_labels,
      mapping = ggplot2::aes(
        x = .data$x, y = .data$y, label = .data$label
      ),
      inherit.aes = FALSE, size = label_size, fontface = "bold",
      lineheight = 0.88, label.padding = grid::unit(0.10, "lines"),
      linewidth = 0, fill = scales::alpha("white", 0.88),
      colour = "grey10"
    ) +
    ggrepel::geom_label_repel(
      data = path_labels,
      mapping = ggplot2::aes(
        x = .data$label_x, y = .data$endpoint_y,
        label = .data$label, colour = .data$response_state
      ),
      inherit.aes = FALSE, size = 3.25, fontface = "bold",
      direction = "y", xlim = c(1.57, 1.94),
      box.padding = 0.12, point.padding = 0.05,
      min.segment.length = 0, segment.colour = "grey62",
      segment.size = 0.28, force = 0.7, max.overlaps = Inf,
      label.padding = grid::unit(0.08, "lines"),
      linewidth = 0, fill = scales::alpha("white", 0.88), seed = 47163
    ) +
    ggplot2::scale_colour_manual(
      values = significant_palette, guide = "none"
    )
  save_plot_atomic(plot, path, width, height)
}

make_climate_response_percentage_bars <- function(
    data, climate_axis_path, response_axis_path) {
  bar_data <- data |>
    dplyr::count(
      .data$climate_context, .data$response_state,
      name = "species_experiment_count"
    ) |>
    tidyr::complete(
      climate_context = factor(climate_levels, levels = climate_levels),
      response_state = factor(significant_levels, levels = significant_levels),
      fill = list(species_experiment_count = 0L)
    ) |>
    dplyr::mutate(
      climate_context = factor(
        as.character(.data$climate_context), levels = climate_levels
      ),
      response_state = factor(
        as.character(.data$response_state), levels = significant_levels
      )
    ) |>
    dplyr::group_by(.data$response_state) |>
    dplyr::mutate(
      response_total = sum(.data$species_experiment_count),
      percentage = .data$species_experiment_count / .data$response_total
    ) |>
    dplyr::ungroup()

  expected_totals <- c(Loser = 132L, Winner = 69L, Colonizer = 34L)
  observed_totals <- bar_data |>
    dplyr::summarise(
      response_total = sum(.data$species_experiment_count),
      percentage_total = sum(.data$percentage),
      .by = "response_state"
    ) |>
    dplyr::arrange(match(as.character(.data$response_state), significant_levels))
  if (nrow(bar_data) != length(climate_levels) * length(significant_levels) ||
      sum(bar_data$species_experiment_count) != 235L ||
      !identical(
        stats::setNames(observed_totals$response_total,
                        as.character(observed_totals$response_state)),
        expected_totals
      ) ||
      any(abs(observed_totals$percentage_total - 1) > 1e-7)) {
    stop("Climate response-percentage denominators changed.", call. = FALSE)
  }

  climate_axis_plot <- ggplot2::ggplot(
    bar_data,
    ggplot2::aes(
      x = .data$climate_context, y = .data$percentage,
      fill = .data$response_state, group = .data$response_state
    )
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.82),
      width = 0.72, colour = "white", linewidth = 0.45
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$species_experiment_count),
      position = ggplot2::position_dodge(width = 0.82),
      vjust = -0.45, size = 4.7, fontface = "bold", colour = "grey15"
    ) +
    ggplot2::scale_x_discrete(
      limits = climate_levels,
      labels = c(
        "Cold/Dry", "Intermediate/\nDry", "Warm/Dry",
        "Cold/\nIntermediate", "Cold/Wet"
      ),
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0, 0.13))
    ) +
    ggplot2::scale_fill_manual(
      values = significant_palette,
      limits = significant_levels, drop = FALSE, name = NULL
    ) +
    transplant_figure_theme(base_size = 17) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 13, lineheight = 0.95),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom"
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, byrow = TRUE)) +
    ggplot2::labs(
      x = "Destination climate",
      y = "Percentage of response category"
    )
  save_plot_atomic(climate_axis_plot, climate_axis_path, 12, 7)

  response_axis_plot <- ggplot2::ggplot(
    bar_data,
    ggplot2::aes(
      x = .data$response_state, y = .data$percentage,
      fill = .data$climate_context, group = .data$climate_context
    )
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.86),
      width = 0.76, colour = "white", linewidth = 0.45
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$species_experiment_count),
      position = ggplot2::position_dodge(width = 0.86),
      vjust = -0.45, size = 4.5, fontface = "bold", colour = "grey15"
    ) +
    ggplot2::scale_x_discrete(
      limits = significant_levels, labels = significant_levels, drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0, 0.13))
    ) +
    ggplot2::scale_fill_manual(
      values = climate_palette,
      limits = climate_levels, drop = FALSE, name = NULL,
      labels = c(
        "Cold/Dry", "Intermediate/Dry", "Warm/Dry",
        "Cold/Intermediate", "Cold/Wet"
      )
    ) +
    transplant_figure_theme(base_size = 17) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 14),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom"
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, byrow = TRUE)) +
    ggplot2::labs(
      x = "Response category",
      y = "Percentage of response category"
    )
  save_plot_atomic(response_axis_plot, response_axis_path, 12, 7)
}

make_functional_family_alluvial <- function(data, path, caption = NULL) {
  active_response_levels <- response_levels[
    response_levels %in% as.character(unique(data$response_state))
  ]
  unavailable_state <- "not_displayed"
  display_data <- data |>
    dplyr::mutate(
      response_state = factor(
        .data$response_state, levels = active_response_levels
      ),
      response_order = match(
        as.character(.data$response_state), active_response_levels
      ),
      pool_state = dplyr::if_else(
        .data$pool == "resident", "Resident pool", "Nonresident pool"
      ),
      functional_state = dplyr::if_else(
        as.character(.data$functional_group_display) == "Other/unresolved",
        .env$unavailable_state,
        as.character(.data$functional_group_display)
      ),
      family_state = dplyr::if_else(
        as.character(.data$functional_group_display) == "Other/unresolved" |
          as.character(.data$family_display) == "Other",
        .env$unavailable_state,
        as.character(.data$family_display)
      )
    )
  omitted_functional <- sum(
    as.character(data$functional_group_display) == "Other/unresolved"
  )
  omitted_family <- sum(as.character(data$family_display) == "Other")
  omitted_union <- sum(
    as.character(data$functional_group_display) == "Other/unresolved" |
      as.character(data$family_display) == "Other"
  )
  if (nrow(display_data) != nrow(data) || anyNA(display_data$response_order) ||
      omitted_functional != 12L || omitted_family != 31L ||
      omitted_union != 38L ||
      sum(display_data$pool_state == "Resident pool") != 201L ||
      sum(display_data$pool_state == "Nonresident pool") != 34L) {
    stop(
      "The hidden functional-group/family snapshot changed.", call. = FALSE
    )
  }
  visible_family_map <- display_data |>
    dplyr::filter(
      .data$family_state != .env$unavailable_state,
      .data$functional_state != .env$unavailable_state
    ) |>
    dplyr::distinct(.data$family_state, .data$functional_state)
  duplicate_family_groups <- visible_family_map |>
    dplyr::count(.data$family_state) |>
    dplyr::filter(.data$n != 1L)
  if (nrow(duplicate_family_groups) ||
      any(!visible_family_map$functional_state %in% c("Forb", "Graminoid"))) {
    stop(
      "Every displayed family must map to exactly one functional group.",
      call. = FALSE
    )
  }
  visible_family_levels <- display_data |>
    dplyr::filter(.data$family_state != .env$unavailable_state) |>
    dplyr::count(
      .data$functional_state, .data$family_state,
      name = "significant_count"
    ) |>
    dplyr::mutate(
      functional_order = match(
        .data$functional_state, c("Forb", "Graminoid")
      )
    ) |>
    dplyr::arrange(
      .data$functional_order, dplyr::desc(.data$significant_count),
      .data$family_state
    ) |>
    dplyr::pull("family_state")
  if (!setequal(visible_family_levels, visible_family_map$family_state)) {
    stop("The ordered family set does not match the visible family set.", call. = FALSE)
  }
  response_legend_labels <- stats::setNames(
    active_response_levels, active_response_levels
  )
  state_levels <- c(
    "Resident pool", "Nonresident pool",
    active_response_levels,
    setdiff(functional_levels, "Other/unresolved"),
    visible_family_levels,
    unavailable_state
  )
  flow_palette <- c(
    `Resident pool` = "#BDBDBD",
    `Nonresident pool` = "#D0D0D0",
    response_palette,
    functional_palette[names(functional_palette) != "Other/unresolved"],
    stats::setNames(
      rep("#F2F2F2", length(visible_family_levels)),
      visible_family_levels
    ),
    stats::setNames("white", unavailable_state)
  )
  stratum_lodes <- dplyr::bind_rows(
    display_data |>
      dplyr::transmute(
        x = 1, state = .data$pool_state,
        weight = 1, stratum_alpha = 1
      ),
    display_data |>
      dplyr::transmute(
        x = 1.38, state = as.character(.data$response_state),
        weight = 1, stratum_alpha = 1
      ),
    display_data |>
      dplyr::transmute(
        x = 2.55, state = .data$functional_state, weight = 1,
        stratum_alpha = dplyr::if_else(
          .data$state == .env$unavailable_state, 0, 1
        )
      ),
    display_data |>
      dplyr::transmute(
        x = 3.70, state = .data$family_state, weight = 1,
        stratum_alpha = dplyr::if_else(
          .data$state == .env$unavailable_state, 0, 1
        )
      )
  ) |>
    dplyr::mutate(state = factor(.data$state, levels = state_levels))

  make_pair_lodes <- function(source_column, destination_column,
                              source_x, destination_x, pair_id) {
    display_data |>
      dplyr::transmute(
        response_state = .data$response_state,
        response_order = .data$response_order,
        source_state = as.character(.data[[source_column]]),
        destination_state = as.character(.data[[destination_column]])
      ) |>
      dplyr::count(
        .data$response_state, .data$response_order,
        .data$source_state, .data$destination_state,
        name = "species_experiment_count"
      ) |>
      dplyr::mutate(
        pair_id = .env$pair_id,
        flow_id = paste(
          .env$pair_id, .data$response_state, .data$source_state,
          .data$destination_state, sep = "::"
        ),
        flow_order = 100L * .data$response_order + dplyr::if_else(
          .data$source_state == .env$unavailable_state |
            .data$destination_state == .env$unavailable_state,
          99L, 0L
        ),
        flow_alpha = dplyr::if_else(
          .data$source_state == .env$unavailable_state |
            .data$destination_state == .env$unavailable_state,
          0, 0.46
        )
      ) |>
      tidyr::pivot_longer(
        c("source_state", "destination_state"),
        names_to = "endpoint", values_to = "state"
      ) |>
      dplyr::mutate(
        x = dplyr::if_else(
          .data$endpoint == "source_state", source_x, destination_x
        ),
        state = factor(.data$state, levels = state_levels)
      )
  }
  pair_lodes <- dplyr::bind_rows(
    make_pair_lodes(
      "pool_state", "response_state", 1, 1.38,
      "pool_to_response"
    ),
    make_pair_lodes(
      "response_state", "functional_state", 1.38, 2.55,
      "response_to_functional"
    ),
    make_pair_lodes(
      "functional_state", "family_state", 2.55, 3.70,
      "functional_to_family"
    )
  )
  flow_layer <- function(pair_id) {
    ggalluvial::geom_flow(
      data = dplyr::filter(pair_lodes, .data$pair_id == .env$pair_id),
      mapping = ggplot2::aes(
        x = .data$x, stratum = .data$state,
        alluvium = .data$flow_id,
        y = .data$species_experiment_count,
        fill = .data$response_state,
        alpha = .data$flow_alpha,
        order = .data$flow_order
      ),
      inherit.aes = FALSE, width = 0.26, knot.pos = 0.42,
      # Do not outline ribbons: transparent spacer lodes for omitted
      # classifications otherwise leave thin white "ghost" paths.
      colour = "transparent", linewidth = 0.1, aes.flow = "forward"
    )
  }
  plot <- ggplot2::ggplot() +
    flow_layer("pool_to_response") +
    flow_layer("response_to_functional") +
    flow_layer("functional_to_family") +
    ggalluvial::geom_stratum(
      data = dplyr::filter(stratum_lodes, .data$x != 1),
      mapping = ggplot2::aes(
        x = .data$x, stratum = .data$state, y = .data$weight,
        fill = .data$state, alpha = .data$stratum_alpha
      ),
      inherit.aes = FALSE,
      width = 0.34, colour = "white", linewidth = 0.8,
      show.legend = FALSE
    ) +
    ggalluvial::geom_stratum(
      data = dplyr::filter(stratum_lodes, .data$x == 1),
      mapping = ggplot2::aes(
        x = .data$x, stratum = .data$state, y = .data$weight,
        fill = .data$state, alpha = .data$stratum_alpha
      ),
      inherit.aes = FALSE,
      width = 0.40, colour = "white", linewidth = 0.8,
      show.legend = FALSE
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(1, 1.38, 2.55, 3.70),
      labels = c("Pool", "Response", "Functional group", "Family"),
      limits = c(0.45, 4.18),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::scale_fill_manual(
      values = flow_palette,
      breaks = active_response_levels,
      labels = response_legend_labels,
      name = NULL
    ) +
    ggplot2::scale_alpha_identity(guide = "none") +
    ggplot2::scale_size_identity() +
    transplant_figure_theme(base_size = 22) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      legend.position = "bottom"
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        nrow = 1, byrow = TRUE, override.aes = list(alpha = 1)
      )
    ) +
    ggplot2::labs(
      x = NULL, y = "Significant species × experiment count"
    )
  built <- ggplot2::ggplot_build(plot)$data

  # The omitted functional-group/family records are needed while calculating
  # node positions, but drawing them with alpha zero creates curved white
  # gaps between otherwise visible ribbons. Repack only the visible source
  # lodes into their existing response-colour blocks, then render those lodes
  # directly. Missing classifications remain as unconnected space beside the
  # source node rather than as a white alluvial path.
  repack_visible_source <- function(flow_data, source_x, target_data,
                                    target_x) {
    visible <- tibble::as_tibble(flow_data) |>
      dplyr::filter(.data$alpha > 0)
    target_ranges <- tibble::as_tibble(target_data) |>
      dplyr::filter(abs(.data$x - .env$target_x) < 1e-7) |>
      dplyr::summarise(
        target_ymin = min(.data$ymin),
        target_ymax = max(.data$ymax),
        .by = c("stratum", "fill")
      )
    source_adjustment <- visible |>
      dplyr::filter(abs(.data$x - .env$source_x) < 1e-7) |>
      dplyr::arrange(.data$stratum, .data$fill, .data$ymin) |>
      dplyr::mutate(flow_height = .data$ymax - .data$ymin) |>
      dplyr::left_join(target_ranges, by = c("stratum", "fill")) |>
      dplyr::group_by(.data$stratum, .data$fill) |>
      dplyr::mutate(
        visible_total = sum(.data$flow_height),
        packed_ymin = .data$target_ymax - .data$visible_total +
          dplyr::lag(cumsum(.data$flow_height), default = 0),
        packed_ymax = .data$packed_ymin + .data$flow_height,
        packed_y = (.data$packed_ymin + .data$packed_ymax) / 2
      ) |>
      dplyr::ungroup() |>
      dplyr::select(
        "group", "packed_ymin", "packed_ymax", "packed_y",
        "target_ymin", "target_ymax"
      )
    if (anyNA(source_adjustment) ||
        any(source_adjustment$packed_ymin <
              source_adjustment$target_ymin - 1e-7) ||
        any(source_adjustment$packed_ymax >
              source_adjustment$target_ymax + 1e-7)) {
      stop("Visible alluvial source repacking failed.", call. = FALSE)
    }
    visible |>
      dplyr::left_join(source_adjustment, by = "group") |>
      dplyr::mutate(
        ymin = dplyr::if_else(
          abs(.data$x - .env$source_x) < 1e-7,
          .data$packed_ymin, .data$ymin
        ),
        ymax = dplyr::if_else(
          abs(.data$x - .env$source_x) < 1e-7,
          .data$packed_ymax, .data$ymax
        ),
        y = dplyr::if_else(
          abs(.data$x - .env$source_x) < 1e-7,
          .data$packed_y, .data$y
        ),
        response_state = names(significant_palette)[
          match(as.character(.data$fill), unname(significant_palette))
        ],
        alluvium = .data$group
      ) |>
      dplyr::select(
        -dplyr::starts_with("packed_"),
        -dplyr::starts_with("target_")
      )
  }
  visible_identity_flows <- list(
    tibble::as_tibble(built[[1]]) |>
      dplyr::mutate(
        response_state = names(significant_palette)[
          match(as.character(.data$fill), unname(significant_palette))
        ],
        alluvium = .data$group
      ),
    repack_visible_source(built[[2]], 1.38, built[[1]], 1.38),
    repack_visible_source(built[[3]], 2.55, built[[2]], 2.55)
  )
  if (any(vapply(
        visible_identity_flows,
        function(layer) anyNA(layer$response_state), logical(1)
      )) ||
      !identical(
        vapply(
          visible_identity_flows,
          function(layer) sum(
            layer$ymax[layer$flow == "from"] -
              layer$ymin[layer$flow == "from"]
          ), numeric(1)
        ),
        c(235, 223, 197)
      )) {
    stop("Visible alluvial flow conservation failed.", call. = FALSE)
  }
  identity_flow_layer <- function(layer_data) {
    suppressWarnings(
      ggalluvial::geom_flow(
        data = layer_data, stat = "identity",
        mapping = ggplot2::aes(
          x = .data$x, y = .data$y,
          ymin = .data$ymin, ymax = .data$ymax,
          alluvium = .data$alluvium,
          fill = .data$response_state,
          colour = .data$response_state,
          alpha = .data$alpha
        ),
        inherit.aes = FALSE, knot.pos = 0.42,
        linewidth = 0.18, outline.type = "full",
        aes.flow = "forward"
      )
    )
  }
  plot$layers <- c(
    lapply(visible_identity_flows, identity_flow_layer),
    plot$layers[4:5]
  )

  pool_bar_labels <- tibble::as_tibble(built[[5]]) |>
    dplyr::filter(.data$x == 1) |>
    dplyr::transmute(
      state = as.character(.data$stratum),
      x = .data$x, y = .data$y, n = .data$ymax - .data$ymin
    ) |>
    dplyr::mutate(
      eligible_pool_n = dplyr::if_else(
        .data$state == "Resident pool", 1625L, 614L
      ),
      significant_percentage = .data$n / .data$eligible_pool_n,
      label = sprintf(
        "%s\n%d/%d (%.1f%%)",
        dplyr::if_else(
          .data$state == "Resident pool", "Residents", "Nonresidents"
        ),
        round(.data$n), .data$eligible_pool_n,
        100 * .data$significant_percentage
      )
    )
  response_bar_labels <- tibble::as_tibble(built[[4]]) |>
    dplyr::filter(abs(.data$x - 1.38) < 1e-7) |>
    dplyr::transmute(
      state = as.character(.data$stratum),
      x = .data$x, y = .data$y,
      n = .data$ymax - .data$ymin,
      label = sprintf("%d", round(.data$n))
    )
  functional_bar_labels <- tibble::as_tibble(built[[4]]) |>
    dplyr::filter(
      abs(.data$x - 2.55) < 1e-7,
      as.character(.data$stratum) != .env$unavailable_state
    ) |>
    dplyr::transmute(
      x = .data$x, y = .data$y,
      label = paste0(
        as.character(.data$stratum), "\n",
        format(
          round(.data$ymax - .data$ymin),
          scientific = FALSE, trim = TRUE
        )
      ),
      n = .data$ymax - .data$ymin
    )
  functional_flow_labels <- tibble::as_tibble(built[[2]]) |>
    dplyr::transmute(
      PANEL = .data$PANEL, group = .data$group,
      x = .data$x, flow = .data$flow,
      state = as.character(.data$stratum),
      n = .data$ymax - .data$ymin, y = .data$y,
      alpha = .data$alpha
    ) |>
    dplyr::group_by(.data$PANEL, .data$group) |>
    dplyr::filter(dplyr::n() == 2L) |>
    dplyr::summarise(
      source_state = .data$state[.data$flow == "from"],
      destination_state = .data$state[.data$flow == "to"],
      source_y = .data$y[.data$flow == "from"],
      destination_y = .data$y[.data$flow == "to"],
      n = dplyr::first(.data$n),
      visible = all(.data$alpha > 0),
      .groups = "drop"
    ) |>
    dplyr::filter(.data$visible) |>
    dplyr::group_by(.data$source_state) |>
    dplyr::mutate(
      percentage = .data$n / sum(.data$n),
      percentage_label = dplyr::case_when(
        .data$percentage == 0 ~ "0%",
        100 * .data$percentage < 1 ~ sprintf(
          "%.1f%%", 100 * .data$percentage
        ),
        TRUE ~ sprintf("%.0f%%", 100 * .data$percentage)
      ),
      label = sprintf("%d (%s)", round(.data$n), .data$percentage_label),
      label_x = mean(c(1.38, 2.55)),
      label_y = (.data$source_y + .data$destination_y) / 2
    ) |>
    dplyr::ungroup()
  response_totals <- display_data |>
    dplyr::count(.data$response_state, name = "response_total")
  family_flow_labels <- tibble::as_tibble(built[[3]]) |>
    dplyr::transmute(
      PANEL = .data$PANEL, group = .data$group,
      flow = .data$flow, state = as.character(.data$stratum),
      n = .data$ymax - .data$ymin, y = .data$y,
      alpha = .data$alpha, fill = as.character(.data$fill)
    ) |>
    dplyr::group_by(.data$PANEL, .data$group) |>
    dplyr::filter(dplyr::n() == 2L) |>
    dplyr::summarise(
      functional_state = .data$state[.data$flow == "from"],
      family_state = .data$state[.data$flow == "to"],
      source_y = .data$y[.data$flow == "from"],
      destination_y = .data$y[.data$flow == "to"],
      n = dplyr::first(.data$n),
      visible = all(.data$alpha > 0),
      fill = dplyr::first(.data$fill),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      response_state = names(significant_palette)[
        match(.data$fill, unname(significant_palette))
      ]
    ) |>
    dplyr::filter(
      .data$visible,
      .data$family_state != .env$unavailable_state,
      !is.na(.data$response_state)
    ) |>
    dplyr::left_join(
      dplyr::mutate(
        response_totals,
        response_state = as.character(.data$response_state)
      ),
      by = "response_state"
    ) |>
    dplyr::mutate(
      percentage = .data$n / .data$response_total,
      label = sprintf(
        "%d (%s)", round(.data$n),
        dplyr::if_else(
          100 * .data$percentage < 1,
          sprintf("%.1f%%", 100 * .data$percentage),
          sprintf("%.0f%%", 100 * .data$percentage)
        )
      ),
      label_x = mean(c(2.55, 3.70)),
      label_y = (.data$source_y + .data$destination_y) / 2
    ) |>
    dplyr::filter(.data$percentage >= 0.05)
  family_centres <- tibble::as_tibble(built[[4]]) |>
    dplyr::filter(
      abs(.data$x - 3.70) < 1e-7,
      as.character(.data$stratum) != .env$unavailable_state
    ) |>
    dplyr::transmute(
      family_state = as.character(.data$stratum), family_y = .data$y,
      family_n = .data$ymax - .data$ymin
    ) |>
    dplyr::arrange(dplyr::desc(.data$family_y))
  compact_family_positions <- family_centres$family_y
  small_family_indices <- which(family_centres$family_n < 10)
  if (length(small_family_indices) > 1L) {
    for (small_index in seq_along(small_family_indices)[-1L]) {
      family_index <- small_family_indices[[small_index]]
      previous_index <- small_family_indices[[small_index - 1L]]
      compact_family_positions[[family_index]] <- min(
        family_centres$family_y[[family_index]],
        compact_family_positions[[previous_index]] - 9.5
      )
    }
  }
  family_centres <- family_centres |>
    dplyr::mutate(
      family_block_y = .env$compact_family_positions,
      small_family = .data$family_n < 10,
      show_leader = .data$small_family
    )
  family_response_labels <- display_data |>
    dplyr::filter(.data$family_state != .env$unavailable_state) |>
    dplyr::count(
      .data$family_state, .data$response_state,
      name = "species_experiment_count"
    ) |>
    tidyr::complete(
      family_state = visible_family_levels,
      response_state = factor(
        active_response_levels, levels = active_response_levels
      ),
      fill = list(species_experiment_count = 0L)
    ) |>
    dplyr::mutate(
      family_state = factor(
        .data$family_state, levels = visible_family_levels
      ),
      response_state = factor(
        .data$response_state, levels = active_response_levels
      )
    ) |>
    dplyr::group_by(.data$family_state) |>
    dplyr::mutate(
      family_total = sum(.data$species_experiment_count),
      percentage = .data$species_experiment_count / .data$family_total,
      percentage_label = dplyr::case_when(
        .data$species_experiment_count == 0L ~ "0%",
        100 * .data$percentage < 1 ~ sprintf(
          "%.1f%%", 100 * .data$percentage
        ),
        TRUE ~ sprintf("%.0f%%", 100 * .data$percentage)
      ),
      label = sprintf("%d", .data$species_experiment_count)
    ) |>
    dplyr::ungroup() |>
    dplyr::left_join(family_centres, by = "family_state") |>
    dplyr::mutate(
      label_x = 4.14,
      line_offset = dplyr::recode(
        as.character(.data$response_state),
        Loser = 3.3, Winner = 0, Colonizer = -3.3
      ),
      label_y = .data$family_block_y + .data$line_offset
    )
  family_name_check <- display_data |>
    dplyr::filter(.data$family_state != .env$unavailable_state) |>
    dplyr::count(.data$family_state, name = "expected_family_n") |>
    dplyr::left_join(family_centres, by = "family_state")
  if (nrow(pool_bar_labels) != 2L ||
      abs(sum(pool_bar_labels$n) - 235) > 1e-7 ||
      sum(pool_bar_labels$eligible_pool_n) != 2239L ||
      nrow(functional_flow_labels) != 6L ||
      nrow(family_flow_labels) != 14L ||
      any(family_flow_labels$percentage < 0.05) ||
      anyNA(family_flow_labels$response_total) ||
      nrow(response_bar_labels) != 3L ||
      abs(sum(response_bar_labels$n) - 235) > 1e-7 ||
      nrow(functional_bar_labels) != 2L ||
      sum(functional_bar_labels$n) != 223L ||
      abs(sum(functional_flow_labels$n) - 223) > 1e-7 ||
      nrow(family_centres) != length(visible_family_levels) ||
      nrow(family_name_check) != length(visible_family_levels) ||
      anyNA(family_name_check$family_y) ||
      any(abs(
        family_name_check$expected_family_n - family_name_check$family_n
      ) > 1e-7) ||
      !setequal(family_name_check$family_state, visible_family_levels) ||
      nrow(family_response_labels) !=
        length(visible_family_levels) * 3L ||
      sum(family_response_labels$species_experiment_count) != 197L ||
      any(abs(
        family_response_labels |>
          dplyr::summarise(total = sum(.data$percentage),
                           .by = "family_state") |>
          dplyr::pull("total") - 1
      ) > 1e-7)) {
    stop("The alluvial annotation denominators changed.", call. = FALSE)
  }
  plot <- plot +
    ggplot2::geom_label(
      data = pool_bar_labels,
      mapping = ggplot2::aes(
        x = .data$x, y = .data$y, label = .data$label
      ),
      inherit.aes = FALSE, size = 4.7, fontface = "bold",
      lineheight = 0.90, label.padding = grid::unit(0.13, "lines"),
      linewidth = 0, fill = scales::alpha("white", 0.92),
      colour = "grey10"
    ) +
    ggplot2::geom_label(
      data = response_bar_labels,
      mapping = ggplot2::aes(
        x = .data$x, y = .data$y, label = .data$label
      ),
      inherit.aes = FALSE, size = 4.7, fontface = "bold",
      lineheight = 0.88, label.padding = grid::unit(0.12, "lines"),
      linewidth = 0, fill = scales::alpha("white", 0.90),
      colour = "grey10"
    ) +
    ggplot2::geom_label(
      data = functional_bar_labels,
      mapping = ggplot2::aes(
        x = .data$x, y = .data$y, label = .data$label
      ),
      inherit.aes = FALSE, size = 4.6, fontface = "bold",
      lineheight = 0.88, label.padding = grid::unit(0.12, "lines"),
      linewidth = 0, fill = scales::alpha("white", 0.90),
      colour = "grey10"
    ) +
    ggplot2::geom_label(
      data = functional_flow_labels,
      mapping = ggplot2::aes(
        x = .data$label_x, y = .data$label_y,
        label = .data$label, colour = .data$source_state
      ),
      inherit.aes = FALSE, size = 4.2, fontface = "bold",
      fill = scales::alpha("white", 0.86),
      linewidth = 0, label.padding = grid::unit(0.08, "lines")
    ) +
    ggrepel::geom_label_repel(
      data = family_flow_labels,
      mapping = ggplot2::aes(
        x = .data$label_x, y = .data$label_y,
        label = .data$label, colour = .data$response_state
      ),
      inherit.aes = FALSE, size = 3.7, fontface = "bold",
      direction = "y", xlim = c(3.00, 3.28),
      box.padding = 0.08, point.padding = 0.03,
      min.segment.length = 0, segment.colour = "grey65",
      segment.size = 0.25, force = 0.35, max.overlaps = Inf,
      label.padding = grid::unit(0.06, "lines"), linewidth = 0,
      fill = scales::alpha("white", 0.84), seed = 47163
    ) +
    ggplot2::geom_segment(
      data = dplyr::filter(family_centres, .data$show_leader),
      mapping = ggplot2::aes(
        x = 3.87, y = .data$family_y,
        xend = 3.895, yend = .data$family_block_y
      ),
      inherit.aes = FALSE, colour = "grey72", linewidth = 0.45
    ) +
    ggplot2::geom_label(
      data = dplyr::filter(family_centres, !.data$small_family),
      mapping = ggplot2::aes(
        x = 3.70, y = .data$family_y,
        label = .data$family_state
      ),
      inherit.aes = FALSE, size = 4.2, fontface = "bold",
      lineheight = 0.88, label.padding = grid::unit(0.10, "lines"),
      linewidth = 0, fill = scales::alpha("white", 0.88),
      colour = "grey12"
    ) +
    ggplot2::geom_label(
      data = dplyr::filter(family_centres, .data$small_family),
      mapping = ggplot2::aes(
        x = 3.90, y = .data$family_block_y,
        label = .data$family_state
      ),
      inherit.aes = FALSE, size = 3.8, fontface = "bold", hjust = 0,
      lineheight = 0.88, label.padding = grid::unit(0.08, "lines"),
      linewidth = 0, fill = scales::alpha("white", 0.92),
      colour = "grey12"
    ) +
    ggplot2::scale_colour_manual(
      values = significant_palette, guide = "none"
    )
  validate_node_alignment <- function(left_layer, right_layer, node_x,
                                      node_label) {
    left_positions <- tibble::as_tibble(left_layer) |>
      dplyr::filter(.data$x == .env$node_x, .data$flow == "to") |>
      dplyr::summarise(
        ymin = min(.data$ymin), ymax = max(.data$ymax),
        .by = c("stratum", "fill")
      )
    right_positions <- tibble::as_tibble(right_layer) |>
      dplyr::filter(.data$x == .env$node_x, .data$flow == "from") |>
      dplyr::summarise(
        ymin = min(.data$ymin), ymax = max(.data$ymax),
        .by = c("stratum", "fill")
      )
    aligned_positions <- dplyr::full_join(
      left_positions, right_positions,
      by = c("stratum", "fill"), suffix = c("_left", "_right")
    )
    if (any(!stats::complete.cases(aligned_positions)) ||
        any(abs(aligned_positions$ymin_left -
                aligned_positions$ymin_right) > 1e-7) ||
        any(abs(aligned_positions$ymax_left -
                aligned_positions$ymax_right) > 1e-7)) {
      stop(
        paste("The", node_label, "alluvial node does not align."),
        call. = FALSE
      )
    }
  }
  validate_node_alignment(built[[1]], built[[2]], 1.38, "response")
  validate_node_alignment(built[[2]], built[[3]], 2.55, "functional-group")
  save_plot_atomic(plot, path, 28, 20)
}

significant_data <- analysis_data |>
  dplyr::filter(.data$category %in% c("loser", "winner", "colonizer"))
if (nrow(significant_data) != 235L) {
  stop("Expected 235 significant species x experiment responses.", call. = FALSE)
}

make_alluvial(
  significant_data, "climate_context", "Destination climate",
  climate_levels, climate_palette,
  file.path(figure_dir, "response_to_climate_alluvial.pdf"),
  17, 11, 4
)

make_climate_response_percentage_bars(
  significant_data,
  file.path(figure_dir, "response_by_climate_percentages.pdf"),
  file.path(
    figure_dir, "response_by_climate_percentages_response_axis.pdf"
  )
)

functional_levels <- c("Forb", "Graminoid", "Other/unresolved")
functional_palette <- c(
  Forb = "#CC79A7", Graminoid = "#0072B2",
  `Other/unresolved` = "#999999"
)
make_functional_family_alluvial(
  significant_data,
  file.path(figure_dir, "response_to_family_alluvial.pdf"),
  paste0(
    "Families shown have at least 50 species × experiment records across ",
    "resident and nonresident pools; blank space preserves omitted ",
    "classifications without changing the 235-record denominator."
  )
)

# Exact experiment values along three continuous climate axes -----------------
axis_configuration <- tibble::tribble(
  ~climate_variable, ~climate_axis,
  "destination_T_JJA", "Destination JJA temperature (°C)",
  "destination_P_JJA", "Destination JJA precipitation (mm)",
  "diff_T_JJA", "JJA warming (°C)"
)
climate_long <- experiment_climate |>
  dplyr::select(
    "experiment_id",
    dplyr::all_of(axis_configuration$climate_variable)
  ) |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(axis_configuration$climate_variable),
    names_to = "climate_variable", values_to = "climate_value"
  ) |>
  dplyr::left_join(axis_configuration, by = "climate_variable") |>
  dplyr::mutate(
    climate_axis = factor(
      .data$climate_axis, levels = axis_configuration$climate_axis
    )
  )
if (nrow(climate_long) != 42L * 3L ||
    anyDuplicated(dplyr::select(
      climate_long, "experiment_id", "climate_variable"
    ))) {
  stop("The continuous experiment climate table is invalid.", call. = FALSE)
}

make_event_rows <- function(data) {
  dplyr::bind_rows(
    data |>
      dplyr::filter(.data$pool == "resident") |>
      dplyr::mutate(outcome = "Loser", event = .data$category == "loser"),
    data |>
      dplyr::filter(.data$pool == "resident") |>
      dplyr::mutate(outcome = "Winner", event = .data$category == "winner"),
    data |>
      dplyr::filter(.data$pool %in% c("colonizer", "outsider")) |>
      dplyr::mutate(
        outcome = "Colonizer", event = .data$category == "colonizer"
      )
  ) |>
    dplyr::mutate(outcome = factor(.data$outcome, levels = significant_levels))
}

summarise_continuous <- function(
    data, entity_column, entity_order, entity_type) {
  entity_symbol <- rlang::sym(entity_column)
  data |>
    dplyr::filter(as.character(!!entity_symbol) %in% entity_order) |>
    dplyr::inner_join(
      climate_long, by = "experiment_id", relationship = "many-to-many"
    ) |>
    make_event_rows() |>
    dplyr::mutate(entity = as.character(!!entity_symbol)) |>
    dplyr::summarise(
      eligible_pool_count = dplyr::n(),
      significant_count = sum(.data$event),
      .by = c(
        "entity", "experiment_id", "climate_variable",
        "climate_axis", "climate_value", "outcome"
      )
    ) |>
    dplyr::mutate(
      entity_type = entity_type,
      entity = factor(.data$entity, levels = entity_order),
      pool_rate = .data$significant_count / .data$eligible_pool_count
    )
}

species_continuous <- analysis_data |>
  dplyr::mutate(entity = "All species") |>
  summarise_continuous("entity", "All species", "all_species")
family_continuous <- summarise_continuous(
  analysis_data, "family_display", family_levels, "family"
)
functional_continuous <- summarise_continuous(
  analysis_data, "functional_group_display",
  functional_levels, "functional_group"
)

continuous_table <- dplyr::bind_rows(
  species_continuous, family_continuous, functional_continuous
) |>
  dplyr::mutate(
    entity = as.character(.data$entity),
    outcome = as.character(.data$outcome),
    climate_axis = as.character(.data$climate_axis),
    pool_percentage = 100 * .data$pool_rate
  )
write_csv_atomic(
  continuous_table,
  file.path(review_output_dir, "response_climate_gradient_counts_rates.csv")
)

# Final validation -------------------------------------------------------------
if (any(continuous_table$significant_count >
        continuous_table$eligible_pool_count) ||
    any(continuous_table$pool_rate < 0 | continuous_table$pool_rate > 1,
        na.rm = TRUE)) {
  stop("A continuous-axis count or pool rate failed validation.", call. = FALSE)
}
expected_events <- c(Loser = 132L, Winner = 69L, Colonizer = 34L)
species_axis_totals <- species_continuous |>
  dplyr::summarise(
    observed = sum(.data$significant_count),
    .by = c("climate_axis", "outcome")
  )
if (any(species_axis_totals$observed !=
        expected_events[as.character(species_axis_totals$outcome)])) {
  stop("The continuous species panels do not conserve events.", call. = FALSE)
}
if (length(unique(analysis_data$family_display)) != 14L ||
    length(unique(analysis_data$functional_group_display)) != 3L) {
  stop("The family or functional-group display mapping is incomplete.", call. = FALSE)
}

cat("Significant-response alluvials and continuous summaries completed.\n")
cat("Alluvials conserve 235 significant species x experiment records.\n")
cat("Family display: 13 common families plus Other.\n")
cat("Continuous summaries use exact experiment climate values; no bins.\n")
