# ==============================================================================
# Script: 07_climate_faceted_descriptive_review.R
#
# Purpose:
#   Compare significant loser, winner, and colonizer trait distributions across
#   an all-climates reference and the five reviewed destination-climate groups.
#
# Inputs:
#   output/tables/analysis_dataset.csv
#   output/tables/climate_context_experiments.csv
#
# Output:
#   plot/trait_categories_with_overall.pdf
# ==============================================================================

source("config/analysis.R")
source("R/pipeline/utils.R")
source("R/pipeline/figures.R")

config <- transplant_config()
analysis_path <- file.path(config$output_dir, "tables", "analysis_dataset.csv")
context_path <- file.path(
  config$output_dir, "tables", "climate_context_experiments.csv"
)
figure_dir <- config$plot_dir

if (!file.exists(analysis_path)) {
  stop("Missing review input: ", analysis_path, call. = FALSE)
}
if (!file.exists(context_path)) {
  stop("Missing review input: ", context_path, call. = FALSE)
}

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

climate_levels <- c(
  "All climates", "Cold/Dry", "Intermediate/Dry", "Warm/Dry",
  "Cold/Intermediate", "Cold/Wet"
)
response_levels <- c("loser", "winner", "colonizer")
response_palette <- category_palette[response_levels]

build_trait_plot_data <- function(data, climate_map) {
  specific <- data |>
    dplyr::mutate(experiment_id = as.character(.data$experiment_id)) |>
    dplyr::left_join(climate_map, by = "experiment_id") |>
    dplyr::filter(.data$category %in% response_levels) |>
    tidyr::pivot_longer(
      c("LDMC", "LMA", "vegetative_height"),
      names_to = "trait", values_to = "raw_value"
    ) |>
    dplyr::filter(is.finite(.data$raw_value), .data$raw_value > 0) |>
    dplyr::mutate(
      value = log(.data$raw_value),
      climate_plot = as.character(.data$climate_context_label)
    )
  if (anyNA(specific$climate_context_label)) {
    stop("A trait observation is missing its climate context.", call. = FALSE)
  }

  plot_data <- dplyr::bind_rows(
    dplyr::mutate(specific, climate_plot = "All climates"),
    specific
  ) |>
    dplyr::mutate(
      category = factor(.data$category, levels = response_levels),
      category_label = factor(
        unname(category_display_labels[as.character(.data$category)]),
        levels = unname(category_display_labels[response_levels])
      ),
      climate_plot = factor(.data$climate_plot, levels = climate_levels),
      trait = factor(
        .data$trait,
        levels = c("LDMC", "LMA", "vegetative_height"),
        labels = c("LDMC", "LMA", "Vegetative height")
      )
    )

  overall_counts <- plot_data |>
    dplyr::filter(.data$climate_plot == "All climates") |>
    dplyr::count(.data$trait, .data$category, name = "overall_n")
  context_counts <- plot_data |>
    dplyr::filter(.data$climate_plot != "All climates") |>
    dplyr::count(.data$trait, .data$category, name = "context_n")
  comparison <- dplyr::full_join(
    overall_counts, context_counts, by = c("trait", "category")
  ) |>
    dplyr::mutate(
      overall_n = tidyr::replace_na(.data$overall_n, 0L),
      context_n = tidyr::replace_na(.data$context_n, 0L)
    )
  if (any(comparison$overall_n != comparison$context_n)) {
    stop("Trait climate groups do not conserve the all-climate sample size.", call. = FALSE)
  }
  list(data = plot_data, counts = comparison)
}

plot_traits_with_overall <- function(trait_data, path) {
  plot <- ggplot2::ggplot(
    trait_data$data,
    ggplot2::aes(
      .data$climate_plot, .data$value,
      fill = .data$category, colour = .data$category
    )
  ) +
    ggplot2::geom_boxplot(
      outlier.shape = NA, alpha = 0.52, linewidth = 0.55, width = 0.76,
      position = ggplot2::position_dodge2(width = 0.84, preserve = "single")
    ) +
    ggplot2::geom_point(
      alpha = 0.22, size = 0.55,
      position = ggplot2::position_jitterdodge(
        jitter.width = 0.055, jitter.height = 0, dodge.width = 0.76
      )
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(.data$trait), scales = "free_y", nrow = 1, drop = FALSE
    ) +
    ggplot2::scale_x_discrete(
      limits = climate_levels,
      labels = c(
        "All\nclimates", "Cold/Dry", "Intermediate/\nDry", "Warm/Dry",
        "Cold/\nIntermediate", "Cold/Wet"
      ),
      drop = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = response_palette, limits = response_levels, drop = FALSE,
      labels = unname(category_display_labels[response_levels]),
      name = "Response"
    ) +
    ggplot2::scale_colour_manual(
      values = response_palette, limits = response_levels, drop = FALSE,
      labels = unname(category_display_labels[response_levels]),
      name = "Response"
    ) +
    transplant_figure_theme(base_size = 16) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 11.5, lineheight = 0.95),
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      panel.spacing.x = grid::unit(1.2, "lines"),
      strip.text = ggplot2::element_text(face = "bold", size = 15)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        nrow = 1, override.aes = list(alpha = 0.8)
      ),
      colour = "none"
    ) +
    ggplot2::labs(x = NULL, y = "Log-transformed trait value")
  save_plot_atomic(plot, path, 23, 8.5)
}

cat("Reading the current documented analysis and climate-context milestones...\n")
analysis_data <- readr::read_csv(analysis_path, show_col_types = FALSE) |>
  dplyr::mutate(
    category = as.character(.data$analysis_category),
    experiment_id = as.character(.data$experiment_id)
  )
climate_map <- readr::read_csv(context_path, show_col_types = FALSE) |>
  dplyr::transmute(
    experiment_id = as.character(.data$experiment_id),
    climate_context_label = as.character(.data$climate_context_label)
  ) |>
  dplyr::distinct()

assert_unique(climate_map, "experiment_id", "reviewed climate-context map")
if (nrow(climate_map) != 42L) {
  stop("Expected 42 reviewed experiment climate profiles; found ", nrow(climate_map), call. = FALSE)
}
unexpected_contexts <- setdiff(
  unique(climate_map$climate_context_label), climate_levels[-1]
)
if (length(unexpected_contexts)) {
  stop(
    "Unexpected reviewed climate contexts: ",
    paste(unexpected_contexts, collapse = ", "), call. = FALSE
  )
}
if (!setequal(unique(analysis_data$experiment_id), climate_map$experiment_id)) {
  stop("Analysis experiments and the climate-context map do not match.", call. = FALSE)
}

cat("Building and validating the climate-faceted trait summary...\n")
trait_data <- build_trait_plot_data(analysis_data, climate_map)

cat("Writing the retained trait figure...\n")
plot_traits_with_overall(
  trait_data,
  file.path(figure_dir, "trait_categories_with_overall.pdf")
)

cat("\nClimate-faceted trait figure completed.\n")
cat("Experiments assigned to one climate context: ", nrow(climate_map), "\n", sep = "")
