# Module: R/figures/climate
# Responsibility: Implement the named scientific workflow stage or shared infrastructure.
# Inputs: Function arguments and explicitly declared target dependencies.
# Outputs: Validated in-memory objects and declared workflow artifacts.
# Side effects: Limited to functions whose contracts explicitly write, copy, or render files.

# Canonical destination-climate display figure functions.

#' Destination Climate Display Context.
#'
#' @param experiment_id A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param grid_context A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
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

#' Plot Destination Climate Warming Gradient.
#'
#' @param data A data frame whose required columns are checked before use.
#' @param climate_contexts A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param path A project-relative or absolute filesystem path.
#' @return A plot object, or the declared figure path when rendering to disk.
plot_destination_climate_warming_gradient <- function(data, climate_contexts, path) {
  plot_data <- data |>
    dplyr::ungroup() |>
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

#' Plot Destination Climate Warming Gradient Continuous.
#'
#' @param data A data frame whose required columns are checked before use.
#' @param path A project-relative or absolute filesystem path.
#' @param base_size A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param experiment_label_size A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param axis_title_size A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param axis_text_size A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param legend_title_size A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param legend_text_size A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param width A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param height A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A plot object, or the declared figure path when rendering to disk.
plot_destination_climate_warming_gradient_continuous <- function(
    data,
    path,
    base_size = 25,
    experiment_label_size = 5.0,
    axis_title_size = 25,
    axis_text_size = 22,
    legend_title_size = 23,
    legend_text_size = 21,
    width = 18,
    height = 12) {
  plot_data <- data |>
    dplyr::ungroup() |>
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
  if (!nrow(plot_data)) {
    stop("The destination-climate figure has no finite experiment rows.", call. = FALSE)
  }
  if (anyDuplicated(plot_data$experiment_id)) {
    stop("The destination-climate figure requires one row per experiment.", call. = FALSE)
  }
  if (!isTRUE(all.equal(
    plot_data$diff_T_JJA,
    plot_data$destination_T_JJA - plot_data$origin_T_JJA,
    tolerance = 1e-10
  ))) {
    stop("diff_T_JJA must equal destination_T_JJA - origin_T_JJA.", call. = FALSE)
  }

  warming_colours <- c("#8C6D00", "#D8A600", "#F16913", "#BD0026", "#67000D")
  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(.data$destination_T_JJA, .data$destination_P_JJA)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = .data$diff_T_JJA),
      shape = 21, size = 5.6, stroke = 0.7, colour = "white"
    ) +
    ggrepel::geom_text_repel(
      ggplot2::aes(
        label = .data$experiment_label,
        colour = .data$diff_T_JJA
      ),
      size = experiment_label_size,
      max.overlaps = Inf,
      min.segment.length = 0,
      box.padding = 0.6,
      point.padding = 0.35,
      seed = 20260827,
      force = 4,
      force_pull = 0.15,
      max.time = 10,
      max.iter = 200000,
      segment.colour = "grey55",
      segment.size = 0.4
    ) +
    ggplot2::scale_fill_gradientn(
      colours = warming_colours,
      name = "JJA warming (destination - origin, °C)",
      guide = ggplot2::guide_colourbar(
        title.position = "top",
        title.hjust = 0.5,
        direction = "horizontal",
        barwidth = grid::unit(16, "cm"),
        barheight = grid::unit(0.75, "cm")
      )
    ) +
    ggplot2::scale_colour_gradientn(colours = warming_colours, guide = "none") +
    ggplot2::scale_x_continuous(
      breaks = scales::breaks_pretty(n = 7),
      labels = scales::label_number(accuracy = 0.1, big.mark = ""),
      expand = ggplot2::expansion(mult = 0.035)
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::breaks_pretty(n = 7),
      labels = scales::label_number(accuracy = 1, big.mark = ""),
      expand = ggplot2::expansion(mult = 0.04)
    ) +
    transplant_figure_theme(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey88", linewidth = 0.4),
      axis.title = ggplot2::element_text(size = axis_title_size),
      axis.text = ggplot2::element_text(size = axis_text_size),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = ggplot2::element_text(size = legend_title_size),
      legend.text = ggplot2::element_text(size = legend_text_size),
      legend.margin = ggplot2::margin(t = 8),
      plot.margin = ggplot2::margin(12, 18, 8, 18)
    ) +
    ggplot2::labs(
      x = "Destination JJA mean temperature (°C)",
      y = "Destination JJA precipitation (mm)"
    )

  save_plot_atomic(plot, path, width, height, bg = "white")
}

# Additional audited data helpers retained for climate-display tests and provenance.

#' Category Density Data Xy.
#'
#' @param data A data frame whose required columns are checked before use.
#' @param x_column A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param y_column A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param grid_size A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
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
#' Destination Climate Experiment Map Data.
#'
#' @param data A data frame whose required columns are checked before use.
#' @return An R object described by the function name; its structure is validated before later stages use it.
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
