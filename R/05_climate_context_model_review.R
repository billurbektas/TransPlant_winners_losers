# ==============================================================================
# Script: 05_climate_context_model_review.R
#
# Purpose:
#   Review a simplified four-context climate-response model without changing the
#   maintained targets pipeline. The two fitted models are:
#
#   resident_outcome ~ climate_context4 * z_warming
#   colonizer_event  ~ climate_context4 * z_warming
#
#   This script deliberately fits fixed-effects models only. Earlier diagnostics
#   showed that species, region, comparison, and combined random-intercept
#   multinomial fits were not numerically reliable. It is a preserved manual
#   trial outside _targets.R. Its reconstructed figure is intentionally not
#   retained in the clean checkpoint or listed in the figure manifest.
#
# Input:
#   output/tables/analysis_dataset.csv
#   output/tables/climate_context_experiments.csv
#
# Output:
#   plot/climate_context_counts_odds_probabilities.pdf
# ==============================================================================

source("config/analysis.R")
source("R/pipeline/utils.R")
source("R/pipeline/figures.R")
source("R/pipeline/climate_response_model.R")

config <- transplant_config()
analysis_path <- file.path(config$output_dir, "tables", "analysis_dataset.csv")
context_path <- file.path(
  config$output_dir, "tables", "climate_context_experiments.csv"
)
figure_path <- file.path(
  config$plot_dir, "climate_context_counts_odds_probabilities.pdf"
)

if (!file.exists(analysis_path)) {
  stop("Missing review-model input: ", analysis_path, call. = FALSE)
}
if (!file.exists(context_path)) {
  stop("Missing climate-context input: ", context_path, call. = FALSE)
}

context_levels <- c(
  "cold_dry", "intermediate_dry", "warm_dry", "cold_non_dry"
)
context_labels <- c(
  cold_dry = "Cold/dry",
  intermediate_dry = "Intermediate/dry",
  warm_dry = "Warm/dry",
  cold_non_dry = "Cold/non-dry"
)
context_colours <- c(
  cold_dry = "#1B9E77",
  intermediate_dry = "#FDAE61",
  warm_dry = "#D7191C",
  cold_non_dry = "#2C7BB6"
)
display_context_colours <- stats::setNames(
  unname(context_colours[context_levels]), unname(context_labels[context_levels])
)
response_levels <- c("loser", "winner", "colonizer")
response_labels <- c(
  loser = "Loser | resident pool",
  winner = "Winner | resident pool",
  colonizer = "Colonizer | nonresident pool"
)

analysis_data <- readr::read_csv(analysis_path, show_col_types = FALSE)
climate_contexts <- readr::read_csv(context_path, show_col_types = FALSE)
prepared <- prepare_pool_climate_response_data(analysis_data)

context4 <- climate_contexts |>
  dplyr::transmute(
    experiment_id = as.character(.data$experiment_id),
    climate_context4 = dplyr::if_else(
      as.character(.data$climate_context) %in%
        c("cold_intermediate", "cold_wet"),
      "cold_non_dry", as.character(.data$climate_context)
    )
  ) |>
  dplyr::mutate(
    climate_context4 = factor(
      .data$climate_context4, levels = .env$context_levels
    )
  )

if (anyNA(context4$climate_context4)) {
  stop("Every experiment must belong to one reviewed four-context group.", call. = FALSE)
}
assert_unique(context4, "experiment_id", "four-context climate table")

add_context <- function(data) {
  result <- data |>
    dplyr::mutate(experiment_id = as.character(.data$experiment_id)) |>
    dplyr::left_join(context4, by = "experiment_id") |>
    dplyr::mutate(
      climate_context4 = factor(
        as.character(.data$climate_context4), levels = context_levels
      )
    )
  if (anyNA(result$climate_context4)) {
    stop("A pool-model row is missing its four-context assignment.", call. = FALSE)
  }
  result
}

resident_data <- add_context(prepared$resident)
nonresident_data <- add_context(prepared$nonresident)
review_formula <- stats::as.formula(
  "response ~ climate_context4 * z_warming"
)

resident_formula <- stats::update.formula(
  review_formula, resident_outcome ~ .
)
colonizer_formula <- stats::update.formula(
  review_formula, colonizer_event ~ .
)

resident_model <- mclogit::mblogit(
  resident_formula, data = resident_data, trace = FALSE,
  control = mclogit::mclogit.control(maxit = 100, trace = FALSE)
)
colonizer_model <- stats::glm(
  colonizer_formula, data = nonresident_data, family = stats::binomial()
)

if (!isTRUE(resident_model$converged) || !isTRUE(colonizer_model$converged) ||
    any(!is.finite(stats::coef(resident_model))) ||
    any(!is.finite(stats::coef(colonizer_model))) ||
    any(!is.finite(stats::vcov(resident_model))) ||
    any(!is.finite(stats::vcov(colonizer_model)))) {
  stop("A four-context review model did not converge cleanly.", call. = FALSE)
}

# Observed event counts retain the corresponding eligible-pool denominator and
# percentage. The bar height is the raw species x experiment event count.
resident_counts <- resident_data |>
  dplyr::summarise(
    pool_count = dplyr::n(),
    loser = sum(.data$resident_outcome == "loser"),
    winner = sum(.data$resident_outcome == "winner"),
    n_experiments = dplyr::n_distinct(.data$experiment_id),
    .by = "climate_context4"
  ) |>
  tidyr::pivot_longer(
    c("loser", "winner"), names_to = "response", values_to = "event_count"
  )

colonizer_counts <- nonresident_data |>
  dplyr::summarise(
    pool_count = dplyr::n(),
    event_count = sum(.data$colonizer_event),
    n_experiments = dplyr::n_distinct(.data$experiment_id),
    .by = "climate_context4"
  ) |>
  dplyr::mutate(response = "colonizer")

observed_counts <- dplyr::bind_rows(resident_counts, colonizer_counts) |>
  dplyr::mutate(
    proportion = .data$event_count / .data$pool_count,
    count_label = sprintf(
      "%d / %d  (%.1f%%)",
      .data$event_count, .data$pool_count, 100 * .data$proportion
    ),
    climate_context4 = factor(
      as.character(.data$climate_context4), levels = context_levels,
      labels = unname(context_labels[context_levels])
    ),
    response = factor(
      .data$response, levels = response_levels,
      labels = unname(response_labels[response_levels])
    )
  )

if (sum(dplyr::filter(observed_counts, response == response_labels[["loser"]])$event_count) != 132L ||
    sum(dplyr::filter(observed_counts, response == response_labels[["winner"]])$event_count) != 69L ||
    sum(dplyr::filter(observed_counts, response == response_labels[["colonizer"]])$event_count) != 34L) {
  stop("Four-context observed counts do not conserve the pool-model events.", call. = FALSE)
}

contrast_warming_or <- function(model, outcomes) {
  coefficients <- stats::coef(model)
  covariance <- stats::vcov(model)
  purrr::map_dfr(outcomes, function(outcome) {
    purrr::map_dfr(context_levels, function(context) {
      contrast <- stats::setNames(rep(0, length(coefficients)), names(coefficients))
      prefix <- if (length(outcomes) > 1L) paste0(outcome, "~") else ""
      contrast[[paste0(prefix, "z_warming")]] <- 1
      if (context != "cold_dry") {
        interaction_name <- paste0(
          prefix, "climate_context4", context, ":z_warming"
        )
        contrast[[interaction_name]] <- 1
      }
      estimate <- sum(contrast * coefficients)
      standard_error <- sqrt(drop(t(contrast) %*% covariance %*% contrast))
      tibble::tibble(
        response = outcome,
        climate_context4 = context,
        estimate = estimate,
        standard_error = standard_error,
        odds_ratio = exp(estimate),
        odds_low = exp(estimate - stats::qnorm(0.975) * standard_error),
        odds_high = exp(estimate + stats::qnorm(0.975) * standard_error)
      )
    })
  })
}

contrast_context_or <- function(model, outcomes) {
  coefficients <- stats::coef(model)
  covariance <- stats::vcov(model)
  purrr::map_dfr(outcomes, function(outcome) {
    purrr::map_dfr(context_levels, function(context) {
      prefix <- if (length(outcomes) > 1L) paste0(outcome, "~") else ""
      if (context == "cold_dry") {
        return(tibble::tibble(
          response = outcome, climate_context4 = context,
          estimate = 0, standard_error = NA_real_,
          odds_ratio = 1, odds_low = 1, odds_high = 1,
          is_reference = TRUE
        ))
      }
      coefficient_name <- paste0(prefix, "climate_context4", context)
      estimate <- unname(coefficients[[coefficient_name]])
      standard_error <- sqrt(covariance[coefficient_name, coefficient_name])
      tibble::tibble(
        response = outcome,
        climate_context4 = context,
        estimate = estimate,
        standard_error = standard_error,
        odds_ratio = exp(estimate),
        odds_low = exp(estimate - stats::qnorm(0.975) * standard_error),
        odds_high = exp(estimate + stats::qnorm(0.975) * standard_error),
        is_reference = FALSE
      )
    })
  })
}

context_odds <- dplyr::bind_rows(
  contrast_context_or(resident_model, c("loser", "winner")),
  contrast_context_or(colonizer_model, "colonizer")
) |>
  dplyr::mutate(
    climate_context4 = factor(
      .data$climate_context4, levels = context_levels,
      labels = unname(context_labels[context_levels])
    ),
    response = factor(
      .data$response, levels = response_levels,
      labels = unname(response_labels[response_levels])
    )
  )

warming_odds <- dplyr::bind_rows(
  contrast_warming_or(resident_model, c("loser", "winner")),
  contrast_warming_or(colonizer_model, "colonizer")
) |>
  dplyr::mutate(
    climate_context4 = factor(
      .data$climate_context4, levels = context_levels,
      labels = unname(context_labels[context_levels])
    ),
    response = factor(
      .data$response, levels = response_levels,
      labels = unname(response_labels[response_levels])
    )
  )

warming_values <- stats::quantile(
  prepared$experiments$warming, probs = c(0.25, 0.5, 0.75), names = FALSE
)
warming_scale <- prepared$scaling |>
  dplyr::filter(.data$variable == "warming")
warming_reference <- tibble::tibble(
  warming_level = factor(
    c("Low", "Median", "High"), levels = c("Low", "Median", "High")
  ),
  warming = warming_values
)
warming_grid <- tidyr::crossing(
  climate_context4 = factor(context_levels, levels = context_levels),
  warming_reference
) |>
  dplyr::mutate(
    z_warming = (.data$warming - warming_scale$center[[1]]) /
      warming_scale$scale[[1]]
  )

simulate_probabilities <- function(
    resident_model, colonizer_model, new_data, n_draws = 2000L,
    seed = 47163L) {
  set.seed(seed)
  resident_draws <- MASS::mvrnorm(
    n_draws, mu = stats::coef(resident_model),
    Sigma = stats::vcov(resident_model)
  )
  colonizer_draws <- MASS::mvrnorm(
    n_draws, mu = stats::coef(colonizer_model),
    Sigma = stats::vcov(colonizer_model)
  )
  colnames(resident_draws) <- names(stats::coef(resident_model))
  colnames(colonizer_draws) <- names(stats::coef(colonizer_model))
  design <- stats::model.matrix(
    ~ climate_context4 * z_warming, data = new_data
  )
  resident_point <- stats::predict(
    resident_model, newdata = new_data, type = "response"
  )
  colonizer_point <- stats::predict(
    colonizer_model, newdata = new_data, type = "response"
  )
  purrr::map_dfr(seq_len(nrow(new_data)), function(index) {
    x <- design[index, ]
    loser_eta <- resident_draws[, paste0("loser~", names(x)), drop = FALSE] %*% x
    winner_eta <- resident_draws[, paste0("winner~", names(x)), drop = FALSE] %*% x
    denominator <- 1 + exp(loser_eta) + exp(winner_eta)
    loser_probability <- drop(exp(loser_eta) / denominator)
    winner_probability <- drop(exp(winner_eta) / denominator)
    colonizer_eta <- colonizer_draws[, names(x), drop = FALSE] %*% x
    colonizer_probability <- stats::plogis(drop(colonizer_eta))
    probabilities <- list(
      loser = loser_probability,
      winner = winner_probability,
      colonizer = colonizer_probability
    )
    points <- c(
      loser = resident_point[index, "loser"],
      winner = resident_point[index, "winner"],
      colonizer = colonizer_point[[index]]
    )
    purrr::imap_dfr(probabilities, function(values, response) {
      interval <- stats::quantile(
        values, probs = c(0.025, 0.975), names = FALSE, na.rm = TRUE
      )
      tibble::tibble(
        grid_id = index, response = response,
        estimate = unname(points[[response]]),
        conf_low = interval[[1]], conf_high = interval[[2]]
      )
    })
  })
}

probabilities <- simulate_probabilities(
  resident_model, colonizer_model, warming_grid
) |>
  dplyr::left_join(
    dplyr::mutate(warming_grid, grid_id = dplyr::row_number()), by = "grid_id"
  ) |>
  dplyr::mutate(
    climate_context4 = factor(
      as.character(.data$climate_context4), levels = context_levels,
      labels = unname(context_labels[context_levels])
    ),
    response = factor(
      .data$response, levels = response_levels,
      labels = unname(response_labels[response_levels])
    ),
    warming_label = factor(
      sprintf("%s\n%.2f°C", .data$warming_level, .data$warming),
      levels = sprintf(
        "%s\n%.2f°C", c("Low", "Median", "High"), warming_values
      )
    )
  )

counts_plot <- ggplot2::ggplot(
  observed_counts,
  ggplot2::aes(
    .data$event_count, .data$climate_context4,
    fill = .data$climate_context4
  )
) +
  ggplot2::geom_col(width = 0.68, alpha = 0.88) +
  ggplot2::geom_text(
    ggplot2::aes(label = .data$count_label),
    hjust = -0.08, size = 3.6, fontface = "bold"
  ) +
  ggplot2::facet_grid(rows = ggplot2::vars(.data$response), scales = "free_x") +
  ggplot2::scale_fill_manual(values = display_context_colours, guide = "none") +
  ggplot2::scale_x_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.55)), breaks = scales::breaks_pretty()
  ) +
  transplant_figure_theme(base_size = 14) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    strip.text.y = ggplot2::element_text(face = "bold", size = 13),
    axis.text.y = ggplot2::element_text(size = 11),
    panel.spacing.y = grid::unit(0.55, "cm")
  ) +
  ggplot2::labs(
    x = "Observed count\nevent / eligible pool  (percentage)", y = NULL
  )

context_odds_plot <- ggplot2::ggplot(
  context_odds,
  ggplot2::aes(
    .data$odds_ratio, .data$climate_context4,
    colour = .data$climate_context4
  )
) +
  ggplot2::geom_vline(xintercept = 1, colour = "grey55", linewidth = 0.55) +
  ggplot2::geom_errorbar(
    data = dplyr::filter(context_odds, !.data$is_reference),
    ggplot2::aes(xmin = .data$odds_low, xmax = .data$odds_high),
    orientation = "y", width = 0, linewidth = 0.85
  ) +
  ggplot2::geom_point(
    ggplot2::aes(shape = .data$is_reference), size = 3.4, stroke = 1.0
  ) +
  ggplot2::facet_grid(rows = ggplot2::vars(.data$response)) +
  ggplot2::scale_x_log10() +
  ggplot2::scale_colour_manual(values = display_context_colours, guide = "none") +
  ggplot2::scale_shape_manual(
    values = c(`FALSE` = 16, `TRUE` = 4), guide = "none"
  ) +
  transplant_figure_theme(base_size = 14) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    strip.text.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    panel.spacing.y = grid::unit(0.55, "cm")
  ) +
  ggplot2::labs(
    x = sprintf(
      "Destination-context odds ratio\nat mean warming (%.2f°C)",
      warming_scale$center[[1]]
    ),
    y = NULL
  )

odds_plot <- ggplot2::ggplot(
  warming_odds,
  ggplot2::aes(
    .data$odds_ratio, .data$climate_context4,
    colour = .data$climate_context4
  )
) +
  ggplot2::geom_vline(xintercept = 1, colour = "grey55", linewidth = 0.55) +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = .data$odds_low, xmax = .data$odds_high),
    orientation = "y", width = 0, linewidth = 0.85
  ) +
  ggplot2::geom_point(size = 3.4) +
  ggplot2::facet_grid(rows = ggplot2::vars(.data$response)) +
  ggplot2::scale_x_log10() +
  ggplot2::scale_colour_manual(values = display_context_colours, guide = "none") +
  transplant_figure_theme(base_size = 14) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    strip.text.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    panel.spacing.y = grid::unit(0.55, "cm")
  ) +
  ggplot2::labs(
    x = sprintf(
      "Warming odds ratio\nper %.2f°C", warming_scale$scale[[1]]
    ),
    y = NULL
  )

probability_limit <- min(1, max(probabilities$conf_high, na.rm = TRUE) * 1.08)
probability_plot <- ggplot2::ggplot(
  probabilities,
  ggplot2::aes(
    .data$warming_label, .data$estimate,
    group = .data$climate_context4, colour = .data$climate_context4
  )
) +
  ggplot2::geom_line(linewidth = 0.9, position = ggplot2::position_dodge(width = 0.10)) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = .data$conf_low, ymax = .data$conf_high),
    width = 0.08, linewidth = 0.65,
    position = ggplot2::position_dodge(width = 0.10)
  ) +
  ggplot2::geom_point(
    size = 2.8, position = ggplot2::position_dodge(width = 0.10)
  ) +
  ggplot2::facet_grid(rows = ggplot2::vars(.data$response)) +
  ggplot2::scale_colour_manual(
    values = display_context_colours, name = "Destination climate"
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, probability_limit), labels = scales::label_percent(),
    expand = ggplot2::expansion(mult = c(0, 0.04))
  ) +
  transplant_figure_theme(base_size = 14) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.text.y = ggplot2::element_blank(),
    panel.spacing.y = grid::unit(0.55, "cm"),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(hjust = 0.5)
  ) +
  ggplot2::labs(x = "Experimental JJA warming", y = "Predicted probability within pool")

combined <- patchwork::wrap_plots(
  counts_plot, context_odds_plot, odds_plot, probability_plot,
  nrow = 1, widths = c(1.35, 0.72, 0.72, 1.55), guides = "collect"
) &
  ggplot2::theme(legend.position = "bottom")

save_plot_atomic(combined, figure_path, width = 29, height = 12)

cat("\nFour-context review models completed. No model or table was saved.\n")
cat("Figure: ", figure_path, "\n", sep = "")
cat("\nObserved counts:\n")
print(observed_counts)
cat("\nDestination-context odds ratios at mean warming:\n")
print(context_odds)
cat("\nContext-specific warming odds ratios:\n")
print(warming_odds)
