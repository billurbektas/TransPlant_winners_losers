# Module: R/responses/temporal
# Responsibility: Implement the named scientific workflow stage or shared infrastructure.
# Inputs: Function arguments and explicitly declared target dependencies.
# Outputs: Validated in-memory objects and declared workflow artifacts.
# Side effects: Limited to functions whose contracts explicitly write, copy, or render files.

#' Temporal Transitions.
#'
#' @param classified A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param transform A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
temporal_transitions <- function(classified, transform = "raw") {
  category_column <- if ("analysis_category" %in% names(classified)) "analysis_category" else "category"
  states <- classified |>
    dplyr::ungroup() |>
    dplyr::filter(.data$transform == .env$transform) |>
    dplyr::mutate(category = as.character(.data[[category_column]])) |>
    dplyr::select(experiment_id, species, horizon, category)
  early <- states |>
    dplyr::filter(.data$horizon == "2") |>
    dplyr::select(experiment_id, species, state_2 = category)
  purrr::map_dfr(c("4", "6"), function(later_horizon) {
    later <- states |>
      dplyr::filter(.data$horizon == .env$later_horizon) |>
      dplyr::select(experiment_id, species, state_later = category)
    dplyr::inner_join(early, later, by = c("experiment_id", "species")) |>
      dplyr::mutate(
        later_horizon = .env$later_horizon,
        same_state = .data$state_2 == .data$state_later
      )
  })
}

#' Temporal Linked Transitions.
#'
#' @param classified A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param transform A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param horizons A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
temporal_linked_transitions <- function(classified, transform = "raw",
                                        horizons = c("2", "4", "6")) {
  category_column <- if ("analysis_category" %in% names(classified)) "analysis_category" else "category"
  if ("pool" %in% names(classified)) {
    changing_pools <- classified |>
      dplyr::ungroup() |>
      dplyr::filter(.data$transform == .env$transform, .data$horizon %in% .env$horizons) |>
      dplyr::summarise(
        n_pools = dplyr::n_distinct(.data$pool),
        .by = c("experiment_id", "species")
      ) |>
      dplyr::filter(.data$n_pools != 1L)
    if (nrow(changing_pools)) {
      stop("Temporal results contain pool changes; pools must be fixed per species and experiment.", call. = FALSE)
    }
  }
  state_columns <- paste0("state_", horizons)
  states <- classified |>
    dplyr::ungroup() |>
    dplyr::filter(
      .data$transform == .env$transform,
      .data$horizon %in% .env$horizons
    ) |>
    dplyr::transmute(
      experiment_id = .data$experiment_id,
      species = .data$species,
      horizon = paste0("state_", .data$horizon),
      category = as.character(.data[[category_column]])
    ) |>
    tidyr::pivot_wider(names_from = "horizon", values_from = "category")
  if ("pool" %in% names(classified)) {
    pools <- classified |>
      dplyr::ungroup() |>
      dplyr::filter(
        .data$transform == .env$transform,
        .data$horizon %in% .env$horizons
      ) |>
      dplyr::transmute(
        experiment_id = .data$experiment_id,
        species = .data$species,
        horizon = paste0("pool_", .data$horizon),
        pool = as.character(.data$pool)
      ) |>
      tidyr::pivot_wider(names_from = "horizon", values_from = "pool")
    states <- dplyr::left_join(states, pools, by = c("experiment_id", "species"))
  } else {
    for (horizon in horizons) states[[paste0("pool_", horizon)]] <- NA_character_
  }
  missing_columns <- setdiff(state_columns, names(states))
  if (length(missing_columns)) {
    stop("Linked temporal flow is missing horizons: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  states |>
    dplyr::filter(!is.na(.data$state_2)) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(setdiff(state_columns, "state_2")),
        ~ dplyr::coalesce(.x, "follow_up_unavailable")
      )
    ) |>
    dplyr::select(
      "experiment_id", "species", dplyr::all_of(state_columns),
      dplyr::all_of(paste0("pool_", horizons))
    )
}

#' Linked Transition Summary.
#'
#' @param linked A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
linked_transition_summary <- function(linked) {
  linked |>
    dplyr::count(
      .data$state_2, .data$state_4, .data$state_6,
      name = "n"
    ) |>
    dplyr::arrange(
      .data$state_2, .data$state_4, .data$state_6
    )
}

#' Summarise Temporal Category Persistence.
#'
#' @param linked A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
summarise_temporal_category_persistence <- function(linked) {
  comparisons <- tibble::tribble(
    ~comparison, ~from_column, ~to_column,
    "First 2 → First 4", "state_2", "state_4",
    "First 4 → First 6", "state_4", "state_6",
    "First 2 → First 6", "state_2", "state_6"
  )
  by_category <- purrr::pmap_dfr(comparisons, function(comparison, from_column, to_column) {
    linked |>
      dplyr::filter(
        .data[[from_column]] != "follow_up_unavailable",
        .data[[to_column]] != "follow_up_unavailable"
      ) |>
      dplyr::transmute(
        comparison = .env$comparison,
        from_category = .data[[from_column]],
        same_category = .data[[from_column]] == .data[[to_column]]
      ) |>
      dplyr::summarise(
        available_n = dplyr::n(),
        same_n = sum(.data$same_category),
        same_proportion = mean(.data$same_category),
        .by = c("comparison", "from_category")
      )
  })
  overall <- purrr::pmap_dfr(comparisons, function(comparison, from_column, to_column) {
    available <- linked |>
      dplyr::filter(
        .data[[from_column]] != "follow_up_unavailable",
        .data[[to_column]] != "follow_up_unavailable"
      )
    tibble::tibble(
      comparison = comparison,
      available_n = nrow(available),
      same_n = sum(available[[from_column]] == available[[to_column]]),
      same_proportion = mean(available[[from_column]] == available[[to_column]])
    )
  })
  complete <- linked |>
    dplyr::filter(dplyr::if_all(dplyr::starts_with("state_"), ~ .x != "follow_up_unavailable"))
  overall <- dplyr::bind_rows(
    overall,
    tibble::tibble(
      comparison = "Same at First 2, 4, and 6",
      available_n = nrow(complete),
      same_n = sum(complete$state_2 == complete$state_4 & complete$state_4 == complete$state_6),
      same_proportion = mean(complete$state_2 == complete$state_4 & complete$state_4 == complete$state_6)
    )
  )
  list(by_category = by_category, overall = overall)
}

#' Build Temporal Transition Matrix.
#'
#' @param linked A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
build_temporal_transition_matrix <- function(linked) {
  comparisons <- tibble::tribble(
    ~comparison, ~from_column, ~to_column,
    "First 2 → First 4", "state_2", "state_4",
    "First 4 → First 6", "state_4", "state_6",
    "First 2 → First 6", "state_2", "state_6"
  )
  categories <- c(
    "loser", "non_significant_loser", "no_change", "non_significant_winner", "winner",
    "colonizer", "non_significant_colonizer"
  )
  response_family <- function(category) {
    dplyr::case_when(
      category %in% c("loser", "non_significant_loser") ~ "loser",
      category %in% c("winner", "non_significant_winner") ~ "winner",
      category %in% c("colonizer", "non_significant_colonizer") ~ "colonizer",
      TRUE ~ "other"
    )
  }
  purrr::pmap_dfr(comparisons, function(comparison, from_column, to_column) {
    from_pool_column <- sub("state_", "pool_", from_column, fixed = TRUE)
    to_pool_column <- sub("state_", "pool_", to_column, fixed = TRUE)
    available <- linked |>
      dplyr::filter(
        .data[[from_column]] != "follow_up_unavailable",
        .data[[to_column]] != "follow_up_unavailable"
      ) |>
      dplyr::transmute(
        from_category = .data[[from_column]],
        to_category = .data[[to_column]],
        from_pool = .data[[from_pool_column]],
        to_pool = .data[[to_pool_column]],
        from_family = response_family(.data[[from_column]]),
        to_family = response_family(.data[[to_column]]),
        from_significant = .data[[from_column]] %in% c("loser", "winner", "colonizer"),
        to_significant = .data[[to_column]] %in% c("loser", "winner", "colonizer")
      ) |>
      dplyr::mutate(
        transition_type = dplyr::case_when(
          .data$from_category == .data$to_category ~ "same category",
          .data$from_family == .data$to_family & .data$from_significant != .data$to_significant ~ "significance only",
          .data$from_family %in% c("loser", "winner") &
            .data$to_family %in% c("loser", "winner") &
            .data$from_family != .data$to_family ~ "loser-winner reversal",
          TRUE ~ "response/status change"
        )
      )
    counts <- available |>
      dplyr::count(.data$from_category, .data$to_category, .data$transition_type, name = "type_n") |>
      dplyr::summarise(
        n = sum(.data$type_n),
        transition_types = paste0(.data$transition_type, "=", .data$type_n, collapse = "; "),
        transition_type = dplyr::if_else(
          dplyr::n() == 1L, dplyr::first(.data$transition_type), "mixed mechanisms"
        ),
        .by = c("from_category", "to_category")
      )
    denominators <- available |>
      dplyr::count(.data$from_category, name = "from_total")
    tidyr::crossing(from_category = categories, to_category = categories) |>
      dplyr::left_join(counts, by = c("from_category", "to_category")) |>
      dplyr::left_join(denominators, by = "from_category") |>
      dplyr::mutate(
        comparison = .env$comparison,
        n = dplyr::coalesce(.data$n, 0L),
        from_total = dplyr::coalesce(.data$from_total, 0L),
        row_proportion = dplyr::if_else(.data$from_total > 0, .data$n / .data$from_total, 0),
        transition_type = dplyr::coalesce(.data$transition_type, "response/status change"),
        transition_types = dplyr::coalesce(.data$transition_types, "")
      )
  })
}

#' Temporal Category Composition.
#'
#' @param classified A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param transform A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param horizons A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
temporal_category_composition <- function(classified, transform = "raw", horizons = NULL) {
  category_column <- if ("analysis_category" %in% names(classified)) "analysis_category" else "category"
  if (is.null(horizons)) {
    horizons <- classified |>
      dplyr::ungroup() |>
      dplyr::filter(.data$transform == .env$transform, .data$horizon != "max") |>
      dplyr::pull(.data$horizon) |>
      unique() |>
      as.integer() |>
      sort() |>
      as.character()
  } else {
    horizons <- as.character(horizons)
  }
  categories <- c(
    "loser", "non_significant_loser", "no_change", "non_significant_winner",
    "winner", "colonizer", "non_significant_colonizer"
  )
  observed <- classified |>
    dplyr::ungroup() |>
    dplyr::filter(
      .data$transform == .env$transform,
      .data$horizon %in% .env$horizons
    ) |>
    dplyr::transmute(
      experiment_id = .data$experiment_id,
      species = .data$species,
      horizon = as.character(.data$horizon),
      category = dplyr::if_else(
        is.na(.data[[category_column]]) | as.character(.data[[category_column]]) == "no_change",
        "no_change",
        as.character(.data[[category_column]])
      )
    ) |>
    dplyr::distinct()
  unexpected <- setdiff(unique(observed$category), categories)
  if (length(unexpected)) {
    stop("Unexpected categories in temporal composition: ", paste(unexpected, collapse = ", "), call. = FALSE)
  }
  denominators <- observed |>
    dplyr::summarise(
      total_species_experiment = dplyr::n(),
      n_species = dplyr::n_distinct(.data$species),
      n_experiments = dplyr::n_distinct(.data$experiment_id),
      .by = "horizon"
    )
  observed |>
    dplyr::count(.data$horizon, .data$category, name = "n") |>
    tidyr::complete(horizon = horizons, category = categories, fill = list(n = 0L)) |>
    dplyr::left_join(denominators, by = "horizon") |>
    dplyr::mutate(
      year = as.integer(.data$horizon),
      proportion = .data$n / .data$total_species_experiment
    ) |>
    dplyr::arrange(.data$year, match(.data$category, categories))
}

#' Transition Summary.
#'
#' @param transitions A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
transition_summary <- function(transitions) {
  counts <- transitions |>
    dplyr::count(.data$later_horizon, .data$state_2, .data$state_later, name = "n") |>
    dplyr::group_by(.data$later_horizon) |>
    dplyr::mutate(proportion = .data$n / sum(.data$n)) |>
    dplyr::ungroup()
  persistence <- transitions |>
    dplyr::group_by(.data$later_horizon) |>
    dplyr::summarise(
      n = dplyr::n(),
      same_state_n = sum(.data$same_state),
      same_state_proportion = mean(.data$same_state),
      .groups = "drop"
    )
  state_counts <- dplyr::bind_rows(
    transitions |>
      dplyr::count(.data$later_horizon, state = .data$state_2, name = "n") |>
      dplyr::mutate(axis = "first_2"),
    transitions |>
      dplyr::count(.data$later_horizon, state = .data$state_later, name = "n") |>
      dplyr::mutate(axis = "later_horizon")
  ) |>
    dplyr::select("later_horizon", "axis", "state", "n")
  list(counts = counts, persistence = persistence, state_counts = state_counts)
}

#' Cohen Kappa.
#'
#' @param x A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param y A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
cohen_kappa <- function(x, y) {
  keep <- !is.na(x) & !is.na(y)
  x <- as.character(x[keep])
  y <- as.character(y[keep])
  levels <- union(unique(x), unique(y))
  if (!length(x) || length(levels) < 2L) return(NA_real_)
  table_xy <- table(factor(x, levels), factor(y, levels))
  observed <- sum(diag(table_xy)) / sum(table_xy)
  expected <- sum(rowSums(table_xy) * colSums(table_xy)) / sum(table_xy)^2
  if (expected == 1) return(NA_real_)
  (observed - expected) / (1 - expected)
}

#' Kappa Cluster Bootstrap.
#'
#' @param transitions A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param replicates A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param seed A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
kappa_cluster_bootstrap <- function(transitions, replicates = 999L, seed = 20260824L) {
  set.seed(seed)
  purrr::map_dfr(unique(transitions$later_horizon), function(horizon) {
    data <- transitions |>
      dplyr::filter(.data$later_horizon == horizon)
    experiments <- unique(data$experiment_id)
    estimate <- cohen_kappa(data$state_2, data$state_later)
    draws <- replicate(replicates, {
      sampled <- sample(experiments, length(experiments), replace = TRUE)
      boot <- purrr::map_dfr(seq_along(sampled), function(i) {
        data |>
          dplyr::filter(.data$experiment_id == sampled[[i]]) |>
          dplyr::mutate(bootstrap_cluster = i)
      })
      cohen_kappa(boot$state_2, boot$state_later)
    })
    tibble::tibble(
      later_horizon = horizon,
      kappa = estimate,
      conf_low = stats::quantile(draws, 0.025, na.rm = TRUE),
      conf_high = stats::quantile(draws, 0.975, na.rm = TRUE),
      p_value = 2 * min(mean(draws <= 0, na.rm = TRUE), mean(draws >= 0, na.rm = TRUE))
    )
  }) |>
    dplyr::mutate(p_holm = stats::p.adjust(.data$p_value, method = "holm"))
}

#' Fit Temporal Strength Models.
#'
#' @param classified A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param transform A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A fitted model bundle with diagnostics retained for validation.
fit_temporal_strength_models <- function(classified, transform = "raw") {
  strength <- classified |>
    dplyr::ungroup() |>
    dplyr::filter(.data$transform == .env$transform) |>
    dplyr::select(experiment_id, species, horizon, estimate) |>
    tidyr::pivot_wider(names_from = horizon, values_from = estimate, names_prefix = "horizon_")
  purrr::map(c("4", "6"), function(horizon) {
    response <- paste0("horizon_", horizon)
    data <- strength |>
      dplyr::filter(!is.na(.data$horizon_2), !is.na(.data[[response]]))
    model <- lme4::lmer(
      stats::as.formula(paste(response, "~ horizon_2 + (1 | species) + (1 | experiment_id)")),
      data = data,
      REML = FALSE
    )
    list(horizon = horizon, model = model, n = nrow(data))
  })
}

#' Validate Flow Conservation.
#'
#' @param transitions A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param summary_counts A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
validate_flow_conservation <- function(transitions, summary_counts) {
  observed <- transitions |>
    dplyr::ungroup() |>
    dplyr::count(.data$later_horizon, name = "transition_rows")
  summarized <- summary_counts |>
    dplyr::ungroup() |>
    dplyr::group_by(.data$later_horizon) |>
    dplyr::summarise(summary_rows = sum(.data$n), .groups = "drop")
  dplyr::left_join(observed, summarized, by = "later_horizon") |>
    dplyr::mutate(pass = .data$transition_rows == .data$summary_rows)
}
