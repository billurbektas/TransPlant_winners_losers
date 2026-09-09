experiment_community <- function(community, experiment) {
  region_value <- experiment$region[[1]]
  origin_value <- experiment$origin_site[[1]]
  destination_value <- experiment$destination_site[[1]]
  dplyr::bind_rows(
    origin_controls = community |>
      dplyr::filter(.data$region == region_value, .data$origin_site == origin_value,
                    .data$destination_site == origin_value, .data$treatment == "LocalControl"),
    destination_controls = community |>
      dplyr::filter(.data$region == region_value, .data$origin_site == destination_value,
                    .data$destination_site == destination_value, .data$treatment == "LocalControl"),
    warmed = community |>
      dplyr::filter(.data$region == region_value, .data$origin_site == origin_value,
                    .data$destination_site == destination_value, .data$treatment == "Warm"),
    .id = "ODT"
  ) |>
    dplyr::mutate(
      ODT = factor(.data$ODT, levels = c("destination_controls", "origin_controls", "warmed")),
      experiment_id = experiment$experiment_id[[1]]
    )
}

apply_legacy_observation_filter <- function(data, minimum_observations = 2L) {
  data |>
    dplyr::group_by(.data$ODT, .data$species) |>
    dplyr::filter(dplyr::n() >= minimum_observations) |>
    dplyr::ungroup()
}

experiment_pool <- function(data) {
  # Match the legacy analysis: pool membership is assigned once from all
  # observations over the complete experiment, after the n > 1 filter.
  present <- data |>
    dplyr::filter(!is.na(.data$species)) |>
    dplyr::distinct(.data$ODT, .data$species)
  origin <- present |>
    dplyr::filter(.data$ODT == "origin_controls") |>
    dplyr::pull(.data$species)
  destination <- present |>
    dplyr::filter(.data$ODT == "destination_controls") |>
    dplyr::pull(.data$species)
  warmed <- present |>
    dplyr::filter(.data$ODT == "warmed") |>
    dplyr::pull(.data$species)
  universe <- union(origin, warmed)
  tibble::tibble(species = universe) |>
    dplyr::mutate(
      in_origin = .data$species %in% origin,
      in_destination = .data$species %in% destination,
      in_warmed = .data$species %in% warmed,
      pool = dplyr::case_when(
        .data$in_origin ~ "resident",
        !.data$in_origin & .data$in_destination & .data$in_warmed ~ "colonizer",
        !.data$in_origin & !.data$in_destination & .data$in_warmed ~ "outsider",
        TRUE ~ "not_analysed"
      )
    )
}

complete_horizon_community <- function(data, pool) {
  species_universe <- pool |>
    dplyr::filter(.data$pool != "not_analysed") |>
    dplyr::pull(.data$species)
  surveyed_plot_years <- data |>
    dplyr::filter(!is.na(.data$relative_year)) |>
    dplyr::distinct(.data$ODT, .data$stable_plot_id, .data$relative_year)
  observed <- data |>
    dplyr::filter(.data$species %in% species_universe, !is.na(.data$relative_year)) |>
    dplyr::select(ODT, stable_plot_id, relative_year, species, relative_cover) |>
    dplyr::group_by(.data$ODT, .data$stable_plot_id, .data$relative_year, .data$species) |>
    dplyr::summarise(relative_cover = sum(.data$relative_cover), .groups = "drop")
  surveyed_plot_years |>
    tidyr::crossing(species = species_universe) |>
    dplyr::left_join(observed, by = c("ODT", "stable_plot_id", "relative_year", "species")) |>
    dplyr::mutate(relative_cover = dplyr::coalesce(.data$relative_cover, 0)) |>
    dplyr::group_by(.data$ODT) |>
    tidyr::complete(
      stable_plot_id, relative_year, species = species_universe
    ) |>
    dplyr::ungroup() |>
    dplyr::group_by(.data$ODT, .data$relative_year, .data$species) |>
    dplyr::mutate(
      relative_cover = dplyr::if_else(
        is.na(.data$relative_cover), mean(.data$relative_cover, na.rm = TRUE), .data$relative_cover
      ),
      relative_cover = dplyr::coalesce(.data$relative_cover, 0)
    ) |>
    dplyr::ungroup()
}

auc_by_plot <- function(completed, horizon, experiment_id) {
  to_value <- if (identical(horizon, "max")) NULL else as.numeric(horizon)
  completed |>
    tidyr::crossing(transform = c("raw", "sqrt")) |>
    dplyr::mutate(value = dplyr::if_else(.data$transform == "sqrt", sqrt(.data$relative_cover), .data$relative_cover)) |>
    dplyr::group_by(.data$ODT, .data$stable_plot_id, .data$species, .data$transform) |>
    dplyr::summarise(
      auc = if (dplyr::n_distinct(.data$relative_year) < 2L) {
        NA_real_
      } else if (is.null(to_value)) {
        DescTools::AUC(.data$relative_year, .data$value, method = "trapezoid")
      } else {
        DescTools::AUC(.data$relative_year, .data$value, method = "trapezoid", to = to_value)
      },
      .groups = "drop"
    ) |>
    dplyr::mutate(experiment_id = experiment_id, horizon = as.character(horizon))
}

calculate_auc_inputs <- function(community, experiments, horizons = 2L:9L) {
  included <- experiments |>
    dplyr::filter(.data$included)
  all_horizons <- c(as.list(as.integer(horizons)), list("max"))
  results <- purrr::map(seq_len(nrow(included)), function(i) {
    experiment <- included[i, ]
    data <- experiment_community(community, experiment) |>
      apply_legacy_observation_filter()
    fixed_pool <- experiment_pool(data)
    completed <- complete_horizon_community(data, fixed_pool)
    purrr::map(all_horizons, function(horizon) {
      eligible <- identical(horizon, "max") ||
        (min(data$relative_year, na.rm = TRUE) < as.integer(horizon) &&
           max(data$relative_year, na.rm = TRUE) >= as.integer(horizon))
      if (!eligible) return(NULL)
      pool <- fixed_pool |>
        dplyr::mutate(experiment_id = experiment$experiment_id[[1]], horizon = as.character(horizon))
      list(
        auc = auc_by_plot(completed, horizon, experiment$experiment_id[[1]]),
        pools = pool
      )
    }) |>
      purrr::compact()
  }) |>
    purrr::flatten()
  list(
    auc = purrr::map_dfr(results, "auc"),
    pools = purrr::map_dfr(results, "pools")
  )
}

stable_contrast_seed <- function(experiment_id, horizon, transform) {
  key <- paste(experiment_id, horizon, transform, sep = "::")
  code_points <- utf8ToInt(key)
  as.integer((sum(as.double(code_points) * (seq_along(code_points) + 17)) %% 2147483646) + 1)
}

fit_one_auc_model <- function(data, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  data <- data |>
    dplyr::filter(!is.na(.data$auc)) |>
    dplyr::mutate(ODT = factor(.data$ODT, levels = c("destination_controls", "origin_controls", "warmed")))
  model <- stats::lm(auc ~ ODT * species, data = data)
  means <- emmeans::emmeans(model, ~ ODT | species)
  contrasts <- emmeans::contrast(
    means,
    method = list(
      between_controls = c(1, -1, 0),
      to_origin = c(0, -1, 1),
      to_destination = c(-1, 0, 1)
    ),
    adjust = "mvt"
  )
  as.data.frame(summary(contrasts, infer = c(TRUE, TRUE))) |>
    tibble::as_tibble()
}

fit_auc_contrasts <- function(auc) {
  auc |>
    dplyr::group_by(.data$experiment_id, .data$horizon, .data$transform) |>
    tidyr::nest() |>
    dplyr::mutate(
      model_seed = purrr::pmap_int(
        list(.data$experiment_id, .data$horizon, .data$transform),
        stable_contrast_seed
      ),
      result = purrr::map2(.data$data, .data$model_seed, ~ tryCatch(
        fit_one_auc_model(.x, seed = .y),
        error = function(e) tibble::tibble(model_error = conditionMessage(e))
      ))
    ) |>
    dplyr::select(-data, -model_seed) |>
    tidyr::unnest(result)
}

classify_winners_losers <- function(contrasts, pools, alpha = 0.05) {
  usable <- contrasts |>
    dplyr::filter(.data$contrast == "to_origin", !is.na(.data$p.value))
  usable |>
    dplyr::left_join(
      pools |>
        dplyr::select(experiment_id, horizon, species, pool),
      by = c("experiment_id", "horizon", "species")
    ) |>
    dplyr::mutate(
      significant = .data$p.value < alpha,
      category = dplyr::case_when(
        .data$pool == "resident" & .data$significant & .data$estimate < 0 ~ "loser",
        .data$pool == "resident" & .data$significant & .data$estimate > 0 ~ "winner",
        .data$pool %in% c("colonizer", "outsider") &
          .data$significant & .data$estimate > 0 ~ "colonizer",
        .data$pool == "resident" & !.data$significant & .data$estimate < 0 ~ "non_significant_loser",
        .data$pool == "resident" & !.data$significant & .data$estimate > 0 ~ "non_significant_winner",
        .data$pool == "resident" & .data$estimate == 0 ~ "no_change",
        .data$pool == "colonizer" & !.data$significant ~ "non_significant_colonizer",
        .data$pool == "outsider" & !.data$significant ~ "non_significant_outsider",
        TRUE ~ NA_character_
      ),
      category = factor(
        .data$category,
        levels = c(
          "loser", "non_significant_loser", "no_change", "non_significant_winner",
          "winner", "colonizer", "non_significant_colonizer", "non_significant_outsider"
        )
      ),
      analysis_category = dplyr::if_else(
        as.character(.data$category) == "non_significant_outsider",
        "non_significant_colonizer",
        as.character(.data$category)
      ),
      analysis_category = factor(
        .data$analysis_category,
        levels = c(
          "loser", "non_significant_loser", "no_change", "non_significant_winner",
          "winner", "colonizer", "non_significant_colonizer"
        )
      )
    )
}

validate_winner_loser_outputs <- function(classified, pools) {
  duplicate_keys <- classified |>
    dplyr::count(.data$experiment_id, .data$horizon, .data$transform, .data$species) |>
    dplyr::filter(.data$n > 1L)
  if (nrow(duplicate_keys)) stop("Winner/loser output has duplicate species-experiment-horizon rows.", call. = FALSE)
  unclassified <- classified |>
    dplyr::filter(is.na(.data$category))
  if (nrow(unclassified)) {
    stop("At least one usable contrast did not satisfy a category definition.", call. = FALSE)
  }
  bad_colonizers <- classified |>
    dplyr::filter(
      .data$category == "colonizer",
      !.data$pool %in% c("colonizer", "outsider") | !.data$significant | .data$estimate <= 0
    )
  if (nrow(bad_colonizers)) stop("Significant colonizer classification failed.", call. = FALSE)
  bad_non_significant <- classified |>
    dplyr::filter(stringr::str_starts(as.character(.data$category), "non_significant_"), .data$significant)
  if (nrow(bad_non_significant)) {
    stop("A significant response was labelled non-significant.", call. = FALSE)
  }
  changing_pools <- pools |>
    dplyr::distinct(.data$experiment_id, .data$species, .data$horizon, .data$pool) |>
    dplyr::summarise(n_pools = dplyr::n_distinct(.data$pool), .by = c("experiment_id", "species")) |>
    dplyr::filter(.data$n_pools != 1L)
  if (nrow(changing_pools)) {
    stop("Pool membership changes across horizons; pools must be fixed per species and experiment.", call. = FALSE)
  }
  resident_categories <- c("loser", "non_significant_loser", "no_change",
                           "non_significant_winner", "winner")
  nonresident_categories <- c("colonizer", "non_significant_colonizer",
                              "non_significant_outsider")
  incompatible_categories <- classified |>
    dplyr::filter(
      (.data$pool == "resident" & !as.character(.data$category) %in% .env$resident_categories) |
        (.data$pool %in% c("colonizer", "outsider") &
           !as.character(.data$category) %in% .env$nonresident_categories)
    )
  if (nrow(incompatible_categories)) {
    stop("A fixed colonizer/outsider pool was classified as a resident response, or vice versa.", call. = FALSE)
  }
  bad_analysis_merge <- classified |>
    dplyr::filter(
      (.data$category == "non_significant_outsider" &
         .data$analysis_category != "non_significant_colonizer") |
        (.data$category != "non_significant_outsider" &
           as.character(.data$analysis_category) != as.character(.data$category))
    )
  if (nrow(bad_analysis_merge)) {
    stop("The reversible outsider-to-colonizer downstream mapping is invalid.", call. = FALSE)
  }
  tibble::tibble(
    check = c(
      "unique result keys", "all usable contrasts classified", "colonizers significant positive",
      "non-significant labels non-significant", "pools fixed across horizons",
      "pool and category definitions compatible", "downstream outsider merge reversible"
    ),
    pass = c(nrow(duplicate_keys) == 0L,
             nrow(unclassified) == 0L,
             nrow(bad_colonizers) == 0L,
             nrow(bad_non_significant) == 0L,
             nrow(changing_pools) == 0L,
             nrow(incompatible_categories) == 0L,
             nrow(bad_analysis_merge) == 0L)
  )
}

compare_legacy_auc <- function(classified, legacy_path, tolerance = 1e-8) {
  if (!file.exists(legacy_path)) {
    return(tibble::tibble(check = "legacy AUC comparison", pass = NA, n = 0L,
                          max_absolute_difference = NA_real_, note = "legacy result file absent"))
  }
  legacy <- readRDS(legacy_path) |>
    dplyr::filter(.data$rel_abund == "auc.trap", as.character(.data$step) == "max", .data$contrast == "to_origin") |>
    dplyr::transmute(
      experiment_id = paste(.data$Region, .data$originSiteID, .data$destSiteID, sep = "::"),
      species = .data$SpeciesName,
      legacy_estimate = .data$estimate
    )
  comparison <- classified |>
    dplyr::filter(.data$transform == "raw", .data$horizon == "max") |>
    dplyr::select(experiment_id, species, revised_estimate = estimate) |>
    dplyr::inner_join(legacy, by = c("experiment_id", "species")) |>
    dplyr::mutate(absolute_difference = abs(.data$revised_estimate - .data$legacy_estimate))
  maximum <- max(comparison$absolute_difference, na.rm = TRUE)
  tibble::tibble(
    check = "legacy AUC comparison",
    pass = is.finite(maximum) && maximum <= tolerance,
    n = nrow(comparison),
    max_absolute_difference = maximum,
    note = paste0("tolerance=", tolerance)
  )
}
