pool_model_fixture <- function() {
  data <- readr::read_csv(
    file.path(PROJECT_ROOT, "output", "tables", "analysis_dataset.csv"),
    show_col_types = FALSE
  )
  prepared <- prepare_pool_climate_response_data(data)
  list(data = data, prepared = prepared)
}

testthat::test_that("pool preparation reproduces reviewed denominators and events", {
  fixture <- pool_model_fixture()
  prepared <- fixture$prepared
  testthat::expect_equal(nrow(prepared$data), 2215L)
  testthat::expect_equal(nrow(prepared$resident), 1604L)
  testthat::expect_equal(nrow(prepared$nonresident), 611L)
  testthat::expect_equal(prepared$snapshot$n_loser, 130L)
  testthat::expect_equal(prepared$snapshot$n_winner, 69L)
  testthat::expect_equal(prepared$snapshot$n_colonizer, 34L)
  testthat::expect_equal(prepared$snapshot$n_significant_outsider, 3L)
  testthat::expect_equal(prepared$snapshot$species_in_multiple_pools, 230L)
  testthat::expect_equal(dplyr::n_distinct(prepared$experiments$experiment_id), 42L)
  testthat::expect_equal(
    prepared$category_counts$n,
    c(130L, 1405L, 69L, 34L, 577L)
  )
  duplicate_keys <- prepared$data |>
    dplyr::count(.data$accepted_species, .data$experiment_id) |>
    dplyr::filter(.data$n != 1L)
  testthat::expect_equal(nrow(duplicate_keys), 0L)
  z_columns <- paste0("z_", prepared$scaling$variable)
  testthat::expect_true(all(vapply(
    prepared$experiments[z_columns],
    function(values) abs(mean(values)) < 1e-12,
    logical(1)
  )))
  testthat::expect_true(all(prepared$scaling$weighting == "42 distinct experiments"))
})

testthat::test_that("pool membership may differ among experiments but remains unique within one", {
  prepared <- pool_model_fixture()$prepared
  changing_species <- prepared$pool_by_species |>
    dplyr::filter(.data$n_pools > 1L)
  testthat::expect_equal(nrow(changing_species), 230L)
  testthat::expect_true(any(stringr::str_detect(
    changing_species$pools, "resident"
  )))
  testthat::expect_false(anyDuplicated(
    paste(prepared$data$accepted_species, prepared$data$experiment_id)
  ) > 0L)
})

testthat::test_that("primary hybrid models use exact fixed structure and converge", {
  prepared <- pool_model_fixture()$prepared
  fit <- fit_pool_climate_response_models(prepared)
  coefficients <- tidy_pool_climate_coefficients(fit)
  expected_terms <- c(
    "(Intercept)", "z_destination_T", "z_log_destination_P", "z_warming",
    "z_destination_T:z_log_destination_P", "z_destination_T:z_warming",
    "z_log_destination_P:z_warming"
  )
  testthat::expect_true(isTRUE(fit$resident$converged))
  testthat::expect_true(isTRUE(fit$colonizer$converged))
  testthat::expect_setequal(unique(coefficients$term), expected_terms)
  testthat::expect_equal(nrow(coefficients), 21L)
  model_counts <- table(coefficients$model)
  testthat::expect_equal(
    as.integer(model_counts[c("nonresident_binomial", "resident_multinomial")]),
    c(7L, 14L)
  )
  testthat::expect_true(all(fit$diagnostics$random_effects == "none"))
  formula_text <- paste(
    deparse(pool_climate_formula("response")), collapse = " "
  )
  testthat::expect_false(stringr::str_detect(formula_text, "diff_P|three"))
  testthat::expect_equal(
    sort(unique(coefficients$reference)),
    sort(c("non_significant_resident", "non_significant_nonresident"))
  )
  non_intercept <- dplyr::filter(coefficients, .data$term != "(Intercept)")
  non_intercept_counts <- table(non_intercept$model)
  testthat::expect_equal(
    as.integer(non_intercept_counts[c(
      "nonresident_binomial", "resident_multinomial"
    )]),
    c(6L, 12L)
  )
  testthat::expect_true(all(non_intercept$p_holm >= non_intercept$p_value))
})

testthat::test_that("pool likelihood-ratio tests have prespecified degrees of freedom", {
  prepared <- pool_model_fixture()$prepared
  fit <- fit_pool_climate_response_models(prepared)
  result <- pool_climate_likelihood_ratio_tests(fit)
  testthat::expect_equal(nrow(result), 4L)
  testthat::expect_equal(result$degrees_freedom, c(12L, 6L, 6L, 3L))
  testthat::expect_equal(
    result$comparison,
    rep(c("full_vs_intercept", "full_pairwise_vs_additive"), 2L)
  )
  testthat::expect_true(all(result$statistic > 0))
  testthat::expect_true(all(result$p_value >= 0 & result$p_value <= 1))
})

testthat::test_that("pool-model performance outputs are finite and calibrated by outcome", {
  prepared <- pool_model_fixture()$prepared
  fit <- fit_pool_climate_response_models(prepared)
  performance <- pool_climate_model_performance(fit)
  testthat::expect_equal(nrow(performance$metrics), 2L)
  testthat::expect_true(all(is.finite(
    unlist(performance$metrics[, setdiff(
      names(performance$metrics), c("model", "limitation")
    )])
  )))
  testthat::expect_equal(nrow(performance$calibration), 40L)
  testthat::expect_setequal(
    unique(performance$calibration$outcome),
    c("non_significant", "loser", "winner", "colonizer")
  )
  testthat::expect_true(all(
    performance$calibration$observed_fraction >= 0 &
      performance$calibration$observed_fraction <= 1
  ))
})

testthat::test_that("binary sensitivity uses species random intercepts only", {
  prepared <- pool_model_fixture()$prepared
  sensitivity <- fit_pool_binary_sensitivity(prepared)
  coefficients <- tidy_pool_binary_sensitivity(sensitivity)
  testthat::expect_equal(sensitivity$diagnostics$response, c(
    "loser", "winner", "colonizer"
  ))
  testthat::expect_false(any(
    sensitivity$diagnostics$experiment_random_intercept
  ))
  testthat::expect_false(any(stringr::str_detect(
    sensitivity$diagnostics$mixed_formula, "experiment_id"
  )))
  testthat::expect_equal(
    sensitivity$diagnostics$selected_model,
    c("species_random_binomial", "species_random_binomial", "fixed_binomial")
  )
  colonizer <- sensitivity$diagnostics |>
    dplyr::filter(.data$response == "colonizer")
  testthat::expect_true(colonizer$mixed_singular)
  testthat::expect_lt(colonizer$species_variance, 1e-6)
  testthat::expect_equal(nrow(coefficients), 21L)
  coefficient_counts <- table(coefficients$response)
  testthat::expect_equal(
    as.integer(coefficient_counts[c("colonizer", "loser", "winner")]),
    c(7L, 7L, 7L)
  )
  non_intercept <- dplyr::filter(coefficients, .data$term != "(Intercept)")
  testthat::expect_true(all(non_intercept$p_holm >= non_intercept$p_value))
})

testthat::test_that("pool probability surfaces are bounded and reproducible", {
  fixture <- pool_model_fixture()
  fit <- fit_pool_climate_response_models(fixture$prepared)
  contexts <- readr::read_csv(
    file.path(PROJECT_ROOT, "output", "tables", "climate_context_experiments.csv"),
    show_col_types = FALSE
  )
  first <- pool_climate_prediction_surfaces(
    fit, contexts, grid_size = 12L, n_draws = 30L,
    seed = 47163L, chunk_size = 50L
  )
  second <- pool_climate_prediction_surfaces(
    fit, contexts, grid_size = 12L, n_draws = 30L,
    seed = 47163L, chunk_size = 50L
  )
  testthat::expect_equal(nrow(first$grid), 12L * 12L * 3L * 4L)
  testthat::expect_equal(first$grid$conf_low, second$grid$conf_low)
  testthat::expect_equal(first$grid$conf_high, second$grid$conf_high)
  testthat::expect_equal(first$simulation$n_draws, 30L)
  testthat::expect_equal(first$simulation$seed, 47163L)
  testthat::expect_lt(
    first$simulation$resident_probability_sum_maximum_error, 1e-12
  )
  testthat::expect_true(all(first$grid$estimate >= 0 & first$grid$estimate <= 1))
  testthat::expect_true(all(first$grid$conf_low >= 0 & first$grid$conf_low <= 1))
  testthat::expect_true(all(first$grid$conf_high >= 0 & first$grid$conf_high <= 1))
  testthat::expect_true(all(first$grid$conf_low <= first$grid$conf_high))
  resident_sum <- first$grid |>
    dplyr::filter(.data$estimand %in% c(
      "non_significant_resident", "loser", "winner"
    )) |>
    dplyr::summarise(total = sum(.data$estimate), .by = "grid_id")
  testthat::expect_equal(resident_sum$total, rep(1, nrow(resident_sum)), tolerance = 1e-12)
  testthat::expect_true(all(
    !first$grid$supported |
      (first$grid$reviewed_climate_domain & first$grid$warming_hull_support)
  ))
  testthat::expect_equal(
    first$warming_levels$warming,
    as.numeric(stats::quantile(
      fixture$prepared$experiments$warming, c(0.25, 0.5, 0.75)
    ))
  )
})

testthat::test_that("surface interaction calculation matches direct model prediction", {
  fixture <- pool_model_fixture()
  fit <- fit_pool_climate_response_models(fixture$prepared)
  contexts <- readr::read_csv(
    file.path(PROJECT_ROOT, "output", "tables", "climate_context_experiments.csv"),
    show_col_types = FALSE
  )
  prediction <- pool_climate_prediction_surfaces(
    fit, contexts, grid_size = 8L, n_draws = 20L,
    seed = 47163L, chunk_size = 32L
  )
  focal <- prediction$grid |>
    dplyr::filter(.data$estimand == "colonizer") |>
    dplyr::slice(11L)
  direct <- stats::predict(
    fit$colonizer, newdata = focal, type = "response"
  )
  testthat::expect_equal(focal$estimate, as.numeric(direct), tolerance = 1e-12)
  resident_focal <- prediction$grid |>
    dplyr::filter(.data$grid_id == focal$grid_id, .data$estimand == "loser")
  resident_direct <- stats::predict(
    fit$resident, newdata = resident_focal, type = "response"
  )
  testthat::expect_equal(
    resident_focal$estimate,
    unname(resident_direct[, "loser"]), tolerance = 1e-12
  )
})

testthat::test_that("leave-one-experiment-out covers both models and all experiments", {
  prepared <- pool_model_fixture()$prepared
  fit <- fit_pool_climate_response_models(prepared)
  result <- pool_climate_leave_one_experiment_out(fit)
  testthat::expect_equal(nrow(result), 84L)
  testthat::expect_equal(dplyr::n_distinct(result$omitted_experiment), 42L)
  model_counts <- table(result$model)
  testthat::expect_equal(
    as.integer(model_counts[c("nonresident_binomial", "resident_multinomial")]),
    c(42L, 42L)
  )
  testthat::expect_true(all(result$converged))
  testthat::expect_true(all(result$warning == ""))
  testthat::expect_true(all(is.finite(
    result$maximum_absolute_coefficient_difference
  )))
  testthat::expect_true(all(is.finite(
    result$maximum_absolute_probability_difference
  )))
})

testthat::test_that("maintained target graph and figure inventory use the explicit allowlist", {
  target_text <- paste(readLines(file.path(PROJECT_ROOT, "_targets.R")), collapse = "\n")
  testthat::expect_false(stringr::str_detect(
    target_text, "climate_response_omnibus_table"
  ))
  testthat::expect_false(stringr::str_detect(
    target_text, "climate_response_multinomial[.]rds"
  ))
  testthat::expect_false(stringr::str_detect(
    target_text, "climate_response_pool_models[.]rds"
  ))
  testthat::expect_false(stringr::str_detect(
    target_text,
    paste0(
      "fit_pool_climate_response_models|fit_pool_binary_sensitivity|",
      "fit_multinomial_model|fit_strength_meta_model"
    )
  ))
  testthat::expect_false(stringr::str_detect(
    target_text, "plot/(questions|review|archive)"
  ))
  figure_target_names <- stringr::str_match_all(
    target_text, "tar_target\\(\\s*(figure_[A-Za-z0-9_]+)"
  )[[1]][, 2]
  testthat::expect_setequal(
    figure_target_names,
    c(
      "figure_transient_strength",
      "figure_temporal_flows",
      "figure_temporal_transition_matrices",
      "figure_destination_climate_warming_gradient",
      "figure_species_consistency",
      "figure_sensitivity"
    )
  )

  retained_root <- c(
    "climate_response_data_model_probability_prototype.pdf",
    "compare_rel.abund_vs.sqrt.rel.abund.pdf",
    "destination_climate_warming_gradient.pdf",
    "eg.sp.1.pdf",
    "eg.sp.2.pdf",
    "repeated_species_consistency.pdf",
    "response_by_climate_percentages.pdf",
    "response_by_climate_percentages_response_axis.pdf",
    "response_to_climate_alluvial.pdf",
    "response_to_family_alluvial.pdf",
    "temporal_category_flows.pdf",
    "temporal_transition_matrices.pdf",
    "trait_categories_with_overall.pdf",
    "transient_winning_losing_strength.pdf",
    "winners_losers.pdf"
  )
  retained_presentation <- file.path(
    "figures_presentation",
    c(
      "focal4_trait_phenology_warming.png",
      "winner_loser_boxplot.pdf",
      "winner_loser_boxplot.png",
      "winner_loser_boxplot_festuca.pdf",
      "winner_loser_boxplot_festuca.png",
      "winner_loser_boxplot_poa.pdf",
      "winner_loser_genus.pdf",
      "winner_loser_genus.png"
    )
  )
  expected_relative <- c(retained_root, retained_presentation)

  automatically_produced <- unique(unlist(stringr::str_extract_all(
    target_text, "[A-Za-z0-9_.-]+[.](?:pdf|png)"
  )))
  testthat::expect_false(any(automatically_produced %in% c(
    "climate_context_counts_odds_probabilities.pdf",
    "climate_response_binary_sensitivity_odds_ratios.pdf"
  )))

  manifest <- figure_manifest(PROJECT_ROOT)
  expected_manifest <- file.path(PROJECT_ROOT, "plot", expected_relative)
  testthat::expect_setequal(manifest$filename, expected_manifest)

  actual_relative <- list.files(
    file.path(PROJECT_ROOT, "plot"),
    pattern = "[.](pdf|png)$", recursive = TRUE, full.names = FALSE,
    ignore.case = TRUE
  )
  testthat::expect_setequal(actual_relative, expected_relative)
  generated_receipts <- list.files(
    file.path(PROJECT_ROOT, "plot"), pattern = "[.]generated$",
    recursive = TRUE, full.names = FALSE
  )
  testthat::expect_setequal(
    generated_receipts,
    paste0(setdiff(retained_root, c("eg.sp.1.pdf", "eg.sp.2.pdf", "winners_losers.pdf")), ".generated")
  )
  testthat::expect_setequal(
    basename(list.dirs(
      file.path(PROJECT_ROOT, "plot"), recursive = FALSE, full.names = TRUE
    )),
    "figures_presentation"
  )
  testthat::expect_false(any(dir.exists(file.path(
    PROJECT_ROOT, "plot", c("questions", "review", "archive")
  ))))
})
