testthat::test_that("review preparation preserves pools and denominator-selected families", {
  path <- project_path("output", "tables", "analysis_dataset.csv")
  testthat::skip_if_not(file.exists(path))
  data <- readr::read_csv(path, show_col_types = FALSE)
  prepared <- prepare_biological_context_review(data, family_minimum = 50L)

  testthat::expect_equal(nrow(prepared$resident), 1604L)
  testthat::expect_equal(nrow(prepared$nonresident), 611L)
  testthat::expect_equal(sum(prepared$resident$loser_event), 130L)
  testthat::expect_equal(sum(prepared$resident$winner_event), 69L)
  testthat::expect_equal(sum(prepared$nonresident$colonizer_event), 34L)
  testthat::expect_equal(prepared$snapshot$n_common_families, 13L)
  testthat::expect_true(all(
    prepared$family_selection$n_species_experiment[
      prepared$family_selection$common_family
    ] >= 50L
  ))
  testthat::expect_true(all(
    prepared$family_selection$n_species_experiment[
      !prepared$family_selection$common_family
    ] < 50L
  ))
  testthat::expect_equal(
    sort(unique(as.character(prepared$data$family_model))),
    sort(c(
      prepared$family_selection$family[
        prepared$family_selection$common_family
      ],
      "Other"
    ))
  )
  keys <- prepared$data |>
    dplyr::count(.data$accepted_species, .data$experiment_id)
  testthat::expect_true(all(keys$n == 1L))
})

testthat::test_that("functional groups and traits use the reviewed inference rules", {
  path <- project_path("output", "tables", "analysis_dataset.csv")
  testthat::skip_if_not(file.exists(path))
  data <- readr::read_csv(path, show_col_types = FALSE)
  prepared <- prepare_biological_context_review(data)
  group_data <- review_subsystem_data(prepared, "functional_group")

  testthat::expect_setequal(
    unique(as.character(group_data$resident$functional_group_model)),
    c("forb", "graminoid")
  )
  testthat::expect_gt(prepared$snapshot$n_functional_group_excluded, 0L)

  purrr::walk(c("LDMC", "LMA", "vegetative_height"), function(trait) {
    distinct_values <- data |>
      dplyr::filter(
        as.character(.data$transform) == "raw",
        as.character(.data$horizon) == "max"
      ) |>
      filter_downstream_model_species("accepted_name") |>
      dplyr::filter(is.finite(.data[[trait]]), .data[[trait]] > 0) |>
      dplyr::distinct(.data$accepted_name, value = .data[[trait]]) |>
      dplyr::pull(.data$value) |>
      log()
    definition <- prepared$trait_scaling |>
      dplyr::filter(.data$trait == .env$trait)
    testthat::expect_equal(definition$center, mean(distinct_values))
    testthat::expect_equal(definition$scale, stats::sd(distinct_values))
    testthat::expect_equal(definition$n_species, length(distinct_values))
    complete <- review_subsystem_data(prepared, "trait", trait)
    testthat::expect_false(anyNA(
      complete$resident[[paste0("z_trait_", trait)]]
    ))
    testthat::expect_false(anyNA(
      complete$nonresident[[paste0("z_trait_", trait)]]
    ))
  })
})

testthat::test_that("review formulas contain only prespecified climate interactions", {
  climate <- biological_context_base_rhs()
  testthat::expect_true(grepl("z_destination_T:z_log_destination_P", climate))
  testthat::expect_true(grepl("z_destination_T:z_warming", climate))
  testthat::expect_true(grepl("z_log_destination_P:z_warming", climate))
  testthat::expect_false(grepl("diff_P", climate))
  testthat::expect_false(grepl(
    "z_destination_T:z_log_destination_P:z_warming", climate, fixed = TRUE
  ))

  binary <- biological_context_formula(
    "binary_event", "functional_group_model", mixed_binomial = TRUE
  )
  formula_text <- paste(deparse(binary), collapse = " ")
  testthat::expect_true(grepl("(1 | accepted_species)", formula_text, fixed = TRUE))
  testthat::expect_true(grepl("(1 | experiment_id)", formula_text, fixed = TRUE))

  family <- biological_context_specs("family")
  testthat::expect_setequal(
    names(family),
    c("climate", "main", "temperature", "precipitation", "warming")
  )
  testthat::expect_true(grepl(
    "family_model:z_destination_T",
    paste(deparse(family$temperature), collapse = " "), fixed = TRUE
  ))
})

testthat::test_that("zero-event diagnostics distinguish supported from empty cells", {
  supported <- tibble::tibble(
    functional_group_model = factor(rep(c("forb", "graminoid"), each = 4)),
    binary_event = rep(c(FALSE, TRUE), 4)
  )
  testthat::expect_length(
    review_zero_event_levels(
      supported, "binomial_glmm", "loser", "functional_group_model"
    ),
    0L
  )
  unsupported <- supported |>
    dplyr::mutate(
      binary_event = .data$functional_group_model == "forb" & .data$binary_event
    )
  testthat::expect_match(
    review_zero_event_levels(
      unsupported, "binomial_glmm", "loser", "functional_group_model"
    ),
    "graminoid::0/4"
  )
})

testthat::test_that("multivariate-normal review draws are reproducible", {
  coefficients <- c(a = 0.2, b = -0.1)
  covariance <- matrix(c(0.2, 0.03, 0.03, 0.1), nrow = 2)
  set.seed(47163)
  first <- pool_climate_mvn_draws(coefficients, covariance, 20L)
  set.seed(47163)
  second <- pool_climate_mvn_draws(coefficients, covariance, 20L)
  testthat::expect_equal(first, second)
})
