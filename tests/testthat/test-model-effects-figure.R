testthat::test_that("fixed-effect labels retain variables and interactions", {
  terms <- c(
    "intrcpt", "Trait_z:T_z:W_z",
    "GraminoidStatusgraminoid:P_z", "PEV2_z:TP_z"
  )
  keys <- c("Climate", "LDMC", "Graminoid", "Phylogeny")
  testthat::expect_identical(
    model_effect_pretty_term(terms, keys),
    c(
      "Intercept",
      "LDMC × Destination temperature × Warming",
      "Graminoid × Destination precipitation",
      "PEV2 × Temperature × precipitation"
    )
  )
})

testthat::test_that("the fixed-effect table retains all seven current models", {
  climate_path <- project_path("output", "tables", "resident_climate_robust_inference.csv")
  predictor_path <- project_path(
    "output", "tables", "predictor_model_robust_inference.csv"
  )
  r2_path <- project_path("output", "tables", "model_fit_r2.csv")
  testthat::skip_if_not(all(file.exists(c(
    climate_path, predictor_path, r2_path
  ))))
  stored_r2 <- readr::read_csv(r2_path, show_col_types = FALSE)
  testthat::skip_if_not(
    "weighting" %in% names(stored_r2),
    "Stored model_fit_r2.csv predates the approved Figure S7 revision."
  )
  fixed <- assemble_model_fixed_effects(
    readr::read_csv(climate_path, show_col_types = FALSE),
    readr::read_csv(predictor_path, show_col_types = FALSE),
    stored_r2
  )
  testthat::expect_equal(nrow(fixed), 111L)
  testthat::expect_identical(
    sort(unique(fixed$model_key)),
    sort(model_performance_key()$model_key)
  )
  testthat::expect_identical(
    anyDuplicated(fixed[c("model_key", "term")]), 0L
  )
  testthat::expect_identical(
    fixed$supported,
    fixed$CR1_CI_low > 0 | fixed$CR1_CI_high < 0
  )
})

testthat::test_that("predictor profile jobs cover five components per model", {
  fake_sequence <- function(key) {
    list(key = key, fits = list(full_REML = list(fit = list())))
  }
  jobs <- build_predictor_variance_profile_jobs(
    fake_sequence("LDMC"),
    fake_sequence("LMA"),
    fake_sequence("Height"),
    fake_sequence("Graminoid"),
    fake_sequence("Woodiness"),
    fake_sequence("Phylogeny")
  )
  testthat::expect_length(jobs, 30L)
  testthat::expect_true(all(vapply(
    jobs,
    function(job) job$specification$component[[1]] %in% c(
      "species_intercept", "species_W", "species_T",
      "species_P", "species_TP"
    ),
    logical(1)
  )))
})

testthat::test_that("variance profiles are displayed on the SD scale", {
  components <- c(
    "species_intercept", "species_W", "species_T", "species_P", "species_TP"
  )
  make_profiles <- function(model_key) {
    tibble::tibble(
      model_key = model_key,
      model = model_performance_display(model_key),
      component = components,
      estimate = seq_along(components)^2 / 10000,
      CI_low = seq_along(components)^2 / 20000,
      CI_high = seq_along(components)^2 / 5000,
      warnings = "",
      error = ""
    )
  }
  climate <- make_profiles("Climate") |>
    dplyr::select(-"model_key", -"model")
  predictors <- lapply(
    model_performance_key()$model_key[-1L], make_profiles
  )
  profiles <- assemble_model_random_effect_profiles(
    climate, predictors
  )
  testthat::expect_equal(nrow(profiles), 35L)
  testthat::expect_equal(profiles$SD_estimate^2, profiles$estimate)
  testthat::expect_equal(profiles$SD_CI_low^2, profiles$CI_low)
  testthat::expect_equal(profiles$SD_CI_high^2, profiles$CI_high)
  testthat::expect_true(all(profiles$supported))
})

testthat::test_that("model-effect targets keep profiling separate from primary fits", {
  target_text <- paste(
    readLines(project_path("R", "targets", "definitions", "models.R")),
    collapse = "\n"
  )
  testthat::expect_match(
    target_text, "pattern = map(model_random_profile_job)", fixed = TRUE
  )
  testthat::expect_match(
    target_text, "model_effects = model_effects_release_validation",
    fixed = TRUE
  )
})

testthat::test_that("the model-effects figure is a reviewable manuscript appendix", {
  path <- "figures/manuscript/appendix_model_fixed_and_random_effects.pdf"
  producer <- known_figure_producer(path)
  testthat::expect_identical(
    producer$generating_script_if_known,
    "R/evaluation/model_effects.R"
  )
  testthat::expect_identical(
    producer$target_if_known,
    "appendix_model_fixed_and_random_effects_pdf"
  )
  testthat::expect_identical(
    default_figure_status(path), "CURRENT_ANALYSIS_FIGURE"
  )
})
