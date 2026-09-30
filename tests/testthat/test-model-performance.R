testthat::test_that("species random covariance preserves the five canonical coefficients", {
  context <- diag(c(0.3, 0.4, 0.5))
  dimnames(context) <- list(c("T_z", "P_z", "TP_z"), c("T_z", "P_z", "TP_z"))
  fit <- list(
    sigma2 = 0.2,
    tau2 = 0.1,
    G = matrix(0.1, 1, 1, dimnames = list("W_z", "W_z")),
    H = context
  )
  data <- data.frame(
    species_id = c("a", "a", "b"),
    W_z = c(0, 1, 2),
    T_z = c(1, 2, 3),
    P_z = c(2, 1, 0),
    TP_z = c(2, 2, 0)
  )
  covariance <- model_performance_random_covariance(fit, data)
  testthat::expect_equal(covariance[1, 2], 3.6)
  testthat::expect_equal(covariance[1, 3], 0)
  testthat::expect_equal(
    diag(covariance), model_performance_random_variance(fit, data)
  )
})

testthat::test_that("equal-gradient and effect-pooled weights are explicit", {
  data <- data.frame(gradient_id = c("a", "a", "a", "b"))
  equal <- model_performance_weights(data, "equal_gradient")
  pooled <- model_performance_weights(data, "effect_pooled")
  testthat::expect_equal(sum(equal), 1)
  testthat::expect_equal(sum(equal[1:3]), 0.5)
  testthat::expect_equal(equal[[4]], 0.5)
  testthat::expect_equal(pooled, rep(0.25, 4))
})

testthat::test_that("hierarchical marginal and conditional R2 use the stated denominator", {
  context <- matrix(0, 3, 3)
  dimnames(context) <- list(
    c("T_z", "P_z", "TP_z"), c("T_z", "P_z", "TP_z")
  )
  fit <- list(
    data = data.frame(
      species_id = paste0("s", 1:4),
      gradient_id = c("a", "a", "b", "b"),
      W_z = 0, T_z = 0, P_z = 0, TP_z = 0
    ),
    k = 4L,
    X = cbind(intrcpt = 1, x = c(-1, 1, -1, 1)),
    coefficients = c(intrcpt = 0, x = 1),
    V = diag(0.5, 4),
    sigma2 = 0.2,
    tau2 = 0,
    G = matrix(0, 1, 1, dimnames = list("W_z", "W_z")),
    H = context
  )
  result <- model_performance_fit_r2(fit, "Climate", "standalone")
  expected_total <- 1 + 0.2 + 0.5
  testthat::expect_equal(result$marginal_R2, rep(1 / expected_total, 2))
  testthat::expect_equal(result$conditional_R2, rep(1.2 / expected_total, 2))
  testthat::expect_true(all(result$conditional_R2 >= result$marginal_R2))
})

testthat::test_that("Gaussian conditioning has the correct unseen-species limit", {
  train_y <- c(1, -1)
  train_mean <- c(0, 0)
  train_covariance <- diag(2)
  train_design <- cbind(intrcpt = 1, x = c(-1, 1))
  test_mean <- c(0.2, -0.2)
  test_design <- cbind(intrcpt = 1, x = c(0.5, -0.5))
  test_marginal <- diag(c(0.7, 0.8))
  beta_covariance <- diag(c(0.1, 0.2))
  conditioned <- condition_model_performance_gaussian(
    train_y, train_mean, train_covariance, train_design,
    test_mean, test_design, test_marginal,
    cross_random_covariance = matrix(0, 2, 2),
    beta_covariance = beta_covariance
  )
  expected_covariance <- test_marginal +
    test_design %*% beta_covariance %*% t(test_design)
  testthat::expect_equal(conditioned$mean, test_mean)
  testthat::expect_equal(conditioned$covariance, expected_covariance)
})

testthat::test_that("Gaussian conditioning uses shared-species training information", {
  conditioned <- condition_model_performance_gaussian(
    train_y = 2,
    train_mean = 0,
    train_covariance = matrix(2, 1, 1),
    train_design = matrix(1, 1, 1),
    test_mean = 0,
    test_design = matrix(1, 1, 1),
    test_marginal_covariance = matrix(3, 1, 1),
    cross_random_covariance = matrix(1, 1, 1),
    beta_covariance = matrix(0.25, 1, 1)
  )
  testthat::expect_equal(conditioned$mean, 1)
  testthat::expect_equal(conditioned$covariance, matrix(2.5625, 1, 1))
})

testthat::test_that("Gaussian log score uses the complete covariance", {
  score <- model_performance_gaussian_score(
    observed = c(0, 0), predicted = c(0, 0), covariance = diag(2)
  )
  testthat::expect_equal(
    unname(score[["negative_log_predictive_density_total"]]),
    log(2 * pi)
  )
  testthat::expect_equal(
    unname(score[["negative_log_predictive_density"]]),
    0.5 * log(2 * pi)
  )
})

testthat::test_that("cross-validated point metrics retain negative predictive R2", {
  statistics <- tibble::tibble(
    heldout_gradient = c("a", "b"),
    n = c(2, 1),
    sum_y = c(2, 3),
    sum_y2 = c(2, 9),
    sum_error = c(2, 2),
    sum_absolute_error = c(2, 2),
    sum_squared_error = c(2, 4),
    sum_covered = c(2, 1),
    sum_interval_width = c(4, 2),
    n_effects = c(2, 1),
    negative_log_predictive_density_total = c(3, 2)
  )
  equal <- aggregate_model_performance_statistics(
    statistics, "equal_gradient"
  )
  pooled <- aggregate_model_performance_statistics(
    statistics, "effect_pooled"
  )
  testthat::expect_lt(equal[["predictive_R2"]], 0)
  testthat::expect_false(isTRUE(all.equal(equal[["RMSE"]], pooled[["RMSE"]])))
  testthat::expect_equal(equal[["coverage_95"]], 1)
})

testthat::test_that("metric and paired-skill bootstraps are reproducible", {
  make_predictions <- function(role, scope) {
    base <- tidyr::expand_grid(
      heldout_gradient = c("g1", "g2"), row = seq_len(4L)
    ) |>
      dplyr::mutate(
        model_key = "LDMC", model = "LDMC", candidate_role = role,
        prediction_scope = scope,
        effect_id = paste(.data$heldout_gradient, .data$row, sep = "_"),
        species_id = paste0("s", .data$row),
        observed_delta_auc = c(-1, 0, 1, 2, -2, -1, 1, 3),
        predicted_delta_auc = if (role == "full") {
          0.8 * .data$observed_delta_auc
        } else {
          0.5 * .data$observed_delta_auc
        },
        residual = .data$predicted_delta_auc - .data$observed_delta_auc,
        predictive_se = 1,
        prediction_interval_low = .data$predicted_delta_auc - 1.96,
        prediction_interval_high = .data$predicted_delta_auc + 1.96,
        interval_covered = .data$observed_delta_auc >=
          .data$prediction_interval_low &
          .data$observed_delta_auc <= .data$prediction_interval_high,
        interval_width = 3.92,
        species_seen_in_training = .data$row <= 2L,
        climate_supported = TRUE
      ) |>
      dplyr::select(-"row")
    if (scope == "repeat_species_conditional") {
      base <- dplyr::filter(base, .data$species_seen_in_training)
    }
    base
  }
  predictions <- dplyr::bind_rows(
    make_predictions("full", "general"),
    make_predictions("matched_baseline", "general"),
    make_predictions("full", "repeat_species_conditional"),
    make_predictions("matched_baseline", "repeat_species_conditional")
  )
  gaussian <- tidyr::expand_grid(
    heldout_gradient = c("g1", "g2"),
    candidate_role = c("full", "matched_baseline"),
    prediction_scope = c("general", "repeat_species_conditional"),
    support_stratum = c("all", "in_support")
  ) |>
    dplyr::mutate(
      model_key = "LDMC", model = "LDMC",
      evaluation_subset = dplyr::if_else(
        .data$prediction_scope == "general", "all_eligible", "repeat_species"
      ),
      n_effects = dplyr::if_else(
        .data$prediction_scope == "general", 4L, 2L
      ),
      negative_log_predictive_density_total = dplyr::if_else(
        .data$candidate_role == "full", 2, 3
      )
    )
  gaussian_repeat_general <- gaussian |>
    dplyr::filter(.data$prediction_scope == "general") |>
    dplyr::mutate(
      evaluation_subset = "repeat_species", n_effects = 2L,
      negative_log_predictive_density_total =
        .data$negative_log_predictive_density_total / 2
    )
  gaussian <- dplyr::bind_rows(gaussian, gaussian_repeat_general)
  first <- summarise_model_cross_validation(
    predictions, gaussian, bootstrap_replicates = 20L, seed = 7L
  )
  second <- summarise_model_cross_validation(
    predictions, gaussian, bootstrap_replicates = 20L, seed = 7L
  )
  skill <- summarise_model_predictive_skill(
    predictions, gaussian, bootstrap_replicates = 20L, seed = 7L
  )
  testthat::expect_identical(first, second)
  testthat::expect_true(all(is.finite(first$estimate)))
  testthat::expect_true(all(is.finite(skill$estimate)))
  testthat::expect_true(all(
    skill$estimate[skill$metric == "RMSE_percent_improvement"] > 0
  ))
})

testthat::test_that("the canonical covariance is block-independent by gradient", {
  data_path <- project_path("output", "data", "resident_effects.csv")
  covariance_path <- project_path("output", "data", "resident_sampling_covariance.rds")
  testthat::skip_if_not(file.exists(data_path) && file.exists(covariance_path))
  data <- readr::read_csv(data_path, show_col_types = FALSE)
  covariance <- readRDS(covariance_path)
  testthat::expect_equal(dplyr::n_distinct(data$gradient_id), 20L)
  for (gradient in unique(as.character(data$gradient_id))) {
    test <- as.character(data$gradient_id) == gradient
    cross <- covariance[!test, test, drop = FALSE]
    testthat::expect_lte(max(abs(as.numeric(cross))), 1e-12)
  }
})

testthat::test_that("constructed species covariance matches canonical fitted M minus V", {
  model_path <- project_path("output", "models", "resident_climate_model.rds")
  testthat::skip_if_not(file.exists(model_path))
  model <- readRDS(model_path)
  index <- seq_len(min(100L, model$k))
  expected <- as.matrix(model$M[index, index] - model$V[index, index])
  reconstructed <- model_performance_random_covariance(
    model, model$data[index, , drop = FALSE]
  )
  testthat::expect_equal(reconstructed, expected, tolerance = 1e-12)
})

testthat::test_that("performance targets preserve validation-only boundaries", {
  target_text <- paste(
    readLines(project_path("R", "targets", "definitions", "models.R")),
    collapse = "\n"
  )
  testthat::expect_match(
    target_text, "pattern = map(model_performance_gradient)", fixed = TRUE
  )
  testthat::expect_match(
    target_text, "bootstrap_replicates = 5000L", fixed = TRUE
  )
  testthat::expect_match(
    target_text, "predictive_performance = model_performance_release_validation",
    fixed = TRUE
  )
})

testthat::test_that("the performance synthesis is a manuscript appendix", {
  path <- "figures/manuscript/appendix_model_predictive_performance.pdf"
  producer <- known_figure_producer(path)
  testthat::expect_identical(
    producer$generating_script_if_known,
    "R/evaluation/model_performance.R"
  )
  testthat::expect_identical(
    producer$target_if_known,
    "appendix_model_predictive_performance_pdf"
  )
  testthat::expect_identical(
    default_figure_status(path), "FINAL_MANUSCRIPT"
  )
})
