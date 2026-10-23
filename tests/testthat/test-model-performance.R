testthat::test_that("species random variance uses the five stated coefficients", {
  context <- diag(c(0.3, 0.4, 0.5))
  dimnames(context) <- list(
    c("T_z", "P_z", "TP_z"), c("T_z", "P_z", "TP_z")
  )
  fit <- list(
    sigma2 = 0.2,
    tau2 = 0.1,
    G = matrix(0.1, 1, 1, dimnames = list("W_z", "W_z")),
    H = context
  )
  data <- data.frame(
    W_z = c(0, 1), T_z = c(1, 2), P_z = c(2, 1), TP_z = c(2, 2)
  )
  testthat::expect_equal(
    model_performance_random_variance(fit, data), c(4.1, 3.9)
  )
})

testthat::test_that("equal-gradient weights give each gradient equal influence", {
  data <- data.frame(gradient_id = c("a", "a", "a", "b"))
  weights <- equal_gradient_weights(data)
  testthat::expect_equal(sum(weights), 1)
  testthat::expect_equal(sum(weights[1:3]), 0.5)
  testthat::expect_equal(weights[[4]], 0.5)
})

testthat::test_that("marginal and conditional R2 use the stated denominator", {
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
  testthat::expect_equal(result$marginal_R2, 1 / expected_total)
  testthat::expect_equal(result$conditional_R2, 1.2 / expected_total)
  testthat::expect_gte(result$conditional_R2, result$marginal_R2)
  testthat::expect_identical(result$weighting, "equal_gradient")
})

testthat::test_that("cross-validation partitions are ordered and disjoint", {
  data <- data.frame(
    effect_id = paste0("e", 1:4),
    gradient_id = c("g1", "g1", "g2", "g2"),
    experiment_id = paste0("x", 1:4),
    species_id = paste0("s", 1:4),
    auc_effect = 1:4,
    T_z = 0, P_z = 0, W_z = 0, TP_z = 0
  )
  variances <- diag(c(0.1, 0.2, 0.3, 0.4))
  dimnames(variances) <- list(data$effect_id, data$effect_id)
  partition <- model_performance_partition(data, variances, "g2")
  testthat::expect_identical(partition$train$effect_id, c("e1", "e2"))
  testthat::expect_identical(partition$test$effect_id, c("e3", "e4"))
  testthat::expect_true(partition$effect_partitions_disjoint)
  testthat::expect_equal(unname(diag(partition$test_V)), c(0.3, 0.4))
})

testthat::test_that("equal-gradient RMSE and predictive R2 retain poor prediction", {
  statistics <- tibble::tibble(
    heldout_gradient = c("a", "b"),
    mean_observed = c(0, 2),
    mean_observed_squared = c(1, 5),
    MSE = c(4, 4)
  )
  result <- aggregate_equal_gradient_performance(statistics)
  testthat::expect_equal(unname(result[["RMSE"]]), 2)
  testthat::expect_lt(unname(result[["predictive_R2"]]), 0)
})

testthat::test_that("complete-gradient bootstrap is reproducible", {
  statistics <- tidyr::expand_grid(
    model_key = c("Climate", "LDMC"),
    candidate_role = c("standalone", "full"),
    heldout_gradient = paste0("g", 1:4)
  ) |>
    dplyr::mutate(
      model = model_performance_display(.data$model_key),
      n_test = 2L,
      mean_observed = rep(c(-1, 0, 1, 2), 4),
      mean_observed_squared = .data$mean_observed^2 + 1,
      MSE = ifelse(.data$model_key == "LDMC", 0.5, 1)
    )
  first <- summarise_model_cross_validation(
    statistics, bootstrap_replicates = 20L, seed = 7L
  )
  second <- summarise_model_cross_validation(
    statistics, bootstrap_replicates = 20L, seed = 7L
  )
  testthat::expect_identical(first, second)
  testthat::expect_setequal(unique(first$metric), c("RMSE", "predictive_R2"))
  testthat::expect_true(all(first$weighting == "equal_gradient"))
  testthat::expect_true(all(is.finite(first$estimate)))
})

testthat::test_that("the active workflow retains only the four S7 measures", {
  workflow_text <- paste(
    readLines(project_path("analysis_workflow.qmd"), warn = FALSE),
    collapse = "\n"
  )
  performance_text <- paste(
    readLines(
      project_path("R", "evaluation", "model_performance.R"), warn = FALSE
    ),
    collapse = "\n"
  )
  testthat::expect_match(workflow_text, "bootstrap_replicates = 5000L", fixed = TRUE)
  testthat::expect_match(workflow_text, "marginal_R2", fixed = TRUE)
  testthat::expect_match(workflow_text, "conditional_R2", fixed = TRUE)
  testthat::expect_match(workflow_text, "predictive_R2", fixed = TRUE)
  testthat::expect_match(workflow_text, "RMSE", fixed = TRUE)
  testthat::expect_false(grepl(
    "summarise_model_predictive_skill|model_performance_gaussian_score|effect_pooled",
    performance_text
  ))
})

testthat::test_that("the performance synthesis remains a manuscript appendix", {
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
  testthat::expect_identical(default_figure_status(path), "FINAL_MANUSCRIPT")
})
