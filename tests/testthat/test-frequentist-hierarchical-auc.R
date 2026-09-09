testthat::test_that("resident AUC formulas preserve hierarchy and random slopes", {
  testthat::expect_setequal(
    attr(stats::terms(auc_core_formula()), "term.labels"),
    c("T_z", "P_z", "W_z", "T_z:P_z", "T_z:W_z", "P_z:W_z", "T_z:P_z:W_z")
  )
  random <- vapply(
    auc_random_structure("species_slope"),
    function(formula) paste(deparse(formula), collapse = " "), character(1)
  )
  testthat::expect_equal(random, c("~1 | species_id", "~0 + W_z | species_id"))
  testthat::expect_equal(
    auc_model_specifications()$model_id,
    c(
      "S1_PRIMARY", "S2_EXP_RE", "S3_GRADIENT_RE", "S4_DESIGN",
      "S5_SPECIES_INTERCEPT_ONLY"
    )
  )
})

testthat::test_that("resident preparation uses continuous effects, not categories", {
  path <- project_path("output", "tables", "winner_loser_results.csv")
  testthat::skip_if_not(file.exists(path), "winner/loser results unavailable")
  data <- read_resident_auc_data(PROJECT_ROOT)
  testthat::expect_setequal(unique(data$horizon), c("max", "3"))
  testthat::expect_true(all(data$auc_effect == data$auc_effect))
  testthat::expect_true(all(is.finite(data$auc_se) & data$auc_se > 0))
  testthat::expect_false(any(c("category", "significant", "p.value") %in% names(data)))
  testthat::expect_false(anyDuplicated(data$effect_id) > 0L)
  testthat::expect_false(any(
    data$species_id %in% resident_downstream_model_exclusions()
  ))
  testthat::expect_equal(
    dplyr::n_distinct(data$experiment_id[data$horizon == "max"]), 42L
  )
})

testthat::test_that("climate scaling gives equal weight to experiments", {
  path <- project_path("output", "tables", "analysis_dataset.csv")
  testthat::skip_if_not(file.exists(path), "analysis dataset unavailable")
  scaled <- scale_resident_auc_data(read_resident_auc_data(PROJECT_ROOT))
  testthat::expect_true(all(scaled$scaling$weighting == "distinct experiments"))
  experiments <- scaled$data |>
    dplyr::distinct(.data$experiment_id, .data$T_z, .data$P_z, .data$W_z)
  testthat::expect_equal(mean(experiments$T_z), 0, tolerance = 1e-12)
  testthat::expect_equal(stats::sd(experiments$P_z), 1, tolerance = 1e-12)
  testthat::expect_equal(stats::sd(experiments$W_z), 1, tolerance = 1e-12)
})

testthat::test_that("covariance validation never repairs an invalid matrix", {
  data <- tibble::tibble(
    effect_id = c("a", "b"), auc_se = c(1, 2)
  )
  valid <- diag(c(1, 4))
  dimnames(valid) <- list(data$effect_id, data$effect_id)
  testthat::expect_true(validate_auc_covariance(valid, data)$valid)

  invalid <- matrix(c(1, 3, 0, 4), 2, 2)
  dimnames(invalid) <- list(data$effect_id, data$effect_id)
  selected <- select_auc_covariance(invalid, data, reconstruction_error = "asymmetric")
  testthat::expect_identical(
    selected$selected_v_type, "diagonal sampling-variance approximation"
  )
  testthat::expect_equal(as.matrix(selected$V), valid)
  testthat::expect_match(selected$fallback_reason, "asymmetric")
  testthat::expect_equal(invalid[1, 2], 0)
})

testthat::test_that("reconstructed covariance matches the saved contrasts", {
  path <- project_path("output", "intermediate", "auc_by_plot.rds")
  testthat::skip_if_not(file.exists(path), "plot AUC data unavailable")
  scaled <- scale_resident_auc_data(read_resident_auc_data(PROJECT_ROOT))$data
  experiment <- levels(scaled$experiment_id)[[1]]
  data <- scaled |>
    dplyr::filter(.data$horizon == "max", .data$experiment_id == experiment) |>
    droplevels()
  reconstructed <- reconstruct_auc_covariance(data, readRDS(path))
  testthat::expect_true(reconstructed$audit$full_covariance_valid)
  testthat::expect_equal(
    unname(sqrt(Matrix::diag(reconstructed$V))), data$auc_se, tolerance = 1e-8
  )
  testthat::expect_identical(rownames(reconstructed$V), data$effect_id)
})

testthat::test_that("trait subsets freeze rows and subset V by effect_id", {
  path <- project_path("output", "tables", "analysis_dataset.csv")
  testthat::skip_if_not(file.exists(path), "analysis dataset unavailable")
  data <- scale_resident_auc_data(read_resident_auc_data(PROJECT_ROOT))$data |>
    dplyr::filter(.data$horizon == "max") |>
    droplevels()
  V <- diag(data$auc_se^2)
  dimnames(V) <- list(data$effect_id, data$effect_id)
  prepared <- trait_complete_data(data, "LDMC")
  subset <- subset_auc_v(V, prepared$data$effect_id)
  testthat::expect_equal(nrow(subset), nrow(prepared$data))
  testthat::expect_identical(rownames(subset), prepared$data$effect_id)
  testthat::expect_true(all(is.finite(prepared$data$trait_z)))
  testthat::expect_equal(prepared$scaling$weighting, "distinct species")
})

testthat::test_that("warming derivative agrees with an analytic linear example", {
  set.seed(10)
  data <- tibble::tibble(
    auc_effect = stats::rnorm(80), auc_se = rep(0.2, 80),
    T_z = stats::rnorm(80), P_z = stats::rnorm(80), W_z = stats::rnorm(80),
    species_id = factor(rep(seq_len(20), each = 4)),
    experiment_id = factor(rep(seq_len(8), each = 10)),
    gradient_id = factor(rep(seq_len(4), each = 20)),
    effect_id = paste0("e", seq_len(80))
  )
  V <- diag(data$auc_se^2)
  model <- fit_auc_rma(data, V, auc_core_formula(), "species_intercept")
  grid <- data[1:3, c("T_z", "P_z", "W_z")]
  derivative <- predict_auc_derivative(model, list(model = NULL), grid)
  beta <- stats::coef(model$fit)
  expected <- beta[["W_z"]] + beta[["T_z:W_z"]] * grid$T_z +
    beta[["P_z:W_z"]] * grid$P_z + beta[["T_z:P_z:W_z"]] * grid$T_z * grid$P_z
  testthat::expect_equal(derivative$estimate, as.numeric(expected), tolerance = 1e-7)
})

testthat::test_that("support mask excludes points beyond observed climate space", {
  observed <- tibble::tibble(
    experiment_id = factor(1:5),
    T_z = c(0, 1, 0, 0, 0.2), P_z = c(0, 0, 1, 0, 0.2),
    W_z = c(0, 0, 0, 1, 0.2)
  )
  grid <- tibble::tibble(
    T_z = c(0.1, 3), P_z = c(0.1, 3), W_z = c(0.1, 3)
  )
  support <- auc_support_mask(observed, grid)$grid
  testthat::expect_true(support$supported[[1]])
  testthat::expect_false(support$supported[[2]])
})

testthat::test_that("pre-model audit scopes every metric to its requested horizon", {
  result_path <- project_path("output", "tables", "winner_loser_results.csv")
  auc_path <- project_path("output", "intermediate", "auc_by_plot.rds")
  testthat::skip_if_not(
    file.exists(result_path) && file.exists(auc_path), "resident AUC inputs unavailable"
  )
  data <- scale_resident_auc_data(read_resident_auc_data(PROJECT_ROOT))$data
  auc_by_plot <- readRDS(auc_path)
  audits <- dplyr::bind_rows(
    auc_data_audit(data, "max", auc_by_plot),
    auc_data_audit(data, "3", auc_by_plot)
  )
  extract_metric <- function(horizon, metric) {
    audits$value[audits$horizon == horizon & audits$metric == metric]
  }
  testthat::expect_equal(extract_metric("max", "total_effect_rows"), 1604)
  testthat::expect_equal(extract_metric("3", "total_effect_rows"), 1604)
  testthat::expect_equal(
    extract_metric("max", "species_experiment_combinations"), 1604
  )
  testthat::expect_equal(
    extract_metric("3", "species_experiment_combinations"), 1604
  )
  for (horizon in c("max", "3")) {
    signs <- sum(vapply(
      c("negative_auc_effects", "zero_auc_effects", "positive_auc_effects"),
      function(metric) extract_metric(horizon, metric), numeric(1)
    ))
    testthat::expect_equal(signs, extract_metric(horizon, "total_effect_rows"))
  }
  distinct_experiments <- data |>
    dplyr::filter(.data$horizon == "max") |>
    dplyr::distinct(.data$experiment_id, .data$T)
  testthat::expect_equal(
    extract_metric("max", "T_mean"), mean(distinct_experiments$T), tolerance = 1e-12
  )
})

testthat::test_that("pre-model validation rejects invalid inferential inputs", {
  path <- project_path("output", "tables", "winner_loser_results.csv")
  testthat::skip_if_not(file.exists(path), "resident AUC inputs unavailable")
  data <- scale_resident_auc_data(read_resident_auc_data(PROJECT_ROOT))$data |>
    dplyr::filter(.data$horizon == "max")
  testthat::expect_true(all(validate_auc_pre_model_data(data, "max")$passed))

  duplicate <- dplyr::bind_rows(data, data[1, ])
  testthat::expect_error(
    validate_auc_pre_model_data(duplicate, "max"), "duplicate effect_id"
  )
  bad_se <- data
  bad_se$auc_se[[1]] <- 0
  testthat::expect_error(validate_auc_pre_model_data(bad_se, "max"), "SEs")
  missing_id <- data
  missing_id$species_id[[1]] <- NA
  testthat::expect_error(
    validate_auc_pre_model_data(missing_id, "max"), "missing species_id"
  )
  inconsistent <- data
  same_experiment <- which(
    inconsistent$experiment_id == inconsistent$experiment_id[[1]]
  )
  inconsistent$T[same_experiment[[1]]] <- inconsistent$T[same_experiment[[1]]] + 1
  testthat::expect_error(
    validate_auc_pre_model_data(inconsistent, "max"),
    "inconsistent experiment-level covariates"
  )
})

testthat::test_that("audit covers traits, plot replication, and all fields", {
  path <- project_path("output", "intermediate", "auc_by_plot.rds")
  testthat::skip_if_not(file.exists(path), "plot AUC data unavailable")
  data <- scale_resident_auc_data(read_resident_auc_data(PROJECT_ROOT))$data
  trait_coverage <- auc_trait_coverage(data, c("max", "3"))
  testthat::expect_equal(nrow(trait_coverage), 6L)
  testthat::expect_setequal(
    trait_coverage$trait, c("LDMC", "LMA", "vegetative_height")
  )
  testthat::expect_true(all(trait_coverage$effect_coverage > 0 &
                              trait_coverage$effect_coverage < 1))
  experiments <- auc_experiment_audit(data, "max", readRDS(path))
  testthat::expect_equal(nrow(experiments), 42L)
  testthat::expect_true(all(experiments$origin_plots > 0))
  testthat::expect_true(all(experiments$warmed_plots > 0))
  missingness <- auc_missingness_audit(data, c("max", "3"))
  testthat::expect_setequal(unique(missingness$horizon), c("max", "3"))
  testthat::expect_true(all(c("auc_effect", "auc_se", "LDMC") %in% missingness$column))
})

testthat::test_that("audit artifacts precede and survive a later failure", {
  path <- project_path("output", "intermediate", "auc_by_plot.rds")
  testthat::skip_if_not(file.exists(path), "plot AUC data unavailable")
  testthat::skip_if_not(rmarkdown::pandoc_available(), "Pandoc unavailable")
  scaled <- scale_resident_auc_data(read_resident_auc_data(PROJECT_ROOT))
  output_dir <- tempfile("resident-auc-audit-")
  dir.create(output_dir)
  write_auc_pre_model_audit(
    scaled$data, c("max", "3"), readRDS(path), scaled,
    output_dir, PROJECT_ROOT, render_html = TRUE
  )
  testthat::expect_error(stop("simulated downstream model failure"), "simulated")
  expected <- c(
    "audit/data_audit.html", "audit/pre_model_validation.csv",
    "audit/experiment_audit.csv", "audit/trait_coverage.csv",
    "audit/figures/effect_and_biological_coverage.png",
    "audit/figures/experiment_and_design_coverage.pdf",
    "audit/figures/inferential_field_missingness.png"
  )
  paths <- file.path(output_dir, expected)
  testthat::expect_true(all(file.exists(paths)))
  testthat::expect_true(all(file.info(paths)$size > 0))
  html <- paste(readLines(file.path(output_dir, "audit", "data_audit.html")), collapse = "\n")
  testthat::expect_match(html, "data:image/png;base64", fixed = TRUE)
  testthat::expect_match(html, ">\\s*526\\s*<")
  testthat::expect_false(grepl("src=[^>]*[.]png", html))
})

testthat::test_that("saved visual report is self-contained and source-driven", {
  output_dir <- project_path("output", "frequentist_hierarchical_auc")
  testthat::skip_if_not(
    file.exists(file.path(output_dir, "tables", "model_robustness_ledger.csv")),
    "resident AUC model outputs unavailable"
  )
  testthat::skip_if_not(rmarkdown::pandoc_available(), "Pandoc unavailable")
  rendered <- render_frequentist_hierarchical_auc_visual_report(
    output_dir, report_file = tempfile(fileext = ".html"), root = PROJECT_ROOT
  )
  html <- paste(readLines(rendered, warn = FALSE), collapse = "\n")
  testthat::expect_match(html, "Data audit and effective replication", fixed = TRUE)
  testthat::expect_match(html, "Trait models: selection versus evidence", fixed = TRUE)
  testthat::expect_match(html, "experiment-clustered CR1", ignore.case = TRUE)
  saved_audit <- readr::read_csv(
    file.path(output_dir, "tables", "data_audit.csv"),
    show_col_types = FALSE
  )
  saved_species <- saved_audit$value[
    saved_audit$horizon == "max" & saved_audit$metric == "unique_species"
  ]
  testthat::expect_match(
    html, paste0(">\\s*", format(saved_species, big.mark = ","), "\\s*<")
  )
  testthat::expect_match(html, "data:image/png;base64", fixed = TRUE)
  testthat::expect_false(grepl("src=[^>]*[.]png", html))
})
