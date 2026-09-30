source(project_path("config", "workflow.R"))
source(project_path("R", "models", "meta_regression.R"))
source(project_path("R", "models", "predictors.R"))
source(project_path("R", "models", "species_context.R"))
source(project_path("R", "evaluation", "heterogeneity.R"))
source(project_path("R", "models", "family.R"))

test_that("the validated formulas preserve the settled resident models", {
  formulas <- predictor_formulas()
  expect_setequal(
    names(formulas),
    c(
      "climate", "continuous_full", "Graminoid_full",
      "Woodiness_full", "Phylogeny_full"
    )
  )
  expect_false(any(grepl("additive", names(formulas), ignore.case = TRUE)))
  expect_identical(
    paste(deparse(formulas$climate), collapse = " "),
    paste(deparse(~ T_z * P_z * W_z), collapse = " ")
  )
  expect_identical(
    paste(deparse(formulas$Graminoid_full), collapse = " "),
    paste(deparse(~ GraminoidStatus * T_z * P_z * W_z), collapse = " ")
  )
  expect_identical(
    paste(deparse(formulas$Woodiness_full), collapse = " "),
    paste(deparse(~ Woodiness * T_z * P_z * W_z), collapse = " ")
  )
  phylogeny_terms <- attr(stats::terms(formulas$Phylogeny_full), "term.labels")
  expect_true(all(c("W_z:PEV1_z", "T_z:PEV2_z", "P_z:PEV3_z", "PEV3_z:TP_z") %in%
                    phylogeny_terms))
  expect_false(any(grepl("PEV4|PEV5", phylogeny_terms)))
})

test_that("unique-species predictor scaling is response-independent", {
  species <- tibble::tibble(
    species_id = letters[1:4], PEV1 = c(4, 2, 8, 6), PEV2 = c(-1, 0, 1, 2)
  )
  scaled <- standardize_unique_species_predictors(species, c("PEV1", "PEV2"))
  expect_equal(mean(scaled$data$PEV1_z), 0, tolerance = 1e-15)
  expect_equal(stats::sd(scaled$data$PEV1_z), 1, tolerance = 1e-15)
  expect_false("auc_effect" %in% names(scaled$data))
  expect_identical(scaled$scaling$species_N, c(4L, 4L))
})

test_that("callout species and their selection rules remain explicit", {
  path <- project_path(
    "output", "tables",
    "species_context_classification.csv"
  )
  expect_true(file.exists(path), info = "Build species_classification_file")
  classification <- readr::read_csv(path, show_col_types = FALSE)
  examples <- select_species_callout_examples(classification)
  expect_identical(examples$species_id, callout_species())
  expect_identical(
    examples$selection_rule,
    c(
      "highest fitted minimum", "lowest fitted maximum",
      "widest fitted min-max span", "widest fitted min-max span"
    )
  )
  expect_identical(examples$selection_rank_within_class, c(1L, 2L, 1L, 2L))
})

test_that("heterogeneity reductions retain all response-scale components", {
  baseline <- c(H_overall = 2, H_W = 4, H_context = 5, H_total = 11)
  predictor <- c(H_overall = 1, H_W = 2, H_context = 4, H_total = 7)
  result <- heterogeneity_resolution(baseline, predictor)
  expect_identical(names(result), names(baseline))
  expect_equal(unname(result), c(50, 50, 20, 100 * 4 / 11))
})

test_that("M_both variance profiles retain the five anchored components", {
  specifications <- variance_profile_specifications()
  expect_identical(
    names(specifications),
    c("species_intercept", "species_W", "species_T", "species_P", "species_TP")
  )
  expect_identical(
    vapply(specifications, `[[`, character(1), "argument"),
    c(
      species_intercept = "sigma2", species_W = "tau2",
      species_T = "gamma2", species_P = "gamma2", species_TP = "gamma2"
    )
  )
  expect_identical(
    unname(vapply(specifications, `[[`, integer(1), "index")),
    c(1L, 1L, 1L, 2L, 3L)
  )
})

test_that("the sole family correction is explicit and documented", {
  ledger <- readr::read_csv(
    testthat::test_path("..", "..", "config", "resident_species_overrides.csv"),
    show_col_types = FALSE
  )
  override <- family_overrides(ledger)
  expect_identical(nrow(override), 1L)
  expect_identical(override$species_id, "Viola palustris")
  expect_identical(override$family_id, "Violaceae")
  expect_match(override$decision_basis, "Step 19, Step 20, and Step 22")
  expect_identical(ledger$graminoid_status, "non_graminoid")
  expect_identical(ledger$woodiness, "non_woody")
})

test_that("family sampling-variance comparison requires exact order and close values", {
  candidate <- diag(c(1, 2))
  dimnames(candidate) <- list(c("a", "b"), c("a", "b"))
  expect_true(all(compare_family_sampling_variances(candidate, candidate)$pass))

  reordered <- candidate[c(2, 1), c(2, 1)]
  expect_error(
    compare_family_sampling_variances(candidate, reordered),
    "family_sampling_variance_order"
  )
})

test_that("sampling variances are squared SEs in exact effect order", {
  effects <- tibble::tibble(
    effect_id = c("effect-b", "effect-a"), auc_se = c(0.2, 0.5)
  )
  sampling_variances <- build_sampling_variance_matrix(effects, "auc_se")
  expect_identical(rownames(sampling_variances), effects$effect_id)
  expect_equal(as.numeric(Matrix::diag(sampling_variances)), c(0.04, 0.25))
  off_diagonal <- as.matrix(sampling_variances)
  diag(off_diagonal) <- 0
  expect_equal(unname(off_diagonal), matrix(0, 2, 2))
  expect_true(all(
    validate_sampling_variance_matrix(
      sampling_variances, effects, "auc_se"
    )$valid
  ))
})

test_that("historical statistics derive the direct climate-full comparison", {
  historical <- tibble::tribble(
    ~predictor, ~model, ~method, ~n_effects, ~n_fixed, ~logLik, ~AIC, ~AICc, ~BIC, ~optimizer_code, ~warnings,
    "test", "M0", "ML", 10L, 2L, -10, 24, 26, 25, 0L, "",
    "test", "M1", "ML", 10L, 3L, -9, 24, 28, 26, 0L, "",
    "test", "M2", "ML", 10L, 5L, -6, 22, 37, 27, 0L, ""
  )
  pair_fit <- legacy_predictor_pair_fit_statistics(historical)
  comparison <- derive_legacy_predictor_pair_comparisons(historical)
  expect_identical(
    pair_fit$model, c("climate only", "full predictor-climate")
  )
  expect_equal(comparison$likelihood_ratio, 8)
  expect_equal(comparison$df, 3)
  expect_equal(comparison$p_value, stats::pchisq(8, 3, lower.tail = FALSE))
})

test_that("compact model signatures preserve every fitted regression check", {
  path <- project_path(
    "output", "models",
    "resident_climate_model.rds"
  )
  expect_true(file.exists(path), info = "Build resident_climate_model_file")
  model <- readRDS(path)
  signature <- model_signature(model)
  expect_s3_class(signature, "model_signature")
  expect_identical(
    names(signature),
    c(
      "fixed_coefficients", "fixed_covariance", "sigma2", "tau2",
      "gamma2", "rho", "phi", "log_likelihood"
    )
  )
  expect_true(all(compare_model_fit(model, signature)$pass))
})
