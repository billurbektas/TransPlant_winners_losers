source(project_path("config", "canonical_analysis.R"))
source(project_path("R", "functions", "model_core.R"))
source(project_path("R", "functions", "predictor_models.R"))
source(project_path("R", "functions", "species_context.R"))
source(project_path("R", "functions", "heterogeneity.R"))
source(project_path("R", "functions", "family.R"))

test_that("the canonical formulas preserve the settled Step 22 models", {
  formulas <- canonical_predictor_formulas()
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
    "species_minmax_classification_526.csv"
  )
  expect_true(file.exists(path), info = "Build canonical_species_classification_file")
  classification <- readr::read_csv(path, show_col_types = FALSE)
  examples <- select_species_callout_examples(classification)
  expect_identical(examples$species_id, canonical_callout_species())
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
  specifications <- canonical_variance_profile_specifications()
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
    testthat::test_path("..", "..", "config", "downstream_species_overrides.csv"),
    show_col_types = FALSE
  )
  override <- canonical_family_overrides(ledger)
  expect_identical(nrow(override), 1L)
  expect_identical(override$species_id, "Viola palustris")
  expect_identical(override$family_id, "Violaceae")
  expect_match(override$decision_basis, "Step 19, Step 20, and Step 22")
  expect_identical(ledger$graminoid_status, "non_graminoid")
  expect_identical(ledger$woodiness, "non_woody")
})

test_that("family covariance comparison requires exact order and close values", {
  candidate <- matrix(c(1, 0.2, 0.2, 2), 2, 2)
  dimnames(candidate) <- list(c("a", "b"), c("a", "b"))
  expect_true(all(compare_family_covariance(candidate, candidate)$pass))

  reordered <- candidate[c(2, 1), c(2, 1)]
  expect_error(
    compare_family_covariance(candidate, reordered),
    "family_covariance_order"
  )
})

test_that("compact model signatures preserve every fitted regression check", {
  path <- project_path(
    "output", "models",
    "M_both_526.rds"
  )
  expect_true(file.exists(path), info = "Build canonical_mboth_model_file")
  model <- readRDS(path)
  signature <- canonical_model_signature(model)
  expect_s3_class(signature, "canonical_model_signature")
  expect_identical(
    names(signature),
    c(
      "fixed_coefficients", "fixed_covariance", "sigma2", "tau2",
      "gamma2", "rho", "phi", "log_likelihood"
    )
  )
  expect_true(all(compare_canonical_model_fit(model, signature)$pass))
})
