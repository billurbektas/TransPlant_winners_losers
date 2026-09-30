source(project_path("R", "validation", "legacy_reference.R"))
source(project_path("R", "infrastructure", "provenance.R"))

test_that("the frozen historical fixture has the required schema and scope", {
  reference <- readRDS(project_path("tests", "fixtures", "step22_reference_results.rds"))
  expect_invisible(validate_legacy_reference_schema(reference))
  expect_identical(reference$data_validation$effects, 1604L)
  expect_identical(reference$data_validation$species, 526L)
  expect_identical(reference$data_validation$experiments, 42L)
  expect_identical(
    reference$metadata$excluded_downstream_species,
    resident_model_exclusions()
  )
  expect_true(all(reference$data_validation$consistency_checks$pass))
})

test_that("the Step 22 data-validation adapter preserves integer check IDs", {
  source(project_path("R", "validation", "legacy_artifacts.R"))
  scope <- list(
    effects = 1L,
    species = 1L,
    experiments = 1L,
    sampling_variance_validation = tibble::tibble(
      dimensions_valid = TRUE,
      order_valid = TRUE,
      maximum_absolute_off_diagonal = 0,
      minimum_sampling_variance = 0.1,
      maximum_diagonal_variance_error = 0,
      valid = TRUE
    ),
    consistency_checks = tibble::tibble(
      check_number = c(1, 2), check = c("a", "b"), pass = TRUE
    )
  )
  adapted <- legacy_data_validation_component(scope)
  expect_type(adapted$consistency_checks$check_number, "integer")
})

test_that("the resident sampling-variance check exposes the release interface", {
  source(project_path("R", "validation", "legacy_artifacts.R"))
  validation <- tibble::tibble(
    dimensions_valid = TRUE,
    order_valid = TRUE,
    maximum_absolute_off_diagonal = 0,
    minimum_sampling_variance = 0.1,
    maximum_diagonal_variance_error = 0,
    valid = TRUE
  )
  gate <- resident_sampling_variance_release_check(validation)
  expect_identical(gate$pass, TRUE)
  expect_match(gate$detail, "minimum variance")
})

test_that("the historical comparator translates retired predictor rows and rejects drift", {
  reference <- readRDS(project_path("tests", "fixtures", "step22_reference_results.rds"))
  candidate <- translate_legacy_reference_for_active_analysis(reference)
  expect_true(all(validate_legacy_reference(candidate, reference)$pass))

  changed <- candidate
  changed$M_both$fixed$estimate[[1]] <-
    changed$M_both$fixed$estimate[[1]] + 1e-4
  validation <- validate_legacy_reference(changed, reference)
  expect_false(validation$pass[[2]])
  expect_match(validation$detail[[2]], "M_both/fixed/estimate")

  changed_count <- candidate
  changed_count$data_validation$effects <-
    changed_count$data_validation$effects + 1L
  count_validation <- validate_legacy_reference(changed_count, reference)
  expect_false(count_validation$pass[[2]])
  expect_match(count_validation$detail[[2]], "data_validation/effects")
})

test_that("empty warning and error columns compare semantically", {
  candidate <- tibble::tibble(
    component = c("a", "b"),
    contrast_warnings = c(NA_character_, ""),
    error = c(NA_character_, NA_character_)
  )
  reference <- tibble::tibble(
    component = c("a", "b"),
    contrast_warnings = c(NA, NA),
    error = c(NA, NA)
  )
  expect_length(compare_reference_node(candidate, reference), 0L)

  candidate$contrast_warnings[[2]] <- "optimizer warning"
  expect_match(
    compare_reference_node(candidate, reference),
    "exact diagnostic values differ"
  )
})

test_that("the tolerance registry matches the frozen release policy", {
  expect_identical(
    reference_tolerance(c("scaling", "climate", "center")),
    c(absolute = 1e-12, relative = 1e-10)
  )
  expect_identical(
    reference_tolerance(c("data_validation", "V_validation", "minimum_eigenvalue")),
    c(absolute = 1e-12, relative = 1e-10)
  )
  expect_identical(
    reference_tolerance(c("M_both", "fixed", "estimate")),
    c(absolute = 1e-8, relative = 1e-6)
  )
  expect_identical(
    reference_tolerance(c("predictor_CR1", "model_p_value")),
    c(absolute = 1e-5, relative = 1e-6)
  )
  expect_identical(
    reference_tolerance(c("heterogeneity_components", "H_total")),
    c(absolute = 2e-7, relative = 1e-6)
  )
  expect_identical(
    reference_tolerance(c(
      "heterogeneity_synthesis", "Total heterogeneity resolved (%)"
    )),
    c(absolute = 2e-3, relative = 1e-6)
  )
  expect_identical(
    reference_tolerance(c("M_both", "variance_profiles", "CI_low")),
    c(absolute = 1e-6, relative = 1e-4)
  )
})

test_that("the immutable Step 22 fixture retains its frozen checksum", {
  fixture <- project_path("tests", "fixtures", "step22_reference_results.rds")
  expected <- "8d5d7544a6e6bcb757827263763ffae95115fa10b52b187669eaa5627af6931f"
  expect_identical(sha256_file(fixture), expected)
})

test_that("release validation fails closed and reports named gates", {
  checks <- list(
    first = tibble::tibble(check = c("a", "b"), pass = c(TRUE, TRUE)),
    second = tibble::tibble(check = "c", pass = TRUE)
  )
  validation <- assemble_release_validation(checks)
  expect_identical(validation$gate, c("first", "second"))
  expect_true(all(assert_release_validation(validation)$pass))

  checks$second$pass <- FALSE
  failed <- assemble_release_validation(checks)
  expect_error(assert_release_validation(failed), "second")
})

test_that("run provenance is separate, complete, and release-gated", {
  gate <- tibble::tibble(gate = "numerical", pass = TRUE)
  provenance <- build_run_provenance(project_path(), gate)
  required <- c(
    "recorded_at_utc", "git_branch", "git_commit", "git_worktree_dirty",
    "r_version", "r_platform", "targets_version", "tarchetypes_version",
    "quarto_r_package_version", "v_phylomaker2_version",
    "renv_lock_sha256", "targets_script_sha256"
  )
  expect_true(all(required %in% provenance$field))
  expect_false(any(grepl("reference", provenance$field, fixed = TRUE)))
  expect_error(
    build_run_provenance(
      project_path(), tibble::tibble(gate = "numerical", pass = FALSE)
    ),
    "numerical"
  )
})
