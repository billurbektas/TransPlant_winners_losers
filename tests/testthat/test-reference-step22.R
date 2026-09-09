source(project_path("R", "functions", "reference_validation.R"))
source(project_path("R", "functions", "provenance.R"))

test_that("the frozen Step 22 fixture has the required schema and scope", {
  reference <- readRDS(project_path("tests", "fixtures", "step22_reference_results.rds"))
  expect_invisible(validate_step22_schema(reference))
  expect_identical(reference$data_validation$effects, 1604L)
  expect_identical(reference$data_validation$species, 526L)
  expect_identical(reference$data_validation$experiments, 42L)
  expect_identical(
    reference$metadata$excluded_downstream_species,
    resident_downstream_model_exclusions()
  )
  expect_true(all(reference$data_validation$consistency_checks$pass))
})

test_that("the Step 22 comparator passes identical objects and rejects drift", {
  reference <- readRDS(project_path("tests", "fixtures", "step22_reference_results.rds"))
  expect_true(all(validate_step22_reference(reference, reference)$pass))

  changed <- reference
  changed$M_both$fixed$estimate[[1]] <-
    changed$M_both$fixed$estimate[[1]] + 1e-4
  validation <- validate_step22_reference(changed, reference)
  expect_false(validation$pass[[2]])
  expect_match(validation$detail[[2]], "M_both/fixed/estimate")

  changed_count <- reference
  changed_count$data_validation$effects <-
    changed_count$data_validation$effects + 1L
  count_validation <- validate_step22_reference(changed_count, reference)
  expect_false(count_validation$pass[[2]])
  expect_match(count_validation$detail[[2]], "data_validation/effects")
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
    reference_tolerance(c("M_both", "variance_profiles", "CI_low")),
    c(absolute = 1e-6, relative = 1e-4)
  )
})

test_that("the immutable project reference and fixture are byte-identical", {
  project_reference <- project_path(
    "output", "step_by_step_modeling", "downstream_526",
    "downstream_526_reference_results.rds"
  )
  fixture <- project_path("tests", "fixtures", "step22_reference_results.rds")
  skip_if_not(file.exists(project_reference))
  expected <- "8d5d7544a6e6bcb757827263763ffae95115fa10b52b187669eaa5627af6931f"
  expect_identical(sha256_file(project_reference), expected)
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
  provenance <- build_canonical_run_provenance(project_path(), gate)
  required <- c(
    "recorded_at_utc", "git_branch", "git_commit", "git_worktree_dirty",
    "r_version", "r_platform", "targets_version", "tarchetypes_version",
    "quarto_r_package_version", "v_phylomaker2_version",
    "renv_lock_sha256", "targets_script_sha256",
    "step22_reference_sha256", "transition_reference_sha256"
  )
  expect_true(all(required %in% provenance$field))
  values <- stats::setNames(provenance$value, provenance$field)
  expect_identical(
    values[["step22_reference_sha256"]],
    "8d5d7544a6e6bcb757827263763ffae95115fa10b52b187669eaa5627af6931f"
  )
  expect_identical(
    values[["transition_reference_sha256"]],
    "6aa52701c604bf23f7366970fbf80e4e3b05f7253ae1961e1f7251ff08c962e6"
  )
  expect_error(
    build_canonical_run_provenance(
      project_path(), tibble::tibble(gate = "numerical", pass = FALSE)
    ),
    "numerical"
  )
})
