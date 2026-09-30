test_that("the frozen release reference is immutable and matches its target", {
  path <- project_path(
    "output", "reference", "canonical_reference_results.rds"
  )
  expect_true(file.exists(path))
  expect_identical(
    sha256_file(path),
    "75f5659904616f3a1d8c2046f356ddbd363506bd21d55398af7b73995eb94792"
  )

  frozen <- readRDS(path)
  expect_invisible(validate_legacy_reference_schema(frozen, "canonical reference"))

  if (targets::tar_exist_objects("legacy_comparison_candidate")) {
    active <- targets::tar_read(legacy_comparison_candidate)
    expect_true(all(validate_legacy_reference(active, frozen)$pass))
  }
})

test_that("canonical reference writing is fail-closed", {
  reference <- readRDS(project_path(
    "output", "reference", "canonical_reference_results.rds"
  ))
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)

  expect_identical(write_release_reference_once(reference, path), path)
  expect_true(all(validate_legacy_reference(reference, readRDS(path))$pass))

  changed <- reference
  changed$M_both$fixed$estimate[[1]] <-
    changed$M_both$fixed$estimate[[1]] + 1e-3
  expect_error(
    write_release_reference_once(changed, path),
    "will not be overwritten"
  )
})
