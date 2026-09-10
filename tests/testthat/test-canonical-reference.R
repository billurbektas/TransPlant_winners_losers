test_that("the active canonical reference is frozen and matches its target", {
  path <- project_path(
    "output", "reference", "canonical_reference_results.rds"
  )
  expect_true(file.exists(path))
  expect_identical(
    sha256_file(path),
    "75f5659904616f3a1d8c2046f356ddbd363506bd21d55398af7b73995eb94792"
  )

  frozen <- readRDS(path)
  expect_invisible(validate_step22_schema(frozen, "canonical reference"))

  if (targets::tar_exist_objects("canonical_reference_results")) {
    active <- targets::tar_read(canonical_reference_results)
    expect_true(all(validate_step22_reference(active, frozen)$pass))
  }
})

test_that("canonical reference writing is fail-closed", {
  reference <- readRDS(project_path(
    "output", "reference", "canonical_reference_results.rds"
  ))
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)

  expect_identical(write_canonical_reference_once(reference, path), path)
  expect_true(all(validate_step22_reference(reference, readRDS(path))$pass))

  changed <- reference
  changed$M_both$fixed$estimate[[1]] <-
    changed$M_both$fixed$estimate[[1]] + 1e-3
  expect_error(
    write_canonical_reference_once(changed, path),
    "will not be overwritten"
  )
})
