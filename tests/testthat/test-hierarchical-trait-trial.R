test_that("hierarchical LDMC trial formulas preserve the requested structure", {
  expect_identical(
    attr(stats::terms(hierarchical_ldmc_fixed_formula()), "term.labels"),
    c(
      "temp_z", "precip_z", "warming_z", "temp_z:LDMC_z",
      "precip_z:LDMC_z", "warming_z:LDMC_z"
    )
  )
  random_text <- vapply(
    hierarchical_ldmc_random_formulas(),
    function(x) paste(deparse(x), collapse = " "), character(1)
  )
  expect_equal(
    random_text,
    c(
      "~1 | species", "~0 + temp_z | species",
      "~0 + precip_z | species", "~0 + warming_z | species",
      "~1 | experiment", "~1 | gradient"
    )
  )
  binomial_text <- paste(
    deparse(hierarchical_ldmc_binomial_formula()), collapse = " "
  )
  expect_match(binomial_text, "\\(1 \\| species\\)")
  expect_match(binomial_text, "\\(1 \\| experiment\\)")
  expect_match(binomial_text, "\\(1 \\| gradient\\)")
})

test_that("hierarchical LDMC preparation uses pools, references, and regions", {
  path <- project_path("output", "tables", "analysis_dataset.csv")
  skip_if_not(file.exists(path), "analysis dataset is not available")
  data <- readr::read_csv(path, show_col_types = FALSE)
  prepared <- prepare_hierarchical_ldmc_trial(data)

  expect_identical(
    levels(prepared$resident$response),
    c("non_significant", "loser", "winner")
  )
  expect_identical(
    levels(prepared$nonresident$response),
    c("non_significant", "colonizer")
  )
  expect_equal(nrow(prepared$experiment_gradient_map), 42L)
  expect_equal(dplyr::n_distinct(prepared$experiment_gradient_map$gradient), 20L)
  expect_true(all(prepared$trait_scaling$weighting == "distinct species"))
  expect_true(all(is.finite(prepared$resident$LDMC_z)))
  expect_true(all(is.finite(prepared$nonresident$LDMC_z)))
  expect_false(anyDuplicated(paste(
    prepared$resident$species, prepared$resident$experiment
  )) > 0L)
  expect_false(anyDuplicated(paste(
    prepared$nonresident$species, prepared$nonresident$experiment
  )) > 0L)
})
