test_that("binomial reduction grid removes experiment and gradient effects", {
  plan <- hierarchical_ldmc_reduction_plan()
  expect_equal(nrow(plan), 9L)
  expect_false(any(plan$includes_experiment_random))
  expect_false(any(plan$includes_gradient_random))
  expect_true(all(grepl("species", plan$formula, fixed = TRUE)))
  expect_setequal(
    plan$specification,
    c(
      "full_species_all_slopes", "full_species_intercept",
      "drop_temp_interaction", "drop_precip_interaction",
      "drop_warming_interaction", "climate_main_only",
      "temperature_axis", "precipitation_axis", "warming_axis"
    )
  )
})

test_that("fixed-effect reductions retain required main effects", {
  specs <- hierarchical_ldmc_binomial_reduction_specs()
  purrr::walk(specs, function(spec) {
    for (interaction in spec$fixed_terms[grepl(":", spec$fixed_terms)]) {
      climate <- sub(":LDMC_z$", "", interaction)
      expect_true(climate %in% spec$fixed_terms)
    }
  })
})
