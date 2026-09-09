test_that("the project dictionary contains only maintained traits", {
  dictionary <- readr::read_csv(project_path("config", "try_trait_map.csv"), show_col_types = FALSE)
  expect_setequal(unique(dictionary$canonical_trait), c("LDMC", "SLA", "vegetative_height"))
  expect_equal(sum(dictionary$canonical_trait == "LDMC"), 2L)
  expect_equal(sum(dictionary$canonical_trait == "SLA"), 7L)
  expect_equal(sum(dictionary$canonical_trait == "vegetative_height"), 2L)
})

test_that("LMA is the reciprocal of positive standardized SLA", {
  observations <- tibble::tibble(
    accepted_name = c("Species one", "Species one", "Species one"),
    TraitName = c("SLA", "LDMC", "vegetative_height"),
    StdValue = c(20, 0.3, 0.5),
    UnitName = c("mm2 mg-1", "g/g", "m"),
    DatasetID = 1:3,
    climate_code = "Cfb"
  )
  result <- summarise_try_traits(observations)$summary
  expect_equal(result$median[result$trait == "LMA"], 0.05)
  observations$StdValue[[1]] <- 0
  expect_error(summarise_try_traits(observations), "SLA contains")
})
