test_that("the project dictionary contains only maintained traits", {
  dictionary <- readr::read_csv(project_path("config", "try_trait_map.csv"), show_col_types = FALSE)
  expect_setequal(unique(dictionary$analysis_trait), c("LDMC", "SLA", "vegetative_height"))
  expect_equal(sum(dictionary$analysis_trait == "LDMC"), 2L)
  expect_equal(sum(dictionary$analysis_trait == "SLA"), 7L)
  expect_equal(sum(dictionary$analysis_trait == "vegetative_height"), 2L)
})

test_that("TRY climate classification uses the pinned bundled raster", {
  raster <- system.file("extdata", "koppen_geiger_0p1.tif", package = "tidyTRY")
  manifest <- readr::read_csv(project_path("data", "source_manifest.csv"),
                              show_col_types = FALSE)
  row <- manifest[manifest$path == "tidyTRY::extdata/koppen_geiger_0p1.tif", ]
  expect_equal(nrow(row), 1L)
  expect_true(nzchar(raster))
  expect_identical(sha256_file(raster), row$sha256[[1]])
  expect_identical(transplant_config()$try_raster_sha256, row$sha256[[1]])
  expect_error(
    filter_try_climate(list(), NULL, character(), "wrong checksum"),
    "does not match"
  )
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
