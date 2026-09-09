test_that("monthly annual, JJA, and July summaries use the requested operations", {
  months <- 1:12
  starts <- as.Date(sprintf("2012-%02d-01", months))
  month_days <- c(31L, 29L, 31L, 30L, 31L, 30L, 31L, 31L, 30L, 31L, 30L, 31L)
  fixture <- tidyr::crossing(site_key = "R::S", variable = c("tas", "pr"), month = months) |>
    dplyr::mutate(
      region = "R", site_id = "S",
      year = 2012L,
      period = starts[.data$month],
      days_in_month = month_days[.data$month],
      value = dplyr::if_else(variable == "tas", as.numeric(month), 1)
    )
  result <- summarise_site_climate(fixture)
  expect_equal(result$P_annual, 12)
  expect_equal(result$P_JJA, 3)
  expect_equal(result$P_July, 1)
  expect_equal(result$T_July, 7)
  expect_equal(result$T_JJA, weighted.mean(6:8, c(30, 31, 31)))
})

test_that("a valid monthly checkpoint is reused without network access", {
  cache <- tempfile("chelsa-cache-")
  sites <- tibble::tibble(region = "R", site_id = "S", site_key = "R::S", longitude = 8, latitude = 47)
  path <- chelsa_checkpoint_path(cache, "tas", 2012, 7)
  checkpoint <- tibble::tibble(
    site_key = "R::S", region = "R", site_id = "S", longitude = 8, latitude = 47,
    variable = "tas", year = 2012L, month = 7L, period = as.Date("2012-07-01"),
    days_in_month = 31L, unit = "degC", value = 12, url = "fixture"
  )
  write_rds_atomic(checkpoint, path)
  expect_true(valid_chelsa_checkpoint(path, sites, "tas", 2012, 7))
  expect_identical(checkpoint_chelsa_month(sites, "tas", 2012, 7, cache), path)
})

test_that("CHELSA URLs and coordinate gates are deterministic", {
  expect_match(chelsa_monthly_url("tas", 2012, 7), "CHELSA_tas_07_2012_V.2.1.tif$")
  expect_error(validate_site_coordinates(tibble::tibble(region = "R", site_id = "S", longitude = 181, latitude = 47)))
})

test_that("coordinate audit reports shared locations and rejects identical experiment endpoints", {
  sites <- tibble::tribble(
    ~region, ~site_id, ~longitude, ~latitude,
    "R1", "O", 8, 47,
    "R1", "D", 9, 46,
    "R2", "same_D", 9, 46
  )
  experiments <- tibble::tribble(
    ~experiment_id, ~region, ~origin_site, ~destination_site, ~included,
    "e1", "R1", "O", "D", TRUE,
    "e2", "R2", "same_D", "same_D", FALSE
  )
  audit <- audit_chelsa_coordinates(sites, experiments)
  expect_equal(nrow(audit$duplicate_coordinates), 1L)
  expect_equal(nrow(audit$identical_experiment_pairs), 0L)
  bad <- dplyr::bind_rows(
    experiments,
    tibble::tibble(
      experiment_id = "e3", region = "R2", origin_site = "same_D",
      destination_site = "same_D", included = TRUE
    )
  )
  expect_error(audit_chelsa_coordinates(sites, bad), "identical origin and destination")
})

test_that("experiment climate flags and calculates partial available-year averages", {
  site_year <- tidyr::crossing(region = "R", site_id = c("O", "D"), year = 2012:2013) |>
    dplyr::mutate(
      site_key = paste(region, site_id, sep = "::"),
      T_annual = dplyr::if_else(site_id == "O", 5, 8), P_annual = 100,
      T_JJA = dplyr::if_else(site_id == "O", 10, 14), P_JJA = dplyr::if_else(site_id == "O", 50, 60),
      T_July = dplyr::if_else(site_id == "O", 12, 16), P_July = 20
    )
  experiment <- tibble::tibble(
    region = "R", origin_site = "O", destination_site = "D", experiment_id = "R::O::D",
    first_year = 2012L, last_year = 2013L, included = TRUE
  )
  result <- summarise_experiment_climate(site_year, experiment)
  expect_equal(result$diff_T_JJA, 4)
  expect_equal(result$destination_P_JJA, 60)
  expect_true(result$climate_complete)
  incomplete <- summarise_experiment_climate(
    dplyr::filter(site_year, !(site_id == "D" & year == 2013)), experiment
  )
  expect_false(incomplete$climate_complete)
  expect_equal(incomplete$climate_summary_status, "partial")
  expect_equal(incomplete$destination_missing_years, "2013")
  expect_equal(incomplete$destination_T_JJA, 14)
  expect_match(incomplete$climate_coverage_note, "1/2 years")
})

test_that("climate joins preserve existing experiment metadata without suffixes", {
  base <- tibble::tibble(experiment_id = "e1", year_range = 4L, species = "sp")
  climate <- tibble::tibble(
    experiment_id = "e1", year_range = 99L,
    destination_T_JJA = 12, climate_complete = TRUE
  )
  joined <- join_experiment_climate(base, climate)
  expect_equal(joined$year_range, 4L)
  expect_equal(joined$destination_T_JJA, 12)
  expect_false(any(grepl("[.][xy]$", names(joined))))
})
