# Module: R/data/climate
# Responsibility: Validate sites and derive CHELSA site and experiment climate.
# Inputs: Site coordinates, monthly checkpoints or CHELSA service, and experiment dates.
# Outputs: Validated monthly, annual, summer, and experiment climate tables.
# Side effects: Online mode may download CHELSA months and write cache or missing-month reports; offline mode only reads.

#' Validate Site Coordinates.
#'
#' @param sites A data frame whose required columns are checked before use.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
validate_site_coordinates <- function(sites) {
  assert_unique(sites, c("region", "site_id"), "CHELSA sites")
  bad <- sites |>
    dplyr::filter(
      is.na(.data$longitude) | is.na(.data$latitude) |
        .data$longitude < -180 | .data$longitude > 180 |
        .data$latitude < -90 | .data$latitude > 90
    )
  if (nrow(bad)) stop("Invalid or missing site coordinates for ", nrow(bad), " sites.", call. = FALSE)
  sites |>
    dplyr::mutate(site_key = paste(.data$region, .data$site_id, sep = "::"))
}

#' Audit Chelsa Coordinates.
#'
#' @param sites A data frame whose required columns are checked before use.
#' @param experiments A data frame whose required columns are checked before use.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
audit_chelsa_coordinates <- function(sites, experiments) {
  sites <- validate_site_coordinates(sites)
  included <- experiments |>
    dplyr::filter(.data$included)
  duplicate_coordinates <- sites |>
    dplyr::summarise(
      n_site_ids = dplyr::n_distinct(.data$site_key),
      site_keys = paste(sort(unique(.data$site_key)), collapse = "; "),
      .by = c("longitude", "latitude")
    ) |>
    dplyr::filter(.data$n_site_ids > 1L) |>
    dplyr::arrange(.data$longitude, .data$latitude)
  usage <- dplyr::bind_rows(
    included |>
      dplyr::transmute(experiment_id, region, role = "origin", site_id = .data$origin_site),
    included |>
      dplyr::transmute(experiment_id, region, role = "destination", site_id = .data$destination_site)
  ) |>
    dplyr::left_join(
      sites |>
        dplyr::select(region, site_id, site_key, longitude, latitude),
      by = c("region", "site_id")
    )
  reused_coordinates <- usage |>
    dplyr::summarise(
      n_experiments = dplyr::n_distinct(.data$experiment_id),
      n_site_ids = dplyr::n_distinct(.data$site_key),
      site_keys = paste(sort(unique(.data$site_key)), collapse = "; "),
      roles = paste(sort(unique(.data$role)), collapse = "; "),
      experiment_ids = paste(sort(unique(.data$experiment_id)), collapse = "; "),
      .by = c("longitude", "latitude")
    ) |>
    dplyr::filter(.data$n_experiments > 1L) |>
    dplyr::arrange(dplyr::desc(.data$n_experiments), .data$longitude, .data$latitude)
  paired <- included |>
    dplyr::left_join(
      sites |>
        dplyr::select(region, origin_site = site_id,
                      origin_longitude = longitude, origin_latitude = latitude),
      by = c("region", "origin_site")
    ) |>
    dplyr::left_join(
      sites |>
        dplyr::select(region, destination_site = site_id,
                      destination_longitude = longitude, destination_latitude = latitude),
      by = c("region", "destination_site")
    ) |>
    dplyr::filter(
      .data$origin_longitude == .data$destination_longitude,
      .data$origin_latitude == .data$destination_latitude
    ) |>
    dplyr::select(experiment_id, region, origin_site, destination_site,
                  origin_longitude, origin_latitude)
  if (nrow(paired)) {
    stop("At least one warmed experiment has identical origin and destination coordinates.", call. = FALSE)
  }
  summary <- tibble::tibble(
    check = c(
      "unique site keys", "coordinate bounds and completeness",
      "unique coordinate locations", "shared coordinate groups",
      "coordinates reused by multiple experiments",
      "identical origin/destination coordinates"
    ),
    value = c(
      nrow(sites), nrow(sites),
      dplyr::n_distinct(paste(sites$longitude, sites$latitude)),
      nrow(duplicate_coordinates), nrow(reused_coordinates), nrow(paired)
    ),
    pass = c(TRUE, TRUE, TRUE, TRUE, TRUE, nrow(paired) == 0L)
  )
  list(
    summary = summary,
    duplicate_coordinates = duplicate_coordinates,
    reused_coordinates = reused_coordinates,
    identical_experiment_pairs = paired
  )
}

#' Chelsa Monthly Url.
#'
#' @param variable A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param year A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param month A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
chelsa_monthly_url <- function(variable, year, month) {
  sprintf(
    "https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/monthly/%s/%04d/CHELSA_%s_%02d_%04d_V.2.1.tif",
    variable, as.integer(year), variable, as.integer(month), as.integer(year)
  )
}

#' Chelsa Jobs.
#'
#' @param years A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param variables A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
chelsa_jobs <- function(years, variables = c("tas", "pr")) {
  jobs <- tidyr::crossing(year = as.integer(years), month = 1:12, variable = variables) |>
    dplyr::mutate(job_id = sprintf("%04d-%02d-%s", .data$year, .data$month, .data$variable))
  split(jobs, seq_len(nrow(jobs)))
}

#' Chelsa Checkpoint Path.
#'
#' @param cache_dir A project-relative or absolute filesystem path.
#' @param variable A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param year A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param month A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
chelsa_checkpoint_path <- function(cache_dir, variable, year, month) {
  file.path(
    cache_dir, sprintf("%04d", as.integer(year)), sprintf("%02d", as.integer(month)),
    sprintf("%s_%04d-%02d.rds", variable, as.integer(year), as.integer(month))
  )
}

#' Check that a saved CHELSA month matches the current sites and coordinates.
#'
#' @param path A project-relative or absolute filesystem path.
#' @param sites A data frame whose required columns are checked before use.
#' @param variable A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param year A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param month A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return TRUE only for a readable month with one finite value per current
#'   site, matching site keys and coordinates; otherwise FALSE.
valid_chelsa_checkpoint <- function(path, sites, variable, year, month) {
  if (!file.exists(path)) return(FALSE)
  value <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(value) || !is.data.frame(value)) return(FALSE)
  needed <- c(
    "site_key", "longitude", "latitude", "variable", "year", "month",
    "days_in_month", "unit", "value"
  )
  if (!all(needed %in% names(value)) || nrow(value) != nrow(sites)) return(FALSE)
  if (!identical(sort(value$site_key), sort(sites$site_key))) return(FALSE)
  site_order <- match(sites$site_key, value$site_key)
  if (anyNA(site_order)) return(FALSE)
  same_coordinates <- isTRUE(all.equal(
    as.numeric(value$longitude[site_order]), as.numeric(sites$longitude),
    tolerance = 0
  )) && isTRUE(all.equal(
    as.numeric(value$latitude[site_order]), as.numeric(sites$latitude),
    tolerance = 0
  ))
  isTRUE(same_coordinates &&
    all(value$variable == variable) &&
    all(value$year == as.integer(year)) &&
    all(value$month == as.integer(month)) &&
    all(value$days_in_month %in% 28:31) &&
    dplyr::n_distinct(value$site_key) == nrow(sites) &&
    all(is.finite(value$value)))
}

#' Extract Chelsa Month.
#'
#' @param sites A data frame whose required columns are checked before use.
#' @param variable A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param year A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param month A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
extract_chelsa_month <- function(sites, variable, year, month) {
  url <- chelsa_monthly_url(variable, year, month)
  raster <- terra::rast(paste0("/vsicurl/", url))
  values <- as.numeric(terra::extract(raster, sites[, c("longitude", "latitude")])[, 2])
  if (variable == "tas") values <- values - 273.15
  start <- as.Date(sprintf("%04d-%02d-01", as.integer(year), as.integer(month)))
  next_month <- seq(start, by = "month", length.out = 2L)[[2]]
  tibble::tibble(
    site_key = sites$site_key,
    region = sites$region,
    site_id = sites$site_id,
    longitude = sites$longitude,
    latitude = sites$latitude,
    variable = variable,
    year = as.integer(year),
    month = as.integer(month),
    period = start,
    days_in_month = as.integer(next_month - start),
    unit = dplyr::if_else(variable == "tas", "degC", "kg m-2 month-1"),
    value = values,
    url = url
  )
}

#' Reuse or download one CHELSA monthly checkpoint.
#'
#' @param sites A data frame whose required columns are checked before use.
#' @param variable A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param year A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param month A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param cache_dir A project-relative or absolute filesystem path.
#' @param retries A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param backoff_seconds A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param offline If TRUE, reject a missing or invalid checkpoint without
#'   modifying the cache or contacting CHELSA.
#' @return Path to a verified RDS checkpoint, or NA with download-error
#'   attributes after bounded retries in online mode.
checkpoint_chelsa_month <- function(sites, variable, year, month, cache_dir, retries = 4L,
                                    backoff_seconds = 2, offline = FALSE) {
  path <- chelsa_checkpoint_path(cache_dir, variable, year, month)
  if (valid_chelsa_checkpoint(path, sites, variable, year, month)) return(path)
  if (isTRUE(offline)) {
    stop(
      "Offline CHELSA checkpoint missing or invalid: ",
      sprintf("%04d-%02d-%s", as.integer(year), as.integer(month), variable),
      " (", path, ")", call. = FALSE
    )
  }
  if (file.exists(path)) {
    quarantine <- paste0(path, ".invalid-", format(Sys.time(), "%Y%m%d%H%M%S"))
    file.rename(path, quarantine)
  }
  last_error <- NULL
  for (attempt in seq_len(retries)) {
    result <- tryCatch(
      extract_chelsa_month(sites, variable, year, month),
      error = function(e) {
        last_error <<- conditionMessage(e)
        NULL
      }
    )
    if (!is.null(result) && nrow(result) == nrow(sites) && all(is.finite(result$value))) {
      write_rds_atomic(result, path)
      if (valid_chelsa_checkpoint(path, sites, variable, year, month)) return(path)
      last_error <- "checkpoint failed validation after writing"
    }
    if (attempt < retries) Sys.sleep(backoff_seconds * 2^(attempt - 1L))
  }
  structure(NA_character_, url = chelsa_monthly_url(variable, year, month), error = last_error)
}

#' Verify or retrieve one CHELSA month for a pipeline branch.
#'
#' @param job A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param sites A data frame whose required columns are checked before use.
#' @param cache_dir A project-relative or absolute filesystem path.
#' @param retries A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param backoff_seconds A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param missing_report_dir A project-relative or absolute filesystem path.
#' @param offline Reuse only validated monthly cache files when TRUE; fail with
#'   the affected year-month-variable and do not download.
#' @return Path to the validated monthly RDS file. Online download failure
#'   writes a missing-month report and stops; offline failure writes nothing.
run_chelsa_job <- function(job, sites, cache_dir, retries = 4L, backoff_seconds = 2,
                           missing_report_dir = file.path(dirname(cache_dir), "missing_monthly"),
                           offline = FALSE) {
  sites <- validate_site_coordinates(sites)
  year <- job$year[[1]]
  month <- job$month[[1]]
  variable <- job$variable[[1]]
  report_path <- file.path(missing_report_dir, sprintf("%04d-%02d-%s.csv", year, month, variable))
  checkpoint <- checkpoint_chelsa_month(
    sites, variable, year, month, cache_dir, retries, backoff_seconds,
    offline = offline
  )
  if (is.na(checkpoint[[1]])) {
    report <- tibble::tibble(
      variable = variable,
      year = year,
      month = month,
      url = attr(checkpoint, "url"),
      error = dplyr::coalesce(attr(checkpoint, "error"), "unknown error")
    )
    write_csv_atomic(report, report_path)
    stop("CHELSA branch failed after bounded retries; see ", report_path, call. = FALSE)
  }
  if (!isTRUE(offline) && file.exists(report_path)) unlink(report_path)
  checkpoint
}

#' Consolidate Chelsa.
#'
#' @param monthly_files A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param sites A data frame whose required columns are checked before use.
#' @param years A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
consolidate_chelsa <- function(monthly_files, sites, years) {
  monthly <- purrr::map_dfr(monthly_files, readRDS) |>
    dplyr::mutate(period = as.Date(.data$period))
  expected_months <- length(years) * 12L
  coverage <- monthly |>
    dplyr::count(.data$site_key, .data$variable, name = "n_months") |>
    dplyr::mutate(expected_months = expected_months, pass = .data$n_months == .data$expected_months)
  if (any(!coverage$pass)) {
    stop("CHELSA monthly coverage is incomplete; annual summaries were not calculated.", call. = FALSE)
  }
  list(monthly = monthly, coverage = coverage)
}

#' Summarise Site Climate.
#'
#' @param monthly A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
summarise_site_climate <- function(monthly) {
  annual_coverage <- monthly |>
    dplyr::count(.data$site_key, .data$variable, .data$year, name = "n_months")
  if (any(annual_coverage$n_months != 12L)) {
    stop("At least one CHELSA site-variable-year does not contain 12 monthly values.", call. = FALSE)
  }
  temperature <- monthly |>
    dplyr::filter(.data$variable == "tas") |>
    dplyr::group_by(.data$site_key, .data$region, .data$site_id, .data$year) |>
    dplyr::summarise(
      T_annual = stats::weighted.mean(.data$value, .data$days_in_month),
      T_JJA = stats::weighted.mean(
        .data$value[.data$month %in% 6:8], .data$days_in_month[.data$month %in% 6:8]
      ),
      T_July = dplyr::first(.data$value[.data$month == 7]),
      n_T_annual = dplyr::n(),
      n_T_JJA = sum(.data$month %in% 6:8),
      n_T_July = sum(.data$month == 7),
      .groups = "drop"
    )
  precipitation <- monthly |>
    dplyr::filter(.data$variable == "pr") |>
    dplyr::group_by(.data$site_key, .data$region, .data$site_id, .data$year) |>
    dplyr::summarise(
      P_annual = sum(.data$value),
      P_JJA = sum(.data$value[.data$month %in% 6:8]),
      P_July = dplyr::first(.data$value[.data$month == 7]),
      n_P_annual = dplyr::n(),
      n_P_JJA = sum(.data$month %in% 6:8),
      n_P_July = sum(.data$month == 7),
      .groups = "drop"
    )
  dplyr::inner_join(temperature, precipitation, by = c("site_key", "region", "site_id", "year"))
}

#' Summarise Experiment Climate.
#'
#' @param site_year A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param experiments A data frame whose required columns are checked before use.
#' @return An R object described by the function name; its structure is validated before later stages use it.
summarise_experiment_climate <- function(site_year, experiments) {
  included <- experiments |>
    dplyr::filter(.data$included) |>
    dplyr::rowwise() |>
    dplyr::mutate(year = list(seq.int(.data$first_year, .data$last_year))) |>
    tidyr::unnest(year) |>
    dplyr::ungroup()
  metrics <- c("T_annual", "P_annual", "T_JJA", "P_JJA", "T_July", "P_July")
  summarize_side <- function(site_column, prefix) {
    result <- included |>
      dplyr::left_join(site_year, by = c("region", setNames("site_id", site_column), "year")) |>
      dplyr::group_by(.data$experiment_id) |>
      dplyr::summarise(
        dplyr::across(
          dplyr::all_of(metrics),
          ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE),
                      .names = paste0(prefix, "_{.col}")),
        years_expected = dplyr::n_distinct(.data$year),
        years_available = sum(stats::complete.cases(dplyr::pick(dplyr::all_of(metrics)))),
        missing_years = {
          missing <- .data$year[
            !stats::complete.cases(dplyr::pick(dplyr::all_of(metrics)))
          ]
          if (length(missing)) paste(missing, collapse = ";") else "none"
        },
        .groups = "drop"
      )
    names(result)[names(result) == "years_expected"] <- paste0(prefix, "_years_expected")
    names(result)[names(result) == "years_available"] <- paste0(prefix, "_years_available")
    names(result)[names(result) == "missing_years"] <- paste0(prefix, "_missing_years")
    result
  }
  origin <- summarize_side("origin_site", "origin")
  destination <- summarize_side("destination_site", "destination")
  out <- experiments |>
    dplyr::filter(.data$included) |>
    dplyr::left_join(origin, by = "experiment_id") |>
    dplyr::left_join(destination, by = "experiment_id")
  for (metric in metrics) {
    out[[paste0("diff_", metric)]] <- out[[paste0("destination_", metric)]] - out[[paste0("origin_", metric)]]
  }
  out <- out |>
    dplyr::mutate(
      years_expected = .data$origin_years_expected,
      climate_complete =
        .data$origin_years_available == .data$origin_years_expected &
        .data$destination_years_available == .data$destination_years_expected,
      climate_summary_status = dplyr::if_else(.data$climate_complete, "complete", "partial"),
      climate_coverage_note = dplyr::if_else(
        .data$climate_complete,
        "complete monthly CHELSA coverage",
        paste0(
          "partial average; origin ", .data$origin_years_available, "/",
          .data$origin_years_expected, " years (missing ", .data$origin_missing_years,
          "); destination ", .data$destination_years_available, "/",
          .data$destination_years_expected, " years (missing ",
          .data$destination_missing_years, ")"
        )
      )
    )
  out
}
