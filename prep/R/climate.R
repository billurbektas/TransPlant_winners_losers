# Experiment climate from the cached CHELSA V2.1 monthly files (P4).
# Logic follows legacy/R/data/climate.R (lines named below); the cache is
# read, never downloaded.

chelsa_cache_files <- function(dir = "data/cache/chelsa_monthly") {
  files <- sort(list.files(dir, pattern = "^(tas|pr)_[0-9]{4}-[0-9]{2}[.]rds$", recursive = TRUE, full.names = TRUE))
  ym <- sub("^(tas|pr)_", "", sub("[.]rds$", "", basename(files)))
  both <- table(ym)
  if (any(both != 2L)) stop("Year-months without both tas and pr: ", paste(names(both)[both != 2L], collapse = ", "))
  files
}

# All months, with the site keys and coordinates checked against the prep
# sites exactly (legacy/R/data/climate.R:154-180)
read_chelsa <- function(files, sites) {
  d <- dplyr::bind_rows(lapply(files, readRDS))
  want <- dplyr::transmute(sites, site_key = paste(region, site_id, sep = "::"), longitude, latitude)
  have <- dplyr::distinct(d, site_key, longitude, latitude)
  if (!setequal(have$site_key, want$site_key) || anyDuplicated(have$site_key)) stop("CHELSA cache sites differ from the prep sites")
  m <- dplyr::inner_join(have, want, by = "site_key", suffix = c("_cache", "_site"))
  if (any(m$longitude_cache != m$longitude_site | m$latitude_cache != m$latitude_site)) stop("CHELSA cache coordinates differ from the prep sites")
  if (any(!is.finite(d$value))) stop("Non-finite CHELSA values")
  d
}

# Per site and year: day-weighted mean temperature, precipitation sums
# (legacy/R/data/climate.R:322-355)
summarise_site_years <- function(chelsa) {
  months <- dplyr::count(chelsa, site_key, variable, year)
  if (any(months$n != 12L)) stop("Site-years without 12 months of CHELSA data")
  tas <- chelsa |>
    dplyr::filter(variable == "tas") |>
    dplyr::summarise(T_annual = stats::weighted.mean(value, days_in_month),
                     T_JJA = stats::weighted.mean(value[month %in% 6:8], days_in_month[month %in% 6:8]),
                     .by = c(site_key, region, site_id, year))
  pr <- chelsa |>
    dplyr::filter(variable == "pr") |>
    dplyr::summarise(P_annual = sum(value), P_JJA = sum(value[month %in% 6:8]),
                     .by = c(site_key, region, site_id, year))
  dplyr::inner_join(tas, pr, by = c("site_key", "region", "site_id", "year"))
}

# Per included experiment: origin and destination means over the experiment
# years, destination minus origin (legacy/R/data/climate.R:362-423)
summarise_experiment_climate <- function(experiments, site_years) {
  side <- function(exp, site_col, prefix) {
    years <- tidyr::uncount(dplyr::transmute(exp, experiment_id, region, site_id = .data[[site_col]],
                                             first_year, n = last_year - first_year + 1L), n, .id = "k") |>
      dplyr::mutate(year = first_year + k - 1L)
    # years are counted before the means replace the yearly columns
    s <- dplyr::left_join(years, site_years, by = c("region", "site_id", "year")) |>
      dplyr::summarise(
        years_expected = dplyr::n(), years_available = sum(!is.na(T_JJA) & !is.na(P_annual)),
        missing_years = if (any(is.na(T_JJA))) paste(year[is.na(T_JJA)], collapse = ";") else "none",
        T_JJA = if (all(is.na(T_JJA))) NA_real_ else mean(T_JJA, na.rm = TRUE),
        P_JJA = if (all(is.na(P_JJA))) NA_real_ else mean(P_JJA, na.rm = TRUE),
        P_annual = if (all(is.na(P_annual))) NA_real_ else mean(P_annual, na.rm = TRUE),
        .by = experiment_id)
    names(s)[-1] <- paste(prefix, names(s)[-1], sep = "_")
    s
  }
  exp <- dplyr::filter(experiments, included)
  out <- exp |>
    dplyr::select(experiment_id, region, origin_site, destination_site, first_year, last_year) |>
    dplyr::left_join(side(exp, "origin_site", "origin"), by = "experiment_id") |>
    dplyr::left_join(side(exp, "destination_site", "destination"), by = "experiment_id") |>
    dplyr::mutate(
      diff_T_JJA = destination_T_JJA - origin_T_JJA,
      T = destination_T_JJA, W = diff_T_JJA,
      P_summer = log(destination_P_JJA), P_annual = log(destination_P_annual),
      climate_complete = origin_years_available == origin_years_expected &
        destination_years_available == destination_years_expected,
      climate_coverage_note = dplyr::if_else(
        climate_complete, "complete monthly CHELSA coverage",
        sprintf("partial average; origin %d/%d years (missing %s); destination %d/%d years (missing %s)",
                origin_years_available, origin_years_expected, origin_missing_years,
                destination_years_available, destination_years_expected, destination_missing_years)))
  if (any(!is.finite(c(out$T, out$W, out$P_summer, out$P_annual)))) stop("Non-finite T, W, P_summer or P_annual")
  dplyr::select(out, experiment_id, T, W, P_summer, P_annual, destination_T_JJA, origin_T_JJA, diff_T_JJA,
                destination_P_JJA, destination_P_annual, origin_years_expected, origin_years_available,
                origin_missing_years, destination_years_expected, destination_years_available,
                destination_missing_years, climate_complete, climate_coverage_note)
}

# CHELSA 1981-2010 summer temperature at the destination site (Lynn's site
# column, the bio10 analogue behind T_opt) against the ladder's JJA T
bio10_check <- function(niche_file, experiment_climate, experiments) {
  lynn_sites <- utils::read.csv(niche_file) |>
    dplyr::distinct(region = Region, destination_site = destSiteID, bio10_site = sum_temp_site) |>
    dplyr::summarise(bio10_site = mean(bio10_site), .by = c(region, destination_site))
  experiments |>
    dplyr::filter(included) |>
    dplyr::select(experiment_id, region, destination_site) |>
    dplyr::left_join(dplyr::select(experiment_climate, experiment_id, T), by = "experiment_id") |>
    dplyr::left_join(lynn_sites, by = c("region", "destination_site")) |>
    dplyr::mutate(difference = T - bio10_site)
}
