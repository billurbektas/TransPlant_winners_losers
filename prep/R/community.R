# Community, sites and experiments from the network SQLite (P1).
# Logic follows legacy/R/data/community.R; each block names its source lines.

read_sqlite_table <- function(path, table) {
  con <- DBI::dbConnect(RSQLite::SQLite(), path, flags = RSQLite::SQLITE_RO)
  on.exit(DBI::dbDisconnect(con))
  tibble::as_tibble(DBI::dbReadTable(con, table))
}

# The SQLite's taxonomy columns (submitted_name, Accepted_name,
# Taxonomic_status, Overall_score) come from TNRS "best" matches on
# unexpanded codes and are never read.
read_community <- function(path) {
  raw <- read_sqlite_table(path, "community")
  out <- tibble::tibble(
    region = raw$Region,
    year = as.integer(raw$Year),
    origin_site = raw$originSiteID,
    destination_site = raw$destSiteID,
    origin_block = raw$originBlockID,
    destination_block = raw$destBlockID,
    # CH_Lavey plot IDs carry the survey year as prefix ("2017_CRE_CRE_T1");
    # without it the ID is the same plot in every year
    destination_plot = ifelse(raw$Region == "CH_Lavey",
                              sub("^[0-9]{4}_", "", raw$destPlotID), raw$destPlotID),
    treatment = raw$Treatment,
    local_species_code = raw$SpeciesName,
    cover = raw$Cover
  ) |>
    dplyr::mutate(plot_id = paste(region, destination_plot, sep = "::"),
                  plot_year_id = paste(plot_id, year, sep = "::"))

  keys <- c("region", "year", "origin_site", "destination_site", "destination_plot", "treatment", "cover")
  na_keys <- vapply(out[keys], function(x) sum(is.na(x)), integer(1))
  if (any(na_keys > 0)) stop("SQLite community has NA in: ", paste(names(na_keys)[na_keys > 0], collapse = ", "))
  if (!all(out$treatment %in% c("LocalControl", "Warm"))) stop("Unexpected treatment values: ",
    paste(setdiff(unique(out$treatment), c("LocalControl", "Warm")), collapse = ", "))
  if (any(out$cover <= 0)) stop("SQLite community has non-positive cover")
  if (anyDuplicated(out[c("plot_year_id", "local_species_code")])) stop("Duplicate plot-year x species rows in the SQLite")
  # a plot keeps one origin, destination and treatment over its years
  per_plot <- dplyr::summarise(out, n = dplyr::n_distinct(origin_site, destination_site, treatment), .by = plot_id)
  if (any(per_plot$n > 1)) stop("Plots that change origin, destination or treatment across years: ",
    paste(head(per_plot$plot_id[per_plot$n > 1]), collapse = ", "))
  out
}

# Site metadata, with the legacy hard-coded corrections applied from
# prep/config/site_edits.csv (legacy/R/data/community.R:76-78)
read_sites <- function(path, site_edits) {
  meta <- read_sqlite_table(path, "meta")
  sites <- tibble::tibble(
    region = meta$Region, gradient = meta$Gradient, site_id = meta$destSiteID,
    longitude = meta$Longitude, latitude = meta$Latitude, elevation = meta$Elevation,
    year_established = as.integer(meta$YearEstablished),
    survey_year_min = as.integer(meta$YearMin), survey_year_max = as.integer(meta$YearMax),
    supplied_year_range = as.integer(meta$YearRange), plot_size = meta$PlotSize_m2,
    country = meta$Country
  )
  if (anyDuplicated(sites[c("region", "site_id")])) stop("Duplicate region x site_id in SQLite meta")
  for (i in seq_len(nrow(site_edits))) {
    e <- site_edits[i, ]
    hit <- sites$region == e$region & (e$site_id == "*" | sites$site_id == e$site_id) &
      as.character(sites[[e$field]]) == e$from
    if (!any(hit)) stop("Site edit ", e$edit_id, " matches no row")
    sites[[e$field]][hit] <- methods::as(e$to, class(sites[[e$field]]))
  }
  sites
}

# Local species codes to names, from the RData dictionary `cleanedspecies`
# plus the reviewed decisions for ambiguous codes
# (legacy/R/data/community.R:89-126; the RData is not used for anything else)
code_dictionary <- function(rdata_path, decisions) {
  env <- new.env()
  load(rdata_path, envir = env)
  cs <- env$cleanedspecies
  map <- tibble::tibble(region = cs$Region, site_id = as.character(cs$destSiteID),
                        local_species_code = as.character(cs$SpeciesName),
                        matched_name = cs$matched_name,
                        match_score = suppressWarnings(as.numeric(cs$score))) |>
    dplyr::filter(!is.na(matched_name), nzchar(matched_name))
  ambiguous <- map |>
    dplyr::summarise(n_matches = dplyr::n_distinct(matched_name),
                     .by = c(region, site_id, local_species_code)) |>
    dplyr::filter(n_matches > 1L) |>
    dplyr::left_join(decisions, by = c("region", "site_id", "local_species_code"))
  if (anyNA(ambiguous$decision)) stop("Ambiguous local codes without a decision in prep/config/species_mapping_decisions.csv")
  map |>
    dplyr::left_join(decisions, by = c("region", "site_id", "local_species_code")) |>
    dplyr::mutate(decision = dplyr::coalesce(decision, "accept"),
                  name_original = dplyr::if_else(!is.na(accepted_name), accepted_name, matched_name)) |>
    dplyr::slice_max(match_score, n = 1L, with_ties = FALSE, by = c(region, site_id, local_species_code)) |>
    dplyr::select(region, site_id, local_species_code, name_original, decision)
}

# Attach names; keys missing from the dictionary keep their recorded name
name_community <- function(community_raw, dictionary) {
  out <- community_raw |>
    dplyr::left_join(dictionary, by = c("region", destination_site = "site_id", "local_species_code"),
                     relationship = "many-to-one") |>
    dplyr::mutate(in_dictionary = !is.na(decision),
                  name_original = dplyr::case_when(
                    decision == "exclude" ~ NA_character_,
                    in_dictionary ~ name_original,
                    TRUE ~ local_species_code),
                  mapping_excluded = decision %in% "exclude")
  stopifnot(nrow(out) == nrow(community_raw))
  out
}

# relative_year from the destination site's establishment year; the
# establishment-year survey (relative_year 0) is dropped, as in legacy
# (legacy/R/data/community.R:229-240). Year edits drop whole surveys.
date_community <- function(community_named, sites, community_edits) {
  out <- community_named |>
    dplyr::left_join(dplyr::select(sites, region, destination_site = site_id, year_established),
                     by = c("region", "destination_site"), relationship = "many-to-one") |>
    dplyr::mutate(relative_year = year - year_established)
  if (anyNA(out$relative_year) || any(out$relative_year < 0)) stop("Surveys before establishment or without an establishment year")
  year_edits <- dplyr::filter(community_edits, field == "year")
  drop_year <- rep(FALSE, nrow(out))
  for (i in seq_len(nrow(year_edits))) {
    hit <- out$region == year_edits$region[i] & out$year == as.integer(year_edits$value[i])
    if (!any(hit)) stop("Community edit ", year_edits$edit_id[i], " matches no row")
    drop_year <- drop_year | hit
  }
  establishment <- out$relative_year == 0L
  dropped <- tibble::tibble(reason = c("establishment year", "year edit"),
                            rows = c(sum(establishment & !drop_year), sum(drop_year)))
  out <- out[!drop_year & !establishment, ]
  attr(out, "dropped") <- dropped
  out
}

# Records that do not enter the analysis keep their row and cover, so that
# relative cover can later be computed over the whole vegetation of a
# plot-year; `drop_reason` says why they are set aside.
finish_community <- function(community_dated, taxonomy, community_edits) {
  d <- community_dated |>
    dplyr::left_join(dplyr::select(taxonomy, name_original, species = accepted_species, identification_level),
                     by = "name_original", relationship = "many-to-one")
  edit_reason <- rep(NA_character_, nrow(d))
  taxon_edits <- dplyr::filter(community_edits, field != "year")
  for (i in seq_len(nrow(taxon_edits))) {
    e <- taxon_edits[i, ]
    origins <- strsplit(e$origin_site, ";", fixed = TRUE)[[1]]
    hit <- d$region == e$region & (e$origin_site == "*" | d$origin_site %in% origins) &
      (!as.logical(e$home_site_only) | d$origin_site == d$destination_site) &
      d[[e$field]] %in% e$value
    if (!any(hit)) stop("Community edit ", e$edit_id, " matches no row")
    edit_reason[hit & is.na(edit_reason)] <- paste0("community edit ", e$edit_id)
  }
  d$drop_reason <- dplyr::case_when(
    d$mapping_excluded ~ "ambiguous code excluded",
    !is.na(edit_reason) ~ edit_reason,
    is.na(d$identification_level) ~ "unidentified",
    d$identification_level != "species_or_lower" ~ d$identification_level,
    is.na(d$species) ~ "excluded by decision",
    TRUE ~ NA_character_)
  d$species[!is.na(d$drop_reason)] <- NA_character_

  # The rule only screens for non-plant rows that a future SQLite release
  # might bring; today every record is a plant
  non_plant <- "\\b(litter|bare|rocks?|stones?|soil|moss(es)?|bryophytes?|lichens?|dung|gravel|dead|total)\\b"
  d$is_vegetation <- !grepl(non_plant, dplyr::coalesce(d$name_original, ""), ignore.case = TRUE) &
    !grepl(non_plant, dplyr::coalesce(d$local_species_code, ""), ignore.case = TRUE)
  if (!all(d$is_vegetation)) stop("Non-vegetation records found: ",
    paste(unique(d$local_species_code[!d$is_vegetation]), collapse = ", "))

  # One row per plot-year x accepted species: names that resolve to the same
  # species are summed (synonyms, subspecies, split codes)
  keep_cols <- c("region", "year", "relative_year", "origin_site", "destination_site", "origin_block",
                 "destination_block", "destination_plot", "treatment", "plot_id", "plot_year_id")
  analysis <- d |>
    dplyr::filter(is.na(drop_reason)) |>
    dplyr::summarise(cover = sum(cover),
                     local_species_code = paste(sort(unique(local_species_code)), collapse = "; "),
                     name_original = paste(sort(unique(name_original)), collapse = "; "),
                     n_records = dplyr::n(),
                     .by = dplyr::all_of(c(keep_cols, "species"))) |>
    dplyr::mutate(drop_reason = NA_character_)
  set_aside <- d |>
    dplyr::filter(!is.na(drop_reason)) |>
    dplyr::mutate(n_records = 1L)
  out <- dplyr::bind_rows(analysis, set_aside[c(keep_cols, "species", "cover", "local_species_code",
                                                "name_original", "n_records", "drop_reason")]) |>
    dplyr::mutate(is_vegetation = TRUE) |>
    dplyr::arrange(plot_year_id, is.na(species), species)
  if (anyDuplicated(out[!is.na(out$species), c("plot_year_id", "species")])) stop("Duplicate plot-year x species after collapsing")
  stopifnot(isTRUE(all.equal(sum(out$cover), sum(d$cover))))
  out
}

# Experiments as legacy/R/data/community.R:345-394: one per warmed
# region x origin x destination; years from the warmed rows; exclusion rules
# in legacy order
build_experiments <- function(community, sites) {
  warmed <- dplyr::filter(community, treatment == "Warm", origin_site != destination_site)
  cells <- unique(paste(community$region, community$origin_site, community$destination_site,
                        community$treatment, sep = "::"))
  exps <- warmed |>
    dplyr::summarise(first_year = min(year), last_year = max(year),
                     .by = c(region, origin_site, destination_site)) |>
    dplyr::mutate(
      experiment_id = paste(region, origin_site, destination_site, sep = "::"),
      year_range = last_year - first_year + 1L,
      has_origin_control = paste(region, origin_site, origin_site, "LocalControl", sep = "::") %in% cells,
      has_destination_control = paste(region, destination_site, destination_site, "LocalControl", sep = "::") %in% cells,
      exclusion_reason = dplyr::case_when(
        # legacy rule by region; no duration threshold was ever stated
        region == "IN_Kashmir" ~ "insufficient duration",
        region == "FR_AlpeHuez" & is.na(first_year) ~ "warmed observations have missing Year",
        !has_origin_control ~ "missing origin control",
        !has_destination_control ~ "missing destination control",
        TRUE ~ NA_character_),
      included = is.na(exclusion_reason)) |>
    dplyr::arrange(region, origin_site, destination_site)
  region_info <- sites |>
    dplyr::summarise(year_established = dplyr::first(stats::na.omit(year_established)),
                     plot_size = dplyr::first(stats::na.omit(plot_size)), .by = region)
  exps <- dplyr::left_join(exps, region_info, by = "region", relationship = "many-to-one")
  if (nrow(exps) != 44L) stop("Expected 44 warmed comparisons, found ", nrow(exps))
  if (!all(exps$has_origin_control & exps$has_destination_control)) stop("Experiments without controls")
  dplyr::relocate(exps, experiment_id, .after = destination_site)
}
