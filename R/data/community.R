# Module: R/data/community
# Responsibility: Implement the named scientific workflow stage or shared infrastructure.
# Inputs: Function arguments and explicitly declared target dependencies.
# Outputs: Validated in-memory objects and declared workflow artifacts.
# Side effects: Limited to functions whose contracts explicitly write, copy, or render files.

#' Read Current Transplant Input.
#'
#' @param path A project-relative or absolute filesystem path.
#' @return An R object described by the function name; its structure is validated before later stages use it.
read_current_transplant_input <- function(path) {
  if (!file.exists(path)) stop("Current TransPlant input is missing: ", path, call. = FALSE)
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  required <- c("dat", "cleanedspecies")
  missing <- setdiff(required, ls(env, all.names = TRUE))
  if (length(missing)) stop("Input RData is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  regional <- mget(ls(env), envir = env) |>
    purrr::keep(~ is.list(.x) && all(c("meta", "community", "cover", "taxa") %in% names(.x)))
  list(dat = env$dat, cleanedspecies = env$cleanedspecies, regional = regional)
}

#' Current Input Inventory.
#'
#' @param input A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
current_input_inventory <- function(input) {
  tibble::tibble(
    object = c("dat", "cleanedspecies", "regional_objects"),
    rows = c(nrow(input$dat), nrow(input$cleanedspecies), length(input$regional)),
    description = c(
      "Long-format vegetation observations",
      "Site-specific local-code to botanical-name mappings",
      "Regional lists containing meta/community/cover/taxa and sometimes trait"
    )
  )
}

#' Site Metadata from Regions.
#'
#' @param regional A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
site_metadata_from_regions <- function(regional) {
  purrr::imap_dfr(regional, function(region_object, region_name) {
    meta <- tibble::as_tibble(region_object$meta)
    names(meta) <- stringr::str_replace_all(names(meta), "[.]", "_")
    candidate <- function(pattern) {
      hit <- names(meta)[stringr::str_detect(tolower(names(meta)), pattern)]
      if (length(hit)) hit[[1]] else NA_character_
    }
    site_col <- candidate("site.*id|siteid|^site$")
    lon_col <- candidate("longitude|^lon$|long$")
    lat_col <- candidate("latitude|^lat$")
    elev_col <- candidate("elevation|^elev$")
    established_col <- candidate("yearestablished|year_established")
    year_min_col <- candidate("yearmin|year_min")
    year_max_col <- candidate("yearmax|year_max")
    year_range_col <- candidate("yearrange|year_range")
    plot_size_col <- candidate("plotsize|plot_size")
    country_col <- candidate("country")
    if (is.na(site_col)) return(tibble::tibble())
    tibble::tibble(
      region = region_name,
      site_id = as.character(meta[[site_col]]),
      longitude = if (!is.na(lon_col)) as.numeric(meta[[lon_col]]) else NA_real_,
      latitude = if (!is.na(lat_col)) as.numeric(meta[[lat_col]]) else NA_real_,
      elevation = if (!is.na(elev_col)) as.numeric(meta[[elev_col]]) else NA_real_,
      year_established = if (!is.na(established_col)) as.integer(meta[[established_col]]) else NA_integer_,
      survey_year_min = if (!is.na(year_min_col)) as.integer(meta[[year_min_col]]) else NA_integer_,
      survey_year_max = if (!is.na(year_max_col)) as.integer(meta[[year_max_col]]) else NA_integer_,
      supplied_year_range = if (!is.na(year_range_col)) as.integer(meta[[year_range_col]]) else NA_integer_,
      plot_size = if (!is.na(plot_size_col)) as.numeric(meta[[plot_size_col]]) else NA_real_,
      country = if (!is.na(country_col)) as.character(meta[[country_col]]) else NA_character_
    ) |>
      dplyr::filter(!is.na(.data$site_id)) |>
      dplyr::mutate(
        year_established = dplyr::if_else(.data$region == "FR_Lautaret", 2016L, .data$year_established),
        elevation = dplyr::if_else(.data$region == "CH_Calanda" & .data$elevation == 2800, 2600, .data$elevation)
      )
  }) |>
    dplyr::distinct()
}

#' Resolve Local Species Mappings.
#'
#' @param cleanedspecies A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param decisions A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
resolve_local_species_mappings <- function(cleanedspecies, decisions) {
  mapping <- cleanedspecies |>
    dplyr::transmute(
      region = .data$Region,
      site_id = as.character(.data$destSiteID),
      local_species_code = as.character(.data$SpeciesName),
      matched_name = as.character(.data$matched_name),
      match_score = suppressWarnings(as.numeric(.data$score))
    ) |>
    dplyr::filter(!is.na(.data$matched_name), nzchar(.data$matched_name))
  ambiguity <- mapping |>
    dplyr::distinct(.data$region, .data$site_id, .data$local_species_code, .data$matched_name) |>
    dplyr::count(.data$region, .data$site_id, .data$local_species_code, name = "n_matches") |>
    dplyr::filter(.data$n_matches > 1L)
  reviewed <- ambiguity |>
    dplyr::left_join(decisions, by = c("region", "site_id", "local_species_code"))
  if (any(is.na(reviewed$decision))) {
    stop("Unreviewed ambiguous local species mappings remain; inspect output/tables/species_mapping_review.csv", call. = FALSE)
  }

  mapping |>
    dplyr::left_join(
      decisions |>
        dplyr::select(region, site_id, local_species_code, decision, accepted_name),
      by = c("region", "site_id", "local_species_code")
    ) |>
    dplyr::mutate(
      decision = dplyr::coalesce(.data$decision, "accept"),
      accepted_name = dplyr::if_else(
        !is.na(.data$accepted_name), .data$accepted_name, .data$matched_name
      )
    ) |>
    dplyr::filter(.data$decision != "exclude") |>
    dplyr::group_by(.data$region, .data$site_id, .data$local_species_code) |>
    dplyr::slice_max(.data$match_score, n = 1L, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::distinct(.data$region, .data$site_id, .data$local_species_code, .keep_all = TRUE)
}

#' Local Species Mapping Review.
#'
#' @param cleanedspecies A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param decisions A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
local_species_mapping_review <- function(cleanedspecies, decisions) {
  cleanedspecies |>
    dplyr::transmute(
      region = .data$Region,
      site_id = as.character(.data$destSiteID),
      local_species_code = as.character(.data$SpeciesName),
      matched_name = as.character(.data$matched_name),
      score = suppressWarnings(as.numeric(.data$score))
    ) |>
    dplyr::filter(!is.na(.data$matched_name), nzchar(.data$matched_name)) |>
    dplyr::distinct() |>
    dplyr::group_by(.data$region, .data$site_id, .data$local_species_code) |>
    dplyr::filter(dplyr::n_distinct(.data$matched_name) > 1L) |>
    dplyr::ungroup() |>
    dplyr::left_join(decisions, by = c("region", "site_id", "local_species_code")) |>
    dplyr::arrange(.data$region, .data$site_id, .data$local_species_code, dplyr::desc(.data$score))
}

#' Standardize Community.
#'
#' @param dat A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param species_mapping A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param site_metadata A value supplied by the calling workflow; admissible values are checked in the function body.
#' @details `dat` must contain region, year, origin/destination site and plot,
#' treatment, species, and cover fields. Reviewed species mappings must be
#' unique by region/site/local code. The function constructs stable identifiers
#' and relative cover but does not apply resident-model exclusions. Missing
#' required columns, ambiguous mappings, or invalid cover totals raise errors.
#' @return An R object described by the function name; its structure is validated before later stages use it.
standardize_community <- function(dat, species_mapping, site_metadata = NULL) {
  assert_columns(dat, c("Region", "Year", "originSiteID", "destSiteID", "destPlotID", "Treatment", "SpeciesName", "Cover"), "dat")
  out <- dat |>
    dplyr::transmute(
      region = as.character(.data$Region),
      year = suppressWarnings(as.integer(.data$Year)),
      origin_site = as.character(.data$originSiteID),
      destination_site = as.character(.data$destSiteID),
      origin_block = as.character(.data$originBlockID),
      destination_block = as.character(.data$destBlockID),
      destination_plot = as.character(.data$destPlotID),
      treatment = as.character(.data$Treatment),
      turf_id = as.character(.data$turfID),
      local_species_code = as.character(.data$SpeciesName),
      cover = as.numeric(.data$Cover),
      supplied_relative_cover = as.numeric(.data$Rel_Cover),
      elevation = as.numeric(.data$Elevation)
    ) |>
    dplyr::mutate(
      stable_plot_id = paste(
        .data$region, .data$origin_site, dplyr::coalesce(.data$origin_block, ""),
        .data$destination_site, dplyr::coalesce(.data$destination_block, ""),
        .data$treatment, dplyr::coalesce(.data$destination_plot, ""),
        dplyr::coalesce(.data$turf_id, ""), sep = "::"
      ),
      plot_year_id = paste(.data$stable_plot_id, .data$year, sep = "::")
    ) |>
    dplyr::filter(
      (.data$cover > 0) | (.data$region == "US_Arizona" & .data$supplied_relative_cover > 0)
    ) |>
    dplyr::group_by(.data$region, .data$plot_year_id) |>
    dplyr::mutate(
      relative_cover = dplyr::if_else(
        .data$region == "US_Arizona" & !is.na(.data$supplied_relative_cover),
        .data$supplied_relative_cover,
        .data$cover / sum(.data$cover, na.rm = TRUE)
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::left_join(
      species_mapping |>
        dplyr::select(region, site_id, local_species_code, species = accepted_name),
      by = c("region", "destination_site" = "site_id", "local_species_code")
    ) |>
    dplyr::mutate(
      species = dplyr::case_when(
        .data$region == "CN_Damxung" & .data$species %in% c("Stipa capillacea", "Stipa capillata") ~ "Stipa",
        TRUE ~ .data$species
      ),
      relative_cover = dplyr::case_when(
        .data$region == "CN_Damxung" & .data$local_species_code == "Polygonaceae" ~ 0,
        .data$region == "DE_Susalps" & .data$origin_site %in% c("FE", "BT") &
          .data$destination_site == .data$origin_site & .data$species == "Festuca pratensis" ~ 0,
        .data$region == "DE_Susalps" & .data$origin_site == "GW" & .data$destination_site == "GW" &
          .data$species == "Elymus repens" ~ 0,
        TRUE ~ .data$relative_cover
      )
    ) |>
    dplyr::filter(!is.na(.data$species)) |>
    dplyr::filter(is.finite(.data$relative_cover), .data$relative_cover > 0) |>
    dplyr::filter(!(.data$region == "DE_Susalps" & .data$year == 2016L),
                  !(.data$region == "CH_Calanda2" & .data$year == 2020L))

  if (!is.null(site_metadata)) {
    established <- site_metadata |>
      dplyr::select(region, destination_site = site_id, year_established) |>
      dplyr::distinct()
    out <- out |>
      dplyr::left_join(established, by = c("region", "destination_site")) |>
      dplyr::mutate(relative_year = .data$year - .data$year_established) |>
      dplyr::filter(is.na(.data$relative_year) | .data$relative_year != 0L)
  } else {
    out <- out |>
      dplyr::mutate(year_established = NA_integer_, relative_year = NA_integer_)
  }

  out |>
    dplyr::group_by(.data$region, .data$origin_site, .data$destination_site, .data$treatment, .data$species) |>
    dplyr::filter(dplyr::n_distinct(.data$plot_year_id) > 1L) |>
    dplyr::ungroup()
}

#' Taxonomic Resolution Review.
#'
#' @param reviewed_taxonomy A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
taxonomic_resolution_review <- function(reviewed_taxonomy) {
  word_count <- function(x) {
    cleaned <- stringr::str_squish(x)
    dplyr::if_else(
      is.na(cleaned) | cleaned == "", NA_integer_,
      stringr::str_count(cleaned, stringr::fixed(" ")) + 1L
    )
  }
  reviewed_taxonomy |>
    dplyr::mutate(
      original_word_count = word_count(.data$name_original),
      accepted_word_count = word_count(.data$accepted_name_reviewed),
      explicit_unspecified_species = stringr::str_detect(
        .data$name_original,
        stringr::regex("(^|[[:space:]_])(sp|spp|spec|indet)[.]?$", ignore_case = TRUE)
      ),
      family_name = stringr::str_detect(
        dplyr::coalesce(.data$accepted_name_reviewed, .data$name_original),
        stringr::regex("aceae$", ignore_case = TRUE)
      ),
      identification_level = dplyr::case_when(
        is.na(.data$accepted_name_reviewed) ~ "excluded_or_unresolved",
        .data$family_name & .data$accepted_word_count == 1L ~ "family_only",
        .data$original_word_count == 1L | .data$accepted_word_count == 1L |
          .data$explicit_unspecified_species ~ "genus_only",
        .data$original_word_count >= 2L & .data$accepted_word_count >= 2L ~ "species_or_lower",
        TRUE ~ "unresolved_rank"
      ),
      include_species_level = .data$identification_level == "species_or_lower",
      resolution_decision = dplyr::if_else(
        .data$include_species_level, "include", "omit_non_species_level"
      )
    )
}

#' Filter Species Level Taxonomy.
#'
#' @param reviewed_taxonomy A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
filter_species_level_taxonomy <- function(reviewed_taxonomy) {
  resolution <- if ("include_species_level" %in% names(reviewed_taxonomy)) {
    reviewed_taxonomy
  } else {
    taxonomic_resolution_review(reviewed_taxonomy)
  }
  resolution |>
    dplyr::filter(.data$include_species_level)
}

#' Filter Community to Reviewed Taxa.
#'
#' @param community A data frame whose required columns are checked before use.
#' @param reviewed_taxonomy A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
filter_community_to_reviewed_taxa <- function(community, reviewed_taxonomy) {
  accepted_originals <- filter_species_level_taxonomy(reviewed_taxonomy) |>
    dplyr::pull(.data$name_original) |>
    unique()
  community |>
    dplyr::filter(.data$species %in% accepted_originals)
}

#' Build Sites.
#'
#' @param community A data frame whose required columns are checked before use.
#' @param regional A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
build_sites <- function(community, regional) {
  regional_sites <- site_metadata_from_regions(regional)
  observed <- community |>
    dplyr::transmute(region = .data$region, site_id = .data$destination_site, elevation_observed = .data$elevation) |>
    dplyr::group_by(.data$region, .data$site_id) |>
    dplyr::summarise(
      elevation_observed = {
        values <- stats::na.omit(.data$elevation_observed)
        if (length(values)) values[[1]] else NA_real_
      },
      .groups = "drop"
    )
  sites <- regional_sites |>
    dplyr::full_join(observed, by = c("region", "site_id")) |>
    dplyr::mutate(elevation = dplyr::coalesce(.data$elevation, .data$elevation_observed)) |>
    dplyr::select(-elevation_observed) |>
    dplyr::distinct()
  assert_unique(sites, c("region", "site_id"), "sites")
  sites
}

#' Build Experiments.
#'
#' @param community A data frame whose required columns are checked before use.
#' @param sites A data frame whose required columns are checked before use.
#' @return An R object described by the function name; its structure is validated before later stages use it.
build_experiments <- function(community, sites = NULL) {
  warmed <- community |>
    dplyr::filter(.data$treatment == "Warm", .data$origin_site != .data$destination_site) |>
    dplyr::distinct(.data$region, .data$origin_site, .data$destination_site) |>
    dplyr::mutate(experiment_id = paste(.data$region, .data$origin_site, .data$destination_site, sep = "::"))
  availability <- community |>
    dplyr::distinct(.data$region, .data$origin_site, .data$destination_site, .data$treatment) |>
    dplyr::mutate(key = paste(.data$region, .data$origin_site, .data$destination_site, .data$treatment, sep = "::")) |>
    dplyr::pull(.data$key)
  year_info <- community |>
    dplyr::filter(.data$treatment == "Warm", .data$origin_site != .data$destination_site, !is.na(.data$year)) |>
    dplyr::group_by(.data$region, .data$origin_site, .data$destination_site) |>
    dplyr::summarise(
      first_year = min(.data$year),
      last_year = max(.data$year),
      year_range = .data$last_year - .data$first_year + 1L,
      plot_size = NA_real_,
      .groups = "drop"
    )
  out <- warmed |>
    dplyr::mutate(
      has_origin_control = paste(.data$region, .data$origin_site, .data$origin_site, "LocalControl", sep = "::") %in% availability,
      has_destination_control = paste(.data$region, .data$destination_site, .data$destination_site, "LocalControl", sep = "::") %in% availability
    ) |>
    dplyr::left_join(year_info, by = c("region", "origin_site", "destination_site")) |>
    dplyr::mutate(
      exclusion_reason = dplyr::case_when(
        .data$region == "IN_Kashmir" ~ "insufficient duration",
        .data$region == "FR_AlpeHuez" & is.na(.data$first_year) ~ "warmed observations have missing Year",
        !.data$has_origin_control ~ "missing origin control",
        !.data$has_destination_control ~ "missing destination control",
        TRUE ~ NA_character_
      ),
      included = is.na(.data$exclusion_reason)
    ) |>
    dplyr::arrange(.data$region, .data$origin_site, .data$destination_site)
  if (!is.null(sites)) {
    region_meta <- sites |>
      dplyr::group_by(.data$region) |>
      dplyr::summarise(
        year_established = dplyr::first(stats::na.omit(.data$year_established)),
        plot_size = dplyr::first(stats::na.omit(.data$plot_size)),
        .groups = "drop"
      )
    out <- out |>
      dplyr::select(-plot_size) |>
      dplyr::left_join(region_meta, by = "region")
  }
  out
}

#' Validate Community Snapshot.
#'
#' @param sites A data frame whose required columns are checked before use.
#' @param community A data frame whose required columns are checked before use.
#' @param experiments A data frame whose required columns are checked before use.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
validate_community_snapshot <- function(sites, community, experiments) {
  checks <- tibble::tribble(
    ~check, ~observed, ~expected,
    "raw warmed comparisons", nrow(experiments), 44,
    "included comparisons", sum(experiments$included), 42,
    "unique sites", nrow(sites), 59,
    "missing origin controls", sum(!experiments$has_origin_control), 0,
    "missing destination controls", sum(!experiments$has_destination_control), 0,
    "stable plot IDs containing terminal year", sum(stringr::str_detect(community$stable_plot_id, "::(19|20)[0-9]{2}$")), 0
  ) |>
    dplyr::mutate(pass = .data$observed == .data$expected)
  if (any(!checks$pass)) {
    stop("Canonical input validation failed: ", paste(checks$check[!checks$pass], collapse = "; "), call. = FALSE)
  }
  checks
}
