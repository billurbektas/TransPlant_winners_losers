# TRY traits (P2): LDMC, LMA and vegetative height per accepted species.
# A TRY record is attached through its own resolved name (the contributor's
# SpeciesName first, TRY's AccSpeciesName only as fallback), never by
# trusting AccSpeciesName, which merges species (Festuca pratensis records
# are filed under Festuca ovina). Cleaning follows
# legacy/R/data/traits.R:103-290.

try_key <- function(x) stringr::str_squish(repair_utf8(x))

# One accepted species per resolved TRY name. A name the rule flagged is used
# only if its accepted species is a study species and the matched name has
# the submitted genus and epithet, so that the doubt lies in the author or
# the rank only (decided 26 Oct). Community name decisions are not carried
# over to TRY names: they are specific to the study records (Thlaspi
# montanum in US_Montana is Noccaea fendleri; TRY's Thlaspi montanum is not).
build_try_name_map <- function(name_resolution, taxonomy) {
  study <- unique(stats::na.omit(taxonomy$accepted_species[taxonomy$identification_level %in% "species_or_lower"]))
  name_resolution |>
    dplyr::transmute(
      name = submitted_name, matched_name, needs_review, review_reason,
      accepted_species = dplyr::coalesce(
        accepted_species,
        dplyr::if_else(word_count(accepted_name) >= 2L, stringr::word(accepted_name, 1, 2), NA_character_))) |>
    dplyr::mutate(
      study = accepted_species %in% study,
      same_binomial = dplyr::coalesce(tolower(stringr::word(name, 1, 2)) == tolower(stringr::word(matched_name, 1, 2)), FALSE),
      used = !is.na(accepted_species) & (!needs_review | (study & same_binomial)))
}

lookup_species <- function(names, try_name_map) {
  used <- try_name_map[try_name_map$used, ]
  unname(stats::setNames(used$accepted_species, used$name)[try_key(names)])
}

# TRY species IDs holding at least one name that resolves to a study species
select_try_ids <- function(try_names, try_name_map) {
  study <- unique(try_name_map$accepted_species[try_name_map$used & try_name_map$study])
  hit <- lookup_species(try_names$SpeciesName, try_name_map) %in% study |
    lookup_species(try_names$AccSpeciesName, try_name_map) %in% study
  ids <- sort(unique(try_names$AccSpeciesID[hit]))
  if (!length(ids)) stop("No TRY species ID resolves to a study species")
  ids
}

read_try_ids <- function(try_file, ids) {
  d <- tidyTRY::read_try(try_file, species_ids = ids, progress = FALSE)
  if (!nrow(d)) stop("Reading TRY by species ID returned no rows")
  # some TRY text fields are not valid UTF-8 (windows-1252); string functions
  # downstream fail on them
  text <- names(d)[vapply(d, is.character, logical(1))]
  n_invalid <- vapply(d[text], function(x) sum(!is.na(x) & !stringi::stri_enc_isutf8(x)), integer(1))
  for (col in text[n_invalid > 0]) d[[col]] <- repair_utf8(d[[col]])
  attr(d, "repaired_utf8") <- n_invalid[n_invalid > 0]
  d$AccSpeciesName <- tidyTRY::fix_species_case(d$AccSpeciesName)
  d
}

# Accepted species per row: the contributor's name first, TRY's accepted name
# as fallback (as tidyTRY::assign_accepted_names(), with names compared after
# the same UTF-8 repair and spacing as in P1). The species that trusting
# TRY's AccSpeciesName would give is kept for the report.
attach_try_names <- function(try_raw, try_name_map) {
  from_species <- lookup_species(try_raw$SpeciesName, try_name_map)
  from_try <- lookup_species(try_raw$AccSpeciesName, try_name_map)
  study <- unique(try_name_map$accepted_species[try_name_map$used & try_name_map$study])
  try_raw$species <- dplyr::coalesce(from_species, from_try)
  try_raw$name_source <- dplyr::case_when(!is.na(from_species) ~ "SpeciesName",
                                          !is.na(from_try) ~ "AccSpeciesName")
  try_raw$species_if_trusting_try <- from_try
  out <- try_raw[try_raw$species %in% study | try_raw$species_if_trusting_try %in% study, ]
  attr(out, "rows_read") <- c(rows = nrow(try_raw), observations = dplyr::n_distinct(try_raw$ObservationID))
  attr(out, "repaired_utf8") <- attr(try_raw, "repaired_utf8")
  out
}

stage_count <- function(d, stage) {
  tibble::tibble(stage = stage, rows = nrow(d), observations = dplyr::n_distinct(d$ObservationID))
}

# legacy/R/data/traits.R:114-127, with experiments removed per observation
clean_try <- function(try_attached, dictionary, natural_values, max_error_risk = 4) {
  warns <- character()
  collect <- function(expr) withCallingHandlers(expr, warning = function(w) {
    warns <<- c(warns, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  counts <- stage_count(try_attached, "rows of study species (by our name or TRY's)")
  no_exp <- collect(tidyTRY::remove_experiments(try_attached, level = "observation",
                                                natural_values = natural_values$value))
  removed_values <- attr(no_exp, "experimental_values")
  counts <- dplyr::bind_rows(counts, stage_count(no_exp, "experimental observations removed"))
  dedup <- collect(tidyTRY::remove_duplicates(no_exp))
  counts <- dplyr::bind_rows(counts, stage_count(dedup, "TRY duplicates removed"))
  coordinates <- collect(tidyTRY::extract_coordinates(dedup))
  split <- collect(tidyTRY::split_traits(dedup, qualitative_ids = integer()))
  quant <- collect(tidyTRY::clean_quantitative(split$quantitative, max_error_risk = max_error_risk))
  counts <- dplyr::bind_rows(counts, stage_count(quant, paste0("quantitative, ErrorRisk < ", max_error_risk)))
  map <- tidyTRY::trait_map(.map = dplyr::transmute(dictionary, original = try_trait_name, short_name = analysis_trait))
  quant <- quant |>
    dplyr::filter(TraitName %in% dictionary$try_trait_name) |>
    tidyTRY::rename_traits(map = map) |>
    dplyr::filter(TraitName %in% c("LDMC", "SLA", "vegetative_height"))
  counts <- dplyr::bind_rows(counts, stage_count(quant, "LDMC, SLA and vegetative height"))
  list(quantitative = quant, coordinates = coordinates, removed_values = removed_values,
       counts = counts, warnings = warns)
}

# Köppen-Geiger zone per observation from tidyTRY's bundled raster; no
# coordinates or an excluded zone removes the observation
# (legacy/R/data/traits.R:183-207; zones and sha256 from legacy/config/analysis.R:24-27)
filter_try_climate <- function(try_clean,
                               excluded = c("Af", "Am", "Aw", "BWh", "BWk", "BSh", "BSk"),
                               raster_sha256 = "59b749fe6711e72c720ff903044eb0efb4373524ee434740ed7384a0284c155b") {
  raster <- system.file("extdata", "koppen_geiger_0p1.tif", package = "tidyTRY")
  if (!nzchar(raster) || !identical(digest::digest(file = raster, algo = "sha256"), raster_sha256)) {
    stop("The installed tidyTRY Koppen-Geiger raster does not match the registered sha256")
  }
  zones <- tidyTRY::extract_climate_zones(try_clean$coordinates)
  code_col <- intersect(c("climate_code", "Climate_code"), names(zones))[[1]]
  zones <- dplyr::distinct(tibble::tibble(ObservationID = zones$ObservationID, climate_code = zones[[code_col]]))
  if (anyDuplicated(zones$ObservationID)) stop("More than one climate zone for an observation")
  q <- dplyr::left_join(try_clean$quantitative, zones, by = "ObservationID", relationship = "many-to-one")
  no_coord <- is.na(q$climate_code)
  kept <- q[!no_coord & !q$climate_code %in% excluded, ]
  counts <- dplyr::bind_rows(
    try_clean$counts,
    stage_count(q[!no_coord, ], "with coordinates and a climate zone"),
    stage_count(kept, paste("outside the excluded zones", paste(excluded, collapse = ", "))))
  list(observations = kept, counts = counts, excluded = excluded)
}

# Units, ranges, LMA = 1 / SLA per observation, then medians per species
# (legacy/R/data/traits.R:213-290)
summarise_traits <- function(try_climate, taxonomy) {
  obs <- try_climate$observations
  # a recorded 0 (or less) is not a measurement, and SLA = 0 has no LMA
  # (decided 26 Oct: drop non-positive values of all three traits)
  non_positive <- obs |>
    dplyr::filter(StdValue <= 0) |>
    dplyr::count(trait = TraitName, DatasetID, Dataset, name = "observations")
  obs <- dplyr::filter(obs, StdValue > 0)
  units <- dplyr::distinct(obs, TraitName, UnitName)
  expected <- c(LDMC = "g/g", SLA = "mm2 mg-1", vegetative_height = "m")
  if (anyDuplicated(units$TraitName) || any(units$UnitName != expected[units$TraitName])) {
    stop("Unexpected TRY units: ", paste(units$TraitName, units$UnitName, sep = "=", collapse = ", "))
  }
  bounds <- list(LDMC = c(0, 1), SLA = c(0, 200), vegetative_height = c(0, 100))
  for (t in names(bounds)) {
    v <- obs$StdValue[obs$TraitName == t]
    if (length(v) && (min(v) < bounds[[t]][1] || max(v) > bounds[[t]][2])) {
      stop("TRY values of ", t, " outside the plausible range ", paste(bounds[[t]], collapse = "-"))
    }
  }
  long <- obs |>
    dplyr::transmute(
      species, trait = dplyr::recode(TraitName, SLA = "LMA", vegetative_height = "height"),
      value = dplyr::if_else(TraitName == "SLA", 1 / StdValue, StdValue),
      dataset_id = DatasetID, observation_id = ObservationID, climate_code,
      SpeciesName, AccSpeciesName, name_source, species_if_trusting_try)
  study <- unique(stats::na.omit(taxonomy$accepted_species[taxonomy$identification_level %in% "species_or_lower"]))
  ours <- dplyr::filter(long, species %in% study)
  traits <- ours |>
    dplyr::summarise(median = stats::median(value), n_obs = dplyr::n(), n_datasets = dplyr::n_distinct(dataset_id),
                     .by = c(species, trait)) |>
    tidyr::pivot_wider(names_from = trait, values_from = c(median, n_obs, n_datasets), names_glue = "{.value}_{trait}") |>
    dplyr::rename(dplyr::any_of(c(LDMC = "median_LDMC", LMA = "median_LMA", height = "median_height"))) |>
    dplyr::arrange(species)
  if (anyDuplicated(traits$species)) stop("Duplicate species in the trait table")
  list(traits = traits, long = long, non_positive = non_positive)
}
