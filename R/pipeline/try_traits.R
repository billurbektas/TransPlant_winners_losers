taxonomy_response <- function(species, method = "tnrs", ...) {
  species <- repair_utf8(species)
  species <- unique(species[!is.na(species) & nzchar(stringr::str_squish(species))])
  response <- tidyTRY::resolve_species(species, method = method, ...)
  response |>
    dplyr::mutate(
      name_original = .data$submitted_name,
      name_normalized = normalize_name(.data$submitted_name),
      accepted_normalized = normalize_name(.data$accepted_name)
    )
}

apply_taxonomy_review <- function(response, overrides, min_score = 0.90) {
  reviewed <- response |>
    dplyr::left_join(overrides, by = "name_original", suffix = c("", "_override")) |>
    dplyr::mutate(
      decision = dplyr::coalesce(.data$decision, "pending"),
      accepted_name_reviewed = dplyr::case_when(
        .data$decision == "exclude" ~ NA_character_,
        .data$decision == "accept" & !is.na(.data$accepted_name_override) ~ .data$accepted_name_override,
        TRUE ~ .data$accepted_name
      ),
      needs_review = .data$decision == "pending" &
        (is.na(.data$accepted_name_reviewed) | is.na(.data$score) | .data$score < min_score)
    )
  reviewed
}

validate_taxonomy_review <- function(reviewed, label = "taxonomy") {
  unresolved <- reviewed |>
    dplyr::filter(.data$needs_review)
  if (nrow(unresolved)) {
    stop(
      label, " has ", nrow(unresolved),
      " unresolved/low-score names. Review the saved response and add accept/exclude rows to config/taxonomy_overrides.csv.",
      call. = FALSE
    )
  }
  reviewed |>
    dplyr::filter(.data$decision != "exclude", !is.na(.data$accepted_name_reviewed))
}

try_trait_map <- function(dictionary) {
  map <- dictionary |>
    dplyr::transmute(original = .data$try_trait_name, short_name = .data$canonical_trait)
  tidyTRY::trait_map(.map = map)
}

try_stage_counts <- function(data, stage) {
  tibble::tibble(
    stage = stage,
    observations = nrow(data),
    species = dplyr::n_distinct(dplyr::coalesce(data$AccSpeciesName, data$SpeciesName), na.rm = TRUE),
    datasets = if ("DatasetID" %in% names(data)) dplyr::n_distinct(data$DatasetID, na.rm = TRUE) else NA_integer_
  )
}

clean_try_export <- function(try_file, species_names, dictionary, max_error_risk = 4L,
                             chunk_size = 100000L) {
  if (!file.exists(try_file)) stop("TRY export is missing: ", try_file, call. = FALSE)
  project_map <- try_trait_map(dictionary)
  raw <- tidyTRY::read_try(try_file, species = unique(species_names), chunk_size = chunk_size)
  audit <- try_stage_counts(raw, "species-filtered export")

  if ("AccSpeciesName" %in% names(raw)) {
    caps <- grepl("^[A-Z][A-Z]", raw$AccSpeciesName)
    raw$AccSpeciesName[caps & !is.na(caps)] <- tidyTRY::fix_species_case(raw$AccSpeciesName[caps & !is.na(caps)])
  }
  no_experiments <- tidyTRY::remove_experiments(raw)
  audit <- dplyr::bind_rows(audit, try_stage_counts(no_experiments, "experimental datasets removed"))
  deduplicated <- tidyTRY::remove_duplicates(no_experiments)
  audit <- dplyr::bind_rows(audit, try_stage_counts(deduplicated, "TRY duplicates removed"))

  coordinates <- tidyTRY::extract_coordinates(deduplicated)
  split <- tidyTRY::split_traits(deduplicated, qualitative_ids = integer())
  quantitative <- tidyTRY::clean_quantitative(split$quantitative, max_error_risk = max_error_risk)
  audit <- dplyr::bind_rows(audit, try_stage_counts(quantitative, paste0("quantitative ErrorRisk < ", max_error_risk)))
  quantitative <- quantitative |>
    dplyr::filter(.data$TraitName %in% dictionary$try_trait_name) |>
    tidyTRY::rename_traits(map = project_map) |>
    dplyr::filter(.data$TraitName %in% c("LDMC", "SLA", "vegetative_height"))
  audit <- dplyr::bind_rows(audit, try_stage_counts(quantitative, "project aliases and focal traits"))

  list(
    quantitative = quantitative,
    coordinates = coordinates,
    units = quantitative |>
      dplyr::distinct(.data$TraitName, .data$UnitName),
    audit = audit
  )
}

try_side_species <- function(cleaned_try) {
  unique(c(cleaned_try$quantitative$SpeciesName, cleaned_try$quantitative$AccSpeciesName)) |>
    stats::na.omit() |>
    as.character() |>
    repair_utf8()
}

join_try_taxonomy <- function(quantitative, try_taxonomy) {
  lookup <- try_taxonomy |>
    dplyr::select(name_original = submitted_name, accepted_name_reviewed) |>
    dplyr::distinct()
  quantitative |>
    dplyr::left_join(lookup, by = c("AccSpeciesName" = "name_original")) |>
    dplyr::left_join(
      lookup |>
        dplyr::rename(accepted_from_species = accepted_name_reviewed),
      by = c("SpeciesName" = "name_original")
    ) |>
    dplyr::mutate(
      accepted_name = dplyr::coalesce(.data$accepted_name_reviewed, .data$accepted_from_species)
    ) |>
    dplyr::select(-accepted_name_reviewed, -accepted_from_species)
}

filter_try_climate <- function(cleaned_try, try_taxonomy, excluded_climates) {
  climate <- tidyTRY::extract_climate_zones(cleaned_try$coordinates)
  climate_column <- intersect(c("climate_code", "Climate_code"), names(climate))[[1]]
  climate <- climate |>
    dplyr::transmute(ObservationID = .data$ObservationID, climate_code = .data[[climate_column]])
  joined <- join_try_taxonomy(cleaned_try$quantitative, try_taxonomy) |>
    dplyr::left_join(climate, by = "ObservationID")
  filtered <- joined |>
    dplyr::filter(!is.na(.data$climate_code), !.data$climate_code %in% excluded_climates)
  list(
    observations = filtered,
    climate = climate,
    audit = dplyr::bind_rows(
      cleaned_try$audit,
      try_stage_counts(joined, "climate zone joined"),
      try_stage_counts(filtered, "arid/tropical/missing climates removed")
    )
  )
}

validate_trait_units <- function(observations) {
  unit_table <- observations |>
    dplyr::distinct(.data$TraitName, .data$UnitName)
  multiple <- unit_table |>
    dplyr::count(.data$TraitName) |>
    dplyr::filter(.data$n != 1L)
  if (nrow(multiple)) {
    stop("Expected exactly one standardized unit per focal TRY trait: ",
         paste(multiple$TraitName, collapse = ", "), call. = FALSE)
  }
  expected_units <- c(LDMC = "g/g", SLA = "mm2 mg-1", vegetative_height = "m")
  unexpected <- unit_table |>
    dplyr::mutate(expected = unname(expected_units[.data$TraitName])) |>
    dplyr::filter(.data$UnitName != .data$expected)
  if (nrow(unexpected)) {
    stop("Unexpected standardized TRY unit(s): ",
         paste(paste(unexpected$TraitName, unexpected$UnitName, sep = "="), collapse = ", "), call. = FALSE)
  }
  sla <- observations |>
    dplyr::filter(.data$TraitName == "SLA")
  if (nrow(sla) && any(!is.finite(sla$StdValue) | sla$StdValue <= 0, na.rm = TRUE)) {
    stop("SLA contains non-positive or non-finite standardized values; LMA was not calculated.", call. = FALSE)
  }
  unit_table
}

validate_trait_ranges <- function(observations) {
  bounds <- tibble::tribble(
    ~TraitName, ~plausible_min, ~plausible_max,
    "LDMC", 0, 1,
    "SLA", 0, 200,
    "vegetative_height", 0, 100
  )
  report <- observations |>
    dplyr::group_by(.data$TraitName) |>
    dplyr::summarise(observed_min = min(.data$StdValue, na.rm = TRUE),
                     observed_max = max(.data$StdValue, na.rm = TRUE), .groups = "drop") |>
    dplyr::left_join(bounds, by = "TraitName") |>
    dplyr::mutate(pass = .data$observed_min >= .data$plausible_min & .data$observed_max <= .data$plausible_max)
  if (any(!report$pass)) stop("Focal TRY trait values exceed the documented plausible-range gate.", call. = FALSE)
  report
}

summarise_try_traits <- function(observations) {
  validate_trait_units(observations)
  long <- observations |>
    dplyr::transmute(
      accepted_name = .data$accepted_name,
      trait = dplyr::if_else(.data$TraitName == "SLA", "LMA", .data$TraitName),
      value = dplyr::if_else(.data$TraitName == "SLA", 1 / .data$StdValue, .data$StdValue),
      source_unit = .data$UnitName,
      unit = dplyr::if_else(.data$TraitName == "SLA", paste0("1/(", .data$UnitName, ")"), .data$UnitName),
      dataset_id = .data$DatasetID,
      climate_code = .data$climate_code
    ) |>
    dplyr::filter(!is.na(.data$accepted_name), is.finite(.data$value))
  summary <- long |>
    dplyr::group_by(.data$accepted_name, .data$trait, .data$unit) |>
    dplyr::summarise(
      median = stats::median(.data$value),
      minimum = min(.data$value),
      maximum = max(.data$value),
      n_observations = dplyr::n(),
      n_datasets = dplyr::n_distinct(.data$dataset_id),
      datasets = paste(sort(unique(.data$dataset_id)), collapse = ";"),
      climate_zones = paste(sort(unique(.data$climate_code)), collapse = ";"),
      .groups = "drop"
    )
  list(long = long, summary = summary)
}

join_community_traits <- function(community_taxonomy, trait_summary) {
  community_taxonomy |>
    dplyr::select(species = name_original, accepted_name = accepted_name_reviewed) |>
    dplyr::distinct() |>
    dplyr::inner_join(trait_summary, by = "accepted_name")
}
