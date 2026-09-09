# Canonical upstream analysis helpers extracted from the reviewed pipeline.
# Functions are object-in/object-out and contain no output-directory discovery.

tnrs_family_lookup <- function(accepted_names) {
  response <- TNRS::TNRS(unique(accepted_names[!is.na(accepted_names)]))
  family_col <- intersect(c("Accepted_family", "Name_matched_accepted_family", "Family_matched"), names(response))[[1]]
  tibble::tibble(
    accepted_name = response$Accepted_name,
    family = response[[family_col]],
    family_match_score = response$Overall_score
  ) |>
    dplyr::distinct(.data$accepted_name, .keep_all = TRUE)
}

derive_functional_groups <- function(taxonomy_family, overrides,
                                     downstream_species_overrides = NULL) {
  graminoid_families <- c("Poaceae", "Cyperaceae", "Juncaceae")
  fern_families <- c(
    "Aspleniaceae", "Athyriaceae", "Blechnaceae", "Cystopteridaceae", "Dennstaedtiaceae",
    "Dryopteridaceae", "Equisetaceae", "Lycopodiaceae", "Ophioglossaceae", "Polypodiaceae",
    "Pteridaceae", "Selaginellaceae", "Thelypteridaceae", "Woodsiaceae"
  )
  woody_genera <- c(
    "Arctostaphylos", "Berberis", "Betula", "Calluna", "Cassiope",
    "Daphne", "Dasiphora", "Dryas", "Empetrum", "Helianthemum",
    "Hippophae", "Juniperus", "Kalmia", "Kleinia", "Loiseleuria",
    "Phyllodoce", "Rhododendron", "Rosa", "Salix", "Sorbus", "Spiraea",
    "Symphoricarpos", "Thymus", "Vaccinium"
  )
  required_override_columns <- c(
    "accepted_name", "graminoid_status", "woody_status", "note"
  )
  missing_override_columns <- setdiff(required_override_columns, names(overrides))
  if (length(missing_override_columns)) {
    stop(
      "Functional-status overrides are missing: ",
      paste(missing_override_columns, collapse = ", "), call. = FALSE
    )
  }
  invalid_graminoid <- setdiff(
    unique(trimws(stats::na.omit(overrides$graminoid_status))),
    c("", "graminoid", "non_graminoid")
  )
  invalid_woody <- setdiff(
    unique(trimws(stats::na.omit(overrides$woody_status))),
    c("", "woody", "non_woody")
  )
  if (length(invalid_graminoid) || length(invalid_woody)) {
    stop(
      "Invalid functional-status override values. graminoid_status: ",
      paste(invalid_graminoid, collapse = ", "), "; woody_status: ",
      paste(invalid_woody, collapse = ", "), call. = FALSE
    )
  }
  taxonomy_augmented <- taxonomy_family
  if (!is.null(downstream_species_overrides)) {
    required_species_override_columns <- c(
      "species_id", "family_id", "graminoid_status", "woodiness",
      "decision_basis"
    )
    missing_species_override_columns <- setdiff(
      required_species_override_columns, names(downstream_species_overrides)
    )
    if (length(missing_species_override_columns)) {
      stop(
        "Downstream species overrides are missing: ",
        paste(missing_species_override_columns, collapse = ", "),
        call. = FALSE
      )
    }
    species_overrides <- downstream_species_overrides |>
      dplyr::transmute(
        accepted_name = as.character(.data$species_id),
        family_species_override = as.character(.data$family_id),
        graminoid_species_override = as.character(.data$graminoid_status),
        woody_species_override = as.character(.data$woodiness),
        species_override_basis = as.character(.data$decision_basis)
      )
    invalid_species_override <-
      anyDuplicated(species_overrides$accepted_name) ||
      anyNA(species_overrides) ||
      any(!nzchar(trimws(species_overrides$accepted_name))) ||
      any(!nzchar(trimws(species_overrides$family_species_override))) ||
      any(!species_overrides$graminoid_species_override %in%
            c("graminoid", "non_graminoid")) ||
      any(!species_overrides$woody_species_override %in%
            c("woody", "non_woody"))
    if (invalid_species_override) {
      stop(
        "Downstream species overrides are incomplete, invalid, or non-unique.",
        call. = FALSE
      )
    }
    taxonomy_augmented <- dplyr::full_join(
      taxonomy_augmented, species_overrides, by = "accepted_name"
    ) |>
      dplyr::mutate(
        family = dplyr::coalesce(
          .data$family_species_override, as.character(.data$family)
        )
      )
  } else {
    taxonomy_augmented <- taxonomy_augmented |>
      dplyr::mutate(
        family_species_override = NA_character_,
        graminoid_species_override = NA_character_,
        woody_species_override = NA_character_,
        species_override_basis = NA_character_
      )
  }
  taxonomy_augmented |>
    dplyr::mutate(genus = stringr::word(.data$accepted_name, 1)) |>
    dplyr::left_join(overrides, by = "accepted_name", suffix = c("", "_override")) |>
    dplyr::mutate(
      graminoid_override_supplied =
        (!is.na(.data$graminoid_species_override) &
           nzchar(trimws(.data$graminoid_species_override))) |
        (!is.na(.data$graminoid_status) &
           nzchar(trimws(.data$graminoid_status))),
      woody_override_supplied =
        (!is.na(.data$woody_species_override) &
           nzchar(trimws(.data$woody_species_override))) |
        (!is.na(.data$woody_status) & nzchar(trimws(.data$woody_status))),
      graminoid_status_auto = dplyr::case_when(
        is.na(.data$family) | is.na(.data$accepted_name) ~ NA_character_,
        .data$family %in% graminoid_families ~ "graminoid",
        TRUE ~ "non_graminoid"
      ),
      woody_status_auto = dplyr::case_when(
        is.na(.data$family) | is.na(.data$accepted_name) ~ NA_character_,
        .data$genus %in% woody_genera ~ "woody",
        TRUE ~ "non_woody"
      ),
      graminoid_status = dplyr::coalesce(
        .data$graminoid_species_override,
        .data$graminoid_status,
        .data$graminoid_status_auto
      ),
      woody_status = dplyr::coalesce(
        .data$woody_species_override,
        .data$woody_status,
        .data$woody_status_auto
      ),
      graminoid_review_basis = dplyr::case_when(
        !is.na(.data$graminoid_species_override) ~
          paste0("downstream species ledger: ", .data$species_override_basis),
        .data$graminoid_override_supplied ~
          paste0("override: ", .data$note),
        .data$family %in% graminoid_families ~
          "reviewed graminoid-family rule",
        !is.na(.data$family) ~ "reviewed non-graminoid family",
        TRUE ~ "unresolved family"
      ),
      woody_review_basis = dplyr::case_when(
        !is.na(.data$woody_species_override) ~
          paste0("downstream species ledger: ", .data$species_override_basis),
        .data$woody_override_supplied ~
          paste0("override: ", .data$note),
        .data$genus %in% woody_genera ~ "reviewed woody-genus rule",
        !is.na(.data$genus) ~ "reviewed non-woody genus",
        TRUE ~ "unresolved genus"
      ),
      # Deprecated compatibility label for old descriptive scripts only. New
      # inference must use graminoid_status and woody_status separately.
      functional_group = dplyr::case_when(
        .data$family %in% fern_families ~ "other_unresolved",
        .data$woody_status == "woody" ~ "other_unresolved",
        .data$graminoid_status == "graminoid" ~ "graminoid",
        is.na(.data$graminoid_status) | is.na(.data$woody_status) ~
          "other_unresolved",
        TRUE ~ "forb"
      ),
      review_basis = dplyr::case_when(
        !is.na(.data$species_override_basis) ~
          paste0("downstream species ledger: ", .data$species_override_basis),
        !is.na(.data$note) ~ paste0("override: ", .data$note),
        .data$family %in% graminoid_families ~ "reviewed graminoid-family rule",
        .data$family %in% fern_families ~ "fern-family exclusion",
        .data$genus %in% woody_genera ~ "woody-genus exclusion",
        is.na(.data$family) ~ "unresolved family",
        TRUE ~ "herbaceous non-graminoid candidate"
      )
    ) |>
    dplyr::select(
      accepted_name, family, genus, graminoid_status, woody_status,
      graminoid_review_basis, woody_review_basis,
      functional_group, review_basis
    )
}

join_experiment_climate <- function(data, experiment_climate) {
  duplicated_metadata <- setdiff(
    intersect(names(data), names(experiment_climate)),
    "experiment_id"
  )
  data |>
    dplyr::left_join(
      experiment_climate |>
        dplyr::select(-dplyr::all_of(duplicated_metadata)),
      by = "experiment_id"
    )
}

reviewed_five_climate_contexts <- function(data) {
  temperature_cutoff_low <- 12.86686795491147
  temperature_cutoff_high <- 16.94267310789054
  precipitation_cutoff_dry <- 550
  precipitation_cutoff_wet <- 768
  context_levels <- c(
    "cold_dry", "intermediate_dry", "warm_dry",
    "cold_intermediate", "cold_wet"
  )
  climate <- data |>
    dplyr::distinct(
      .data$experiment_id, .data$destination_T_JJA, .data$destination_P_JJA,
      .data$climate_summary_status
    ) |>
    dplyr::filter(
      is.finite(.data$destination_T_JJA),
      is.finite(.data$destination_P_JJA)
    )
  if (nrow(climate) == 0L) stop("No finite destination climates are available.", call. = FALSE)
  duplicates <- climate |>
    dplyr::count(.data$experiment_id) |>
    dplyr::filter(.data$n != 1L)
  if (nrow(duplicates)) {
    stop("Every experiment must have one unique destination climate.", call. = FALSE)
  }
  classified <- climate |>
    dplyr::mutate(
      climate_context = dplyr::case_when(
        .data$destination_P_JJA < .env$precipitation_cutoff_dry &
          .data$destination_T_JJA < .env$temperature_cutoff_low ~ "cold_dry",
        .data$destination_P_JJA < .env$precipitation_cutoff_dry &
          .data$destination_T_JJA < .env$temperature_cutoff_high ~ "intermediate_dry",
        .data$destination_P_JJA < .env$precipitation_cutoff_dry ~ "warm_dry",
        .data$destination_P_JJA < .env$precipitation_cutoff_wet &
          .data$destination_T_JJA < .env$temperature_cutoff_high ~ "cold_intermediate",
        .data$destination_P_JJA >= .env$precipitation_cutoff_wet &
          .data$destination_T_JJA < .env$temperature_cutoff_high ~ "cold_wet",
        TRUE ~ NA_character_
      ),
      climate_context = factor(.data$climate_context, levels = .env$context_levels),
      climate_context_label = factor(
        stringr::str_to_title(stringr::str_replace(.data$climate_context, "_", "/")),
        levels = c(
          "Cold/Dry", "Intermediate/Dry", "Warm/Dry",
          "Cold/Intermediate", "Cold/Wet"
        )
      ),
      temperature_cutoff_low = .env$temperature_cutoff_low,
      temperature_cutoff_high = .env$temperature_cutoff_high,
      precipitation_cutoff_dry = .env$precipitation_cutoff_dry,
      precipitation_cutoff_wet = .env$precipitation_cutoff_wet,
      context_review_status = dplyr::if_else(
        is.na(.data$climate_context), "review_required", "reviewed"
      ),
      context_review_reason = dplyr::if_else(
        is.na(.data$climate_context),
        "Warm non-dry destination climate is outside the reviewed five-group scheme.",
        NA_character_
      )
    ) |>
    dplyr::arrange(.data$climate_context, .data$experiment_id)
  list(
    assigned = dplyr::filter(classified, !is.na(.data$climate_context)),
    review = dplyr::filter(classified, is.na(.data$climate_context))
  )
}

validate_five_climate_contexts <- function(context_result) {
  if (nrow(context_result$review)) {
    stop(
      "Destination climate contexts require review; see climate_context_review.csv.",
      call. = FALSE
    )
  }
  context_result$assigned
}
analysis_dataset <- function(classified, community_taxonomy, trait_summary,
                             functional_groups, experiments, experiment_climate = NULL) {
  taxonomy_lookup <- community_taxonomy |>
    dplyr::select(species = name_original, accepted_name = accepted_name_reviewed) |>
    dplyr::distinct()
  traits <- trait_summary |>
    dplyr::select(accepted_name, trait, value = median) |>
    tidyr::pivot_wider(names_from = trait, values_from = value)
  result <- classified |>
    dplyr::filter(.data$transform == "raw", .data$horizon == "max") |>
    dplyr::mutate(
      calculation_category = .data$category,
      category = droplevels(.data$analysis_category),
      model_category = dplyr::case_when(
        .data$category == "loser" ~ "loser",
        .data$category == "winner" ~ "winner",
        .data$category == "colonizer" ~ "colonizer",
        stringr::str_starts(as.character(.data$category), "non_significant_") |
          .data$category == "no_change" ~ "no_change",
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::left_join(taxonomy_lookup, by = "species") |>
    dplyr::left_join(traits, by = "accepted_name") |>
    dplyr::left_join(functional_groups, by = "accepted_name") |>
    dplyr::left_join(
      experiments |>
        dplyr::select(experiment_id, plot_size, year_range),
      by = "experiment_id"
    )
  if (!is.null(experiment_climate)) {
    result <- join_experiment_climate(result, experiment_climate)
  }
  result
}

transform_model_predictors <- function(data) {
  positive <- c("LDMC", "LMA", "vegetative_height", "destination_P_JJA", "plot_size")
  missing <- setdiff(c(positive, "diff_T_JJA", "destination_T_JJA", "year_range"), names(data))
  if (length(missing)) stop("Model data are missing predictors: ", paste(missing, collapse = ", "), call. = FALSE)
  data |>
    dplyr::mutate(
      dplyr::across(dplyr::all_of(positive), ~ dplyr::if_else(.x > 0, log(.x), NA_real_), .names = "log_{.col}"),
      z_LDMC = safe_scale(.data$log_LDMC),
      z_LMA = safe_scale(.data$log_LMA),
      z_vegetative_height = safe_scale(.data$log_vegetative_height),
      z_diff_T_JJA = safe_scale(.data$diff_T_JJA),
      z_destination_T_JJA = safe_scale(.data$destination_T_JJA),
      z_destination_P_JJA = safe_scale(.data$log_destination_P_JJA),
      z_plot_size = safe_scale(.data$log_plot_size),
      z_year_range = safe_scale(.data$year_range)
    )
}
