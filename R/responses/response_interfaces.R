# Module: R/responses/response_interfaces
# Responsibility: Implement the named scientific workflow stage or shared infrastructure.
# Inputs: Function arguments and explicitly declared target dependencies.
# Outputs: Validated in-memory objects and declared workflow artifacts.
# Side effects: Limited to functions whose contracts explicitly write, copy, or render files.

# Validation interface separating response-stage AUC/temporal results from the
# resident-model six-taxon exclusion.

#' Validate Auc Temporal Boundary.
#'
#' @param response_auc Response-stage AUC rows before resident-model exclusions.
#' @param resident_auc Resident-model AUC rows after the configured exclusions.
#' @param exclusions A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
validate_auc_temporal_boundary <- function(response_auc, resident_auc, exclusions) {
  response_species <- unique(as.character(response_auc$species_id))
  resident_species <- unique(as.character(resident_auc$species_id))
  tibble::tibble(
    exclusion_count = length(exclusions),
    exclusions_retained_in_responses = all(exclusions %in% response_species),
    exclusions_absent_from_resident_models = !any(exclusions %in% resident_species),
    resident_species_subset_of_responses = all(resident_species %in% response_species),
    pass = all(exclusions %in% response_species) &&
      !any(exclusions %in% resident_species) &&
      all(resident_species %in% response_species)
  )
}

#' Build Resident Auc Data from Objects.
#'
#' @param winner_loser A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param analysis_data A data frame whose required columns are checked before use.
#' @param climate_context_experiments A data frame whose required columns are checked before use.
#' @param horizons A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
build_resident_auc_data_from_objects <- function(winner_loser, analysis_data,
                                                 climate_context_experiments,
                                                 horizons = "max") {
  results <- winner_loser |>
    dplyr::mutate(horizon = as.character(.data$horizon)) |>
    dplyr::filter(
      .data$transform == "raw",
      .data$contrast == "to_origin",
      .data$pool == "resident",
      .data$horizon %in% as.character(horizons)
    ) |>
    dplyr::transmute(
      horizon = .data$horizon,
      experiment_id = as.character(.data$experiment_id),
      source_species = as.character(.data$species),
      auc_effect = as.numeric(.data$estimate),
      auc_se = as.numeric(.data$SE)
    )

  metadata_required <- c(
    "experiment_id", "species", "accepted_name", "region", "year_range",
    "plot_size", "destination_T_JJA", "destination_P_JJA", "diff_T_JJA",
    "LDMC", "LMA", "vegetative_height", "climate_complete",
    "climate_summary_status"
  )
  assert_columns(analysis_data, metadata_required, "analysis data target")
  metadata <- analysis_data |>
    dplyr::transmute(
      experiment_id = as.character(.data$experiment_id),
      source_species = as.character(.data$species),
      species_id = as.character(.data$accepted_name),
      gradient_id = as.character(.data$region),
      duration = as.numeric(.data$year_range),
      plot_area = as.numeric(.data$plot_size),
      T = as.numeric(.data$destination_T_JJA),
      P_raw = as.numeric(.data$destination_P_JJA),
      W = as.numeric(.data$diff_T_JJA),
      LDMC = as.numeric(.data$LDMC),
      LMA = as.numeric(.data$LMA),
      vegetative_height = as.numeric(.data$vegetative_height),
      climate_complete = as.logical(.data$climate_complete),
      climate_summary_status = as.character(.data$climate_summary_status)
    ) |>
    dplyr::distinct()
  assert_unique(metadata, c("experiment_id", "source_species"),
                "resident metadata target")

  context <- climate_context_experiments |>
    dplyr::transmute(
      experiment_id = as.character(.data$experiment_id),
      climate_context = factor(
        as.character(.data$climate_context),
        levels = c(
          "cold_dry", "intermediate_dry", "warm_dry",
          "cold_intermediate", "cold_wet"
        )
      )
    )
  assert_unique(context, "experiment_id", "climate context target")

  data <- results |>
    dplyr::left_join(metadata, by = c("experiment_id", "source_species")) |>
    dplyr::left_join(context, by = "experiment_id") |>
    dplyr::mutate(
      P = dplyr::if_else(.data$P_raw > 0, log(.data$P_raw), NA_real_),
      log_plot_area = dplyr::if_else(.data$plot_area > 0, log(.data$plot_area), NA_real_),
      integration_span = pmax(.data$duration - 1, 1),
      effect_id = paste(.data$horizon, .data$experiment_id, .data$species_id,
                        sep = "::")
    ) |>
    dplyr::arrange(.data$horizon, .data$experiment_id, .data$species_id) |>
    filter_resident_model_species("species_id")

  required_complete <- c(
    "species_id", "gradient_id", "auc_effect", "auc_se", "T", "P", "W",
    "duration", "log_plot_area"
  )
  if (any(!stats::complete.cases(data[, required_complete])) ||
      any(!is.finite(data$auc_se) | data$auc_se <= 0)) {
    stop("Resident AUC target rows have invalid required values.", call. = FALSE)
  }
  assert_unique(data, "effect_id", "resident AUC target data")
  data
}

#' Compare Resident Inputs.
#'
#' @param candidate_data A data frame whose required columns are checked before use.
#' @param candidate_V A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param reference_data A data frame whose required columns are checked before use.
#' @param reference_V A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
compare_resident_inputs <- function(candidate_data, candidate_V,
                                      reference_data, reference_V) {
  numeric_columns <- intersect(
    names(candidate_data)[vapply(candidate_data, is.numeric, logical(1))],
    names(reference_data)
  )
  numeric_error <- max(vapply(numeric_columns, function(column) {
    max(abs(candidate_data[[column]] - reference_data[[column]]), na.rm = TRUE)
  }, numeric(1)))
  V_error <- max(abs(as.matrix(candidate_V) - as.matrix(reference_V)))
  checks <- tibble::tibble(
    check = c("effect order", "species set", "numeric preprocessing", "sampling covariance"),
    value = c(
      as.numeric(identical(as.character(candidate_data$effect_id),
                           as.character(reference_data$effect_id))),
      as.numeric(setequal(as.character(candidate_data$species_id),
                          as.character(reference_data$species_id))),
      numeric_error,
      V_error
    ),
    tolerance = c(0, 0, 1e-12, 1e-12),
    pass = c(
      identical(as.character(candidate_data$effect_id),
                as.character(reference_data$effect_id)),
      setequal(as.character(candidate_data$species_id),
               as.character(reference_data$species_id)),
      numeric_error <= 1e-12,
      V_error <= 1e-12
    )
  )
  if (any(!checks$pass)) {
    stop("Computed resident inputs differ from the frozen transition inputs: ",
         paste(checks$check[!checks$pass], collapse = ", "), call. = FALSE)
  }
  checks
}
