# Module: R/validation/resident_scope
# Responsibility: Check resident-model cohort counts, exclusions, identifiers,
# and diagonal sampling variances.
# Inputs: Resident effect rows, their variance matrix, and the analysis settings.
# Outputs: Validated cohort counts and sampling-variance checks.
# Side effects: None; invalid inputs stop the affected target.

#' Validate the resident-model cohort and its sampling variances
#'
#' @param auc_data Resident effects with `effect_id`, `species_id`,
#'   `experiment_id`, and `auc_se` columns.
#' @param sampling_variance_matrix Diagonal matrix in exact effect-row order.
#' @param config Analysis settings with expected counts and excluded taxa.
#' @details Stops if a count, exclusion, identifier, or variance check fails.
#' @return Cohort counts, identifiers, variance validation, and six checks.
validate_resident_model_scope <- function(
    auc_data, sampling_variance_matrix, config) {
  ids <- as.character(auc_data$effect_id)
  species <- unique(as.character(auc_data$species_id))
  experiments <- unique(as.character(auc_data$experiment_id))
  expected <- config$expected
  validation <- validate_sampling_variance_matrix(
    sampling_variance_matrix, auc_data, "auc_se"
  )
  checks <- tibble::tibble(
    check_number = seq_len(6L),
    check = c(
      "Expected resident effect count",
      "Expected resident species count",
      "Expected resident experiment count",
      "Excluded taxa absent from resident rows",
      "Unique effect identifiers",
      "Valid diagonal sampling variances"
    ),
    pass = c(
      nrow(auc_data) == expected$effects,
      length(species) == expected$species,
      length(experiments) == expected$experiments,
      !any(species %in% config$exclusions),
      !anyDuplicated(ids),
      isTRUE(validation$valid[[1]])
    ),
    detail = c(
      paste(nrow(auc_data), "effects"),
      paste(length(species), "species"),
      paste(length(experiments), "experiments"),
      paste(length(intersect(species, config$exclusions)), "excluded taxa found"),
      paste(anyDuplicated(ids), "first duplicate index"),
      paste(
        "dimensions", validation$dimensions_valid,
        "order", validation$order_valid,
        "maximum off-diagonal", validation$maximum_absolute_off_diagonal
      )
    )
  )
  gates <- c(
    nrow(auc_data) == expected$effects,
    length(species) == expected$species,
    length(experiments) == expected$experiments,
    !any(species %in% config$exclusions),
    !anyDuplicated(ids),
    isTRUE(validation$valid[[1]]),
    all(checks$pass)
  )
  if (!all(gates)) stop("Validated resident-model scope failed.", call. = FALSE)
  list(
    effects = as.integer(nrow(auc_data)),
    species = as.integer(length(species)),
    experiments = as.integer(length(experiments)),
    effect_ids = ids,
    species_ids = sort(species),
    experiment_ids = sort(experiments),
    sampling_variance_validation = validation,
    consistency_checks = checks
  )
}

#' Summarize the resident sampling-variance check for release
#'
#' @param validation One-row result from `validate_sampling_variance_matrix()`.
#' @details Errors on a malformed input.
#' @return One check row with ordering, variance range, and off-diagonal size.
resident_sampling_variance_release_check <- function(validation) {
  if (!is.data.frame(validation) || !"valid" %in% names(validation) ||
      nrow(validation) != 1L) {
    stop("Resident sampling-variance validation must contain one valid row.",
         call. = FALSE)
  }
  tibble::tibble(
    check = "validated resident diagonal sampling variances",
    pass = as.logical(validation$valid),
    detail = paste(
      "dimensions", validation$dimensions_valid,
      "order", validation$order_valid,
      "minimum variance",
      format(validation$minimum_sampling_variance, digits = 8),
      "maximum off-diagonal",
      format(validation$maximum_absolute_off_diagonal, digits = 8)
    )
  )
}
