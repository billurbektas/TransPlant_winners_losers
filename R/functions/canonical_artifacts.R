# Explicit adapters for the validated Step 19-23 artifacts during migration.
# These functions receive paths; they never search an output directory.

read_csv_canonical <- function(path) {
  result <- readr::read_csv(path, show_col_types = FALSE)
  empty_text_fields <- intersect(
    names(result), c("warnings", "error", "CR1_warnings", "contrast_warnings")
  )
  for (field in empty_text_fields) {
    result[[field]][is.na(result[[field]])] <- ""
  }
  result
}

classification_summary_canonical <- function(classification, resolution) {
  classes <- c("Consistent winner", "Consistent loser", "Context-dependent")
  counts <- table(factor(classification$fitted_class, levels = classes))
  total <- nrow(classification)
  winner_n <- as.integer(counts[["Consistent winner"]])
  loser_n <- as.integer(counts[["Consistent loser"]])
  context_n <- as.integer(counts[["Context-dependent"]])
  tibble::tibble(
    Resolution = resolution,
    `N eligible units` = as.integer(total),
    `Winner N` = winner_n,
    `Winner %` = 100 * winner_n / total,
    `Loser N` = loser_n,
    `Loser %` = 100 * loser_n / total,
    `Context-dependent N` = context_n,
    `Context-dependent %` = 100 * context_n / total
  )
}

validate_downstream_scope_artifacts <- function(auc_data, V_auc, checks, config) {
  ids <- as.character(auc_data$effect_id)
  species <- unique(as.character(auc_data$species_id))
  experiments <- unique(as.character(auc_data$experiment_id))
  expected <- config$expected
  dimensions_valid <- nrow(V_auc) == nrow(auc_data) && ncol(V_auc) == nrow(auc_data)
  order_valid <- dimensions_valid && identical(rownames(V_auc), ids) &&
    identical(colnames(V_auc), ids)
  V <- as.matrix(V_auc)
  validation <- tibble::tibble(
    dimensions_valid = dimensions_valid,
    order_valid = order_valid,
    maximum_symmetry_error = max(abs(V - t(V))),
    minimum_eigenvalue = min(eigen(V, symmetric = TRUE, only.values = TRUE)$values),
    maximum_diagonal_se_error = max(abs(diag(V) - auc_data$auc_se^2)),
    valid = dimensions_valid && order_valid &&
      max(abs(V - t(V))) <= 1e-10 &&
      min(eigen(V, symmetric = TRUE, only.values = TRUE)$values) >= -1e-10 &&
      max(abs(diag(V) - auc_data$auc_se^2)) <= 1e-8
  )
  gates <- c(
    nrow(auc_data) == expected$effects,
    length(species) == expected$species,
    length(experiments) == expected$experiments,
    !any(species %in% config$exclusions),
    !anyDuplicated(ids),
    isTRUE(validation$valid[[1]]),
    nrow(checks) == 24L && all(checks$pass)
  )
  if (!all(gates)) stop("Validated downstream 526 scope artifacts failed.", call. = FALSE)
  list(
    effects = as.integer(nrow(auc_data)),
    species = as.integer(length(species)),
    experiments = as.integer(length(experiments)),
    effect_ids = ids,
    species_ids = sort(species),
    experiment_ids = sort(experiments),
    V_validation = validation,
    consistency_checks = checks
  )
}

downstream_covariance_release_check <- function(validation) {
  if (!is.data.frame(validation) || !"valid" %in% names(validation) ||
      nrow(validation) != 1L) {
    stop("Downstream covariance validation must contain one valid row.",
         call. = FALSE)
  }
  tibble::tibble(
    check = "validated downstream sampling covariance",
    pass = as.logical(validation$valid),
    detail = paste(
      "dimensions", validation$dimensions_valid,
      "order", validation$order_valid,
      "minimum eigenvalue", format(validation$minimum_eigenvalue, digits = 8)
    )
  )
}

step22_data_validation_component <- function(scope) {
  result <- scope[
    c("effects", "species", "experiments", "V_validation", "consistency_checks")
  ]
  result$consistency_checks <- result$consistency_checks |>
    dplyr::mutate(check_number = as.integer(.data$check_number))
  result
}

validated_step22_metadata <- function(config) {
  list(
    step = as.integer(config$step),
    script = "R/34_step_22_downstream_526_models.R",
    response = config$response,
    excluded_downstream_species = config$exclusions,
    upstream_modified = FALSE
  )
}

validated_phylogenetic_variation <- function(eigen_summary) {
  value <- eigen_summary$cumulative_percent[eigen_summary$PEV == "PEV3"]
  if (length(value) != 1L || !is.finite(value)) {
    stop("PEV1-PEV3 cumulative positive variation is unavailable.", call. = FALSE)
  }
  as.numeric(value)
}
