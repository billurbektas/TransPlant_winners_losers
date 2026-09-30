# Module: R/responses/auc_model_helpers
# Responsibility: Implement the named scientific workflow stage or shared infrastructure.
# Inputs: Function arguments and explicitly declared target dependencies.
# Outputs: Validated in-memory objects and declared workflow artifacts.
# Side effects: Limited to functions whose contracts explicitly write, copy, or render files.

# Current AUC model helpers extracted from the reviewed pipeline.
# These functions preserve the validated response-model implementations.

#' Auc Capture.
#'
#' @param expression A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
auc_capture <- function(expression) {
  warnings <- character()
  error <- ""
  fit <- tryCatch(
    withCallingHandlers(
      expression,
      warning = function(condition) {
        warnings <<- c(warnings, conditionMessage(condition))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(condition) {
      error <<- conditionMessage(condition)
      NULL
    }
  )
  list(fit = fit, warnings = unique(warnings), error = error)
}

#' Auc Formula Text.
#'
#' @param formula A model formula or random-effects specification.
#' @return An R object described by the function name; its structure is validated before later stages use it.
auc_formula_text <- function(formula) paste(deparse(formula), collapse = " ")

#' Auc Value or.
#'
#' @param x A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param fallback A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
auc_value_or <- function(x, fallback) {
  if (is.null(x) || !length(x)) fallback else x
}
#' Scale Resident Auc Data.
#'
#' @param data A data frame whose required columns are checked before use.
#' @return An R object described by the function name; its structure is validated before later stages use it.
scale_resident_auc_data <- function(data) {
  experiments <- data |>
    dplyr::select(dplyr::all_of(c(
      "experiment_id", "gradient_id", "T", "P", "W", "duration",
      "log_plot_area"
    ))) |>
    dplyr::distinct()
  assert_unique(experiments, "experiment_id", "experiment-level predictors")
  variables <- c(
    T_z = "T", P_z = "P", W_z = "W", duration_z = "duration",
    log_plot_area_z = "log_plot_area"
  )
  scaling <- purrr::imap_dfr(variables, function(source, scaled) {
    values <- experiments[[source]]
    tibble::tibble(
      scaled_field = scaled, source_field = source,
      center = mean(values), scale = stats::sd(values),
      weighting = "distinct experiments", n = length(values)
    )
  })
  if (any(!is.finite(scaling$scale) | scaling$scale <= 0)) {
    stop("All experiment-level variables need a positive finite scale.", call. = FALSE)
  }
  for (i in seq_len(nrow(scaling))) {
    data[[scaling$scaled_field[[i]]]] <-
      (data[[scaling$source_field[[i]]]] - scaling$center[[i]]) /
      scaling$scale[[i]]
  }
  data$species_id <- factor(data$species_id)
  data$experiment_id <- factor(data$experiment_id)
  data$gradient_id <- factor(data$gradient_id)
  list(data = data, scaling = scaling)
}
#' Build a diagonal sampling-variance matrix.
#'
#' The AUC contrast model supplies one standard error per effect. The resident
#' and family meta-regressions treat sampling errors for different effects as
#' independent, so the model input contains squared standard errors on the
#' diagonal and exact zeros elsewhere.
#'
#' @param data Data frame containing unique effect identifiers and standard
#'   errors.
#' @param se_column Name of the standard-error column.
#' @param id_column Name of the effect-identifier column.
#' @return A sparse diagonal matrix ordered exactly like `data`.
build_sampling_variance_matrix <- function(
    data, se_column, id_column = "effect_id") {
  required <- c(id_column, se_column)
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(
      "Sampling-variance data lack: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  effect_ids <- as.character(data[[id_column]])
  standard_errors <- as.numeric(data[[se_column]])
  if (anyNA(effect_ids) || any(!nzchar(effect_ids)) || anyDuplicated(effect_ids)) {
    stop("Sampling-variance effect identifiers must be complete and unique.",
         call. = FALSE)
  }
  if (any(!is.finite(standard_errors) | standard_errors <= 0)) {
    stop("Sampling standard errors must be positive and finite.", call. = FALSE)
  }
  result <- Matrix::Diagonal(x = standard_errors^2)
  dimnames(result) <- list(effect_ids, effect_ids)
  result
}

#' Validate a diagonal sampling-variance matrix.
#'
#' @param sampling_variance_matrix Matrix to validate.
#' @param data Data frame containing the matching effects and standard errors.
#' @param se_column Name of the standard-error column.
#' @param id_column Name of the effect-identifier column.
#' @param tolerance Numerical comparison tolerance.
#' @return One-row validation table.
validate_sampling_variance_matrix <- function(
    sampling_variance_matrix, data, se_column, id_column = "effect_id",
    tolerance = 1e-12) {
  expected_dimension <- c(nrow(data), nrow(data))
  dimensions_valid <- identical(dim(sampling_variance_matrix), expected_dimension)
  if (!dimensions_valid) {
    return(tibble::tibble(
      dimensions_valid = FALSE,
      order_valid = FALSE,
      maximum_diagonal_variance_error = NA_real_,
      maximum_absolute_off_diagonal = NA_real_,
      minimum_sampling_variance = NA_real_,
      valid = FALSE
    ))
  }
  effect_ids <- as.character(data[[id_column]])
  diagonal <- as.numeric(Matrix::diag(sampling_variance_matrix))
  expected_diagonal <- as.numeric(data[[se_column]])^2
  off_diagonal <- sampling_variance_matrix - Matrix::Diagonal(x = diagonal)
  maximum_off_diagonal <- if (Matrix::nnzero(off_diagonal) == 0L) {
    0
  } else {
    max(abs(as.numeric(off_diagonal)))
  }
  order_valid <- identical(rownames(sampling_variance_matrix), effect_ids) &&
    identical(colnames(sampling_variance_matrix), effect_ids)
  diagonal_error <- max(abs(diagonal - expected_diagonal))
  minimum_variance <- min(diagonal)
  tibble::tibble(
    dimensions_valid = dimensions_valid,
    order_valid = order_valid,
    maximum_diagonal_variance_error = diagonal_error,
    maximum_absolute_off_diagonal = maximum_off_diagonal,
    minimum_sampling_variance = minimum_variance,
    valid = dimensions_valid && order_valid &&
      is.finite(diagonal_error) && diagonal_error <= tolerance &&
      is.finite(maximum_off_diagonal) && maximum_off_diagonal <= tolerance &&
      is.finite(minimum_variance) && minimum_variance > 0
  )
}

#' Subset a sampling-variance matrix by effect identifier.
#'
#' @param sampling_variance_matrix Effect-ordered sampling-variance matrix.
#' @param effect_ids Identifiers to retain, in the required order.
#' @return The matching square submatrix.
subset_sampling_variance_matrix <- function(
    sampling_variance_matrix, effect_ids) {
  if (is.null(rownames(sampling_variance_matrix)) ||
      any(!effect_ids %in% rownames(sampling_variance_matrix))) {
    stop("Cannot align sampling variances to requested effect identifiers.",
         call. = FALSE)
  }
  sampling_variance_matrix[effect_ids, effect_ids, drop = FALSE]
}
#' Auc Fit Statistics.
#'
#' @param fit A fitted model object with the structure checked by this function.
#' @return An R object described by the function name; its structure is validated before later stages use it.
auc_fit_statistics <- function(fit) {
  if (is.null(fit)) {
    return(c(logLik = NA_real_, AIC = NA_real_, AICc = NA_real_, BIC = NA_real_))
  }
  values <- tryCatch(metafor::fitstats(fit), error = function(condition) NULL)
  if (is.null(values)) {
    return(c(logLik = NA_real_, AIC = NA_real_, AICc = NA_real_, BIC = NA_real_))
  }
  as.numeric(values[c("logLik", "AIC", "AICc", "BIC"), 1]) |>
    stats::setNames(c("logLik", "AIC", "AICc", "BIC"))
}
#' Auc Design Matrix.
#'
#' @param formula A model formula or random-effects specification.
#' @param newdata A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param coefficient_names A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
auc_design_matrix <- function(formula, newdata, coefficient_names) {
  design <- stats::model.matrix(formula, data = newdata)
  if ("intrcpt" %in% coefficient_names && "(Intercept)" %in% colnames(design)) {
    colnames(design)[colnames(design) == "(Intercept)"] <- "intrcpt"
  }
  missing <- setdiff(coefficient_names, colnames(design))
  if (length(missing)) {
    stop("Prediction matrix is missing coefficients: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  design[, coefficient_names, drop = FALSE]
}
#' Auc Support Mask.
#'
#' @param experiments A data frame whose required columns are checked before use.
#' @param grid A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
auc_support_mask <- function(experiments, grid) {
  observed <- as.matrix(experiments[, c("T_z", "P_z", "W_z")])
  candidate <- as.matrix(grid[, c("T_z", "P_z", "W_z")])
  captured <- auc_capture({
    hull <- geometry::convhulln(observed, options = "Qt")
    geometry::inhulln(hull, candidate)
  })
  if (is.null(captured$fit)) {
    return(list(
      grid = experiments |>
        dplyr::distinct(.data$experiment_id, .data$T_z, .data$P_z, .data$W_z) |>
        dplyr::mutate(supported = TRUE, support_method = "observed experiment points"),
      diagnostic = captured$error
    ))
  }
  list(
    grid = grid |>
      dplyr::mutate(
        supported = as.logical(captured$fit),
        support_method = "three-dimensional convex hull"
      ),
    diagnostic = paste(captured$warnings, collapse = " | ")
  )
}

#' Build Auc Climate Grid.
#'
#' @param data A data frame whose required columns are checked before use.
#' @param points A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
build_auc_climate_grid <- function(data, points = 60L) {
  experiments <- data |>
    dplyr::distinct(.data$experiment_id, .data$T_z, .data$P_z, .data$W_z)
  warming <- stats::quantile(
    experiments$W_z, probs = c(0.25, 0.50, 0.75), names = FALSE
  )
  warming_labels <- c("low (25th percentile)", "median", "high (75th percentile)")
  grid <- tidyr::expand_grid(
    T_z = seq(min(experiments$T_z), max(experiments$T_z), length.out = points),
    P_z = seq(min(experiments$P_z), max(experiments$P_z), length.out = points),
    W_z = warming
  ) |>
    dplyr::mutate(
      warming_level = factor(
        .env$warming_labels[match(.data$W_z, .env$warming)],
        levels = .env$warming_labels
      ),
      duration_z = 0, log_plot_area_z = 0,
      climate_context = factor("cold_dry", levels = levels(data$climate_context))
    )
  auc_support_mask(experiments, grid)
}
