# Module: R/responses/auc_model_helpers
# Responsibility: Implement the named scientific workflow stage or shared infrastructure.
# Inputs: Function arguments and explicitly declared target dependencies.
# Outputs: Validated in-memory objects and declared workflow artifacts.
# Side effects: Limited to functions whose contracts explicitly write, copy, or render files.

# Canonical AUC model helpers extracted from the reviewed pipeline.
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
#' Validate Auc Covariance.
#'
#' @param V A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param data A data frame whose required columns are checked before use.
#' @param tolerance A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
validate_auc_covariance <- function(V, data, tolerance = 1e-8) {
  dimensions_valid <- is.matrix(V) && identical(dim(V), c(nrow(data), nrow(data)))
  if (!dimensions_valid) {
    return(tibble::tibble(
      dimensions_valid = FALSE, order_valid = FALSE,
      maximum_symmetry_error = NA_real_, minimum_eigenvalue = NA_real_,
      maximum_diagonal_se_error = NA_real_, valid = FALSE
    ))
  }
  order_valid <- !is.null(rownames(V)) && !is.null(colnames(V)) &&
    identical(rownames(V), data$effect_id) && identical(colnames(V), data$effect_id)
  symmetry_error <- max(abs(V - t(V)))
  eigen_min <- if (symmetry_error <= tolerance) {
    min(eigen(V, symmetric = TRUE, only.values = TRUE)$values)
  } else {
    NA_real_
  }
  diagonal_error <- if (all(diag(V) >= 0)) {
    max(abs(sqrt(diag(V)) - data$auc_se))
  } else {
    Inf
  }
  tibble::tibble(
    dimensions_valid = dimensions_valid, order_valid = order_valid,
    maximum_symmetry_error = symmetry_error, minimum_eigenvalue = eigen_min,
    maximum_diagonal_se_error = diagonal_error,
    valid = dimensions_valid && order_valid && symmetry_error <= tolerance &&
      is.finite(eigen_min) && eigen_min >= -tolerance && diagonal_error <= tolerance
  )
}

#' Select Auc Covariance.
#'
#' @param full_V A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param data A data frame whose required columns are checked before use.
#' @param tolerance A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param reconstruction_error A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
select_auc_covariance <- function(full_V, data, tolerance = 1e-8,
                                  reconstruction_error = "") {
  validation <- validate_auc_covariance(full_V, data, tolerance)
  if (isTRUE(validation$valid[[1]])) {
    computational_V <- Matrix::drop0(
      Matrix::Matrix(full_V, sparse = TRUE), tol = tolerance
    )
    return(list(
      V = computational_V,
      selected_v_type = paste0(
        "full emmeans covariance reconstructed from per-experiment AUC models; ",
        "entries below ", format(tolerance, scientific = TRUE),
        " stored as exact sparse zeros"
      ),
      fallback_reason = "", validation = validation
    ))
  }
  if (any(!is.finite(data$auc_se) | data$auc_se <= 0)) {
    stop("Full V failed and diagonal sampling variances are invalid.", call. = FALSE)
  }
  diagonal <- diag(data$auc_se^2)
  dimnames(diagonal) <- list(data$effect_id, data$effect_id)
  list(
    V = diagonal, selected_v_type = "diagonal sampling-variance approximation",
    fallback_reason = paste(
      c(reconstruction_error, "full covariance validation failed")[
        nzchar(c(reconstruction_error, "full covariance validation failed"))
      ],
      collapse = " | "
    ),
    validation = validation
  )
}

#' Reconstruct Auc Covariance.
#'
#' @param data A data frame whose required columns are checked before use.
#' @param auc_by_plot A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param tolerance A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
reconstruct_auc_covariance <- function(data, auc_by_plot, tolerance = 1e-8) {
  horizon_value <- unique(as.character(data$horizon))
  if (length(horizon_value) != 1L) stop("Reconstruct V one horizon at a time.", call. = FALSE)
  blocks <- list()
  ordered_ids <- character()
  block_audit <- list()
  errors <- character()

  for (experiment in levels(droplevels(data$experiment_id))) {
    target <- data |>
      dplyr::filter(as.character(.data$experiment_id) == .env$experiment) |>
      dplyr::arrange(.data$source_species)
    raw <- auc_by_plot |>
      dplyr::filter(
        as.character(.data$experiment_id) == .env$experiment,
        as.character(.data$horizon) == .env$horizon_value,
        as.character(.data$transform) == "raw", !is.na(.data$auc)
      ) |>
      dplyr::mutate(
        ODT = factor(
          as.character(.data$ODT),
          levels = c("destination_controls", "origin_controls", "warmed")
        )
      )
    captured <- auc_capture({
      model <- stats::lm(auc ~ ODT * species, data = raw)
      means <- emmeans::emmeans(model, ~ ODT | species)
      contrast <- emmeans::contrast(
        means, method = list(to_origin = c(0, -1, 1)), adjust = "none"
      )
      summary <- tibble::as_tibble(as.data.frame(summary(contrast)))
      covariance <- as.matrix(stats::vcov(contrast))
      index <- match(target$source_species, as.character(summary$species))
      if (anyNA(index)) stop("Resident species missing from reconstructed contrasts.")
      list(summary = summary[index, , drop = FALSE], covariance = covariance[index, index, drop = FALSE])
    })
    if (is.null(captured$fit)) {
      errors <- c(errors, paste(experiment, captured$error, sep = ": "))
      next
    }
    reconstructed <- captured$fit
    estimate_error <- max(abs(reconstructed$summary$estimate - target$auc_effect))
    se_error <- max(abs(reconstructed$summary$SE - target$auc_se))
    symmetry_error <- max(abs(reconstructed$covariance - t(reconstructed$covariance)))
    eigen_min <- min(eigen(reconstructed$covariance, symmetric = TRUE, only.values = TRUE)$values)
    blocks[[experiment]] <- reconstructed$covariance
    ordered_ids <- c(ordered_ids, target$effect_id)
    block_audit[[experiment]] <- tibble::tibble(
      horizon = horizon_value, experiment_id = experiment, n_effects = nrow(target),
      maximum_estimate_error = estimate_error, maximum_se_error = se_error,
      maximum_symmetry_error = symmetry_error, minimum_eigenvalue = eigen_min,
      valid = estimate_error <= tolerance && se_error <= tolerance &&
        symmetry_error <= tolerance && eigen_min >= -tolerance,
      warning = paste(captured$warnings, collapse = " | "), error = ""
    )
  }

  full_valid <- !length(errors) && length(blocks) == dplyr::n_distinct(data$experiment_id)
  V_full <- NULL
  if (full_valid) {
    assembled <- as.matrix(Matrix::bdiag(blocks))
    index <- match(data$effect_id, ordered_ids)
    full_valid <- !anyNA(index)
    if (full_valid) V_full <- assembled[index, index, drop = FALSE]
  }
  audit_blocks <- dplyr::bind_rows(block_audit)
  full_valid <- full_valid && nrow(audit_blocks) && all(audit_blocks$valid)
  if (full_valid) dimnames(V_full) <- list(data$effect_id, data$effect_id)
  selected <- select_auc_covariance(
    if (full_valid) V_full else matrix(numeric(), 0, 0), data, tolerance,
    paste(errors, collapse = " | ")
  )
  V <- selected$V
  v_type <- selected$selected_v_type
  fallback_reason <- selected$fallback_reason
  validation <- selected$validation
  full_valid <- isTRUE(validation$valid[[1]])
  dim_valid <- validation$dimensions_valid[[1]]
  symmetry_error <- validation$maximum_symmetry_error[[1]]
  eigen_min <- validation$minimum_eigenvalue[[1]]
  diagonal_error <- validation$maximum_diagonal_se_error[[1]]
  off_diagonal <- if (!is.null(V_full) && length(V_full)) {
    V_full[row(V_full) != col(V_full)]
  } else {
    V[row(V) != col(V)]
  }
  audit <- tibble::tibble(
    horizon = horizon_value, n_effects = nrow(data), n_experiments = length(blocks),
    dimensions_valid = dim_valid, maximum_symmetry_error = symmetry_error,
    minimum_eigenvalue = eigen_min, maximum_diagonal_se_error = diagonal_error,
    maximum_absolute_off_diagonal = if (length(off_diagonal)) max(abs(off_diagonal)) else 0,
    full_covariance_valid = full_valid, selected_v_type = v_type,
    fallback_reason = fallback_reason
  )
  list(V = V, audit = audit, block_audit = audit_blocks, selected_v_type = v_type)
}

#' Subset Auc v.
#'
#' @param V A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param effect_ids A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
subset_auc_v <- function(V, effect_ids) {
  if (is.null(rownames(V)) || any(!effect_ids %in% rownames(V))) {
    stop("Cannot align V to requested effect_id values.", call. = FALSE)
  }
  V[effect_ids, effect_ids, drop = FALSE]
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
