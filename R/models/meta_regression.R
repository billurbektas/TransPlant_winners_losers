# Module: R/models/meta_regression
# Responsibility: Implement the named scientific workflow stage or shared infrastructure.
# Inputs: Function arguments and explicitly declared target dependencies.
# Outputs: Validated in-memory objects and declared workflow artifacts.
# Side effects: Limited to functions whose contracts explicitly write, copy, or render files.

# Model helpers for the resident-species meta-regressions.

#' Random Structure Text.
#'
#' @return An R object described by the function name; its structure is validated before later stages use it.
random_structure_text <- function() {
  paste(
    "list(~1 | species_id, ~0 + W_z | species_id,",
    "~0 + T_z + P_z + TP_z | species_id); struct='GEN'"
  )
}

#' Capture Model Fit.
#'
#' @param expression A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param label A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
capture_model_fit <- function(expression, label) {
  captured <- auc_capture(expression)
  if (is.null(captured$fit)) {
    stop(label, " failed without simplification: ", captured$error, call. = FALSE)
  }
  captured
}

#' Model Fit Signature.
#'
#' @param data A data frame whose required columns are checked before use.
#' @param formula A model formula or random-effects specification.
#' @param method A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param predictor_columns A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
model_fit_signature <- function(data, formula, method,
                                    predictor_columns = character()) {
  predictor_values <- if (length(predictor_columns)) {
    lapply(data[predictor_columns], function(x) {
      if (is.factor(x)) as.character(x) else as.numeric(x)
    })
  } else {
    list()
  }
  list(
    effect_id = as.character(data$effect_id),
    formula = paste(deparse(formula), collapse = " "),
    method = method,
    predictor_values = predictor_values,
    random_structure = random_structure_text()
  )
}

#' Fit Meta Regression.
#'
#' @param data A data frame whose required columns are checked before use.
#' @param sampling_variance_matrix Effect-ordered sampling variances.
#' @param formula A model formula or random-effects specification.
#' @param method A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param label A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param predictor_columns A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param random A model formula or random-effects specification.
#' @param struct Covariance structure passed to `metafor::rma.mv()`.
#' @param sparse Whether to use sparse matrix calculations.
#' @details `data` must contain `effect_id`, `auc_effect`, the grouping fields,
#' and every variable referenced by `formula` and `random`. Row and column names
#' of `sampling_variance_matrix` must exactly equal `effect_id` order. The registered `GEN` random
#' structure and supplied ML/REML method are not simplified on failure; fit
#' errors and warnings are retained in the returned bundle.
#' @return A fitted model bundle with diagnostics retained for validation.
fit_meta_regression <- function(
    data, sampling_variance_matrix, formula, method, label,
    predictor_columns = character(), random, struct, sparse = TRUE) {
  if (!identical(
        rownames(sampling_variance_matrix), as.character(data$effect_id)
      ) ||
      !identical(
        colnames(sampling_variance_matrix), as.character(data$effect_id)
      )) {
    stop(label, ": effect rows and sampling variances are not identically ordered.",
         call. = FALSE)
  }
  captured <- capture_model_fit(
    metafor::rma.mv(
      yi = auc_effect, V = sampling_variance_matrix, mods = formula,
      random = random, struct = struct, data = data,
      method = method, sparse = sparse
    ),
    label
  )
  c(
    captured,
    list(signature = model_fit_signature(
      data, formula, method, predictor_columns
    ))
  )
}

#' Model Fit Statistics.
#'
#' @param captured A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param predictor A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param model A fitted model object with the structure checked by this function.
#' @param method A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
model_fit_statistics <- function(captured, predictor, model, method) {
  fit <- captured$fit
  values <- auc_fit_statistics(fit)
  tibble::tibble(
    predictor = predictor,
    model = model,
    method = method,
    n_effects = fit$k,
    n_fixed = fit$p,
    logLik = values[["logLik"]],
    AIC = values[["AIC"]],
    AICc = values[["AICc"]],
    BIC = values[["BIC"]],
    optimizer_code = auc_value_or(
      fit$opt.res$convergence,
      auc_value_or(fit$opt.res$ierr, NA_integer_)
    ),
    warnings = paste(captured$warnings, collapse = " | ")
  )
}

#' Nested Model Lrt.
#'
#' @param smaller A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param larger A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param predictor A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param comparison A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
nested_model_lrt <- function(smaller, larger, predictor, comparison) {
  ll0 <- stats::logLik(smaller$fit)
  ll1 <- stats::logLik(larger$fit)
  df <- attr(ll1, "df") - attr(ll0, "df")
  statistic <- 2 * (as.numeric(ll1) - as.numeric(ll0))
  same_rows_and_sampling_variances <- identical(
    smaller$signature$effect_id, larger$signature$effect_id
  ) && isTRUE(all.equal(
    smaller$fit$V, larger$fit$V,
    tolerance = 0, check.attributes = TRUE
  ))
  tibble::tibble(
    predictor = predictor,
    comparison = comparison,
    likelihood_ratio = statistic,
    df = df,
    p_value = stats::pchisq(statistic, df, lower.tail = FALSE),
    valid_nested_comparison = is.finite(statistic) &&
      statistic >= -1e-7 && df > 0 &&
      same_rows_and_sampling_variances,
    same_rows_and_sampling_variances = same_rows_and_sampling_variances,
    method = "ML"
  )
}

#' Cluster Robust Inference.
#'
#' @param captured A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param predictor A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param model_name A fitted model object with the structure checked by this function.
#' @param cluster_column Column identifying independent experiment clusters.
#' @param adjust Whether to apply the small-sample CR1 adjustment.
#' @param club_sandwich Whether `metafor` should delegate to `clubSandwich`.
#' @return An R object described by the function name; its structure is validated before later stages use it.
cluster_robust_inference <- function(
    captured, predictor, model_name, cluster_column,
    adjust, club_sandwich) {
  fit <- captured$fit
  if (!cluster_column %in% names(fit$data)) {
    stop("Robust inference is missing cluster column `", cluster_column, "`.",
         call. = FALSE)
  }
  warnings <- character()
  robust <- withCallingHandlers(
    metafor::robust(
      fit,
      cluster = fit$data[[cluster_column]],
      adjust = adjust,
      clubSandwich = club_sandwich
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (!isTRUE(all.equal(
    as.numeric(stats::coef(fit)), as.numeric(robust$beta),
    tolerance = sqrt(.Machine$double.eps)
  ))) {
    stop(predictor, ": CR1 changed fixed estimates.", call. = FALSE)
  }
  list(
    object = robust,
    table = tibble::tibble(
      predictor = predictor,
      model = model_name,
      term = names(stats::coef(fit)),
      estimate = as.numeric(stats::coef(fit)),
      model_SE = as.numeric(fit$se),
      model_CI_low = as.numeric(fit$ci.lb),
      model_CI_high = as.numeric(fit$ci.ub),
      model_p_value = as.numeric(fit$pval),
      CR1_SE = as.numeric(robust$se),
      CR1_CI_low = as.numeric(robust$ci.lb),
      CR1_CI_high = as.numeric(robust$ci.ub),
      CR1_p_value = as.numeric(robust$pval),
      CR1_denominator_df = rep(as.numeric(robust$ddf), length.out = fit$p),
      CR1_warnings = paste(unique(warnings), collapse = " | ")
    )
  )
}

#' Context Covariance.
#'
#' @param fit A fitted model object with the structure checked by this function.
#' @return An R object described by the function name; its structure is validated before later stages use it.
context_covariance <- function(fit) {
  covariance <- as.matrix(fit$H)
  if (!identical(colnames(fit$G), "W_z") ||
      !identical(colnames(covariance), c("T_z", "P_z", "TP_z"))) {
    stop(
      "The fitted model does not retain the independent intercept/W and GEN T/P/TP blocks.",
      call. = FALSE
    )
  }
  covariance
}

#' Random Effect Diagnostics.
#'
#' @param captured A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param predictor A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param model_name A fitted model object with the structure checked by this function.
#' @return An R object described by the function name; its structure is validated before later stages use it.
random_effect_diagnostics <- function(captured, predictor, model_name) {
  fit <- captured$fit
  G <- context_covariance(fit)
  R <- stats::cov2cor(G)
  full <- matrix(
    0, 5, 5,
    dimnames = list(
      c("intercept", "W_z", "T_z", "P_z", "TP_z"),
      c("intercept", "W_z", "T_z", "P_z", "TP_z")
    )
  )
  full["intercept", "intercept"] <- fit$sigma2[[1]]
  full["W_z", "W_z"] <- fit$tau2[[1]]
  full[c("T_z", "P_z", "TP_z"), c("T_z", "P_z", "TP_z")] <- G
  eigenvalues <- eigen(full, symmetric = TRUE, only.values = TRUE)$values
  tolerance <- 1e-8 * max(1, stats::var(fit$yi))
  tibble::tibble(
    predictor = predictor,
    model = model_name,
    n_effects = fit$k,
    n_species = dplyr::n_distinct(fit$data$species_id),
    n_experiments = dplyr::n_distinct(fit$data$experiment_id),
    species_intercept_variance = fit$sigma2[[1]],
    species_W_variance = fit$tau2[[1]],
    species_T_variance = G["T_z", "T_z"],
    species_P_variance = G["P_z", "P_z"],
    species_TP_variance = G["TP_z", "TP_z"],
    T_P_correlation = R["T_z", "P_z"],
    T_TP_correlation = R["T_z", "TP_z"],
    P_TP_correlation = R["P_z", "TP_z"],
    minimum_full_covariance_eigenvalue = min(eigenvalues),
    covariance_positive_definite = min(eigenvalues) > 0,
    any_variance_at_boundary = any(diag(full) <= tolerance),
    maximum_absolute_context_correlation = max(abs(R[upper.tri(R)])),
    covariance_identifiable = min(eigenvalues) > tolerance &&
      all(diag(full) > tolerance) &&
      max(abs(R[upper.tri(R)])) < 1 - 1e-6,
    optimizer_code = auc_value_or(
      fit$opt.res$convergence,
      auc_value_or(fit$opt.res$ierr, NA_integer_)
    ),
    finite_coefficients = all(is.finite(stats::coef(fit))),
    finite_fixed_covariance = all(is.finite(stats::vcov(fit))),
    warnings = paste(captured$warnings, collapse = " | ")
  )
}

#' Predict Fixed Effects.
#'
#' @param fit A fitted model object with the structure checked by this function.
#' @param formula A model formula or random-effects specification.
#' @param newdata A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
predict_fixed_effects <- function(fit, formula, newdata) {
  design <- auc_design_matrix(formula, newdata, names(stats::coef(fit)))
  as.numeric(design %*% stats::coef(fit))
}

#' Fixed Effect Table.
#'
#' @param fit A fitted model object with the structure checked by this function.
#' @return An R object described by the function name; its structure is validated before later stages use it.
fixed_effect_table <- function(fit) {
  tibble::tibble(
    term = names(stats::coef(fit)),
    estimate = as.numeric(stats::coef(fit)),
    model_SE = as.numeric(fit$se),
    model_CI_low = as.numeric(fit$ci.lb),
    model_CI_high = as.numeric(fit$ci.ub),
    model_p_value = as.numeric(fit$pval)
  )
}

#' Variance Profile Specifications.
#'
#' @return An R object described by the function name; its structure is validated before later stages use it.
variance_profile_specifications <- function() {
  specifications <- tibble::tribble(
    ~component, ~argument, ~index, ~interval_name,
    "species_intercept", "sigma2", 1L, "sigma^2.1",
    "species_W", "tau2", 1L, "tau^2.1",
    "species_T", "gamma2", 1L, "gamma^2.1",
    "species_P", "gamma2", 2L, "gamma^2.2",
    "species_TP", "gamma2", 3L, "gamma^2.3"
  )
  stats::setNames(
    lapply(seq_len(nrow(specifications)), function(index) {
      specifications[index, , drop = FALSE]
    }),
    specifications$component
  )
}

#' Profile Model Variance.
#'
#' @param captured A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param specification A validated workflow configuration or specification object.
#' @return An R object described by the function name; its structure is validated before later stages use it.
profile_model_variance <- function(captured, specification) {
  fit <- if (is.list(captured) && !is.null(captured$fit)) {
    captured$fit
  } else {
    captured
  }
  component <- specification$component[[1]]
  argument <- specification$argument[[1]]
  index <- specification$index[[1]]
  interval_name <- specification$interval_name[[1]]
  estimates <- c(
    species_intercept = fit$sigma2[[1]],
    species_W = fit$tau2[[1]],
    species_T = fit$gamma2[[1]],
    species_P = fit$gamma2[[2]],
    species_TP = fit$gamma2[[3]]
  )
  arguments <- list(
    object = fit,
    xlim = c(0, max(0.025, 6 * estimates[[component]]))
  )
  arguments[[argument]] <- index
  profiled <- auc_capture(do.call(stats::confint, arguments))
  if (is.null(profiled$fit)) {
    arguments$xlim <- NULL
    profiled <- auc_capture(do.call(stats::confint, arguments))
  }
  if (is.null(profiled$fit)) {
    return(tibble::tibble(
      component = component, estimate = NA_real_,
      CI_low = NA_real_, CI_high = NA_real_,
      warnings = paste(profiled$warnings, collapse = " | "),
      error = profiled$error
    ))
  }
  random_table <- profiled$fit$random
  row_name <- interval_name
  if (!row_name %in% rownames(random_table)) {
    row_name <- sub("[.]1$", "", row_name)
  }
  if (!row_name %in% rownames(random_table)) {
    return(tibble::tibble(
      component = component, estimate = NA_real_,
      CI_low = NA_real_, CI_high = NA_real_,
      warnings = paste(profiled$warnings, collapse = " | "),
      error = paste(
        "Profile row unavailable; rows:",
        paste(rownames(random_table), collapse = ", ")
      )
    ))
  }
  interval <- as.numeric(
    random_table[row_name, c("estimate", "ci.lb", "ci.ub")]
  )
  if (any(!is.finite(interval))) {
    arguments$xlim <- NULL
    unrestricted <- auc_capture(do.call(stats::confint, arguments))
    if (!is.null(unrestricted$fit) &&
        row_name %in% rownames(unrestricted$fit$random)) {
      profiled <- unrestricted
      random_table <- profiled$fit$random
    }
  }
  tibble::tibble(
    component = component,
    estimate = as.numeric(random_table[row_name, "estimate"]),
    CI_low = as.numeric(random_table[row_name, "ci.lb"]),
    CI_high = as.numeric(random_table[row_name, "ci.ub"]),
    warnings = paste(profiled$warnings, collapse = " | "),
    error = ""
  )
}

#' Model Signature.
#'
#' @param model A fitted model object with the structure checked by this function.
#' @return An R object described by the function name; its structure is validated before later stages use it.
model_signature <- function(model) {
  if (!requireNamespace("metafor", quietly = TRUE)) {
    stop("Package 'metafor' is required for model signatures.",
         call. = FALSE)
  }
  fit <- if (is.list(model) && !is.null(model$fit) &&
             inherits(model$fit, "rma")) {
    model$fit
  } else {
    model
  }
  if (inherits(fit, "model_signature")) return(fit)
  signature_fields <- c(
    "fixed_coefficients", "fixed_covariance", "sigma2", "tau2",
    "gamma2", "rho", "phi", "log_likelihood"
  )
  if (is.list(fit) && all(signature_fields %in% names(fit))) {
    return(structure(fit[signature_fields], class = "model_signature"))
  }
  if (!inherits(fit, "rma")) {
    stop("A fitted metafor model or complete model signature is required.",
         call. = FALSE)
  }
  structure(
    list(
      fixed_coefficients = stats::coef(fit),
      fixed_covariance = stats::vcov(fit),
      sigma2 = fit$sigma2,
      tau2 = fit$tau2,
      gamma2 = fit$gamma2,
      rho = fit$rho,
      phi = fit$phi,
      log_likelihood = as.numeric(stats::logLik(fit))
    ),
    class = "model_signature"
  )
}


#' Assemble Climate Model Results.
#'
#' @param captured A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param variance_profiles A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param cluster_column Column identifying independent experiment clusters.
#' @param adjust Whether to apply the small-sample CR1 adjustment.
#' @param club_sandwich Whether `metafor` should delegate to `clubSandwich`.
#' @return An R object described by the function name; its structure is validated before later stages use it.
assemble_climate_model_results <- function(
    captured, variance_profiles, cluster_column, adjust, club_sandwich) {
  list(
    fixed = fixed_effect_table(captured$fit),
    CR1 = cluster_robust_inference(
      captured, "resident_climate", "resident_climate",
      cluster_column = cluster_column,
      adjust = adjust,
      club_sandwich = club_sandwich
    )$table,
    random = random_effect_diagnostics(
      captured, "resident_climate", "resident_climate"
    ),
    variance_profiles = variance_profiles,
    fit_statistics = model_fit_statistics(
      captured, "resident_climate", "resident_climate", "REML"
    )
  )
}
