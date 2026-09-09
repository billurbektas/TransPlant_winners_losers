# Pure model helpers for the anchored 526-species meta-regression lineage.

canonical_random_structure_text <- function() {
  paste(
    "list(~1 | species_id, ~0 + W_z | species_id,",
    "~0 + T_z + P_z + TP_z | species_id); struct='GEN'"
  )
}

capture_canonical_fit <- function(expression, label) {
  captured <- auc_capture(expression)
  if (is.null(captured$fit)) {
    stop(label, " failed without simplification: ", captured$error, call. = FALSE)
  }
  captured
}

canonical_fit_signature <- function(data, formula, method,
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
    random_structure = canonical_random_structure_text()
  )
}

fit_canonical_model <- function(data, V, formula, method, label,
                                predictor_columns = character(),
                                random = canonical_downstream_config()$species_random) {
  if (!identical(rownames(V), as.character(data$effect_id)) ||
      !identical(colnames(V), as.character(data$effect_id))) {
    stop(label, ": effect rows and sampling covariance are not identically ordered.",
         call. = FALSE)
  }
  captured <- capture_canonical_fit(
    metafor::rma.mv(
      yi = auc_effect, V = V, mods = formula,
      random = random, struct = "GEN", data = data,
      method = method, sparse = TRUE
    ),
    label
  )
  c(
    captured,
    list(signature = canonical_fit_signature(
      data, formula, method, predictor_columns
    ))
  )
}

canonical_fit_statistics <- function(captured, predictor, model, method) {
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

canonical_nested_lrt <- function(smaller, larger, predictor, comparison) {
  ll0 <- stats::logLik(smaller$fit)
  ll1 <- stats::logLik(larger$fit)
  df <- attr(ll1, "df") - attr(ll0, "df")
  statistic <- 2 * (as.numeric(ll1) - as.numeric(ll0))
  tibble::tibble(
    predictor = predictor,
    comparison = comparison,
    likelihood_ratio = statistic,
    df = df,
    p_value = stats::pchisq(statistic, df, lower.tail = FALSE),
    valid_nested_comparison = is.finite(statistic) && statistic >= -1e-7 && df > 0,
    same_rows_and_V = TRUE,
    method = "ML"
  )
}

canonical_cr1 <- function(captured, predictor, model_name) {
  fit <- captured$fit
  warnings <- character()
  robust <- withCallingHandlers(
    metafor::robust(
      fit,
      cluster = fit$data$experiment_id,
      adjust = TRUE,
      clubSandwich = FALSE
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

canonical_context_covariance <- function(fit) {
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

canonical_random_diagnostics <- function(captured, predictor, model_name) {
  fit <- captured$fit
  G <- canonical_context_covariance(fit)
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

canonical_fixed_predict <- function(fit, formula, newdata) {
  design <- auc_design_matrix(formula, newdata, names(stats::coef(fit)))
  as.numeric(design %*% stats::coef(fit))
}

canonical_fixed_effect_table <- function(fit) {
  tibble::tibble(
    term = names(stats::coef(fit)),
    estimate = as.numeric(stats::coef(fit)),
    model_SE = as.numeric(fit$se),
    model_CI_low = as.numeric(fit$ci.lb),
    model_CI_high = as.numeric(fit$ci.ub),
    model_p_value = as.numeric(fit$pval)
  )
}

canonical_variance_profile_specifications <- function() {
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

profile_canonical_variance <- function(captured, specification) {
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

canonical_model_signature <- function(model) {
  if (!requireNamespace("metafor", quietly = TRUE)) {
    stop("Package 'metafor' is required for canonical model signatures.",
         call. = FALSE)
  }
  fit <- if (is.list(model) && !is.null(model$fit) &&
             inherits(model$fit, "rma")) {
    model$fit
  } else {
    model
  }
  if (inherits(fit, "canonical_model_signature")) return(fit)
  if (!inherits(fit, "rma")) {
    stop("A metafor model is required to build a canonical model signature.",
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
    class = "canonical_model_signature"
  )
}

compare_canonical_model_fit <- function(candidate, reference) {
  candidate_signature <- canonical_model_signature(candidate)
  reference_signature <- canonical_model_signature(reference)
  components <- list(
    fixed_coefficients = cbind(
      candidate = as.numeric(candidate_signature$fixed_coefficients),
      reference = as.numeric(reference_signature$fixed_coefficients)
    ),
    fixed_covariance = cbind(
      candidate = as.numeric(candidate_signature$fixed_covariance),
      reference = as.numeric(reference_signature$fixed_covariance)
    ),
    species_intercept_variance = cbind(
      candidate = as.numeric(candidate_signature$sigma2),
      reference = as.numeric(reference_signature$sigma2)
    ),
    species_W_variance = cbind(
      candidate = as.numeric(candidate_signature$tau2),
      reference = as.numeric(reference_signature$tau2)
    ),
    context_variances = cbind(
      candidate = as.numeric(candidate_signature$gamma2),
      reference = as.numeric(reference_signature$gamma2)
    ),
    context_correlations = cbind(
      candidate = c(
        as.numeric(candidate_signature$rho),
        as.numeric(candidate_signature$phi)
      ),
      reference = c(
        as.numeric(reference_signature$rho),
        as.numeric(reference_signature$phi)
      )
    ),
    log_likelihood = cbind(
      candidate = as.numeric(candidate_signature$log_likelihood),
      reference = as.numeric(reference_signature$log_likelihood)
    )
  )
  rows <- lapply(names(components), function(name) {
    values <- components[[name]]
    difference <- abs(values[, "candidate"] - values[, "reference"])
    optimizer_derived <- grepl("variance|correlation", name)
    absolute_tolerance <- if (optimizer_derived) 1e-6 else 1e-8
    relative_tolerance <- if (optimizer_derived) 1e-4 else 1e-6
    allowed <- absolute_tolerance +
      relative_tolerance * abs(values[, "reference"])
    tibble::tibble(
      check = name,
      maximum_absolute_difference = max(difference),
      maximum_allowed_difference = max(allowed),
      absolute_tolerance = absolute_tolerance,
      relative_tolerance = relative_tolerance,
      pass = all(is.finite(difference)) && all(difference <= allowed)
    )
  })
  result <- dplyr::bind_rows(rows)
  if (any(!result$pass)) {
    stop("Canonical M_both refit differs from the frozen fit: ",
         paste(result$check[!result$pass], collapse = ", "), call. = FALSE)
  }
  result
}

assemble_mboth_reference_component <- function(captured, variance_profiles) {
  list(
    fixed = canonical_fixed_effect_table(captured$fit),
    CR1 = canonical_cr1(captured, "Updated M_both", "M_both_526")$table,
    random = canonical_random_diagnostics(
      captured, "Updated M_both", "M_both_526"
    ),
    variance_profiles = variance_profiles,
    fit_statistics = canonical_fit_statistics(
      captured, "Updated M_both", "M_both_526", "REML"
    )
  )
}
