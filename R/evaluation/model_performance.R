# Module: R/evaluation/model_performance
# Responsibility: Calculate the four model-performance summaries reported in
# Appendix Figure S7.
# Inputs: Final REML fits, their analysis data, diagonal sampling variances, and
# leave-one-gradient-out REML fits.
# Outputs: Equal-gradient marginal R2, conditional R2, predictive R2, and RMSE.
# Side effects: None. Model fitting occurs only when the corresponding targets
# are explicitly run.

#' Model display order
#'
#' @return A table defining the seven models shown in Appendix Figure S7.
model_performance_key <- function() {
  tibble::tribble(
    ~model_key,   ~model_display,                   ~model_order,
    "Climate",   "Climate-only",                            1L,
    "LDMC",      "LDMC",                                    2L,
    "LMA",       "LMA",                                     3L,
    "Height",    "Vegetative height",                       4L,
    "Graminoid", "Graminoid status",                        5L,
    "Woodiness", "Woodiness",                               6L,
    "Phylogeny", "Phylogeny (PEV1-PEV3)",                  7L
  )
}

#' Convert model keys to publication labels
#'
#' @param model_key Character model keys.
#' @return Character display labels in the same order.
model_performance_display <- function(model_key) {
  key <- model_performance_key()
  display <- key$model_display[match(model_key, key$model_key)]
  if (anyNA(display)) {
    stop(
      "Unknown model-performance key: ",
      paste(unique(model_key[is.na(display)]), collapse = ", "),
      call. = FALSE
    )
  }
  display
}

#' Collapse a formula to one line for audits
#'
#' @param formula An R formula.
#' @return A one-line character representation.
model_performance_formula_text <- function(formula) {
  paste(deparse(formula), collapse = " ")
}

#' Assemble the fitted species random-effect covariance
#'
#' @param fit A fitted `metafor::rma.mv()` model.
#' @return A 5 by 5 covariance matrix ordered as intercept, W, T, P, and T:P.
model_performance_random_block <- function(fit) {
  context <- context_covariance(fit)
  block <- matrix(
    0, nrow = 5L, ncol = 5L,
    dimnames = list(
      c("intercept", "W_z", "T_z", "P_z", "TP_z"),
      c("intercept", "W_z", "T_z", "P_z", "TP_z")
    )
  )
  block["intercept", "intercept"] <- as.numeric(fit$sigma2[[1]])
  block["W_z", "W_z"] <- as.numeric(fit$tau2[[1]])
  block[c("T_z", "P_z", "TP_z"), c("T_z", "P_z", "TP_z")] <- context
  block
}

#' Random-effect design matrix on the response scale
#'
#' @param data Model data containing W, T, P, and T:P.
#' @return A numeric matrix with one row per effect.
model_performance_random_design <- function(data) {
  required <- c("W_z", "T_z", "P_z", "TP_z")
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(
      "Random-effect data lack: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  cbind(
    intercept = 1,
    W_z = as.numeric(data$W_z),
    T_z = as.numeric(data$T_z),
    P_z = as.numeric(data$P_z),
    TP_z = as.numeric(data$TP_z)
  )
}

#' Expected among-species variance for each effect row
#'
#' @param fit A fitted `metafor::rma.mv()` model.
#' @param data Model data in the same variable parameterisation as `fit`.
#' @return A numeric vector of response-scale random-effect variances.
model_performance_random_variance <- function(fit, data) {
  design <- model_performance_random_design(data)
  covariance <- model_performance_random_block(fit)
  rowSums((design %*% covariance) * design)
}

#' Give every geographic gradient equal influence
#'
#' @param data Model data containing `gradient_id`.
#' @return Row weights that sum to one and give each gradient equal total weight.
equal_gradient_weights <- function(data) {
  if (!"gradient_id" %in% names(data)) {
    stop("Performance data require gradient_id.", call. = FALSE)
  }
  gradient <- as.character(data$gradient_id)
  counts <- table(gradient)
  1 / (length(counts) * as.numeric(counts[gradient]))
}

#' Calculate marginal and conditional R2 for a final REML fit
#'
#' @param fit A final REML `metafor::rma.mv()` fit.
#' @param model_key Model key used in Appendix Figure S7.
#' @param candidate_role `standalone`, `matched_baseline`, or `full`.
#' @return One row containing equal-gradient marginal and conditional R2.
model_performance_fit_r2 <- function(fit, model_key, candidate_role) {
  data <- fit$data
  if (!is.data.frame(data) || nrow(data) != fit$k) {
    stop("A fitted rma.mv model with retained data is required.", call. = FALSE)
  }
  weights <- equal_gradient_weights(data)
  fixed <- as.numeric(fit$X %*% stats::coef(fit))
  fixed_mean <- sum(weights * fixed)
  fixed_variance <- sum(weights * (fixed - fixed_mean)^2)
  species_variance <- sum(
    weights * model_performance_random_variance(fit, data)
  )
  sampling_variance <- sum(weights * as.numeric(Matrix::diag(fit$V)))
  total <- fixed_variance + species_variance + sampling_variance
  if (!is.finite(total) || total <= 0) {
    stop("The response-scale variance partition is invalid.", call. = FALSE)
  }
  tibble::tibble(
    model_key = model_key,
    model = model_performance_display(model_key),
    candidate_role = candidate_role,
    weighting = "equal_gradient",
    n_effects = nrow(data),
    n_species = dplyr::n_distinct(data$species_id),
    n_gradients = dplyr::n_distinct(data$gradient_id),
    fixed_effect_variance = fixed_variance,
    species_random_variance = species_variance,
    sampling_variance = sampling_variance,
    total_response_variance = total,
    marginal_R2 = fixed_variance / total,
    conditional_R2 = (fixed_variance + species_variance) / total,
    R2_definition = paste(
      "equal-gradient response-scale variance partition; denominator is",
      "fixed-effect + species-random + sampling variance"
    )
  )
}

#' Assemble final-fit R2 for climate and predictor models
#'
#' @param climate_model Final resident climate REML fit bundle.
#' @param ... Predictor model-pair bundles.
#' @return R2 rows for the standalone climate model and matched/full predictors.
assemble_model_fit_r2 <- function(climate_model, ...) {
  model_pairs <- list(...)
  rows <- list(
    model_performance_fit_r2(climate_model$fit, "Climate", "standalone")
  )
  for (model_pair in model_pairs) {
    rows[[length(rows) + 1L]] <- model_performance_fit_r2(
      model_pair$fits$climate_REML$fit,
      model_pair$key,
      "matched_baseline"
    )
    rows[[length(rows) + 1L]] <- model_performance_fit_r2(
      model_pair$fits$full_REML$fit,
      model_pair$key,
      "full"
    )
  }
  dplyr::bind_rows(rows)
}

#' Split one model dataset into training and held-out gradients
#'
#' @param data Model data.
#' @param sampling_variance_matrix Diagonal, effect-ordered sampling variances.
#' @param heldout_gradient Gradient omitted from fitting.
#' @return Training/test data and matching sampling-variance submatrices.
model_performance_partition <- function(
    data, sampling_variance_matrix, heldout_gradient) {
  required <- c(
    "effect_id", "gradient_id", "experiment_id", "species_id",
    "auc_effect", "T_z", "P_z", "W_z", "TP_z"
  )
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop("Performance data lack: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  effect_ids <- as.character(data$effect_id)
  if (!identical(rownames(sampling_variance_matrix), effect_ids) ||
      !identical(colnames(sampling_variance_matrix), effect_ids)) {
    stop("Performance data and sampling variances are not identically ordered.",
         call. = FALSE)
  }
  test_index <- as.character(data$gradient_id) == heldout_gradient
  if (!any(test_index) || all(test_index)) {
    stop("Each held-out gradient needs non-empty train and test data.", call. = FALSE)
  }
  train <- data[!test_index, , drop = FALSE]
  test <- data[test_index, , drop = FALSE]
  train_ids <- as.character(train$effect_id)
  test_ids <- as.character(test$effect_id)
  disjoint <- !any(train_ids %in% test_ids)
  if (!disjoint) stop("Training and held-out effects overlap.", call. = FALSE)
  list(
    train = train,
    test = test,
    train_V = subset_sampling_variance_matrix(
      sampling_variance_matrix, train_ids
    ),
    test_V = subset_sampling_variance_matrix(
      sampling_variance_matrix, test_ids
    ),
    effect_partitions_disjoint = disjoint
  )
}

#' Fit one leave-gradient-out REML model
#'
#' @param partition Training/test partition.
#' @param formula Fixed-effect formula.
#' @param label Human-readable fit label.
#' @param predictor_columns Predictor columns used in the full model.
#' @param random Random intercept and slope formulas.
#' @param struct Random-slope covariance structure.
#' @param method Estimation method; Appendix Figure S7 requires `"REML"`.
#' @param sparse Whether `metafor::rma.mv()` should use sparse matrices.
#' @return A fitted model bundle with rank, convergence, and finite-value checks.
fit_model_performance_fold <- function(
    partition, formula, label, predictor_columns = character(), random, struct,
    method = "REML", sparse = TRUE) {
  if (!identical(method, "REML")) {
    stop("Appendix Figure S7 fold models require REML.", call. = FALSE)
  }
  design <- stats::model.matrix(formula, data = partition$train)
  full_rank <- qr(design)$rank == ncol(design)
  if (!full_rank) {
    stop(label, " has a rank-deficient training design.", call. = FALSE)
  }
  captured <- fit_meta_regression(
    partition$train,
    partition$train_V,
    formula,
    method = method,
    label = label,
    predictor_columns = predictor_columns,
    random = random,
    struct = struct,
    sparse = sparse
  )
  fit <- captured$fit
  convergence <- auc_value_or(
    fit$opt.res$convergence,
    auc_value_or(fit$opt.res$ierr, NA_integer_)
  )
  finite <- all(is.finite(stats::coef(fit))) &&
    all(is.finite(stats::vcov(fit))) &&
    all(is.finite(c(fit$sigma2, fit$tau2, fit$gamma2)))
  captured$performance_audit <- tibble::tibble(
    full_rank = full_rank,
    convergence_code = as.integer(convergence),
    finite_parameters = finite,
    warnings = paste(captured$warnings, collapse = " | "),
    pass = full_rank && identical(as.integer(convergence), 0L) && finite
  )
  captured
}

#' Audit one cross-validation fit
#'
#' @param captured Fitted fold model.
#' @param partition Training/test partition.
#' @param heldout_gradient Omitted gradient.
#' @param model_key Model key.
#' @param candidate_role Candidate role.
#' @param fit_id Repeated fit identifier.
#' @param formula Fixed-effect formula.
#' @return One audit row.
model_performance_fit_audit_row <- function(
    captured, partition, heldout_gradient, model_key, candidate_role,
    fit_id, formula) {
  audit <- captured$performance_audit
  tibble::tibble(
    heldout_gradient = heldout_gradient,
    model_key = model_key,
    model = model_performance_display(model_key),
    candidate_role = candidate_role,
    fit_id = fit_id,
    formula = model_performance_formula_text(formula),
    n_train = nrow(partition$train),
    n_test = nrow(partition$test),
    effect_partitions_disjoint = partition$effect_partitions_disjoint,
    full_rank = audit$full_rank,
    convergence_code = audit$convergence_code,
    finite_parameters = audit$finite_parameters,
    warnings = audit$warnings,
    pass = audit$pass && partition$effect_partitions_disjoint
  )
}

#' Summarise fixed-effect predictions for one held-out gradient
#'
#' @param captured Fitted fold model.
#' @param partition Training/test partition.
#' @param formula Fixed-effect formula.
#' @param heldout_gradient Omitted gradient.
#' @param model_key Model key.
#' @param candidate_role Candidate role.
#' @return One row containing held-out response moments and mean squared error.
model_performance_fold_statistics <- function(
    captured, partition, formula, heldout_gradient, model_key,
    candidate_role) {
  test_design <- auc_design_matrix(
    formula, partition$test, names(stats::coef(captured$fit))
  )
  observed <- as.numeric(partition$test$auc_effect)
  predicted <- as.numeric(test_design %*% stats::coef(captured$fit))
  error <- predicted - observed
  if (!all(is.finite(c(observed, predicted, error)))) {
    stop(model_key, " produced non-finite held-out predictions.", call. = FALSE)
  }
  tibble::tibble(
    heldout_gradient = heldout_gradient,
    model_key = model_key,
    model = model_performance_display(model_key),
    candidate_role = candidate_role,
    n_test = length(observed),
    n_species_test = dplyr::n_distinct(partition$test$species_id),
    mean_observed = mean(observed),
    mean_observed_squared = mean(observed^2),
    MSE = mean(error^2)
  )
}

#' Fit and score all models for one omitted geographic gradient
#'
#' @param heldout_gradient Omitted gradient.
#' @param predictor_specifications Six predictor specifications.
#' @param full_data Complete resident model data.
#' @param full_sampling_variance_matrix Complete diagonal sampling variances.
#' @param climate_formula Climate-only fixed-effect formula.
#' @param random Random intercept and slope formulas.
#' @param struct Random-slope covariance structure.
#' @param method Estimation method; Appendix Figure S7 requires `"REML"`.
#' @param sparse Whether `metafor::rma.mv()` should use sparse matrices.
#' @return Fold statistics and fit audits.
run_model_performance_gradient_fold <- function(
    heldout_gradient, predictor_specifications, full_data,
    full_sampling_variance_matrix, climate_formula, random, struct,
    method = "REML", sparse = TRUE) {
  full_partition <- model_performance_partition(
    full_data, full_sampling_variance_matrix, heldout_gradient
  )
  full_climate <- fit_model_performance_fold(
    full_partition,
    climate_formula,
    paste(heldout_gradient, "complete-data climate-only REML"),
    random = random,
    struct = struct,
    method = method,
    sparse = sparse
  )
  statistics <- list(model_performance_fold_statistics(
    full_climate, full_partition, climate_formula, heldout_gradient,
    "Climate", "standalone"
  ))
  audits <- list(model_performance_fit_audit_row(
    full_climate, full_partition, heldout_gradient,
    "Climate", "standalone", "Climate_standalone", climate_formula
  ))

  for (spec in predictor_specifications) {
    partition <- model_performance_partition(
      spec$data, spec$sampling_variance_matrix, heldout_gradient
    )
    if (isTRUE(spec$full_data)) {
      baseline <- full_climate
      baseline_partition <- full_partition
    } else {
      baseline <- fit_model_performance_fold(
        partition,
        climate_formula,
        paste(heldout_gradient, spec$label, "matched climate-only REML"),
        random = random,
        struct = struct,
        method = method,
        sparse = sparse
      )
      audits[[length(audits) + 1L]] <- model_performance_fit_audit_row(
        baseline, partition, heldout_gradient, spec$key,
        "matched_baseline", paste0(spec$key, "_matched_climate"),
        climate_formula
      )
      baseline_partition <- partition
    }
    statistics[[length(statistics) + 1L]] <- model_performance_fold_statistics(
      baseline, baseline_partition, climate_formula, heldout_gradient,
      spec$key, "matched_baseline"
    )

    full_fit <- fit_model_performance_fold(
      partition,
      spec$full_formula,
      paste(heldout_gradient, spec$label, "full predictor-climate REML"),
      predictor_columns = spec$predictor_columns,
      random = random,
      struct = struct,
      method = method,
      sparse = sparse
    )
    audits[[length(audits) + 1L]] <- model_performance_fit_audit_row(
      full_fit, partition, heldout_gradient, spec$key, "full",
      paste0(spec$key, "_full"), spec$full_formula
    )
    statistics[[length(statistics) + 1L]] <- model_performance_fold_statistics(
      full_fit, partition, spec$full_formula, heldout_gradient,
      spec$key, "full"
    )
  }

  list(
    heldout_gradient = heldout_gradient,
    fold_statistics = dplyr::bind_rows(statistics),
    fit_audit = dplyr::bind_rows(audits)
  )
}

#' Combine all leave-gradient-out folds
#'
#' @param folds List of fold results.
#' @return Combined fold statistics and fit audits.
combine_model_performance_folds <- function(folds) {
  list(
    fold_statistics = dplyr::bind_rows(lapply(folds, `[[`, "fold_statistics")),
    fit_audit = dplyr::bind_rows(lapply(folds, `[[`, "fit_audit"))
  )
}

#' Aggregate RMSE and predictive R2 across equally weighted gradients
#'
#' @param fold_statistics One model-role series with one row per gradient.
#' @param sample_index Optional bootstrap row indices.
#' @return Named RMSE and predictive R2 values.
aggregate_equal_gradient_performance <- function(
    fold_statistics, sample_index = NULL) {
  if (is.null(sample_index)) sample_index <- seq_len(nrow(fold_statistics))
  x <- fold_statistics[sample_index, , drop = FALSE]
  mean_observed <- mean(x$mean_observed)
  observed_variance <- mean(x$mean_observed_squared) - mean_observed^2
  if (!is.finite(observed_variance) || observed_variance <= 0) {
    stop("Held-out response variance is not positive.", call. = FALSE)
  }
  mean_mse <- mean(x$MSE)
  c(
    RMSE = sqrt(mean_mse),
    predictive_R2 = 1 - mean_mse / observed_variance
  )
}

#' Summarise leave-gradient-out RMSE and predictive R2
#'
#' @param fold_statistics Combined fold statistics.
#' @param bootstrap_replicates Number of complete-gradient bootstrap samples.
#' @param seed Fixed bootstrap seed.
#' @return Estimates and gradient-bootstrap 95 percent intervals.
summarise_model_cross_validation <- function(
    fold_statistics, bootstrap_replicates = 5000L, seed = 20260914L) {
  group_columns <- c("model_key", "model", "candidate_role")
  group_key <- interaction(
    lapply(fold_statistics[group_columns], as.character),
    drop = TRUE, lex.order = TRUE
  )
  groups <- split(seq_len(nrow(fold_statistics)), group_key)
  set.seed(seed)
  rows <- lapply(groups, function(indices) {
    data <- fold_statistics[indices, , drop = FALSE]
    estimate <- aggregate_equal_gradient_performance(data)
    draws <- replicate(
      bootstrap_replicates,
      aggregate_equal_gradient_performance(
        data,
        sample.int(nrow(data), nrow(data), replace = TRUE)
      )
    )
    intervals <- t(apply(
      draws, 1L, stats::quantile,
      probs = c(0.025, 0.975), names = FALSE, type = 7
    ))
    tibble::tibble(
      model_key = data$model_key[[1]],
      model = data$model[[1]],
      candidate_role = data$candidate_role[[1]],
      weighting = "equal_gradient",
      metric = names(estimate),
      estimate = unname(estimate),
      CI_low = intervals[, 1L],
      CI_high = intervals[, 2L],
      n_effects = sum(data$n_test),
      n_gradients = nrow(data),
      bootstrap_replicates = as.integer(bootstrap_replicates),
      bootstrap_seed = as.integer(seed)
    )
  })
  dplyr::bind_rows(rows)
}

#' Validate the four Appendix Figure S7 measures
#'
#' @param folds Individual fold results.
#' @param combined Combined fold results.
#' @param fit_r2 Marginal and conditional R2 table.
#' @param metrics RMSE and predictive R2 table.
#' @param expected_gradients Expected number of geographic gradients.
#' @param expected_fits Expected number of validation-only REML fits.
#' @return Validation rows followed by detailed fit-audit rows.
validate_model_performance <- function(
    folds, combined, fit_r2, metrics,
    expected_gradients = 20L, expected_fits = 200L) {
  fold_statistics <- combined$fold_statistics
  fit_audit <- combined$fit_audit
  expected_series <- fit_r2 |>
    dplyr::select("model_key", "candidate_role", expected_effects = "n_effects")
  series_counts <- fold_statistics |>
    dplyr::summarise(
      heldout_effects = sum(.data$n_test),
      heldout_folds = dplyr::n_distinct(.data$heldout_gradient),
      .by = c("model_key", "candidate_role")
    ) |>
    dplyr::left_join(
      expected_series, by = c("model_key", "candidate_role")
    )
  expected_metric_names <- c("RMSE", "predictive_R2")
  metric_counts <- metrics |>
    dplyr::summarise(
      metrics = paste(sort(.data$metric), collapse = ";"),
      n_metrics = dplyr::n(),
      .by = c("model_key", "candidate_role")
    )
  validation <- tibble::tribble(
    ~check, ~pass, ~detail,
    "Twenty leave-gradient-out folds",
    length(folds) == expected_gradients &&
      dplyr::n_distinct(fold_statistics$heldout_gradient) == expected_gradients,
    paste(dplyr::n_distinct(fold_statistics$heldout_gradient),
          "distinct held-out gradients"),
    "Expected validation-only REML fits",
    nrow(fit_audit) == expected_fits,
    paste(nrow(fit_audit), "REML fold fits; no final fit replaced"),
    "Every eligible effect held out once",
    nrow(series_counts) == nrow(expected_series) &&
      all(!is.na(series_counts$expected_effects)) &&
      all(series_counts$heldout_effects == series_counts$expected_effects) &&
      all(series_counts$heldout_folds == expected_gradients),
    paste(nrow(series_counts), "model-role series checked"),
    "Fold fit convergence and design rank",
    nrow(fit_audit) == expected_fits && all(fit_audit$pass),
    paste(sum(fit_audit$pass), "of", nrow(fit_audit), "fit audits pass"),
    "Disjoint training and held-out effects",
    all(fit_audit$effect_partitions_disjoint),
    paste(sum(fit_audit$effect_partitions_disjoint), "fit partitions checked"),
    "Finite hierarchical R2",
    all(is.finite(fit_r2$marginal_R2)) &&
      all(is.finite(fit_r2$conditional_R2)) &&
      all(fit_r2$marginal_R2 >= 0 & fit_r2$marginal_R2 <= 1) &&
      all(fit_r2$conditional_R2 >= 0 & fit_r2$conditional_R2 <= 1),
    paste(nrow(fit_r2), "equal-gradient final-fit rows"),
    "Only RMSE and predictive R2 exported",
    all(metric_counts$n_metrics == 2L) &&
      all(metric_counts$metrics == paste(sort(expected_metric_names), collapse = ";")),
    paste(nrow(metrics), "cross-validation metric rows"),
    "Finite cross-validation estimates and intervals",
    all(is.finite(metrics$estimate)) && all(is.finite(metrics$CI_low)) &&
      all(is.finite(metrics$CI_high)) &&
      all(metrics$bootstrap_replicates == 5000L) &&
      all(metrics$bootstrap_seed == 20260914L),
    "All estimates use the registered complete-gradient bootstrap"
  )
  detailed_fits <- fit_audit |>
    dplyr::transmute(
      check = paste0("Fold fit: ", .data$fit_id),
      pass = .data$pass,
      detail = paste0(
        "held out ", .data$heldout_gradient,
        "; train/test=", .data$n_train, "/", .data$n_test,
        "; convergence=", .data$convergence_code,
        ifelse(nzchar(.data$warnings), paste0("; warnings=", .data$warnings), "")
      ),
      heldout_gradient = .data$heldout_gradient,
      model_key = .data$model_key,
      candidate_role = .data$candidate_role,
      fit_id = .data$fit_id,
      warnings = .data$warnings
    )
  dplyr::bind_rows(
    validation |>
      dplyr::mutate(
        heldout_gradient = NA_character_, model_key = NA_character_,
        candidate_role = NA_character_, fit_id = NA_character_,
        warnings = NA_character_
      ),
    detailed_fits
  )
}

#' Stop when Appendix Figure S7 validation fails
#'
#' @param validation Output from `validate_model_performance()`.
#' @return The validation table when all checks pass.
assert_model_performance_validation <- function(validation) {
  failed <- is.na(validation$pass) | !validation$pass
  if (any(failed)) {
    stop(
      "Model-performance validation failed: ",
      paste(unique(validation$check[failed]), collapse = ", "),
      call. = FALSE
    )
  }
  validation
}

#' Build Appendix Figure S7
#'
#' @param fit_r2 Marginal and conditional R2 table.
#' @param metrics Leave-gradient-out RMSE and predictive R2 table.
#' @return A four-panel ggplot object.
build_model_performance_figure <- function(fit_r2, metrics) {
  key <- model_performance_key()
  selected_roles <- c("standalone", "full")
  sample_key <- fit_r2 |>
    dplyr::filter(.data$candidate_role %in% .env$selected_roles) |>
    dplyr::distinct(.data$model, .data$n_effects, .data$n_species)
  if (nrow(sample_key) != nrow(key)) {
    stop("The performance figure needs one sample size per model.", call. = FALSE)
  }
  model_labels <- stats::setNames(
    paste0(
      sample_key$model, "\n(n = ", sample_key$n_effects,
      " effects; ", sample_key$n_species, " species)"
    ),
    sample_key$model
  )
  model_levels <- rev(unname(model_labels[key$model_display]))

  fit_data <- fit_r2 |>
    dplyr::filter(.data$candidate_role %in% .env$selected_roles) |>
    tidyr::pivot_longer(
      c("marginal_R2", "conditional_R2"),
      names_to = "metric", values_to = "estimate"
    ) |>
    dplyr::transmute(
      model = .data$model,
      metric = .data$metric,
      estimate = .data$estimate,
      CI_low = NA_real_, CI_high = NA_real_
    )
  predictive_data <- metrics |>
    dplyr::filter(
      .data$candidate_role %in% .env$selected_roles,
      .data$metric %in% c("predictive_R2", "RMSE")
    ) |>
    dplyr::transmute(
      model = .data$model,
      metric = .data$metric,
      estimate = .data$estimate,
      CI_low = .data$CI_low,
      CI_high = .data$CI_high
    )
  metric_labels <- c(
    marginal_R2 = "'Marginal'~R^2~(R[m]^2)",
    conditional_R2 = "'Conditional'~R^2~(R[c]^2)",
    predictive_R2 = "'Leave-gradient-out predictive'~R^2~(R[pred]^2)",
    RMSE = "'Leave-gradient-out RMSE'"
  )
  plot_data <- dplyr::bind_rows(fit_data, predictive_data) |>
    dplyr::mutate(
      model = factor(
        unname(.env$model_labels[as.character(.data$model)]),
        levels = .env$model_levels
      ),
      metric_label = factor(
        unname(.env$metric_labels[.data$metric]),
        levels = unname(.env$metric_labels)
      ),
      value_label = sprintf("%.3f", .data$estimate)
    )
  expected_rows <- nrow(key) * length(metric_labels)
  if (nrow(plot_data) != expected_rows || anyNA(plot_data$metric_label)) {
    stop("The performance figure needs one estimate per model and metric.",
         call. = FALSE)
  }
  metric_colours <- stats::setNames(
    c("#3B4A6B", "#1B9E77", "#D95F02", "#6A3D9A"),
    unname(metric_labels)
  )
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$estimate, y = .data$model, colour = .data$metric_label)
  ) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.7) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$CI_low, xmax = .data$CI_high),
      width = 0.18, linewidth = 0.8, orientation = "y", na.rm = TRUE
    ) +
    ggplot2::geom_point(size = 4.5) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$value_label),
      vjust = -1.05, size = 5.4, fontface = "bold", show.legend = FALSE
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(.data$metric_label), scales = "free_x", ncol = 2,
      labeller = ggplot2::label_parsed
    ) +
    ggplot2::scale_colour_manual(values = metric_colours, guide = "none") +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.08, 0.15))
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      x = "Estimate", y = NULL,
      subtitle = paste(
        "Every gradient has equal influence; whiskers are",
        "complete-gradient bootstrap 95% intervals"
      )
    ) +
    ggplot2::theme_bw(base_size = 21) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.spacing = grid::unit(1.6, "lines"),
      strip.text = ggplot2::element_text(size = 22, face = "bold"),
      axis.title.x = ggplot2::element_text(size = 23),
      axis.text.x = ggplot2::element_text(size = 18, colour = "black"),
      axis.text.y = ggplot2::element_text(size = 16.5, colour = "black"),
      plot.subtitle = ggplot2::element_text(size = 20, hjust = 0.5),
      plot.margin = ggplot2::margin(18, 28, 18, 18)
    )
}
