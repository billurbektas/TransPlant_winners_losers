# Module: R/evaluation/model_performance
# Responsibility: Implement the named scientific workflow stage or shared infrastructure.
# Inputs: Function arguments and explicitly declared target dependencies.
# Outputs: Validated in-memory objects and declared workflow artifacts.
# Side effects: Limited to functions whose contracts explicitly write, copy, or render files.

# Global model-fit and held-out predictive-performance helpers for the
# canonical 526-species meta-regression lineage.

#' Model Performance Key.
#'
#' @return An R object described by the function name; its structure is validated before later stages use it.
model_performance_key <- function() {
  tibble::tribble(
    ~model_key,  ~model_display,          ~model_order,
    "Climate",  "Climate-only",                   1L,
    "LDMC",     "LDMC",                           2L,
    "LMA",      "LMA",                            3L,
    "Height",   "Vegetative height",              4L,
    "Graminoid","Graminoid status",               5L,
    "Woodiness","Woodiness",                      6L,
    "Phylogeny", "Phylogeny (PEV1-PEV3)",         7L
  )
}

#' Model Performance Display.
#'
#' @param model_key A fitted model object with the structure checked by this function.
#' @return An R object described by the function name; its structure is validated before later stages use it.
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

#' Model Performance Formula Text.
#'
#' @param formula A model formula or random-effects specification.
#' @return An R object described by the function name; its structure is validated before later stages use it.
model_performance_formula_text <- function(formula) {
  paste(deparse(formula), collapse = " ")
}

#' Model Performance Random Block.
#'
#' @param fit A fitted model object with the structure checked by this function.
#' @return An R object described by the function name; its structure is validated before later stages use it.
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

#' Model Performance Random Design.
#'
#' @param data A data frame whose required columns are checked before use.
#' @return An R object described by the function name; its structure is validated before later stages use it.
model_performance_random_design <- function(data) {
  required <- c("W_z", "T_z", "P_z", "TP_z")
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(
      "Random-effect prediction data lack: ",
      paste(missing, collapse = ", "), call. = FALSE
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

#' Model Performance Random Covariance.
#'
#' @param fit A fitted model object with the structure checked by this function.
#' @param left A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param right A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
model_performance_random_covariance <- function(fit, left, right = left) {
  if (!"species_id" %in% names(left) || !"species_id" %in% names(right)) {
    stop("Random-effect prediction data require species_id.", call. = FALSE)
  }
  left_design <- model_performance_random_design(left)
  right_design <- model_performance_random_design(right)
  same_species <- outer(
    as.character(left$species_id), as.character(right$species_id), `==`
  )
  covariance <- left_design %*% model_performance_random_block(fit) %*%
    t(right_design)
  covariance * same_species
}

#' Model Performance Random Variance.
#'
#' @param fit A fitted model object with the structure checked by this function.
#' @param data A data frame whose required columns are checked before use.
#' @return An R object described by the function name; its structure is validated before later stages use it.
model_performance_random_variance <- function(fit, data) {
  design <- model_performance_random_design(data)
  block <- model_performance_random_block(fit)
  rowSums((design %*% block) * design)
}

#' Model Performance Weights.
#'
#' @param data A data frame whose required columns are checked before use.
#' @param aggregation A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
model_performance_weights <- function(data, aggregation) {
  if (!"gradient_id" %in% names(data)) {
    stop("Performance data require gradient_id.", call. = FALSE)
  }
  if (identical(aggregation, "effect_pooled")) {
    return(rep(1 / nrow(data), nrow(data)))
  }
  if (!identical(aggregation, "equal_gradient")) {
    stop("Unknown performance aggregation: ", aggregation, call. = FALSE)
  }
  gradient <- as.character(data$gradient_id)
  counts <- table(gradient)
  1 / (length(counts) * as.numeric(counts[gradient]))
}

#' Model Performance Fit R2.
#'
#' @param fit A fitted model object with the structure checked by this function.
#' @param model_key A fitted model object with the structure checked by this function.
#' @param candidate_role A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
model_performance_fit_r2 <- function(fit, model_key, candidate_role) {
  data <- fit$data
  if (!is.data.frame(data) || nrow(data) != fit$k) {
    stop("A fitted rma.mv model with retained data is required.", call. = FALSE)
  }
  fixed <- as.numeric(fit$X %*% stats::coef(fit))
  random <- model_performance_random_variance(fit, data)
  sampling <- as.numeric(Matrix::diag(fit$V))
  rows <- lapply(c("equal_gradient", "effect_pooled"), function(aggregation) {
    weight <- model_performance_weights(data, aggregation)
    fixed_mean <- sum(weight * fixed)
    fixed_variance <- sum(weight * (fixed - fixed_mean)^2)
    random_variance <- sum(weight * random)
    sampling_variance <- sum(weight * sampling)
    total <- fixed_variance + random_variance + sampling_variance
    if (!is.finite(total) || total <= 0) {
      stop("The response-scale variance partition is invalid.", call. = FALSE)
    }
    tibble::tibble(
      model_key = model_key,
      model = model_performance_display(model_key),
      candidate_role = candidate_role,
      aggregation = aggregation,
      n_effects = nrow(data),
      n_species = dplyr::n_distinct(data$species_id),
      n_gradients = dplyr::n_distinct(data$gradient_id),
      fixed_effect_variance = fixed_variance,
      species_random_variance = random_variance,
      sampling_variance = sampling_variance,
      total_response_variance = total,
      marginal_R2 = fixed_variance / total,
      conditional_R2 = (fixed_variance + random_variance) / total,
      R2_definition = paste(
        "response-scale hierarchical variance partition; denominator is",
        "fixed-effect variance + mean species-random variance + mean sampling variance"
      )
    )
  })
  dplyr::bind_rows(rows)
}

#' Assemble Model Fit R2.
#'
#' @param climate_model A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param ... Additional arguments forwarded to the delegated operation.
#' @return An R object described by the function name; its structure is validated before later stages use it.
assemble_model_fit_r2 <- function(climate_model, ...) {
  model_pairs <- list(...)
  rows <- list(
    model_performance_fit_r2(
      climate_model$fit, "Climate", "standalone"
    )
  )
  for (model_pair in model_pairs) {
    rows[[length(rows) + 1L]] <- model_performance_fit_r2(
      model_pair$fits$climate_REML$fit,
      model_pair$key, "matched_baseline"
    )
    rows[[length(rows) + 1L]] <- model_performance_fit_r2(
      model_pair$fits$full_REML$fit, model_pair$key, "full"
    )
  }
  dplyr::bind_rows(rows)
}

#' Model Performance Partition.
#'
#' @param data A data frame whose required columns are checked before use.
#' @param sampling_variance_matrix Diagonal, effect-ordered sampling variances.
#' @param heldout_gradient A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
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
  train_V <- subset_sampling_variance_matrix(
    sampling_variance_matrix, train_ids
  )
  test_V <- subset_sampling_variance_matrix(
    sampling_variance_matrix, test_ids
  )
  effect_partitions_disjoint <- !any(train_ids %in% test_ids)
  if (!effect_partitions_disjoint) {
    stop("Training and held-out effects overlap.", call. = FALSE)
  }
  training_experiments <- train |>
    dplyr::distinct(.data$experiment_id, .data$T_z, .data$P_z, .data$W_z)
  support <- auc_support_mask(
    training_experiments,
    test |> dplyr::select("T_z", "P_z", "W_z")
  )
  if (nrow(support$grid) != nrow(test)) {
    stop("The held-out climate-support calculation failed.", call. = FALSE)
  }
  list(
    train = train,
    test = test,
    train_V = train_V,
    test_V = test_V,
    climate_supported = as.logical(support$grid$supported),
    support_diagnostic = support$diagnostic,
    effect_partitions_disjoint = effect_partitions_disjoint
  )
}

#' Fit Model Performance Fold.
#'
#' @param partition A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param formula A model formula or random-effects specification.
#' @param label A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param predictor_columns A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A fitted model bundle with diagnostics retained for validation.
fit_model_performance_fold <- function(partition, formula, label,
                                       predictor_columns = character()) {
  design <- stats::model.matrix(formula, data = partition$train)
  full_rank <- qr(design)$rank == ncol(design)
  if (!full_rank) {
    stop(label, " has a rank-deficient training design.", call. = FALSE)
  }
  captured <- fit_meta_regression(
    partition$train, partition$train_V, formula, "REML", label,
    predictor_columns = predictor_columns
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

#' Model Performance Covariance Audit.
#'
#' @param covariance A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param label A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
model_performance_covariance_audit <- function(covariance, label) {
  covariance <- (as.matrix(covariance) + t(as.matrix(covariance))) / 2
  eigenvalues <- eigen(covariance, symmetric = TRUE, only.values = TRUE)$values
  minimum <- min(eigenvalues)
  cholesky_ok <- !inherits(try(chol(covariance), silent = TRUE), "try-error")
  if (!all(is.finite(covariance)) || !is.finite(minimum) || !cholesky_ok) {
    stop(label, " predictive covariance is not positive definite.", call. = FALSE)
  }
  list(covariance = covariance, minimum_eigenvalue = minimum)
}

#' Condition Model Performance Gaussian.
#'
#' @param train_y A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param train_mean A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param train_covariance A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param train_design A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param test_mean A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param test_design A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param test_marginal_covariance A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param cross_random_covariance A model formula or random-effects specification.
#' @param beta_covariance A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
condition_model_performance_gaussian <- function(
    train_y, train_mean, train_covariance, train_design,
    test_mean, test_design, test_marginal_covariance,
    cross_random_covariance, beta_covariance) {
  marginal_solution <- Matrix::solve(
    train_covariance,
    cbind(train_y - train_mean, t(cross_random_covariance))
  )
  alpha <- as.numeric(marginal_solution[, 1L])
  cross_times_inverse <- t(as.matrix(marginal_solution[, -1L, drop = FALSE]))
  conditional_mean <- as.numeric(
    test_mean + cross_random_covariance %*% alpha
  )
  adjusted_design <- test_design - cross_times_inverse %*% train_design
  conditional_covariance <- test_marginal_covariance -
    cross_times_inverse %*% t(cross_random_covariance) +
    adjusted_design %*% beta_covariance %*% t(adjusted_design)
  list(
    mean = conditional_mean,
    covariance = (conditional_covariance + t(conditional_covariance)) / 2,
    adjusted_design = adjusted_design
  )
}

#' Predict Model Performance Fold.
#'
#' @param captured A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param partition A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param formula A model formula or random-effects specification.
#' @param label A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
predict_model_performance_fold <- function(captured, partition, formula, label) {
  fit <- captured$fit
  train <- partition$train
  test <- partition$test
  if (!identical(as.character(fit$data$effect_id), as.character(train$effect_id))) {
    stop(label, ": fitted data and training partition differ.", call. = FALSE)
  }
  train_design <- as.matrix(fit$X)
  test_design <- auc_design_matrix(formula, test, names(stats::coef(fit)))
  beta <- as.numeric(stats::coef(fit))
  train_mean <- as.numeric(train_design %*% beta)
  test_mean <- as.numeric(test_design %*% beta)
  test_random <- model_performance_random_covariance(fit, test)
  test_marginal <- as.matrix(partition$test_V) + test_random
  general_covariance <- test_marginal +
    test_design %*% as.matrix(fit$vb) %*% t(test_design)
  general_audit <- model_performance_covariance_audit(
    general_covariance, paste(label, "general")
  )

  cross_random <- model_performance_random_covariance(fit, test, train)
  conditioned <- condition_model_performance_gaussian(
    train_y = as.numeric(fit$yi),
    train_mean = train_mean,
    train_covariance = fit$M,
    train_design = train_design,
    test_mean = test_mean,
    test_design = test_design,
    test_marginal_covariance = test_marginal,
    cross_random_covariance = cross_random,
    beta_covariance = as.matrix(fit$vb)
  )
  conditional_audit <- model_performance_covariance_audit(
    conditioned$covariance, paste(label, "conditional")
  )
  seen <- as.character(test$species_id) %in% as.character(train$species_id)
  unseen <- !seen
  unseen_mean_error <- if (any(unseen)) {
    max(abs(conditioned$mean[unseen] - test_mean[unseen]))
  } else {
    0
  }
  unseen_variance_error <- if (any(unseen)) {
    max(abs(
      diag(conditioned$covariance)[unseen] -
        diag(general_audit$covariance)[unseen]
    ))
  } else {
    0
  }
  make_prediction_rows <- function(index, prediction_scope, predicted,
                                   covariance) {
    predictive_se <- sqrt(diag(covariance))
    critical <- stats::qnorm(0.975)
    test[index, , drop = FALSE] |>
      dplyr::transmute(
        effect_id = as.character(.data$effect_id),
        gradient_id = as.character(.data$gradient_id),
        experiment_id = as.character(.data$experiment_id),
        species_id = as.character(.data$species_id),
        observed_delta_auc = as.numeric(.data$auc_effect),
        predicted_delta_auc = as.numeric(.env$predicted),
        residual = .data$predicted_delta_auc - .data$observed_delta_auc,
        predictive_se = as.numeric(.env$predictive_se),
        prediction_interval_low = .data$predicted_delta_auc -
          .env$critical * .data$predictive_se,
        prediction_interval_high = .data$predicted_delta_auc +
          .env$critical * .data$predictive_se,
        interval_covered = .data$observed_delta_auc >=
          .data$prediction_interval_low &
          .data$observed_delta_auc <= .data$prediction_interval_high,
        interval_width = .data$prediction_interval_high -
          .data$prediction_interval_low,
        prediction_scope = prediction_scope,
        species_seen_in_training = .env$seen[index],
        climate_supported = .env$partition$climate_supported[index]
      )
  }
  general_rows <- make_prediction_rows(
    seq_len(nrow(test)), "general", test_mean,
    general_audit$covariance
  )
  repeat_index <- which(seen)
  repeat_rows <- if (length(repeat_index)) {
    make_prediction_rows(
      repeat_index, "repeat_species_conditional",
      conditioned$mean[repeat_index],
      conditional_audit$covariance[repeat_index, repeat_index, drop = FALSE]
    )
  } else {
    general_rows[0, , drop = FALSE]
  }
  list(
    predictions = dplyr::bind_rows(general_rows, repeat_rows),
    distributions = list(
      general = list(
        rows = seq_len(nrow(test)),
        mean = test_mean,
        covariance = general_audit$covariance
      ),
      repeat_species_conditional = list(
        rows = repeat_index,
        mean = conditioned$mean[repeat_index],
        covariance = conditional_audit$covariance[
          repeat_index, repeat_index, drop = FALSE
        ]
      )
    ),
    audit = tibble::tibble(
      n_test = nrow(test),
      n_repeat_species_effects = sum(seen),
      n_unseen_species_effects = sum(unseen),
      general_minimum_eigenvalue = general_audit$minimum_eigenvalue,
      conditional_minimum_eigenvalue = conditional_audit$minimum_eigenvalue,
      maximum_unseen_mean_difference = unseen_mean_error,
      maximum_unseen_variance_difference = unseen_variance_error,
      pass = general_audit$minimum_eigenvalue > 0 &&
        conditional_audit$minimum_eigenvalue > 0 &&
        unseen_mean_error <= 1e-10 && unseen_variance_error <= 1e-10
    )
  )
}

#' Model Performance Gaussian Score.
#'
#' @param observed A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param predicted A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param covariance A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
model_performance_gaussian_score <- function(observed, predicted, covariance) {
  covariance <- (as.matrix(covariance) + t(as.matrix(covariance))) / 2
  factor <- chol(covariance)
  error <- as.numeric(observed - predicted)
  standardized <- forwardsolve(t(factor), error)
  log_determinant <- 2 * sum(log(diag(factor)))
  total <- 0.5 * (
    length(error) * log(2 * pi) + log_determinant + sum(standardized^2)
  )
  c(
    negative_log_predictive_density_total = total,
    negative_log_predictive_density = total / length(error),
    log_determinant = log_determinant,
    quadratic_form = sum(standardized^2)
  )
}

#' Decorate Model Performance Prediction.
#'
#' @param bundle A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param partition A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param heldout_gradient A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param model_key A fitted model object with the structure checked by this function.
#' @param candidate_role A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param fit_id A fitted model object with the structure checked by this function.
#' @param formula A model formula or random-effects specification.
#' @return An R object described by the function name; its structure is validated before later stages use it.
decorate_model_performance_prediction <- function(
    bundle, partition, heldout_gradient, model_key, candidate_role,
    fit_id, formula) {
  display <- model_performance_display(model_key)
  predictions <- bundle$predictions |>
    dplyr::mutate(
      heldout_gradient = heldout_gradient,
      model_key = model_key,
      model = display,
      candidate_role = candidate_role,
      matched_baseline = ifelse(
        model_key == "Climate", NA_character_, "Climate-only"
      ),
      fit_id = fit_id,
      formula = model_performance_formula_text(formula),
      .before = 1L
    )
  gaussian_rows <- list()
  add_scores <- function(distribution, prediction_scope, evaluation_subset) {
    if (!length(distribution$rows)) return(NULL)
    subset_rows <- distribution$rows
    supported <- partition$climate_supported[subset_rows]
    combinations <- list(
      list(
        keep = rep(TRUE, length(subset_rows)),
        evaluation_subset = evaluation_subset,
        support_stratum = "all"
      ),
      list(
        keep = supported,
        evaluation_subset = evaluation_subset,
        support_stratum = "in_support"
      )
    )
    if (identical(prediction_scope, "general")) {
      repeated <- as.character(partition$test$species_id[subset_rows]) %in%
        as.character(partition$train$species_id)
      combinations <- c(
        combinations,
        list(
          list(
            keep = repeated,
            evaluation_subset = "repeat_species",
            support_stratum = "all"
          ),
          list(
            keep = repeated & supported,
            evaluation_subset = "repeat_species",
            support_stratum = "in_support"
          )
        )
      )
    }
    rows <- lapply(combinations, function(combination) {
      keep <- which(combination$keep)
      if (!length(keep)) return(NULL)
      covariance <- distribution$covariance[keep, keep, drop = FALSE]
      scores <- model_performance_gaussian_score(
        partition$test$auc_effect[subset_rows[keep]],
        distribution$mean[keep], covariance
      )
      tibble::tibble(
        heldout_gradient = heldout_gradient,
        model_key = model_key,
        model = display,
        candidate_role = candidate_role,
        fit_id = fit_id,
        prediction_scope = prediction_scope,
        evaluation_subset = combination$evaluation_subset,
        support_stratum = combination$support_stratum,
        n_effects = length(keep),
        negative_log_predictive_density_total = unname(scores[[1]]),
        negative_log_predictive_density = unname(scores[[2]]),
        log_determinant = unname(scores[[3]]),
        quadratic_form = unname(scores[[4]])
      )
    })
    dplyr::bind_rows(rows)
  }
  gaussian_rows[[1L]] <- add_scores(
    bundle$distributions$general, "general", "all_eligible"
  )
  gaussian_rows[[2L]] <- add_scores(
    bundle$distributions$repeat_species_conditional,
    "repeat_species_conditional", "repeat_species"
  )
  list(
    predictions = predictions,
    gaussian_scores = dplyr::bind_rows(gaussian_rows),
    audit = bundle$audit |>
      dplyr::mutate(
        heldout_gradient = heldout_gradient,
        model_key = model_key,
        model = display,
        candidate_role = candidate_role,
        fit_id = fit_id,
        .before = 1L
      )
  )
}

#' Model Performance Fit Audit Row.
#'
#' @param captured A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param partition A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param heldout_gradient A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param model_key A fitted model object with the structure checked by this function.
#' @param candidate_role A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param fit_id A fitted model object with the structure checked by this function.
#' @param formula A model formula or random-effects specification.
#' @return An R object described by the function name; its structure is validated before later stages use it.
model_performance_fit_audit_row <- function(
    captured, partition, heldout_gradient, model_key, candidate_role,
    fit_id, formula) {
  captured$performance_audit |>
    dplyr::mutate(
      heldout_gradient = heldout_gradient,
      model_key = model_key,
      model = model_performance_display(model_key),
      candidate_role = candidate_role,
      fit_id = fit_id,
      formula = model_performance_formula_text(formula),
      n_train = nrow(partition$train),
      n_test = nrow(partition$test),
      n_train_species = dplyr::n_distinct(partition$train$species_id),
      n_test_species = dplyr::n_distinct(partition$test$species_id),
      effect_partitions_disjoint = partition$effect_partitions_disjoint,
      support_diagnostic = partition$support_diagnostic,
      .before = 1L
    )
}

#' Run Model Performance Gradient Fold.
#'
#' @param heldout_gradient A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param predictor_specifications A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param full_data A data frame whose required columns are checked before use.
#' @param full_sampling_variance_matrix Complete diagonal sampling-variance matrix.
#' @param climate_formula A model formula or random-effects specification.
#' @details Holds out one complete gradient, refits the registered climate and
#' predictor formulas on all remaining gradients, and returns general and
#' repeat-species predictions with covariance-aware Gaussian scores. Sampling
#' variance order, support diagnostics, and fit warnings are retained; no
#' fitted release model is replaced by these evaluation-only fits.
#' @return An R object described by the function name; its structure is validated before later stages use it.
run_model_performance_gradient_fold <- function(
    heldout_gradient, predictor_specifications, full_data,
    full_sampling_variance_matrix,
    climate_formula) {
  heldout_gradient <- as.character(heldout_gradient[[1]])
  predictions <- list()
  gaussian_scores <- list()
  fit_audits <- list()
  prediction_audits <- list()

  full_partition <- model_performance_partition(
    full_data, full_sampling_variance_matrix, heldout_gradient
  )
  baseline_id <- "Climate_full_data"
  full_baseline <- fit_model_performance_fold(
    full_partition, climate_formula,
    paste(heldout_gradient, "full-data climate-only REML")
  )
  fit_audits[[length(fit_audits) + 1L]] <-
    model_performance_fit_audit_row(
      full_baseline, full_partition, heldout_gradient,
      "Climate", "standalone", baseline_id, climate_formula
    )
  full_baseline_bundle <- predict_model_performance_fold(
    full_baseline, full_partition, climate_formula,
    paste(heldout_gradient, "full-data climate-only")
  )
  decorated <- decorate_model_performance_prediction(
    full_baseline_bundle, full_partition, heldout_gradient,
    "Climate", "standalone", baseline_id, climate_formula
  )
  predictions[[length(predictions) + 1L]] <- decorated$predictions
  gaussian_scores[[length(gaussian_scores) + 1L]] <- decorated$gaussian_scores
  prediction_audits[[length(prediction_audits) + 1L]] <- decorated$audit

  for (spec in predictor_specifications) {
    partition <- model_performance_partition(
      spec$data, spec$sampling_variance_matrix, heldout_gradient
    )
    if (isTRUE(spec$full_data)) {
      if (!identical(
        as.character(partition$train$effect_id),
        as.character(full_partition$train$effect_id)
      ) || !identical(
        as.character(partition$test$effect_id),
        as.character(full_partition$test$effect_id)
      )) {
        stop(spec$label, " no longer matches the full-data fold.", call. = FALSE)
      }
      baseline <- full_baseline
      baseline_bundle <- full_baseline_bundle
      baseline_fit_id <- baseline_id
    } else {
      baseline_fit_id <- paste0(spec$key, "_matched_climate")
      baseline <- fit_model_performance_fold(
        partition, climate_formula,
        paste(heldout_gradient, spec$label, "matched climate-only REML")
      )
      fit_audits[[length(fit_audits) + 1L]] <-
        model_performance_fit_audit_row(
          baseline, partition, heldout_gradient, spec$key,
          "matched_baseline", baseline_fit_id, climate_formula
        )
      baseline_bundle <- predict_model_performance_fold(
        baseline, partition, climate_formula,
        paste(heldout_gradient, spec$label, "matched climate-only")
      )
    }
    decorated_baseline <- decorate_model_performance_prediction(
      baseline_bundle, partition, heldout_gradient, spec$key,
      "matched_baseline", baseline_fit_id, climate_formula
    )
    predictions[[length(predictions) + 1L]] <-
      decorated_baseline$predictions
    gaussian_scores[[length(gaussian_scores) + 1L]] <-
      decorated_baseline$gaussian_scores
    prediction_audits[[length(prediction_audits) + 1L]] <-
      decorated_baseline$audit

    full_fit_id <- paste0(spec$key, "_full")
    full_fit <- fit_model_performance_fold(
      partition, spec$full_formula,
      paste(heldout_gradient, spec$label, "full predictor-climate REML"),
      spec$predictor_columns
    )
    fit_audits[[length(fit_audits) + 1L]] <-
      model_performance_fit_audit_row(
        full_fit, partition, heldout_gradient, spec$key,
        "full", full_fit_id, spec$full_formula
      )
    full_bundle <- predict_model_performance_fold(
      full_fit, partition, spec$full_formula,
      paste(heldout_gradient, spec$label, "full predictor-climate")
    )
    decorated_full <- decorate_model_performance_prediction(
      full_bundle, partition, heldout_gradient, spec$key,
      "full", full_fit_id, spec$full_formula
    )
    predictions[[length(predictions) + 1L]] <- decorated_full$predictions
    gaussian_scores[[length(gaussian_scores) + 1L]] <-
      decorated_full$gaussian_scores
    prediction_audits[[length(prediction_audits) + 1L]] <- decorated_full$audit
  }

  list(
    heldout_gradient = heldout_gradient,
    predictions = dplyr::bind_rows(predictions),
    gaussian_scores = dplyr::bind_rows(gaussian_scores),
    fit_audit = dplyr::bind_rows(fit_audits),
    prediction_audit = dplyr::bind_rows(prediction_audits)
  )
}

#' Combine Model Performance Folds.
#'
#' @param folds A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
combine_model_performance_folds <- function(folds) {
  list(
    predictions = dplyr::bind_rows(lapply(folds, `[[`, "predictions")),
    gaussian_scores = dplyr::bind_rows(
      lapply(folds, `[[`, "gaussian_scores")
    ),
    fit_audit = dplyr::bind_rows(lapply(folds, `[[`, "fit_audit")),
    prediction_audit = dplyr::bind_rows(
      lapply(folds, `[[`, "prediction_audit")
    )
  )
}

#' Model Performance Evaluation Rows.
#'
#' @param predictions A data frame whose required columns are checked before use.
#' @return An R object described by the function name; its structure is validated before later stages use it.
model_performance_evaluation_rows <- function(predictions) {
  general <- predictions |>
    dplyr::filter(.data$prediction_scope == "general")
  evaluation <- dplyr::bind_rows(
    general |> dplyr::mutate(evaluation_subset = "all_eligible"),
    general |>
      dplyr::filter(.data$species_seen_in_training) |>
      dplyr::mutate(evaluation_subset = "repeat_species"),
    predictions |>
      dplyr::filter(
        .data$prediction_scope == "repeat_species_conditional"
      ) |>
      dplyr::mutate(evaluation_subset = "repeat_species")
  )
  dplyr::bind_rows(
    evaluation |>
      dplyr::filter(.data$climate_supported) |>
      dplyr::mutate(support_stratum = "in_support"),
    evaluation |> dplyr::mutate(support_stratum = "all")
  )
}

#' Model Performance Sufficient Statistics.
#'
#' @param data A data frame whose required columns are checked before use.
#' @param gaussian_scores A data frame whose required columns are checked before use.
#' @return An R object described by the function name; its structure is validated before later stages use it.
model_performance_sufficient_statistics <- function(data, gaussian_scores) {
  gradient_rows <- lapply(
    split(seq_len(nrow(data)), as.character(data$heldout_gradient)),
    function(index) {
      x <- data[index, , drop = FALSE]
      tibble::tibble(
        heldout_gradient = as.character(x$heldout_gradient[[1]]),
        n = nrow(x),
        sum_y = sum(x$observed_delta_auc),
        sum_y2 = sum(x$observed_delta_auc^2),
        sum_error = sum(x$residual),
        sum_absolute_error = sum(abs(x$residual)),
        sum_squared_error = sum(x$residual^2),
        sum_covered = sum(x$interval_covered),
        sum_interval_width = sum(x$interval_width)
      )
    }
  )
  statistics <- dplyr::bind_rows(gradient_rows) |>
    dplyr::left_join(
      gaussian_scores |>
        dplyr::select(
          "heldout_gradient", "n_effects",
          "negative_log_predictive_density_total"
        ),
      by = "heldout_gradient"
    )
  if (anyNA(statistics) || any(statistics$n != statistics$n_effects)) {
    stop("Point and Gaussian fold summaries do not align.", call. = FALSE)
  }
  statistics
}

#' Aggregate Model Performance Statistics.
#'
#' @param statistics A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param aggregation A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param sample_index A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
aggregate_model_performance_statistics <- function(statistics, aggregation,
                                                   sample_index = NULL) {
  if (is.null(sample_index)) sample_index <- seq_len(nrow(statistics))
  x <- statistics[sample_index, , drop = FALSE]
  if (identical(aggregation, "equal_gradient")) {
    mean_y <- mean(x$sum_y / x$n)
    mean_y2 <- mean(x$sum_y2 / x$n)
    mse <- mean(x$sum_squared_error / x$n)
    mae <- mean(x$sum_absolute_error / x$n)
    bias <- mean(x$sum_error / x$n)
    coverage <- mean(x$sum_covered / x$n)
    interval_width <- mean(x$sum_interval_width / x$n)
    nlpd <- mean(x$negative_log_predictive_density_total / x$n_effects)
  } else if (identical(aggregation, "effect_pooled")) {
    n <- sum(x$n)
    mean_y <- sum(x$sum_y) / n
    mean_y2 <- sum(x$sum_y2) / n
    mse <- sum(x$sum_squared_error) / n
    mae <- sum(x$sum_absolute_error) / n
    bias <- sum(x$sum_error) / n
    coverage <- sum(x$sum_covered) / n
    interval_width <- sum(x$sum_interval_width) / n
    nlpd <- sum(x$negative_log_predictive_density_total) /
      sum(x$n_effects)
  } else {
    stop("Unknown performance aggregation: ", aggregation, call. = FALSE)
  }
  response_variance <- mean_y2 - mean_y^2
  if (!is.finite(response_variance) || response_variance <= 0) {
    stop("Held-out response variance is not positive.", call. = FALSE)
  }
  c(
    RMSE = sqrt(mse),
    MAE = mae,
    bias = bias,
    predictive_R2 = 1 - mse / response_variance,
    gaussian_NLPD = nlpd,
    coverage_95 = coverage,
    mean_interval_width = interval_width,
    MSE = mse
  )
}

#' Model Performance Group Indices.
#'
#' @param data A data frame whose required columns are checked before use.
#' @param columns A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
model_performance_group_indices <- function(data, columns) {
  key <- do.call(
    interaction,
    c(lapply(data[columns], as.character), list(drop = TRUE, lex.order = TRUE))
  )
  split(seq_len(nrow(data)), key)
}

#' Summarise Model Cross Validation.
#'
#' @param predictions A data frame whose required columns are checked before use.
#' @param gaussian_scores A data frame whose required columns are checked before use.
#' @param bootstrap_replicates A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param seed A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
summarise_model_cross_validation <- function(
    predictions, gaussian_scores, bootstrap_replicates = 5000L,
    seed = 20260914L) {
  evaluation <- model_performance_evaluation_rows(predictions)
  group_columns <- c(
    "model_key", "model", "candidate_role", "prediction_scope",
    "evaluation_subset", "support_stratum"
  )
  groups <- model_performance_group_indices(evaluation, group_columns)
  set.seed(seed)
  rows <- list()
  for (indices in groups) {
    data <- evaluation[indices, , drop = FALSE]
    group <- data[1L, group_columns, drop = FALSE]
    gaussian <- gaussian_scores
    for (column in group_columns) {
      gaussian <- gaussian[
        as.character(gaussian[[column]]) == as.character(group[[column]][[1]]),
        , drop = FALSE
      ]
    }
    statistics <- model_performance_sufficient_statistics(data, gaussian)
    for (aggregation in c("equal_gradient", "effect_pooled")) {
      estimate <- aggregate_model_performance_statistics(
        statistics, aggregation
      )
      draws <- replicate(
        bootstrap_replicates,
        aggregate_model_performance_statistics(
          statistics, aggregation,
          sample.int(nrow(statistics), nrow(statistics), replace = TRUE)
        )
      )
      intervals <- t(apply(
        draws, 1L, stats::quantile,
        probs = c(0.025, 0.975), names = FALSE, type = 7
      ))
      output_metrics <- setdiff(names(estimate), "MSE")
      metric_index <- match(output_metrics, names(estimate))
      row <- group[rep(1L, length(output_metrics)), , drop = FALSE]
      row$aggregation <- aggregation
      row$metric <- output_metrics
      row$estimate <- unname(estimate[metric_index])
      row$CI_low <- intervals[metric_index, 1L]
      row$CI_high <- intervals[metric_index, 2L]
      row$n_effects <- nrow(data)
      row$n_species <- dplyr::n_distinct(data$species_id)
      row$n_gradients <- nrow(statistics)
      row$bootstrap_replicates <- as.integer(bootstrap_replicates)
      row$bootstrap_seed <- as.integer(seed)
      rows[[length(rows) + 1L]] <- row
    }
  }
  dplyr::bind_rows(rows)
}

#' Model Performance Skill Values.
#'
#' @param full A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param baseline A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
model_performance_skill_values <- function(full, baseline) {
  c(
    RMSE_percent_improvement = 100 * (1 - full[["RMSE"]] / baseline[["RMSE"]]),
    MAE_percent_improvement = 100 * (1 - full[["MAE"]] / baseline[["MAE"]]),
    incremental_predictive_R2 = 1 - full[["MSE"]] / baseline[["MSE"]],
    absolute_bias_improvement = abs(baseline[["bias"]]) - abs(full[["bias"]]),
    gaussian_log_score_improvement =
      baseline[["gaussian_NLPD"]] - full[["gaussian_NLPD"]]
  )
}

#' Summarise Model Predictive Skill.
#'
#' @param predictions A data frame whose required columns are checked before use.
#' @param gaussian_scores A data frame whose required columns are checked before use.
#' @param bootstrap_replicates A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param seed A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
summarise_model_predictive_skill <- function(
    predictions, gaussian_scores, bootstrap_replicates = 5000L,
    seed = 20260914L) {
  evaluation <- model_performance_evaluation_rows(predictions) |>
    dplyr::filter(.data$model_key != "Climate")
  group_columns <- c(
    "model_key", "model", "prediction_scope", "evaluation_subset",
    "support_stratum"
  )
  groups <- model_performance_group_indices(evaluation, group_columns)
  set.seed(seed)
  rows <- list()
  for (indices in groups) {
    data <- evaluation[indices, , drop = FALSE]
    group <- data[1L, group_columns, drop = FALSE]
    full_data <- data[data$candidate_role == "full", , drop = FALSE]
    baseline_data <- data[
      data$candidate_role == "matched_baseline", , drop = FALSE
    ]
    if (!identical(
      sort(as.character(full_data$effect_id)),
      sort(as.character(baseline_data$effect_id))
    )) {
      stop(group$model[[1]], " skill rows are not matched.", call. = FALSE)
    }
    gaussian_group <- gaussian_scores
    for (column in group_columns) {
      gaussian_group <- gaussian_group[
        as.character(gaussian_group[[column]]) ==
          as.character(group[[column]][[1]]), , drop = FALSE
      ]
    }
    full_statistics <- model_performance_sufficient_statistics(
      full_data,
      gaussian_group[gaussian_group$candidate_role == "full", , drop = FALSE]
    )
    baseline_statistics <- model_performance_sufficient_statistics(
      baseline_data,
      gaussian_group[
        gaussian_group$candidate_role == "matched_baseline", , drop = FALSE
      ]
    )
    if (!identical(
      full_statistics$heldout_gradient,
      baseline_statistics$heldout_gradient
    )) {
      stop(group$model[[1]], " skill folds are not matched.", call. = FALSE)
    }
    for (aggregation in c("equal_gradient", "effect_pooled")) {
      full_estimate <- aggregate_model_performance_statistics(
        full_statistics, aggregation
      )
      baseline_estimate <- aggregate_model_performance_statistics(
        baseline_statistics, aggregation
      )
      estimate <- model_performance_skill_values(
        full_estimate, baseline_estimate
      )
      draws <- replicate(bootstrap_replicates, {
        sample_index <- sample.int(
          nrow(full_statistics), nrow(full_statistics), replace = TRUE
        )
        model_performance_skill_values(
          aggregate_model_performance_statistics(
            full_statistics, aggregation, sample_index
          ),
          aggregate_model_performance_statistics(
            baseline_statistics, aggregation, sample_index
          )
        )
      })
      intervals <- t(apply(
        draws, 1L, stats::quantile,
        probs = c(0.025, 0.975), names = FALSE, type = 7
      ))
      row <- group[rep(1L, length(estimate)), , drop = FALSE]
      row$aggregation <- aggregation
      row$metric <- names(estimate)
      row$estimate <- unname(estimate)
      row$CI_low <- intervals[, 1L]
      row$CI_high <- intervals[, 2L]
      row$n_effects <- nrow(full_data)
      row$n_species <- dplyr::n_distinct(full_data$species_id)
      row$n_gradients <- nrow(full_statistics)
      row$bootstrap_replicates <- as.integer(bootstrap_replicates)
      row$bootstrap_seed <- as.integer(seed)
      rows[[length(rows) + 1L]] <- row
    }
  }
  dplyr::bind_rows(rows)
}

#' Validate Model Performance.
#'
#' @param folds A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param combined A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param fit_r2 A fitted model object with the structure checked by this function.
#' @param metrics A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param skills A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param expected_gradients A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param expected_fits A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
validate_model_performance <- function(
    folds, combined, fit_r2, metrics, skills,
    expected_gradients = 20L, expected_fits = 200L) {
  fold_names <- sort(unique(as.character(combined$predictions$heldout_gradient)))
  fit_audit <- combined$fit_audit
  prediction_audit <- combined$prediction_audit
  expected_counts <- fit_r2 |>
    dplyr::filter(.data$aggregation == "equal_gradient") |>
    dplyr::select(
      "model_key", "candidate_role", expected_effects = "n_effects"
    )
  general_counts <- combined$predictions |>
    dplyr::filter(.data$prediction_scope == "general") |>
    dplyr::summarise(
      prediction_rows = dplyr::n(),
      unique_effects = dplyr::n_distinct(.data$effect_id),
      heldout_folds = dplyr::n_distinct(.data$heldout_gradient),
      .by = c("model_key", "candidate_role")
    ) |>
    dplyr::left_join(
      expected_counts, by = c("model_key", "candidate_role")
    )
  fit_partition_audit <- fit_audit |>
    dplyr::left_join(
      expected_counts, by = c("model_key", "candidate_role")
    )
  fit_id_counts <- dplyr::count(fit_audit, .data$fit_id)
  repeat_key <- c(
    "heldout_gradient", "model_key", "candidate_role", "effect_id"
  )
  general_repeat <- combined$predictions |>
    dplyr::filter(
      .data$prediction_scope == "general",
      .data$species_seen_in_training
    ) |>
    dplyr::select(dplyr::all_of(repeat_key))
  conditional_repeat <- combined$predictions |>
    dplyr::filter(
      .data$prediction_scope == "repeat_species_conditional"
    ) |>
    dplyr::select(dplyr::all_of(repeat_key))
  unmatched_repeat_rows <- dplyr::bind_rows(
    dplyr::anti_join(general_repeat, conditional_repeat, by = repeat_key),
    dplyr::anti_join(conditional_repeat, general_repeat, by = repeat_key)
  )
  validation <- tibble::tribble(
    ~check, ~pass, ~detail,
    "Twenty leave-gradient-out folds",
    length(fold_names) == expected_gradients,
    paste(length(fold_names), "distinct held-out gradients"),
    "Expected validation-only fits",
    nrow(fit_audit) == expected_fits,
    paste(nrow(fit_audit), "REML fold fits; no canonical fit replaced"),
    "Every eligible effect held out once",
    nrow(general_counts) == nrow(expected_counts) &&
      all(!is.na(general_counts$expected_effects)) &&
      all(general_counts$prediction_rows == general_counts$expected_effects) &&
      all(general_counts$unique_effects == general_counts$expected_effects) &&
      all(general_counts$heldout_folds == expected_gradients),
    paste(nrow(general_counts), "model-role prediction series checked"),
    "Exact train/test partitions",
    all(!is.na(fit_partition_audit$expected_effects)) &&
      all(fit_partition_audit$n_train + fit_partition_audit$n_test ==
            fit_partition_audit$expected_effects) &&
      all(fit_id_counts$n == expected_gradients),
    paste(nrow(fit_id_counts), "fit specifications each span 20 folds"),
    "Fold fit convergence and rank",
    nrow(fit_audit) == expected_fits && all(fit_audit$pass),
    paste(sum(fit_audit$pass), "of", nrow(fit_audit), "fit audits pass"),
    "Disjoint training and held-out effects",
    all(fit_audit$effect_partitions_disjoint),
    paste(sum(fit_audit$effect_partitions_disjoint), "fit partitions checked"),
    "Predictive covariance and unseen-species identity",
    all(prediction_audit$pass),
    paste(sum(prediction_audit$pass), "of", nrow(prediction_audit),
          "prediction audits pass"),
    "Matched repeat-species evaluation rows",
    nrow(general_repeat) == nrow(conditional_repeat) &&
      nrow(unmatched_repeat_rows) == 0L,
    paste(
      nrow(conditional_repeat),
      "conditional predictions have general-model predictions on identical rows"
    ),
    "Finite exported predictions",
    all(is.finite(combined$predictions$observed_delta_auc)) &&
      all(is.finite(combined$predictions$predicted_delta_auc)) &&
      all(is.finite(combined$predictions$predictive_se)) &&
      all(combined$predictions$predictive_se > 0),
    paste(nrow(combined$predictions), "general and repeat-species predictions"),
    "Valid interval fields",
    all(combined$predictions$interval_width > 0) &&
      all(!is.na(combined$predictions$interval_covered)) &&
      all(!is.na(combined$predictions$climate_supported)),
    "All predictive intervals are finite and classifiable",
    "Finite hierarchical R2",
    all(is.finite(fit_r2$marginal_R2)) &&
      all(is.finite(fit_r2$conditional_R2)) &&
      all(fit_r2$marginal_R2 >= 0 & fit_r2$marginal_R2 <= 1) &&
      all(fit_r2$conditional_R2 >= 0 & fit_r2$conditional_R2 <= 1),
    paste(nrow(fit_r2), "model-role and weighting rows"),
    "Finite cross-validated metrics",
    all(is.finite(metrics$estimate)) && all(is.finite(metrics$CI_low)) &&
      all(is.finite(metrics$CI_high)),
    paste(nrow(metrics), "metric rows"),
    "Finite matched-baseline skill",
    all(is.finite(skills$estimate)) && all(is.finite(skills$CI_low)) &&
      all(is.finite(skills$CI_high)),
    paste(nrow(skills), "paired skill rows"),
    "Coverage bounds",
    all(metrics$estimate[metrics$metric == "coverage_95"] >= 0 &
          metrics$estimate[metrics$metric == "coverage_95"] <= 1),
    "Every estimated 95% coverage is between zero and one"
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

#' Assert Model Performance Validation.
#'
#' @param validation A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
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

#' Build Model Performance Figure.
#'
#' @param fit_r2 A fitted model object with the structure checked by this function.
#' @param metrics A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A plot object, or the declared figure path when rendering to disk.
build_model_performance_figure <- function(fit_r2, metrics) {
  key <- model_performance_key()
  selected_roles <- c("standalone", "full")
  sample_key <- fit_r2 |>
    dplyr::filter(
      .data$aggregation == "equal_gradient",
      .data$candidate_role %in% .env$selected_roles
    ) |>
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
  if (anyNA(model_levels)) {
    stop("The performance figure model labels are incomplete.", call. = FALSE)
  }

  fit_data <- fit_r2 |>
    dplyr::filter(
      .data$aggregation == "equal_gradient",
      .data$candidate_role %in% .env$selected_roles
    ) |>
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
      .data$aggregation == "equal_gradient",
      .data$support_stratum == "all",
      .data$prediction_scope == "general",
      .data$evaluation_subset == "all_eligible",
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
  if (nrow(plot_data) != expected_rows || anyNA(plot_data$metric_label) ||
      any(table(plot_data$metric_label) != nrow(key))) {
    stop("The performance figure needs one estimate per model and metric.",
         call. = FALSE)
  }

  metric_colours <- stats::setNames(
    c("#3B4A6B", "#1B9E77", "#D95F02", "#6A3D9A"),
    unname(metric_labels)
  )
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data$estimate, y = .data$model,
      colour = .data$metric_label
    )
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
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.15))) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      x = "Estimate", y = NULL,
      subtitle = paste(
        "Primary estimates give each of 20 gradients equal influence;",
        "whiskers are gradient-bootstrap 95% intervals"
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
