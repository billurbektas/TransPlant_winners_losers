climate_response_formula <- function() {
  stats::as.formula(
    paste(
      "climate_category ~",
      "z_destination_T + z_log_destination_P + z_warming +",
      "z_destination_T:z_log_destination_P +",
      "z_destination_T:z_warming +",
      "z_log_destination_P:z_warming"
    )
  )
}

climate_response_additive_formula <- function() {
  climate_category ~ z_destination_T + z_log_destination_P + z_warming
}

climate_response_term_names <- function() {
  c(
    "z_destination_T", "z_log_destination_P", "z_warming",
    "z_destination_T:z_log_destination_P",
    "z_destination_T:z_warming",
    "z_log_destination_P:z_warming"
  )
}

climate_response_incremental_formulas <- function() {
  list(
    intercept = climate_category ~ 1,
    warming = climate_category ~ z_warming,
    context = climate_category ~ z_destination_T + z_log_destination_P,
    additive = climate_response_additive_formula(),
    context_interaction = climate_category ~ z_destination_T + z_log_destination_P +
      z_warming + z_destination_T:z_log_destination_P,
    full_pairwise = climate_response_formula()
  )
}

prepare_climate_response_data <- function(data) {
  required <- c(
    "model_category", "accepted_name", "experiment_id",
    "destination_T_JJA", "destination_P_JJA", "diff_T_JJA"
  )
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(
      "Climate-response data are missing columns: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  data <- filter_downstream_model_species(data, "accepted_name")
  prepared <- data |>
    dplyr::transmute(
      experiment_id = factor(.data$experiment_id),
      accepted_species = factor(.data$accepted_name),
      climate_category = dplyr::case_when(
        .data$model_category %in% c("loser", "winner", "colonizer") ~
          as.character(.data$model_category),
        .data$model_category == "no_change" ~ "non_significant",
        TRUE ~ NA_character_
      ),
      destination_T_JJA = .data$destination_T_JJA,
      destination_P_JJA = .data$destination_P_JJA,
      warming = .data$diff_T_JJA
    ) |>
    dplyr::distinct() |>
    dplyr::mutate(
      climate_category = factor(
        .data$climate_category,
        levels = c("non_significant", "loser", "winner", "colonizer")
      ),
      log_destination_P = dplyr::if_else(
        .data$destination_P_JJA > 0, log(.data$destination_P_JJA), NA_real_
      )
    ) |>
    dplyr::filter(stats::complete.cases(dplyr::pick(dplyr::everything())))
  duplicates <- prepared |>
    dplyr::count(.data$accepted_species, .data$experiment_id) |>
    dplyr::filter(.data$n != 1L)
  if (nrow(duplicates)) {
    stop("Climate-response data must contain one row per species and experiment.", call. = FALSE)
  }
  expected_categories <- c("non_significant", "loser", "winner", "colonizer")
  absent <- setdiff(expected_categories, unique(as.character(prepared$climate_category)))
  if (length(absent)) {
    stop("Climate-response categories are absent: ", paste(absent, collapse = ", "), call. = FALSE)
  }
  scaling_columns <- c(
    destination_T = "destination_T_JJA",
    log_destination_P = "log_destination_P",
    warming = "warming"
  )
  scaling <- purrr::imap_dfr(scaling_columns, function(column, variable) {
    values <- prepared[[column]]
    tibble::tibble(
      variable = variable,
      source_column = column,
      center = mean(values),
      scale = stats::sd(values)
    )
  })
  if (any(!is.finite(scaling$scale)) || any(scaling$scale <= 0)) {
    stop("Every climate-response predictor must have a positive finite scale.", call. = FALSE)
  }
  for (index in seq_len(nrow(scaling))) {
    prepared[[paste0("z_", scaling$variable[[index]])]] <-
      (prepared[[scaling$source_column[[index]]]] - scaling$center[[index]]) /
      scaling$scale[[index]]
  }
  category_counts <- prepared |>
    dplyr::count(.data$climate_category, name = "n") |>
    dplyr::mutate(proportion = .data$n / sum(.data$n))
  species_replication <- prepared |>
    dplyr::summarise(
      n_experiments = dplyr::n_distinct(.data$experiment_id),
      n_significant = sum(.data$climate_category != "non_significant"),
      n_outcome_categories = dplyr::n_distinct(.data$climate_category),
      .by = "accepted_species"
    )
  snapshot <- tibble::tibble(
    n_species_experiment = nrow(prepared),
    n_experiments = dplyr::n_distinct(prepared$experiment_id),
    n_species = dplyr::n_distinct(prepared$accepted_species),
    n_significant = sum(prepared$climate_category != "non_significant"),
    significant_proportion = mean(prepared$climate_category != "non_significant"),
    species_single_experiment = sum(species_replication$n_experiments == 1L),
    species_at_most_two_experiments = sum(species_replication$n_experiments <= 2L),
    species_one_outcome = sum(species_replication$n_outcome_categories == 1L),
    species_without_significant_response = sum(species_replication$n_significant == 0L)
  )
  list(
    data = prepared,
    scaling = scaling,
    category_counts = category_counts,
    species_replication = species_replication,
    snapshot = snapshot
  )
}

capture_climate_fit <- function(expression) {
  warnings <- character()
  error <- NULL
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

fit_climate_response_model <- function(prepared) {
  captured <- capture_climate_fit(
    mclogit::mblogit(
      climate_response_formula(),
      data = prepared$data,
      trace = FALSE,
      control = mclogit::mclogit.control(maxit = 100, trace = FALSE)
    )
  )
  if (is.null(captured$fit)) {
    stop("The fixed climate-response model failed: ", captured$error, call. = FALSE)
  }
  model <- captured$fit
  covariance <- stats::vcov(model)
  if (!isTRUE(model$converged) || length(captured$warnings) ||
      any(!is.finite(stats::coef(model))) || any(!is.finite(covariance))) {
    stop(
      "The fixed climate-response model did not converge cleanly: ",
      paste(c(captured$warnings, captured$error), collapse = " | "),
      call. = FALSE
    )
  }
  covariance_eigenvalues <- eigen(covariance, symmetric = TRUE, only.values = TRUE)$values
  diagnostics <- dplyr::bind_cols(
    prepared$snapshot,
    tibble::tibble(
      selected_model = "fixed-effects multinomial logit",
      reference_category = "non_significant",
      formula = paste(deparse(climate_response_formula()), collapse = " "),
      converged = model$converged,
      iterations = model$iter,
      deviance = model$deviance,
      coefficient_covariance_minimum_eigenvalue = min(covariance_eigenvalues),
      warning = "",
      limitation = paste(
        "No random effect or cluster correction is used; residual dependence among",
        "repeated observations of the same species is not modelled."
      )
    )
  )
  list(
    model = model,
    data = prepared$data,
    scaling = prepared$scaling,
    category_counts = prepared$category_counts,
    diagnostics = diagnostics,
    formula = climate_response_formula()
  )
}

mixed_fit_summary <- function(captured, fit_name, formula_stage, fixed_model = NULL) {
  if (is.null(captured$fit)) {
    return(tibble::tibble(
      fit_name = fit_name, formula_stage = formula_stage,
      outer_converged = FALSE, iterations = NA_integer_, deviance = NA_real_,
      coefficients_finite = FALSE, covariance_finite = FALSE,
      random_covariance_minimum_eigenvalue = NA_real_,
      maximum_fixed_coefficient_difference = NA_real_, deviance_difference = NA_real_,
      warning = paste(c(captured$warnings, captured$error), collapse = " | ")
    ))
  }
  fit <- captured$fit
  covariance <- tryCatch(stats::vcov(fit), error = function(condition) NULL)
  random_eigenvalues <- unlist(lapply(
    fit$VarCov,
    function(matrix) eigen(matrix, symmetric = TRUE, only.values = TRUE)$values
  ))
  tibble::tibble(
    fit_name = fit_name,
    formula_stage = formula_stage,
    outer_converged = isTRUE(fit$converged),
    iterations = fit$iter,
    deviance = fit$deviance,
    coefficients_finite = all(is.finite(stats::coef(fit))),
    covariance_finite = !is.null(covariance) && all(is.finite(covariance)),
    random_covariance_minimum_eigenvalue = min(random_eigenvalues),
    maximum_fixed_coefficient_difference = if (is.null(fixed_model)) NA_real_ else
      max(abs(stats::coef(fit) - stats::coef(fixed_model$model))),
    deviance_difference = if (is.null(fixed_model)) NA_real_ else
      fit$deviance - fixed_model$model$deviance,
    warning = paste(c(captured$warnings, captured$error), collapse = " | ")
  )
}

fit_climate_response_mixed_diagnostics <- function(prepared, fixed_model) {
  fit_mixed <- function(formula, inner_optimizer = "nlminb", maxit_inner = 100L) {
    capture_climate_fit(
      mclogit::mblogit(
        formula,
        random = ~ 1 | accepted_species,
        data = prepared$data,
        estimator = "ML",
        trace = FALSE,
        control = mclogit::mmclogit.control(
          maxit = 100, maxit.inner = maxit_inner,
          inner.optimizer = inner_optimizer, trace = FALSE
        )
      )
    )
  }
  incremental <- purrr::imap_dfr(
    climate_response_incremental_formulas(),
    function(formula, stage) {
      captured <- fit_mixed(formula, maxit_inner = 1000L)
      mixed_fit_summary(captured, "nlminb_1000", stage)
    }
  )
  nlminb_100 <- fit_mixed(climate_response_formula(), maxit_inner = 100L)
  nlminb_1000 <- fit_mixed(climate_response_formula(), maxit_inner = 1000L)
  bfgs_1000 <- fit_mixed(
    climate_response_formula(), inner_optimizer = "BFGS", maxit_inner = 1000L
  )
  optimizer_comparison <- dplyr::bind_rows(
    mixed_fit_summary(nlminb_100, "nlminb_100", "full_pairwise", fixed_model),
    mixed_fit_summary(nlminb_1000, "nlminb_1000", "full_pairwise", fixed_model),
    mixed_fit_summary(bfgs_1000, "BFGS_1000", "full_pairwise", fixed_model)
  ) |>
    dplyr::mutate(
      acceptable_for_inference =
        .data$outer_converged & .data$coefficients_finite & .data$covariance_finite &
        .data$warning == "" & .data$random_covariance_minimum_eigenvalue > 1e-6
    )
  decision <- tibble::tibble(
    selected_model = "fixed-effects multinomial logit",
    mixed_model_accepted = FALSE,
    reason = paste(
      "The default nlminb inner optimizer reports false convergence at both 100 and",
      "1000 iterations. BFGS reproduces the fixed model while the free species",
      "random covariance collapses to the numerical boundary. Sparse within-species",
      "replication does not support category-specific random intercept covariance."
    )
  )
  list(
    nlminb_100 = nlminb_100$fit,
    nlminb_1000 = nlminb_1000$fit,
    bfgs_1000 = bfgs_1000$fit,
    incremental = incremental,
    optimizer_comparison = optimizer_comparison,
    species_replication = prepared$species_replication,
    decision = decision
  )
}

tidy_climate_response_coefficients <- function(fit) {
  coefficients <- stats::coef(fit$model)
  standard_errors <- sqrt(diag(stats::vcov(fit$model)))
  result <- tibble::tibble(
    coefficient = names(coefficients),
    estimate = unname(coefficients),
    standard_error = unname(standard_errors)
  ) |>
    tidyr::separate_wider_delim(
      "coefficient", delim = "~", names = c("category_vs_non_significant", "term")
    ) |>
    dplyr::mutate(
      z = .data$estimate / .data$standard_error,
      p_value = 2 * stats::pnorm(abs(.data$z), lower.tail = FALSE),
      conf_low = .data$estimate - stats::qnorm(0.975) * .data$standard_error,
      conf_high = .data$estimate + stats::qnorm(0.975) * .data$standard_error,
      odds_ratio = exp(.data$estimate),
      odds_low = exp(.data$conf_low),
      odds_high = exp(.data$conf_high),
      raw_significant = .data$p_value < 0.05
    )
  non_intercept <- which(result$term != "(Intercept)")
  result$p_holm <- NA_real_
  result$p_holm[non_intercept] <- stats::p.adjust(
    result$p_value[non_intercept], method = "holm"
  )
  result
}

climate_response_wald_test <- function(coefficients, covariance, contrast) {
  estimate <- as.numeric(contrast %*% coefficients)
  contrast_covariance <- contrast %*% covariance %*% t(contrast)
  if (any(!is.finite(contrast_covariance)) ||
      qr(contrast_covariance)$rank != nrow(contrast)) {
    stop("Climate-response Wald contrast covariance is singular.", call. = FALSE)
  }
  statistic <- drop(t(estimate) %*% solve(contrast_covariance, estimate))
  degrees_freedom <- nrow(contrast)
  tibble::tibble(
    statistic = statistic,
    degrees_freedom = degrees_freedom,
    p_value = stats::pchisq(statistic, df = degrees_freedom, lower.tail = FALSE)
  )
}

climate_response_omnibus_tests <- function(fit) {
  coefficients <- stats::coef(fit$model)
  covariance <- stats::vcov(fit$model)
  coefficient_names <- names(coefficients)
  categories <- c("loser", "winner", "colonizer")
  terms <- climate_response_term_names()

  outcome_distribution <- purrr::map_dfr(terms, function(term) {
    selected <- match(paste0(categories, "~", term), coefficient_names)
    if (anyNA(selected)) stop("Missing multinomial coefficients for ", term, call. = FALSE)
    contrast <- matrix(0, nrow = 3L, ncol = length(coefficients))
    contrast[cbind(seq_len(3L), selected)] <- 1
    climate_response_wald_test(coefficients, covariance, contrast) |>
      dplyr::mutate(
        test_family = "outcome_distribution",
        term = .env$term,
        null_hypothesis = "loser = winner = colonizer = 0",
        .before = 1
      )
  }) |>
    dplyr::mutate(p_holm = stats::p.adjust(.data$p_value, method = "holm"))

  significant_composition <- purrr::map_dfr(terms, function(term) {
    selected <- match(paste0(categories, "~", term), coefficient_names)
    if (anyNA(selected)) stop("Missing multinomial coefficients for ", term, call. = FALSE)
    contrast <- matrix(0, nrow = 2L, ncol = length(coefficients))
    contrast[1, selected[c(1, 2)]] <- c(1, -1)
    contrast[2, selected[c(1, 3)]] <- c(1, -1)
    climate_response_wald_test(coefficients, covariance, contrast) |>
      dplyr::mutate(
        test_family = "significant_composition",
        term = .env$term,
        null_hypothesis = "loser = winner = colonizer",
        .before = 1
      )
  }) |>
    dplyr::mutate(p_holm = stats::p.adjust(.data$p_value, method = "holm"))

  dplyr::bind_rows(outcome_distribution, significant_composition)
}

fit_climate_response_nested <- function(formula, data, model_name) {
  captured <- capture_climate_fit(
    mclogit::mblogit(
      formula,
      data = data,
      trace = FALSE,
      control = mclogit::mclogit.control(maxit = 100, trace = FALSE)
    )
  )
  if (is.null(captured$fit) || !isTRUE(captured$fit$converged) ||
      length(captured$warnings)) {
    stop(
      "Nested climate-response model failed (", model_name, "): ",
      paste(c(captured$warnings, captured$error), collapse = " | "),
      call. = FALSE
    )
  }
  captured$fit
}

climate_response_likelihood_ratio_tests <- function(fit) {
  intercept <- fit_climate_response_nested(
    climate_category ~ 1, fit$data, "intercept_only"
  )
  additive <- fit_climate_response_nested(
    climate_response_additive_formula(), fit$data, "additive"
  )
  full <- fit$model
  comparison <- function(reduced, full, comparison_name, null_hypothesis) {
    statistic <- reduced$deviance - full$deviance
    degrees_freedom <- length(stats::coef(full)) - length(stats::coef(reduced))
    tibble::tibble(
      comparison = comparison_name,
      null_hypothesis = null_hypothesis,
      reduced_deviance = reduced$deviance,
      full_deviance = full$deviance,
      statistic = statistic,
      degrees_freedom = degrees_freedom,
      p_value = stats::pchisq(
        statistic, df = degrees_freedom, lower.tail = FALSE
      ),
      reduced_converged = isTRUE(reduced$converged),
      full_converged = isTRUE(full$converged)
    )
  }
  dplyr::bind_rows(
    comparison(
      intercept, full, "full_vs_intercept",
      "All 18 non-intercept coefficients are zero"
    ),
    comparison(
      additive, full, "full_pairwise_vs_additive",
      "All nine pairwise-interaction coefficients are zero"
    )
  )
}

climate_response_model_performance <- function(fit) {
  probabilities <- stats::predict(fit$model, type = "response")
  observed <- as.character(fit$data$climate_category)
  observed_index <- match(observed, colnames(probabilities))
  observed_probability <- probabilities[
    cbind(seq_len(nrow(probabilities)), observed_index)
  ]
  one_hot <- matrix(
    0, nrow = nrow(probabilities), ncol = ncol(probabilities),
    dimnames = list(NULL, colnames(probabilities))
  )
  one_hot[cbind(seq_len(nrow(one_hot)), observed_index)] <- 1
  rhs <- stats::delete.response(stats::terms(fit$model))
  design <- stats::model.matrix(
    rhs, data = fit$data, contrasts.arg = fit$model$contrasts,
    xlev = fit$model$xlevels
  )
  information <- as.matrix(fit$model$information.matrix)
  information_eigenvalues <- eigen(
    information, symmetric = TRUE, only.values = TRUE
  )$values
  covariance_eigenvalues <- eigen(
    stats::vcov(fit$model), symmetric = TRUE, only.values = TRUE
  )$values
  metrics <- tibble::tibble(
    n_species_experiment = nrow(probabilities),
    n_experiments = dplyr::n_distinct(fit$data$experiment_id),
    n_species = dplyr::n_distinct(fit$data$accepted_species),
    multiclass_brier_score = mean(rowSums((one_hot - probabilities)^2)),
    multinomial_log_loss = -mean(log(pmax(observed_probability, .Machine$double.eps))),
    model_matrix_condition_number = kappa(design, exact = TRUE),
    information_condition_number = kappa(information, exact = TRUE),
    information_minimum_eigenvalue = min(information_eigenvalues),
    covariance_minimum_eigenvalue = min(covariance_eigenvalues),
    minimum_fitted_probability = min(probabilities),
    maximum_fitted_probability = max(probabilities),
    limitation = paste(
      "Wald inference treats species-experiment rows as independent; repeated",
      "species and shared experiments are not represented in the covariance."
    )
  )
  calibration <- purrr::map_dfr(colnames(probabilities), function(category) {
    tibble::tibble(
      category = category,
      predicted_probability = probabilities[, category],
      observed = as.integer(observed == category)
    ) |>
      dplyr::mutate(probability_decile = dplyr::ntile(.data$predicted_probability, 10L)) |>
      dplyr::summarise(
        n = dplyr::n(),
        mean_predicted_probability = mean(.data$predicted_probability),
        observed_proportion = mean(.data$observed),
        .by = c("category", "probability_decile")
      )
  })
  list(metrics = metrics, calibration = calibration)
}

climate_response_leave_one_experiment_out <- function(fit) {
  full_coefficients <- stats::coef(fit$model)
  non_intercept <- !stringr::str_detect(names(full_coefficients), "[(]Intercept[)]")
  experiment_profiles <- fit$data |>
    dplyr::distinct(
      .data$experiment_id, .data$z_destination_T,
      .data$z_log_destination_P, .data$z_warming
    )
  full_predictions <- stats::predict(
    fit$model, newdata = experiment_profiles, type = "response"
  )
  purrr::map_dfr(levels(fit$data$experiment_id), function(omitted_experiment) {
    refit_data <- fit$data |>
      dplyr::filter(as.character(.data$experiment_id) != .env$omitted_experiment) |>
      droplevels()
    captured <- capture_climate_fit(
      mclogit::mblogit(
        climate_response_formula(), data = refit_data, trace = FALSE,
        control = mclogit::mclogit.control(maxit = 100, trace = FALSE)
      )
    )
    if (is.null(captured$fit)) {
      return(tibble::tibble(
        omitted_experiment = omitted_experiment,
        converged = FALSE,
        warning = paste(c(captured$warnings, captured$error), collapse = " | "),
        maximum_absolute_coefficient_difference = NA_real_,
        coefficient_sign_changes = NA_integer_,
        maximum_absolute_probability_difference = NA_real_,
        probability_rmse = NA_real_
      ))
    }
    refit_coefficients <- stats::coef(captured$fit)[names(full_coefficients)]
    refit_predictions <- stats::predict(
      captured$fit, newdata = experiment_profiles, type = "response"
    )
    probability_difference <- refit_predictions - full_predictions
    tibble::tibble(
      omitted_experiment = omitted_experiment,
      converged = isTRUE(captured$fit$converged),
      warning = paste(c(captured$warnings, captured$error), collapse = " | "),
      maximum_absolute_coefficient_difference = max(abs(
        refit_coefficients - full_coefficients
      )),
      coefficient_sign_changes = sum(
        sign(refit_coefficients[non_intercept]) !=
          sign(full_coefficients[non_intercept])
      ),
      maximum_absolute_probability_difference = max(abs(probability_difference)),
      probability_rmse = sqrt(mean(probability_difference^2))
    )
  })
}

climate_response_scale_value <- function(fit, variable, value) {
  definition <- fit$scaling |>
    dplyr::filter(.data$variable == .env$variable)
  if (nrow(definition) != 1L) {
    stop("Missing scaling definition for ", variable, call. = FALSE)
  }
  (value - definition$center[[1]]) / definition$scale[[1]]
}

climate_response_support_mask <- function(grid, experiments, climate_contexts) {
  cutoffs <- climate_contexts |>
    dplyr::distinct(
      .data$temperature_cutoff_low, .data$temperature_cutoff_high,
      .data$precipitation_cutoff_dry, .data$precipitation_cutoff_wet
    )
  if (nrow(cutoffs) != 1L) {
    stop("Climate-response prediction cutoffs must be unique.", call. = FALSE)
  }
  grid <- grid |>
    dplyr::mutate(
      climate_context = dplyr::case_when(
        .data$destination_P_JJA < cutoffs$precipitation_cutoff_dry &
          .data$destination_T_JJA < cutoffs$temperature_cutoff_low ~ "cold_dry",
        .data$destination_P_JJA < cutoffs$precipitation_cutoff_dry &
          .data$destination_T_JJA < cutoffs$temperature_cutoff_high ~ "intermediate_dry",
        .data$destination_P_JJA < cutoffs$precipitation_cutoff_dry ~ "warm_dry",
        .data$destination_P_JJA < cutoffs$precipitation_cutoff_wet &
          .data$destination_T_JJA < cutoffs$temperature_cutoff_high ~ "cold_intermediate",
        .data$destination_P_JJA >= cutoffs$precipitation_cutoff_wet &
          .data$destination_T_JJA < cutoffs$temperature_cutoff_high ~ "cold_wet",
        TRUE ~ NA_character_
      ),
      reviewed_climate_domain = !is.na(.data$climate_context)
    )
  hull_support <- purrr::map_dfr(unique(grid$warming_level), function(level) {
    focal_grid <- dplyr::filter(grid, .data$warming_level == .env$level)
    focal_experiments <- dplyr::filter(
      experiments, .data$warming_level == .env$level
    )
    if (nrow(focal_experiments) < 3L) {
      stop("At least three experiments are required for each warming support hull.", call. = FALSE)
    }
    hull <- focal_experiments |>
      sf::st_as_sf(coords = c("destination_T_JJA", "destination_P_JJA")) |>
      sf::st_union() |>
      sf::st_convex_hull()
    grid_geometry <- focal_grid |>
      sf::st_as_sf(coords = c("destination_T_JJA", "destination_P_JJA"))
    focal_grid |>
      dplyr::transmute(
        grid_id = .data$grid_id,
        warming_hull_support = lengths(sf::st_intersects(grid_geometry, hull)) > 0L
      )
  })
  grid |>
    dplyr::left_join(hull_support, by = "grid_id") |>
    dplyr::mutate(
      supported = .data$reviewed_climate_domain & .data$warming_hull_support,
      temperature_cutoff_low = cutoffs$temperature_cutoff_low,
      temperature_cutoff_high = cutoffs$temperature_cutoff_high,
      precipitation_cutoff_dry = cutoffs$precipitation_cutoff_dry,
      precipitation_cutoff_wet = cutoffs$precipitation_cutoff_wet
    )
}

climate_response_probability_matrices <- function(design, coefficient_draws, categories) {
  category_eta <- purrr::map(categories, function(category) {
    selected <- stringr::str_starts(colnames(coefficient_draws), paste0(category, "~"))
    design %*% t(coefficient_draws[, selected, drop = FALSE])
  })
  names(category_eta) <- categories
  maximum_eta <- pmax(
    0, category_eta$loser, category_eta$winner, category_eta$colonizer
  )
  reference_exp <- exp(-maximum_eta)
  category_exp <- purrr::map(category_eta, ~ exp(.x - maximum_eta))
  denominator <- reference_exp + category_exp$loser +
    category_exp$winner + category_exp$colonizer
  probabilities <- list(
    non_significant = reference_exp / denominator,
    loser = category_exp$loser / denominator,
    winner = category_exp$winner / denominator,
    colonizer = category_exp$colonizer / denominator
  )
  probabilities$any_significant <- 1 - probabilities$non_significant
  probabilities$loser_given_significant <-
    probabilities$loser / probabilities$any_significant
  probabilities$winner_given_significant <-
    probabilities$winner / probabilities$any_significant
  probabilities$colonizer_given_significant <-
    probabilities$colonizer / probabilities$any_significant
  probabilities
}

climate_response_prediction_surfaces <- function(
    fit, climate_contexts, grid_size = 60L, n_draws = 2000L, seed = 47163L,
    chunk_size = 180L) {
  experiments <- fit$data |>
    dplyr::transmute(
      experiment_id = as.character(.data$experiment_id),
      destination_T_JJA = .data$destination_T_JJA,
      destination_P_JJA = .data$destination_P_JJA,
      warming = .data$warming
    ) |>
    dplyr::distinct()
  experiment_duplicates <- experiments |>
    dplyr::count(.data$experiment_id) |>
    dplyr::filter(.data$n != 1L)
  if (nrow(experiment_duplicates)) {
    stop("Each experiment must have one climate profile for prediction.", call. = FALSE)
  }
  warming_quantiles <- stats::quantile(
    experiments$warming, probs = c(0.25, 0.50, 0.75), names = FALSE
  )
  warming_tertiles <- stats::quantile(
    experiments$warming, probs = c(1 / 3, 2 / 3), names = FALSE
  )
  warming_levels <- tibble::tibble(
    warming_level = c("Low warming", "Median warming", "High warming"),
    warming_quantile = c(0.25, 0.50, 0.75),
    warming = warming_quantiles,
    support_lower = c(-Inf, warming_tertiles[[1]], warming_tertiles[[2]]),
    support_upper = c(warming_tertiles[[1]], warming_tertiles[[2]], Inf)
  ) |>
    dplyr::mutate(
      z_warming = climate_response_scale_value(fit, "warming", .data$warming),
      warming_label = factor(
        sprintf(
          "%s\n%.2f °C", .data$warming_level, .data$warming
        ),
        levels = sprintf(
          "%s\n%.2f °C", .data$warming_level, .data$warming
        )
      )
    )
  experiments <- experiments |>
    dplyr::mutate(
      warming_level = dplyr::case_when(
        .data$warming <= warming_tertiles[[1]] ~ "Low warming",
        .data$warming <= warming_tertiles[[2]] ~ "Median warming",
        TRUE ~ "High warming"
      )
    ) |>
    dplyr::left_join(
      dplyr::select(warming_levels, "warming_level", "warming_label"),
      by = "warming_level"
    )

  grid <- tidyr::expand_grid(
    destination_T_JJA = seq(
      min(experiments$destination_T_JJA), max(experiments$destination_T_JJA),
      length.out = grid_size
    ),
    destination_P_JJA = seq(
      min(experiments$destination_P_JJA), max(experiments$destination_P_JJA),
      length.out = grid_size
    ),
    warming_level = warming_levels$warming_level
  ) |>
    dplyr::left_join(warming_levels, by = "warming_level") |>
    dplyr::mutate(
      grid_id = dplyr::row_number(),
      log_destination_P = log(.data$destination_P_JJA),
      z_destination_T = climate_response_scale_value(
        fit, "destination_T", .data$destination_T_JJA
      ),
      z_log_destination_P = climate_response_scale_value(
        fit, "log_destination_P", .data$log_destination_P
      )
    )
  grid <- climate_response_support_mask(grid, experiments, climate_contexts)

  coefficients <- stats::coef(fit$model)
  covariance <- stats::vcov(fit$model)
  set.seed(seed)
  coefficient_draws <- matrix(
    stats::rnorm(n_draws * length(coefficients)), nrow = n_draws
  ) %*% chol(covariance)
  coefficient_draws <- sweep(coefficient_draws, 2, coefficients, "+")
  colnames(coefficient_draws) <- names(coefficients)
  central_coefficients <- matrix(coefficients, nrow = 1)
  colnames(central_coefficients) <- names(coefficients)
  categories <- c("loser", "winner", "colonizer")
  rhs <- stats::delete.response(stats::terms(fit$model))
  chunks <- split(seq_len(nrow(grid)), ceiling(seq_len(nrow(grid)) / chunk_size))
  prediction_grid <- purrr::map_dfr(chunks, function(rows) {
    focal <- grid[rows, , drop = FALSE]
    design <- stats::model.matrix(
      rhs, data = focal, contrasts.arg = fit$model$contrasts,
      xlev = fit$model$xlevels
    )
    draw_probabilities <- climate_response_probability_matrices(
      design, coefficient_draws, categories
    )
    central_probabilities <- climate_response_probability_matrices(
      design, central_coefficients, categories
    )
    purrr::imap_dfr(draw_probabilities, function(draw_matrix, estimand) {
      intervals <- matrixStats::rowQuantiles(
        draw_matrix, probs = c(0.025, 0.975), drop = FALSE
      )
      focal |>
        dplyr::mutate(
          estimand = .env$estimand,
          estimate = as.numeric(.env$central_probabilities[[.env$estimand]][, 1]),
          conf_low = intervals[, 1],
          conf_high = intervals[, 2],
          interval_width = .data$conf_high - .data$conf_low
        )
    })
  })
  list(
    grid = prediction_grid,
    warming_levels = warming_levels,
    experiment_support = experiments,
    simulation = tibble::tibble(
      n_draws = as.integer(n_draws), seed = as.integer(seed),
      grid_size = as.integer(grid_size), chunk_size = as.integer(chunk_size),
      covariance = "independence-based model Wald covariance"
    )
  )
}

climate_response_odds_ratio_plot <- function(coefficients, facet_ncol = 1L) {
  term_labels <- c(
    z_destination_T = "Destination temperature",
    z_log_destination_P = "Destination precipitation",
    z_warming = "Warming",
    `z_destination_T:z_log_destination_P` = "Temperature × precipitation",
    `z_destination_T:z_warming` = "Temperature × warming",
    `z_log_destination_P:z_warming` = "Precipitation × warming"
  )
  category_labels <- c(loser = "Loser", winner = "Winner", colonizer = "Colonizer")
  colours <- c(loser = "#E66101", winner = "#1B9E77", colonizer = "#5E3C99")
  plot_data <- coefficients |>
    dplyr::filter(.data$term != "(Intercept)") |>
    dplyr::mutate(
      term_label = factor(
        unname(term_labels[.data$term]), levels = rev(unname(term_labels))
      ),
      category_label = factor(
        unname(category_labels[.data$category_vs_non_significant]),
        levels = unname(category_labels)
      )
    )
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      .data$odds_ratio, .data$term_label,
      colour = .data$category_vs_non_significant
    )
  ) +
    ggplot2::geom_vline(xintercept = 1, colour = "grey55", linewidth = 0.65) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$odds_low, xmax = .data$odds_high),
      orientation = "y", width = 0, linewidth = 0.8
    ) +
    ggplot2::geom_point(
      ggplot2::aes(shape = .data$raw_significant), size = 3.5, stroke = 1.1
    ) +
    ggplot2::facet_wrap(~ category_label, ncol = facet_ncol) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_colour_manual(values = colours, guide = "none") +
    ggplot2::scale_shape_manual(
      values = c(`FALSE` = 1, `TRUE` = 16),
      labels = c(`FALSE` = "p >= 0.05", `TRUE` = "p < 0.05"),
      name = NULL
    ) +
    ggplot2::theme_bw(base_size = 15) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 14),
      legend.position = "bottom"
    ) +
    ggplot2::labs(x = "Odds ratio versus non-significant response", y = NULL)
}

plot_climate_response_odds_ratios <- function(coefficients, path) {
  plot <- climate_response_odds_ratio_plot(coefficients, facet_ncol = 1L)
  save_plot_atomic(plot, path, 10, 10)
}

# =============================================================================
# Pool-specific climate-response analysis
#
# The functions above are retained as the archived single-multinomial
# implementation. Maintained targets use the pool-specific functions below.
# =============================================================================

pool_climate_formula <- function(response) {
  stats::as.formula(paste(
    response, "~ z_destination_T + z_log_destination_P + z_warming +",
    "z_destination_T:z_log_destination_P +",
    "z_destination_T:z_warming + z_log_destination_P:z_warming"
  ))
}

pool_climate_additive_formula <- function(response) {
  stats::as.formula(paste(
    response, "~ z_destination_T + z_log_destination_P + z_warming"
  ))
}

pool_climate_intercept_formula <- function(response) {
  stats::as.formula(paste(response, "~ 1"))
}

prepare_pool_climate_response_data <- function(data) {
  required <- c(
    "accepted_name", "experiment_id", "pool", "category", "analysis_category",
    "destination_T_JJA", "destination_P_JJA", "diff_T_JJA"
  )
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(
      "Pool climate-response data are missing columns: ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  data <- filter_downstream_model_species(data, "accepted_name")
  base <- data |>
    dplyr::transmute(
      experiment_id = as.character(.data$experiment_id),
      accepted_species = as.character(.data$accepted_name),
      pool = as.character(.data$pool),
      category = as.character(.data$category),
      analysis_category = as.character(.data$analysis_category),
      destination_T_JJA = as.numeric(.data$destination_T_JJA),
      destination_P_JJA = as.numeric(.data$destination_P_JJA),
      warming = as.numeric(.data$diff_T_JJA)
    ) |>
    dplyr::distinct()
  duplicate_keys <- base |>
    dplyr::count(.data$accepted_species, .data$experiment_id) |>
    dplyr::filter(.data$n != 1L)
  if (nrow(duplicate_keys)) {
    stop("Pool models require one row per species and experiment.", call. = FALSE)
  }
  if (any(!base$pool %in% c("resident", "colonizer", "outsider"))) {
    stop("Unexpected pool in climate-response data.", call. = FALSE)
  }
  experiments <- base |>
    dplyr::select(
      "experiment_id", "destination_T_JJA", "destination_P_JJA", "warming"
    ) |>
    dplyr::distinct()
  experiment_duplicates <- experiments |>
    dplyr::count(.data$experiment_id) |>
    dplyr::filter(.data$n != 1L)
  if (nrow(experiment_duplicates)) {
    stop("Each experiment must have exactly one climate profile.", call. = FALSE)
  }
  experiments <- experiments |>
    dplyr::mutate(
      log_destination_P = dplyr::if_else(
        .data$destination_P_JJA > 0,
        log(.data$destination_P_JJA), NA_real_
      )
    )
  if (any(!stats::complete.cases(experiments))) {
    stop("Experiment climate profiles must be complete and positive for precipitation.", call. = FALSE)
  }
  scaling_columns <- c(
    destination_T = "destination_T_JJA",
    log_destination_P = "log_destination_P",
    warming = "warming"
  )
  scaling <- purrr::imap_dfr(scaling_columns, function(column, variable) {
    tibble::tibble(
      variable = variable,
      source_column = column,
      center = mean(experiments[[column]]),
      scale = stats::sd(experiments[[column]]),
      weighting = "42 distinct experiments"
    )
  })
  if (any(!is.finite(scaling$scale)) || any(scaling$scale <= 0)) {
    stop("Every pool-model predictor must have a positive finite scale.", call. = FALSE)
  }
  add_scaled_predictors <- function(frame) {
    frame <- frame |>
      dplyr::mutate(log_destination_P = log(.data$destination_P_JJA))
    for (index in seq_len(nrow(scaling))) {
      frame[[paste0("z_", scaling$variable[[index]])]] <-
        (frame[[scaling$source_column[[index]]]] - scaling$center[[index]]) /
        scaling$scale[[index]]
    }
    frame
  }
  base <- add_scaled_predictors(base) |>
    dplyr::mutate(
      experiment_id = factor(.data$experiment_id),
      accepted_species = factor(.data$accepted_species)
    )
  experiments <- add_scaled_predictors(experiments)
  resident <- base |>
    dplyr::filter(.data$pool == "resident") |>
    dplyr::mutate(
      resident_outcome = factor(
        dplyr::case_when(
          .data$category == "loser" ~ "loser",
          .data$category == "winner" ~ "winner",
          TRUE ~ "non_significant"
        ),
        levels = c("non_significant", "loser", "winner")
      )
    )
  nonresident <- base |>
    dplyr::filter(.data$pool %in% c("colonizer", "outsider")) |>
    dplyr::mutate(
      colonizer_event = .data$analysis_category == "colonizer"
    )
  pool_by_species <- base |>
    dplyr::distinct(.data$accepted_species, .data$pool) |>
    dplyr::summarise(
      n_pools = dplyr::n_distinct(.data$pool),
      pools = paste(sort(unique(.data$pool)), collapse = " + "),
      .by = "accepted_species"
    )
  category_counts <- dplyr::bind_rows(
    resident |>
      dplyr::count(outcome = as.character(.data$resident_outcome), name = "n") |>
      dplyr::mutate(model = "resident_multinomial", pool_denominator = nrow(resident)),
    nonresident |>
      dplyr::count(
        outcome = dplyr::if_else(
          .data$colonizer_event, "colonizer", "non_significant"
        ),
        name = "n"
      ) |>
      dplyr::mutate(model = "nonresident_binomial", pool_denominator = nrow(nonresident))
  ) |>
    dplyr::mutate(proportion = .data$n / .data$pool_denominator) |>
    dplyr::select("model", "outcome", "n", "pool_denominator", "proportion")
  snapshot <- tibble::tibble(
    n_species_experiment = nrow(base),
    n_experiments = dplyr::n_distinct(base$experiment_id),
    n_species = dplyr::n_distinct(base$accepted_species),
    n_resident = nrow(resident),
    n_nonresident = nrow(nonresident),
    n_loser = sum(resident$resident_outcome == "loser"),
    n_winner = sum(resident$resident_outcome == "winner"),
    n_colonizer = sum(nonresident$colonizer_event),
    n_significant_outsider = sum(
      nonresident$pool == "outsider" & nonresident$colonizer_event
    ),
    species_in_multiple_pools = sum(pool_by_species$n_pools > 1L)
  )
  list(
    data = base, resident = resident, nonresident = nonresident,
    experiments = experiments, scaling = scaling,
    category_counts = category_counts, pool_by_species = pool_by_species,
    snapshot = snapshot
  )
}

fit_pool_primary_components <- function(resident, nonresident) {
  resident_captured <- capture_climate_fit(
    mclogit::mblogit(
      pool_climate_formula("resident_outcome"), data = resident,
      trace = FALSE,
      control = mclogit::mclogit.control(maxit = 100, trace = FALSE)
    )
  )
  colonizer_captured <- capture_climate_fit(
    stats::glm(
      pool_climate_formula("colonizer_event"), data = nonresident,
      family = stats::binomial()
    )
  )
  if (is.null(resident_captured$fit) || is.null(colonizer_captured$fit)) {
    stop(
      "A primary pool model failed: ",
      paste(
        c(resident_captured$error, colonizer_captured$error),
        collapse = " | "
      ), call. = FALSE
    )
  }
  resident_model <- resident_captured$fit
  colonizer_model <- colonizer_captured$fit
  warnings <- c(resident_captured$warnings, colonizer_captured$warnings)
  if (!isTRUE(resident_model$converged) || !isTRUE(colonizer_model$converged) ||
      length(warnings) ||
      any(!is.finite(stats::coef(resident_model))) ||
      any(!is.finite(stats::coef(colonizer_model))) ||
      any(!is.finite(stats::vcov(resident_model))) ||
      any(!is.finite(stats::vcov(colonizer_model)))) {
    stop(
      "Primary pool models did not converge cleanly: ",
      paste(unique(warnings), collapse = " | "), call. = FALSE
    )
  }
  list(resident = resident_model, colonizer = colonizer_model)
}

fit_pool_climate_response_models <- function(prepared) {
  models <- fit_pool_primary_components(prepared$resident, prepared$nonresident)
  diagnostics <- dplyr::bind_rows(
    tibble::tibble(
      model = "resident_multinomial",
      n = nrow(prepared$resident),
      events = sum(prepared$resident$resident_outcome != "non_significant"),
      converged = isTRUE(models$resident$converged),
      iterations = models$resident$iter,
      deviance = models$resident$deviance,
      covariance_minimum_eigenvalue = min(eigen(
        stats::vcov(models$resident), symmetric = TRUE, only.values = TRUE
      )$values)
    ),
    tibble::tibble(
      model = "nonresident_binomial",
      n = nrow(prepared$nonresident),
      events = sum(prepared$nonresident$colonizer_event),
      converged = isTRUE(models$colonizer$converged),
      iterations = models$colonizer$iter,
      deviance = stats::deviance(models$colonizer),
      covariance_minimum_eigenvalue = min(eigen(
        stats::vcov(models$colonizer), symmetric = TRUE, only.values = TRUE
      )$values)
    )
  ) |>
    dplyr::mutate(
      random_effects = "none",
      limitation = paste(
        "Independence-based Wald inference does not account for repeated",
        "species or species sharing an experiment."
      )
    )
  list(
    resident = models$resident, colonizer = models$colonizer,
    resident_data = prepared$resident,
    nonresident_data = prepared$nonresident,
    experiments = prepared$experiments,
    scaling = prepared$scaling,
    diagnostics = diagnostics,
    formulas = list(
      resident = pool_climate_formula("resident_outcome"),
      colonizer = pool_climate_formula("colonizer_event")
    )
  )
}

tidy_pool_climate_coefficients <- function(fit) {
  tidy_component <- function(coefficients, covariance, categories, model, reference) {
    standard_errors <- sqrt(diag(covariance))
    result <- tibble::tibble(
      coefficient = names(coefficients), estimate = unname(coefficients),
      standard_error = unname(standard_errors)
    )
    if (length(categories) > 1L) {
      result <- result |>
        tidyr::separate_wider_delim(
          "coefficient", delim = "~",
          names = c("category_vs_non_significant", "term")
        )
    } else {
      result <- result |>
        dplyr::transmute(
          category_vs_non_significant = categories[[1]],
          term = .data$coefficient,
          estimate = .data$estimate,
          standard_error = .data$standard_error
        )
    }
    result |>
      dplyr::mutate(
        model = .env$model, reference = .env$reference,
        z = .data$estimate / .data$standard_error,
        p_value = 2 * stats::pnorm(abs(.data$z), lower.tail = FALSE),
        conf_low = .data$estimate - stats::qnorm(0.975) * .data$standard_error,
        conf_high = .data$estimate + stats::qnorm(0.975) * .data$standard_error,
        odds_ratio = exp(.data$estimate),
        odds_low = exp(.data$conf_low), odds_high = exp(.data$conf_high),
        raw_significant = .data$p_value < 0.05
      )
  }
  resident <- tidy_component(
    stats::coef(fit$resident), stats::vcov(fit$resident),
    c("loser", "winner"), "resident_multinomial", "non_significant_resident"
  )
  colonizer <- tidy_component(
    stats::coef(fit$colonizer), stats::vcov(fit$colonizer),
    "colonizer", "nonresident_binomial", "non_significant_nonresident"
  )
  dplyr::bind_rows(resident, colonizer) |>
    dplyr::group_by(.data$model) |>
    dplyr::mutate(
      p_holm = dplyr::if_else(
        .data$term == "(Intercept)", NA_real_,
        stats::p.adjust(
          dplyr::if_else(.data$term == "(Intercept)", NA_real_, .data$p_value),
          method = "holm"
        )
      )
    ) |>
    dplyr::ungroup()
}

pool_climate_likelihood_ratio_tests <- function(fit) {
  fit_resident <- function(formula) {
    mclogit::mblogit(
      formula, data = fit$resident_data, trace = FALSE,
      control = mclogit::mclogit.control(maxit = 100, trace = FALSE)
    )
  }
  fit_colonizer <- function(formula) {
    stats::glm(formula, data = fit$nonresident_data, family = stats::binomial())
  }
  nested <- list(
    resident_intercept = fit_resident(pool_climate_intercept_formula("resident_outcome")),
    resident_additive = fit_resident(pool_climate_additive_formula("resident_outcome")),
    colonizer_intercept = fit_colonizer(pool_climate_intercept_formula("colonizer_event")),
    colonizer_additive = fit_colonizer(pool_climate_additive_formula("colonizer_event"))
  )
  comparison <- function(reduced, full, model, name) {
    statistic <- stats::deviance(reduced) - stats::deviance(full)
    degrees_freedom <- length(stats::coef(full)) - length(stats::coef(reduced))
    tibble::tibble(
      model = model, comparison = name,
      reduced_deviance = stats::deviance(reduced),
      full_deviance = stats::deviance(full), statistic = statistic,
      degrees_freedom = degrees_freedom,
      p_value = stats::pchisq(
        statistic, df = degrees_freedom, lower.tail = FALSE
      )
    )
  }
  dplyr::bind_rows(
    comparison(
      nested$resident_intercept, fit$resident,
      "resident_multinomial", "full_vs_intercept"
    ),
    comparison(
      nested$resident_additive, fit$resident,
      "resident_multinomial", "full_pairwise_vs_additive"
    ),
    comparison(
      nested$colonizer_intercept, fit$colonizer,
      "nonresident_binomial", "full_vs_intercept"
    ),
    comparison(
      nested$colonizer_additive, fit$colonizer,
      "nonresident_binomial", "full_pairwise_vs_additive"
    )
  )
}

fit_pool_binary_sensitivity <- function(prepared, variance_tolerance = 1e-6) {
  sensitivity_data <- list(
    loser = prepared$resident |>
      dplyr::mutate(binary_event = .data$resident_outcome == "loser"),
    winner = prepared$resident |>
      dplyr::mutate(binary_event = .data$resident_outcome == "winner"),
    colonizer = prepared$nonresident |>
      dplyr::mutate(binary_event = .data$colonizer_event)
  )
  fixed_formula <- pool_climate_formula("binary_event")
  mixed_formula <- stats::update.formula(
    fixed_formula, . ~ . + (1 | accepted_species)
  )
  results <- purrr::imap(sensitivity_data, function(frame, response) {
    fixed_capture <- capture_climate_fit(
      stats::glm(fixed_formula, data = frame, family = stats::binomial())
    )
    mixed_capture <- capture_climate_fit(
      lme4::glmer(
        mixed_formula, data = frame, family = stats::binomial(),
        control = lme4::glmerControl(
          optimizer = "bobyqa", optCtrl = list(maxfun = 200000)
        )
      )
    )
    if (is.null(fixed_capture$fit)) {
      stop("Fixed binary sensitivity failed for ", response, call. = FALSE)
    }
    mixed <- mixed_capture$fit
    optimizer_message <- if (is.null(mixed)) {
      mixed_capture$error
    } else {
      paste(
        c(
          mixed_capture$warnings,
          unlist(mixed@optinfo$conv$lme4$messages, use.names = FALSE)
        ),
        collapse = " | "
      )
    }
    optimizer_message <- paste(unique(optimizer_message[nzchar(optimizer_message)]), collapse = " | ")
    species_variance <- if (is.null(mixed)) {
      NA_real_
    } else {
      as.data.frame(lme4::VarCorr(mixed))$vcov[[1]]
    }
    mixed_converged <- !is.null(mixed) &&
      isTRUE(mixed@optinfo$conv$opt == 0) &&
      !nzchar(optimizer_message) &&
      all(is.finite(lme4::fixef(mixed))) &&
      all(is.finite(stats::vcov(mixed)))
    singular <- is.null(mixed) || lme4::isSingular(mixed, tol = variance_tolerance)
    use_mixed <- mixed_converged && !singular &&
      is.finite(species_variance) && species_variance >= variance_tolerance
    selected <- if (use_mixed) mixed else fixed_capture$fit
    replication <- frame |>
      dplyr::summarise(
        n_records = dplyr::n(), events = sum(.data$binary_event),
        .by = "accepted_species"
      ) |>
      dplyr::mutate(
        response = .env$response,
        pattern = dplyr::case_when(
          .data$events == 0 ~ "never",
          .data$events == .data$n_records ~ "always",
          TRUE ~ "varies"
        )
      )
    diagnostics <- tibble::tibble(
      response = response,
      n = nrow(frame), events = sum(frame$binary_event),
      n_species = dplyr::n_distinct(frame$accepted_species),
      mixed_formula = paste(deparse(mixed_formula), collapse = " "),
      experiment_random_intercept = FALSE,
      mixed_converged = mixed_converged,
      mixed_singular = singular,
      species_variance = species_variance,
      optimizer_message = optimizer_message,
      selected_model = if (use_mixed) "species_random_binomial" else "fixed_binomial",
      selection_reason = if (use_mixed) {
        "Species-random fit converged with supported nonzero variance."
      } else {
        "Species-random fit failed, was singular, or had variance below 1e-6; fixed fallback used."
      }
    )
    list(
      fixed = fixed_capture$fit, mixed_attempt = mixed,
      selected = selected, selected_is_mixed = use_mixed,
      data = frame, diagnostics = diagnostics, replication = replication
    )
  })
  list(
    models = purrr::map(results, "selected"),
    fixed_models = purrr::map(results, "fixed"),
    mixed_attempts = purrr::map(results, "mixed_attempt"),
    data = sensitivity_data,
    diagnostics = purrr::map_dfr(results, "diagnostics"),
    replication = purrr::map_dfr(results, "replication"),
    selected_is_mixed = purrr::map_lgl(results, "selected_is_mixed"),
    fixed_formula = fixed_formula,
    mixed_formula = mixed_formula
  )
}

tidy_pool_binary_sensitivity <- function(sensitivity) {
  purrr::imap_dfr(sensitivity$models, function(model, response) {
    is_mixed <- sensitivity$selected_is_mixed[[response]]
    coefficients <- if (is_mixed) lme4::fixef(model) else stats::coef(model)
    covariance <- stats::vcov(model)
    standard_errors <- sqrt(diag(as.matrix(covariance)))
    result <- tibble::tibble(
      response = response, term = names(coefficients),
      estimate = unname(coefficients),
      standard_error = unname(standard_errors),
      selected_model = if (is_mixed) "species_random_binomial" else "fixed_binomial"
    ) |>
      dplyr::mutate(
        reference = dplyr::if_else(
          .data$response %in% c("loser", "winner"),
          "all_other_resident_outcomes", "non_significant_nonresident"
        ),
        z = .data$estimate / .data$standard_error,
        p_value = 2 * stats::pnorm(abs(.data$z), lower.tail = FALSE),
        conf_low = .data$estimate - stats::qnorm(0.975) * .data$standard_error,
        conf_high = .data$estimate + stats::qnorm(0.975) * .data$standard_error,
        odds_ratio = exp(.data$estimate),
        odds_low = exp(.data$conf_low), odds_high = exp(.data$conf_high),
        raw_significant = .data$p_value < 0.05
      )
    non_intercept <- result$term != "(Intercept)"
    result$p_holm <- NA_real_
    result$p_holm[non_intercept] <- stats::p.adjust(
      result$p_value[non_intercept], method = "holm"
    )
    result
  })
}

pool_climate_model_performance <- function(fit) {
  resident_probabilities <- stats::predict(fit$resident, type = "response")
  resident_observed <- as.character(fit$resident_data$resident_outcome)
  resident_index <- match(resident_observed, colnames(resident_probabilities))
  resident_observed_probability <- resident_probabilities[
    cbind(seq_len(nrow(resident_probabilities)), resident_index)
  ]
  resident_one_hot <- matrix(
    0, nrow = nrow(resident_probabilities), ncol = ncol(resident_probabilities),
    dimnames = list(NULL, colnames(resident_probabilities))
  )
  resident_one_hot[cbind(seq_len(nrow(resident_one_hot)), resident_index)] <- 1
  colonizer_probability <- stats::predict(fit$colonizer, type = "response")
  colonizer_observed <- as.numeric(fit$nonresident_data$colonizer_event)
  metrics <- dplyr::bind_rows(
    tibble::tibble(
      model = "resident_multinomial", n = nrow(fit$resident_data),
      brier_score = mean(rowSums((resident_probabilities - resident_one_hot)^2)),
      log_loss = -mean(log(pmax(resident_observed_probability, 1e-15))),
      minimum_fitted_probability = min(resident_probabilities),
      maximum_fitted_probability = max(resident_probabilities)
    ),
    tibble::tibble(
      model = "nonresident_binomial", n = nrow(fit$nonresident_data),
      brier_score = mean((colonizer_probability - colonizer_observed)^2),
      log_loss = -mean(
        colonizer_observed * log(pmax(colonizer_probability, 1e-15)) +
          (1 - colonizer_observed) * log(pmax(1 - colonizer_probability, 1e-15))
      ),
      minimum_fitted_probability = min(colonizer_probability),
      maximum_fitted_probability = max(colonizer_probability)
    )
  ) |>
    dplyr::mutate(
      limitation = paste(
        "Apparent performance; observations are reused for estimation and",
        "shared species/experiments are not independent."
      )
    )
  resident_calibration <- purrr::map_dfr(
    colnames(resident_probabilities), function(outcome) {
      tibble::tibble(
        observed = as.numeric(resident_observed == .env$outcome),
        probability = resident_probabilities[, .env$outcome]
      ) |>
        dplyr::mutate(bin = dplyr::ntile(.data$probability, 10L)) |>
        dplyr::summarise(
          n = dplyr::n(), mean_probability = mean(.data$probability),
          observed_fraction = mean(.data$observed), .by = "bin"
        ) |>
        dplyr::mutate(model = "resident_multinomial", outcome = .env$outcome)
    }
  )
  colonizer_calibration <- tibble::tibble(
    observed = colonizer_observed, probability = colonizer_probability
  ) |>
    dplyr::mutate(bin = dplyr::ntile(.data$probability, 10L)) |>
    dplyr::summarise(
      n = dplyr::n(), mean_probability = mean(.data$probability),
      observed_fraction = mean(.data$observed), .by = "bin"
    ) |>
    dplyr::mutate(model = "nonresident_binomial", outcome = "colonizer")
  list(
    metrics = metrics,
    calibration = dplyr::bind_rows(resident_calibration, colonizer_calibration)
  )
}

pool_climate_leave_one_experiment_out <- function(fit) {
  experiments <- sort(unique(c(
    as.character(fit$resident_data$experiment_id),
    as.character(fit$nonresident_data$experiment_id)
  )))
  full_resident_coefficients <- stats::coef(fit$resident)
  full_colonizer_coefficients <- stats::coef(fit$colonizer)
  full_resident_probabilities <- stats::predict(fit$resident, type = "response")
  full_colonizer_probabilities <- stats::predict(fit$colonizer, type = "response")
  purrr::map_dfr(experiments, function(omitted) {
    resident_subset <- dplyr::filter(
      fit$resident_data, as.character(.data$experiment_id) != .env$omitted
    ) |>
      droplevels()
    nonresident_subset <- dplyr::filter(
      fit$nonresident_data, as.character(.data$experiment_id) != .env$omitted
    ) |>
      droplevels()
    captured <- capture_climate_fit(
      fit_pool_primary_components(resident_subset, nonresident_subset)
    )
    if (is.null(captured$fit)) {
      return(tibble::tibble(
        omitted_experiment = omitted,
        model = c("resident_multinomial", "nonresident_binomial"),
        converged = FALSE, warning = captured$error,
        maximum_absolute_coefficient_difference = NA_real_,
        coefficient_sign_changes = NA_integer_,
        maximum_absolute_probability_difference = NA_real_,
        probability_rmse = NA_real_
      ))
    }
    refit <- captured$fit
    resident_coefficients <- stats::coef(refit$resident)
    colonizer_coefficients <- stats::coef(refit$colonizer)
    resident_probabilities <- stats::predict(
      refit$resident, newdata = fit$resident_data, type = "response"
    )
    colonizer_probabilities <- stats::predict(
      refit$colonizer, newdata = fit$nonresident_data, type = "response"
    )
    dplyr::bind_rows(
      tibble::tibble(
        omitted_experiment = omitted, model = "resident_multinomial",
        converged = isTRUE(refit$resident$converged),
        warning = paste(captured$warnings, collapse = " | "),
        maximum_absolute_coefficient_difference = max(abs(
          resident_coefficients - full_resident_coefficients
        )),
        coefficient_sign_changes = sum(
          sign(resident_coefficients) != sign(full_resident_coefficients)
        ),
        maximum_absolute_probability_difference = max(abs(
          resident_probabilities - full_resident_probabilities
        )),
        probability_rmse = sqrt(mean(
          (resident_probabilities - full_resident_probabilities)^2
        ))
      ),
      tibble::tibble(
        omitted_experiment = omitted, model = "nonresident_binomial",
        converged = isTRUE(refit$colonizer$converged),
        warning = paste(captured$warnings, collapse = " | "),
        maximum_absolute_coefficient_difference = max(abs(
          colonizer_coefficients - full_colonizer_coefficients
        )),
        coefficient_sign_changes = sum(
          sign(colonizer_coefficients) != sign(full_colonizer_coefficients)
        ),
        maximum_absolute_probability_difference = max(abs(
          colonizer_probabilities - full_colonizer_probabilities
        )),
        probability_rmse = sqrt(mean(
          (colonizer_probabilities - full_colonizer_probabilities)^2
        ))
      )
    )
  })
}

pool_climate_mvn_draws <- function(coefficients, covariance, n_draws) {
  draws <- matrix(
    stats::rnorm(n_draws * length(coefficients)), nrow = n_draws
  ) %*% chol(covariance)
  draws <- sweep(draws, 2, coefficients, "+")
  colnames(draws) <- names(coefficients)
  draws
}

pool_climate_prediction_surfaces <- function(
    fit, climate_contexts, grid_size = 60L, n_draws = 2000L,
    seed = 47163L, chunk_size = 180L) {
  experiments <- fit$experiments |>
    dplyr::select(
      "experiment_id", "destination_T_JJA", "destination_P_JJA", "warming"
    ) |>
    dplyr::distinct()
  warming_quantiles <- stats::quantile(
    experiments$warming, probs = c(0.25, 0.50, 0.75), names = FALSE
  )
  warming_tertiles <- stats::quantile(
    experiments$warming, probs = c(1 / 3, 2 / 3), names = FALSE
  )
  warming_levels <- tibble::tibble(
    warming_level = c("Low warming", "Median warming", "High warming"),
    warming_quantile = c(0.25, 0.50, 0.75),
    warming = warming_quantiles,
    support_lower = c(-Inf, warming_tertiles[[1]], warming_tertiles[[2]]),
    support_upper = c(warming_tertiles[[1]], warming_tertiles[[2]], Inf)
  ) |>
    dplyr::mutate(
      z_warming = climate_response_scale_value(fit, "warming", .data$warming),
      warming_label = factor(
        sprintf("%s\n%.2f °C", .data$warming_level, .data$warming),
        levels = sprintf("%s\n%.2f °C", .data$warming_level, .data$warming)
      )
    )
  experiments <- experiments |>
    dplyr::mutate(
      warming_level = dplyr::case_when(
        .data$warming <= warming_tertiles[[1]] ~ "Low warming",
        .data$warming <= warming_tertiles[[2]] ~ "Median warming",
        TRUE ~ "High warming"
      )
    ) |>
    dplyr::left_join(
      dplyr::select(warming_levels, "warming_level", "warming_label"),
      by = "warming_level"
    )
  grid <- tidyr::expand_grid(
    destination_T_JJA = seq(
      min(experiments$destination_T_JJA), max(experiments$destination_T_JJA),
      length.out = grid_size
    ),
    destination_P_JJA = seq(
      min(experiments$destination_P_JJA), max(experiments$destination_P_JJA),
      length.out = grid_size
    ),
    warming_level = warming_levels$warming_level
  ) |>
    dplyr::left_join(warming_levels, by = "warming_level") |>
    dplyr::mutate(
      grid_id = dplyr::row_number(),
      log_destination_P = log(.data$destination_P_JJA),
      z_destination_T = climate_response_scale_value(
        fit, "destination_T", .data$destination_T_JJA
      ),
      z_log_destination_P = climate_response_scale_value(
        fit, "log_destination_P", .data$log_destination_P
      )
    )
  grid <- climate_response_support_mask(grid, experiments, climate_contexts)
  resident_coefficients <- stats::coef(fit$resident)
  resident_covariance <- stats::vcov(fit$resident)
  colonizer_coefficients <- stats::coef(fit$colonizer)
  colonizer_covariance <- stats::vcov(fit$colonizer)
  set.seed(seed)
  resident_draws <- pool_climate_mvn_draws(
    resident_coefficients, resident_covariance, n_draws
  )
  colonizer_draws <- pool_climate_mvn_draws(
    colonizer_coefficients, colonizer_covariance, n_draws
  )
  resident_central <- matrix(resident_coefficients, nrow = 1)
  colnames(resident_central) <- names(resident_coefficients)
  colonizer_central <- matrix(colonizer_coefficients, nrow = 1)
  colnames(colonizer_central) <- names(colonizer_coefficients)
  resident_rhs <- stats::delete.response(stats::terms(fit$resident))
  colonizer_rhs <- stats::delete.response(stats::terms(fit$colonizer))
  resident_probabilities <- function(design, draws) {
    loser_columns <- match(
      paste0("loser~", colnames(design)), colnames(draws)
    )
    winner_columns <- match(
      paste0("winner~", colnames(design)), colnames(draws)
    )
    if (anyNA(c(loser_columns, winner_columns))) {
      stop("Resident prediction coefficients do not match the design matrix.", call. = FALSE)
    }
    loser_eta <- design %*% t(draws[, loser_columns, drop = FALSE])
    winner_eta <- design %*% t(draws[, winner_columns, drop = FALSE])
    maximum_eta <- pmax(0, loser_eta, winner_eta)
    reference_exp <- exp(-maximum_eta)
    loser_exp <- exp(loser_eta - maximum_eta)
    winner_exp <- exp(winner_eta - maximum_eta)
    denominator <- reference_exp + loser_exp + winner_exp
    list(
      non_significant_resident = reference_exp / denominator,
      loser = loser_exp / denominator,
      winner = winner_exp / denominator
    )
  }
  colonizer_probabilities <- function(design, draws) {
    selected <- match(colnames(design), colnames(draws))
    if (anyNA(selected)) {
      stop("Colonizer prediction coefficients do not match the design matrix.", call. = FALSE)
    }
    stats::plogis(design %*% t(draws[, selected, drop = FALSE]))
  }
  chunks <- split(seq_len(nrow(grid)), ceiling(seq_len(nrow(grid)) / chunk_size))
  prediction_grid <- purrr::map_dfr(chunks, function(rows) {
    focal <- grid[rows, , drop = FALSE]
    resident_design <- stats::model.matrix(
      resident_rhs, data = focal, contrasts.arg = fit$resident$contrasts,
      xlev = fit$resident$xlevels
    )
    colonizer_design <- stats::model.matrix(
      colonizer_rhs, data = focal, contrasts.arg = fit$colonizer$contrasts,
      xlev = fit$colonizer$xlevels
    )
    resident_draw_probability <- resident_probabilities(
      resident_design, resident_draws
    )
    resident_central_probability <- resident_probabilities(
      resident_design, resident_central
    )
    colonizer_draw_probability <- colonizer_probabilities(
      colonizer_design, colonizer_draws
    )
    colonizer_central_probability <- colonizer_probabilities(
      colonizer_design, colonizer_central
    )
    probabilities <- c(
      resident_draw_probability,
      list(colonizer = colonizer_draw_probability)
    )
    central <- c(
      resident_central_probability,
      list(colonizer = colonizer_central_probability)
    )
    purrr::imap_dfr(probabilities, function(draw_matrix, estimand) {
      intervals <- matrixStats::rowQuantiles(
        draw_matrix, probs = c(0.025, 0.975), drop = FALSE
      )
      focal |>
        dplyr::mutate(
          estimand = .env$estimand,
          pool_condition = dplyr::if_else(
            .env$estimand == "colonizer", "nonresident", "resident"
          ),
          estimate = as.numeric(.env$central[[.env$estimand]][, 1]),
          conf_low = intervals[, 1], conf_high = intervals[, 2],
          interval_width = .data$conf_high - .data$conf_low
        )
    })
  })
  resident_sum_error <- prediction_grid |>
    dplyr::filter(.data$estimand %in% c(
      "non_significant_resident", "loser", "winner"
    )) |>
    dplyr::summarise(total = sum(.data$estimate), .by = "grid_id") |>
    dplyr::summarise(error = max(abs(.data$total - 1))) |>
    dplyr::pull(.data$error)
  list(
    grid = prediction_grid,
    warming_levels = warming_levels,
    experiment_support = experiments,
    simulation = tibble::tibble(
      n_draws = as.integer(n_draws), seed = as.integer(seed),
      grid_size = as.integer(grid_size), chunk_size = as.integer(chunk_size),
      covariance = "independence-based model Wald covariance",
      resident_probability_sum_maximum_error = resident_sum_error
    )
  )
}
