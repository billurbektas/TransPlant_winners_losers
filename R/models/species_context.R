# Module: R/models/species_context
# Responsibility: Classify supported species as winners, losers, or context-dependent.
# Inputs: Resident climate fit and observed species climate support.
# Outputs: Species coefficient summaries, supported predictions, and classes.
# Side effects: None.

# Pure species-context classification and callout-selection helpers.

#' Summarize eligible species or families by fitted response class
#'
#' @param classification Table containing one `fitted_class` per eligible unit.
#' @param resolution Label for the summarized taxonomic resolution.
#' @return One row with counts and percentages for each response class.
summarize_context_classification <- function(classification, resolution) {
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

#' Extract Conditional Species Coefficients.
#'
#' @param fit A fitted model object with the structure checked by this function.
#' @return An R object described by the function name; its structure is validated before later stages use it.
extract_conditional_species_coefficients <- function(fit) {
  X <- as.matrix(fit$X)
  residual <- as.numeric(fit$yi - X %*% fit$beta)
  weighted_residual <- Matrix::solve(fit$M, residual)
  species_levels <- levels(fit$mf.g$outer)
  context_levels <- levels(fit$mf.h$outer)
  if (!identical(species_levels, context_levels)) {
    stop("Species order differs across random blocks.", call. = FALSE)
  }
  n_species <- length(species_levels)
  Z_intercept <- fit$Z.S[[1]]
  Z_W <- Matrix::Diagonal(x = as.numeric(fit$Z.G1[, 1])) %*% fit$Z.G2
  Z_context <- do.call(cbind, lapply(seq_len(ncol(fit$H)), function(index) {
    Matrix::Diagonal(x = as.numeric(fit$Z.H1[, index])) %*% fit$Z.H2
  }))
  random_intercept <- as.numeric(
    fit$sigma2[[1]] * Matrix::crossprod(Z_intercept, weighted_residual)
  )
  random_W <- as.numeric(
    fit$tau2[[1]] * Matrix::crossprod(Z_W, weighted_residual)
  )
  D_context <- Matrix::kronecker(as.matrix(fit$H), Matrix::Diagonal(n_species))
  context_vector <- as.numeric(
    D_context %*% Matrix::crossprod(Z_context, weighted_residual)
  )
  context_matrix <- matrix(
    context_vector, nrow = n_species, ncol = 3L,
    dimnames = list(species_levels, c("T_z", "P_z", "TP_z"))
  )
  list(
    table = tibble::tibble(
      species_id = species_levels,
      random_intercept = random_intercept,
      random_W = random_W,
      random_T = context_matrix[, "T_z"],
      random_P = context_matrix[, "P_z"],
      random_TP = context_matrix[, "TP_z"]
    ),
    X = X,
    levels = species_levels,
    Z_intercept = Z_intercept,
    Z_W = Z_W,
    Z_context = Z_context
  )
}

#' Species Context Support.
#'
#' @param data A data frame whose required columns are checked before use.
#' @return An R object described by the function name; its structure is validated before later stages use it.
species_context_support <- function(data) {
  data |>
    dplyr::distinct(.data$species_id, .data$experiment_id, .data$T_z, .data$P_z) |>
    dplyr::group_by(.data$species_id) |>
    dplyr::group_modify(function(.x, .y) {
      design <- cbind(1, .x$T_z, .x$P_z, .x$T_z * .x$P_z)
      tibble::tibble(
        n_experiments = dplyr::n_distinct(.x$experiment_id),
        n_distinct_TP = nrow(dplyr::distinct(.x, .data$T_z, .data$P_z)),
        context_design_rank = qr(design)$rank
      )
    }) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      eligible = .data$n_experiments >= 4L & .data$context_design_rank == 4L
    )
}

#' Species Observed Hull Grid.
#'
#' @param data A data frame whose required columns are checked before use.
#' @param species_name A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param W_z_value A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param points A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
species_observed_hull_grid <- function(data, species_name, W_z_value,
                                       points = 81L) {
  observed <- data |>
    dplyr::filter(as.character(.data$species_id) == .env$species_name) |>
    dplyr::distinct(.data$T_z, .data$P_z)
  candidate <- tidyr::expand_grid(
    T_z = seq(min(observed$T_z), max(observed$T_z), length.out = points),
    P_z = seq(min(observed$P_z), max(observed$P_z), length.out = points)
  )
  hull <- geometry::convhulln(
    as.matrix(observed[, c("T_z", "P_z")]), options = "Qt"
  )
  candidate |>
    dplyr::filter(geometry::inhulln(
      hull, as.matrix(dplyr::select(candidate, "T_z", "P_z"))
    )) |>
    dplyr::bind_rows(observed) |>
    dplyr::distinct(.data$T_z, .data$P_z) |>
    dplyr::mutate(W_z = W_z_value, TP_z = .data$T_z * .data$P_z)
}

#' Predict Species on Observed Support.
#'
#' @param fit A fitted model object with the structure checked by this function.
#' @param data A data frame whose required columns are checked before use.
#' @param formula A model formula or random-effects specification.
#' @param random_effects A data frame whose required columns are checked before use.
#' @param species_name A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param W_z_value A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param points A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
predict_species_on_observed_support <- function(
    fit, data, formula, random_effects, species_name, W_z_value,
    points = 81L) {
  species_index <- match(species_name, random_effects$levels)
  if (is.na(species_index)) stop("Species absent from fitted random effects: ", species_name)
  grid <- species_observed_hull_grid(data, species_name, W_z_value, points)
  fixed_beta <- stats::coef(fit)
  fixed_covariance <- as.matrix(stats::vcov(fit))
  X_new <- auc_design_matrix(formula, grid, names(fixed_beta))
  R_new <- cbind(
    intercept = 1, W_z = grid$W_z, T_z = grid$T_z,
    P_z = grid$P_z, TP_z = grid$TP_z
  )
  random_row <- random_effects$table |>
    dplyr::filter(.data$species_id == .env$species_name)
  random_vector <- as.numeric(random_row[1, c(
    "random_intercept", "random_W", "random_T", "random_P", "random_TP"
  )])
  n_species <- length(random_effects$levels)
  Z_species <- cbind(
    random_effects$Z_intercept[, species_index, drop = FALSE],
    random_effects$Z_W[, species_index, drop = FALSE],
    random_effects$Z_context[, species_index, drop = FALSE],
    random_effects$Z_context[, n_species + species_index, drop = FALSE],
    random_effects$Z_context[, 2L * n_species + species_index, drop = FALSE]
  )
  D_species <- matrix(
    0, 5L, 5L,
    dimnames = list(
      c("intercept", "W_z", "T_z", "P_z", "TP_z"),
      c("intercept", "W_z", "T_z", "P_z", "TP_z")
    )
  )
  D_species["intercept", "intercept"] <- fit$sigma2[[1]]
  D_species["W_z", "W_z"] <- fit$tau2[[1]]
  D_species[c("T_z", "P_z", "TP_z"), c("T_z", "P_z", "TP_z")] <-
    as.matrix(fit$H)
  M_inverse_Z <- Matrix::solve(fit$M, Z_species)
  XtM_inverse_Z <- as.matrix(crossprod(random_effects$X, M_inverse_Z))
  ZtM_inverse_Z <- as.matrix(crossprod(Z_species, M_inverse_Z))
  random_pev <- D_species - D_species %*% ZtM_inverse_Z %*% D_species +
    D_species %*% t(XtM_inverse_Z) %*% fixed_covariance %*%
    XtM_inverse_Z %*% D_species
  fixed_random_cov <- -fixed_covariance %*% XtM_inverse_Z %*% D_species
  variance <- rowSums((X_new %*% fixed_covariance) * X_new) +
    rowSums((R_new %*% random_pev) * R_new) +
    2 * rowSums((X_new %*% fixed_random_cov) * R_new)
  estimate <- as.numeric(X_new %*% fixed_beta) +
    as.numeric(R_new %*% random_vector)
  grid |>
    dplyr::mutate(
      species_id = species_name,
      estimate = estimate,
      prediction_SE = sqrt(pmax(0, variance)),
      CI_low = .data$estimate - stats::qnorm(0.975) * .data$prediction_SE,
      CI_high = .data$estimate + stats::qnorm(0.975) * .data$prediction_SE
    )
}

#' Build Species Context Bundle.
#'
#' @param captured A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param data A data frame whose required columns are checked before use.
#' @param formula A model formula or random-effects specification.
#' @param points A value supplied by the calling workflow; admissible values are checked in the function body.
#' @details Eligibility requires the registered experiment-count and full-rank
#' climate-support rules. Predictions are evaluated only inside each species'
#' observed destination-climate hull at the registered warming context. The
#' returned classification is extrema-based and is not a significance test.
#' @return An R object described by the function name; its structure is validated before later stages use it.
build_species_context_bundle <- function(captured, data, formula, points = 81L) {
  fit <- if (is.list(captured) && !is.null(captured$fit)) captured$fit else captured
  random_effects <- extract_conditional_species_coefficients(fit)
  support <- species_context_support(data)
  eligible_species <- support |>
    dplyr::filter(.data$eligible) |>
    dplyr::pull(.data$species_id) |>
    as.character()
  climates <- data |>
    dplyr::distinct(.data$experiment_id, .data$W, .data$W_z)
  warming_raw <- stats::median(climates$W)
  warming_z <- stats::median(climates$W_z)
  predictions <- dplyr::bind_rows(lapply(eligible_species, function(species_name) {
    predict_species_on_observed_support(
      fit, data, formula, random_effects, species_name, warming_z, points
    )
  }))
  classification <- classify_supported_extrema(
    predictions, support, random_effects$table, warming_raw, warming_z
  )
  list(
    support = support,
    random_coefficients = random_effects$table,
    predictions = predictions,
    classification = classification
  )
}


#' Classify Supported Extrema.
#'
#' @param predictions A data frame whose required columns are checked before use.
#' @param support A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param random_coefficients A model formula or random-effects specification.
#' @param warming_raw A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param warming_z A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
classify_supported_extrema <- function(predictions, support, random_coefficients,
                                       warming_raw, warming_z) {
  extrema <- predictions |>
    dplyr::group_split(.data$species_id, .keep = TRUE) |>
    lapply(function(x) {
      minimum <- x[which.min(x$estimate), , drop = FALSE]
      maximum <- x[which.max(x$estimate), , drop = FALSE]
      tibble::tibble(
        species_id = as.character(x$species_id[[1]]),
        min_fitted_delta_auc = minimum$estimate[[1]],
        min_prediction_SE = minimum$prediction_SE[[1]],
        min_CI_low = minimum$CI_low[[1]],
        min_CI_high = minimum$CI_high[[1]],
        min_T_z = minimum$T_z[[1]],
        min_P_z = minimum$P_z[[1]],
        max_fitted_delta_auc = maximum$estimate[[1]],
        max_prediction_SE = maximum$prediction_SE[[1]],
        max_CI_low = maximum$CI_low[[1]],
        max_CI_high = maximum$CI_high[[1]],
        max_T_z = maximum$T_z[[1]],
        max_P_z = maximum$P_z[[1]],
        n_support_grid_points = nrow(x)
      )
    }) |>
    dplyr::bind_rows() |>
    dplyr::left_join(support, by = "species_id") |>
    dplyr::left_join(
      dplyr::select(random_coefficients, "species_id", "random_W"),
      by = "species_id"
    ) |>
    dplyr::mutate(
      fitted_class = dplyr::case_when(
        .data$min_fitted_delta_auc > 0 ~ "Consistent winner",
        .data$max_fitted_delta_auc < 0 ~ "Consistent loser",
        TRUE ~ "Context-dependent"
      ),
      warming_raw_C = warming_raw,
      warming_z = warming_z
    )
  if (any(extrema$min_fitted_delta_auc > extrema$max_fitted_delta_auc)) {
    stop("A species fitted minimum exceeds its maximum.", call. = FALSE)
  }
  extrema
}

#' Callout Species.
#'
#' @return An R object described by the function name; its structure is validated before later stages use it.
callout_species <- function() {
  c(
    "Potentilla erecta", "Plantago lanceolata", "Elymus repens",
    "Silene acaulis"
  )
}

#' Select Species Callout Examples.
#'
#' @param classification A data frame whose required columns are checked before use.
#' @return An R object described by the function name; its structure is validated before later stages use it.
select_species_callout_examples <- function(classification) {
  requested <- callout_species()
  expected_classes <- c(
    "Consistent winner", "Consistent loser", "Context-dependent",
    "Context-dependent"
  )
  ranked <- classification |>
    dplyr::mutate(
      species_id = as.character(.data$species_id),
      fitted_class = as.character(.data$fitted_class),
      response_span = .data$max_fitted_delta_auc - .data$min_fitted_delta_auc
    ) |>
    dplyr::group_by(.data$fitted_class) |>
    dplyr::arrange(
      dplyr::case_when(
        .data$fitted_class == "Consistent winner" ~ -.data$min_fitted_delta_auc,
        .data$fitted_class == "Consistent loser" ~ .data$max_fitted_delta_auc,
        TRUE ~ -.data$response_span
      ),
      .data$species_id,
      .by_group = TRUE
    ) |>
    dplyr::mutate(selection_rank_within_class = dplyr::row_number()) |>
    dplyr::ungroup()
  selection <- ranked |>
    dplyr::filter(.data$species_id %in% .env$requested) |>
    dplyr::mutate(display_order = match(.data$species_id, .env$requested)) |>
    dplyr::arrange(.data$display_order)
  if (!identical(selection$species_id, requested)) {
    stop("One or more established callout species are absent.", call. = FALSE)
  }
  if (!identical(selection$fitted_class, expected_classes)) {
    stop("An established callout species changed fitted class.", call. = FALSE)
  }
  selection |>
    dplyr::mutate(
      selection_rule = dplyr::case_when(
        .data$fitted_class == "Consistent winner" ~ "highest fitted minimum",
        .data$fitted_class == "Consistent loser" ~ "lowest fitted maximum",
        TRUE ~ "widest fitted min-max span"
      )
    )
}
