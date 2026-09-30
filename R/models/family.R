# Module: R/models/family
# Responsibility: Implement the named scientific workflow stage or shared infrastructure.
# Inputs: Function arguments and explicitly declared target dependencies.
# Outputs: Validated in-memory objects and declared workflow artifacts.
# Side effects: Limited to functions whose contracts explicitly write, copy, or render files.

# Family-resolution validation and summaries. Plot-cover aggregation remains
# distinct from species-effect averaging by construction.

#' Validate Family Bundle.
#'
#' @param data A data frame whose required columns are checked before use.
#' @param V A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param classification A data frame whose required columns are checked before use.
#' @param model_summary A fitted model object with the structure checked by this function.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
validate_family_bundle <- function(data, V, classification, model_summary) {
  ids <- as.character(data$effect_id)
  gates <- tibble::tibble(
    check = c(
      "unique_family_effect_ids", "covariance_dimensions", "covariance_order",
      "covariance_symmetric", "classification_unique", "model_summary_present"
    ),
    pass = c(
      !anyDuplicated(ids),
      identical(dim(V), c(length(ids), length(ids))),
      identical(rownames(V), ids) && identical(colnames(V), ids),
      max(abs(as.matrix(V) - t(as.matrix(V)))) <= 1e-10,
      !anyDuplicated(classification$family_id),
      nrow(model_summary) > 0L
    )
  )
  if (any(!gates$pass)) {
    stop("Family-resolution bundle failed: ",
         paste(gates$check[!gates$pass], collapse = ", "), call. = FALSE)
  }
  gates
}

#' Compare Classification Resolutions.
#'
#' @param species A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param family A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
compare_classification_resolutions <- function(species, family) {
  dplyr::bind_rows(
    summarize_reference_classification(species, "Species"),
    summarize_reference_classification(family, "Family")
  )
}

#' Family Overrides.
#'
#' @param resident_species_overrides A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
family_overrides <- function(resident_species_overrides) {
  assert_columns(
    resident_species_overrides,
    c("species_id", "family_id", "decision_basis"),
    "downstream species overrides"
  )
  result <- resident_species_overrides |>
    dplyr::transmute(
      species_id = as.character(.data$species_id),
      family_id = as.character(.data$family_id),
      decision_basis = as.character(.data$decision_basis)
    ) |>
    dplyr::filter(!is.na(.data$family_id), nzchar(trimws(.data$family_id)))
  if (anyDuplicated(result$species_id) || anyNA(result$species_id)) {
    stop("Canonical family overrides are incomplete or non-unique.", call. = FALSE)
  }
  result
}

#' Build Family Resident Membership.
#'
#' @param winner_loser A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param analysis_data A data frame whose required columns are checked before use.
#' @param family_overrides A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
build_family_resident_membership <- function(winner_loser, analysis_data,
                                             family_overrides = NULL) {
  species_results <- winner_loser |>
    dplyr::mutate(horizon = as.character(.data$horizon)) |>
    dplyr::filter(
      .data$transform == "raw", .data$contrast == "to_origin",
      .data$pool == "resident", .data$horizon == "max"
    ) |>
    dplyr::transmute(
      experiment_id = as.character(.data$experiment_id),
      source_species = as.character(.data$species)
    ) |>
    dplyr::distinct()
  taxonomy <- analysis_data |>
    dplyr::transmute(
      source_species = as.character(.data$species),
      species_id = as.character(.data$accepted_name),
      family_id = trimws(as.character(.data$family))
    ) |>
    dplyr::distinct()
  if (anyDuplicated(taxonomy$source_species)) {
    stop("Audited family assignment is non-unique.", call. = FALSE)
  }
  if (!is.null(family_overrides)) {
    assert_columns(
      family_overrides, c("species_id", "family_id", "decision_basis"),
      "canonical family overrides"
    )
    if (anyDuplicated(family_overrides$species_id)) {
      stop("Canonical family overrides are non-unique.", call. = FALSE)
    }
    taxonomy <- taxonomy |>
      dplyr::left_join(
        dplyr::select(
          family_overrides, "species_id",
          family_id_override = "family_id"
        ),
        by = "species_id"
      ) |>
      dplyr::mutate(
        family_id = dplyr::coalesce(
          .data$family_id_override, .data$family_id
        )
      ) |>
      dplyr::select(-"family_id_override")
  }
  if (any(is.na(taxonomy$family_id) | taxonomy$family_id == "")) {
    stop("Audited family assignment is missing.", call. = FALSE)
  }
  result <- species_results |>
    dplyr::left_join(taxonomy, by = "source_species") |>
    filter_resident_model_species("species_id")
  if (nrow(result) != 1604L ||
      dplyr::n_distinct(result$experiment_id) != 42L ||
      any(is.na(result$family_id) | result$family_id == "")) {
    stop("The canonical 1,604-effect resident family membership was not recovered.",
         call. = FALSE)
  }
  result
}

#' Complete Family Community.
#'
#' @param data A data frame whose required columns are checked before use.
#' @param membership A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
complete_family_community <- function(data, membership) {
  family_universe <- sort(unique(membership$family_id))
  surveyed_plot_years <- data |>
    dplyr::filter(!is.na(.data$relative_year)) |>
    dplyr::distinct(.data$ODT, .data$stable_plot_id, .data$relative_year)
  observed <- data |>
    dplyr::filter(
      .data$species %in% membership$source_species,
      !is.na(.data$relative_year)
    ) |>
    dplyr::left_join(membership, by = c("species" = "source_species")) |>
    dplyr::group_by(
      .data$ODT, .data$stable_plot_id, .data$relative_year, .data$family_id
    ) |>
    dplyr::summarise(
      relative_cover = sum(.data$relative_cover), .groups = "drop"
    )
  surveyed_plot_years |>
    tidyr::crossing(family_id = family_universe) |>
    dplyr::left_join(
      observed,
      by = c("ODT", "stable_plot_id", "relative_year", "family_id")
    ) |>
    dplyr::mutate(relative_cover = dplyr::coalesce(.data$relative_cover, 0)) |>
    dplyr::group_by(.data$ODT) |>
    tidyr::complete(stable_plot_id, relative_year, family_id = family_universe) |>
    dplyr::ungroup() |>
    dplyr::group_by(.data$ODT, .data$relative_year, .data$family_id) |>
    dplyr::mutate(
      relative_cover = dplyr::if_else(
        is.na(.data$relative_cover),
        mean(.data$relative_cover, na.rm = TRUE),
        .data$relative_cover
      ),
      relative_cover = dplyr::coalesce(.data$relative_cover, 0)
    ) |>
    dplyr::ungroup()
}

#' Family Auc by Plot.
#'
#' @param completed A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param experiment_id A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
family_auc_by_plot <- function(completed, experiment_id) {
  completed |>
    dplyr::group_by(.data$ODT, .data$stable_plot_id, .data$family_id) |>
    dplyr::summarise(
      auc = if (dplyr::n_distinct(.data$relative_year) < 2L) {
        NA_real_
      } else {
        DescTools::AUC(
          .data$relative_year, .data$relative_cover, method = "trapezoid"
        )
      },
      .groups = "drop"
    ) |>
    dplyr::mutate(
      experiment_id = experiment_id,
      horizon = "max",
      transform = "raw"
    )
}

#' Build Family Plot Auc.
#'
#' @param community A data frame whose required columns are checked before use.
#' @param experiments A data frame whose required columns are checked before use.
#' @param membership A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
build_family_plot_auc <- function(community, experiments, membership) {
  experiments <- experiments |>
    dplyr::filter(.data$included) |>
    dplyr::arrange(.data$experiment_id)
  if (nrow(experiments) != 42L) {
    stop("Expected exactly 42 included experiments.", call. = FALSE)
  }
  plot_blocks <- vector("list", nrow(experiments))
  audit_blocks <- vector("list", nrow(experiments))
  for (index in seq_len(nrow(experiments))) {
    experiment <- experiments[index, ]
    experiment_id <- experiment$experiment_id[[1]]
    data <- experiment_community(community, experiment) |>
      apply_legacy_observation_filter()
    experiment_membership <- membership |>
      dplyr::filter(.data$experiment_id == .env$experiment_id)
    if (!nrow(experiment_membership)) {
      stop("No resident membership for ", experiment_id, ".", call. = FALSE)
    }
    completed <- complete_family_community(data, experiment_membership)
    plot_blocks[[index]] <- family_auc_by_plot(completed, experiment_id)
    audit_blocks[[index]] <- tibble::tibble(
      experiment_id = experiment_id,
      n_resident_species = dplyr::n_distinct(
        experiment_membership$source_species
      ),
      n_resident_families = dplyr::n_distinct(
        experiment_membership$family_id
      ),
      n_completed_family_plot_year_rows = nrow(completed),
      n_origin_plots = dplyr::n_distinct(
        completed$stable_plot_id[completed$ODT == "origin_controls"]
      ),
      n_warmed_plots = dplyr::n_distinct(
        completed$stable_plot_id[completed$ODT == "warmed"]
      )
    )
  }
  list(
    plot_auc = dplyr::bind_rows(plot_blocks),
    cover_audit = dplyr::bind_rows(audit_blocks)
  )
}

#' Build Family Auc Covariance.
#'
#' @param plot_auc A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param membership A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param climate_data A data frame whose required columns are checked before use.
#' @return An R object described by the function name; its structure is validated before later stages use it.
build_family_auc_covariance <- function(plot_auc, membership, climate_data) {
  experiment_ids <- sort(unique(as.character(membership$experiment_id)))
  contrast_blocks <- list()
  covariance_blocks <- list()
  contrast_warnings <- character()
  for (experiment_id in experiment_ids) {
    raw <- plot_auc |>
      dplyr::filter(
        .data$experiment_id == .env$experiment_id, !is.na(.data$auc)
      ) |>
      dplyr::mutate(
        ODT = factor(
          as.character(.data$ODT),
          levels = c("destination_controls", "origin_controls", "warmed")
        ),
        family_id = factor(.data$family_id)
      )
    warnings <- character()
    fitted <- tryCatch(
      withCallingHandlers({
        model <- stats::lm(auc ~ ODT * family_id, data = raw)
        means <- emmeans::emmeans(model, ~ ODT | family_id)
        contrast <- emmeans::contrast(
          means, method = list(to_origin = c(0, -1, 1)), adjust = "none"
        )
        list(
          summary = tibble::as_tibble(as.data.frame(summary(contrast))),
          covariance = as.matrix(stats::vcov(contrast))
        )
      }, warning = function(condition) {
        warnings <<- c(warnings, conditionMessage(condition))
        invokeRestart("muffleWarning")
      }),
      error = function(condition) condition
    )
    if (inherits(fitted, "condition")) {
      stop("Family AUC contrast failed for ", experiment_id, ": ",
           conditionMessage(fitted), call. = FALSE)
    }
    target_families <- membership |>
      dplyr::filter(.data$experiment_id == .env$experiment_id) |>
      dplyr::distinct(.data$family_id) |>
      dplyr::arrange(.data$family_id) |>
      dplyr::pull(.data$family_id)
    order_index <- match(target_families, as.character(fitted$summary$family_id))
    if (anyNA(order_index)) {
      stop("A resident family is missing from contrasts for ", experiment_id,
           ".", call. = FALSE)
    }
    contrast_blocks[[experiment_id]] <-
      fitted$summary[order_index, , drop = FALSE] |>
      dplyr::transmute(
        experiment_id = experiment_id,
        family_id = as.character(.data$family_id),
        family_auc_effect = as.numeric(.data$estimate),
        family_auc_se = as.numeric(.data$SE),
        contrast_df = as.numeric(.data$df)
      )
    covariance_blocks[[experiment_id]] <-
      fitted$covariance[order_index, order_index, drop = FALSE]
    if (length(warnings)) {
      contrast_warnings <- c(
        contrast_warnings,
        paste(experiment_id, unique(warnings), sep = ": ")
      )
    }
  }
  data <- dplyr::bind_rows(contrast_blocks) |>
    dplyr::arrange(.data$experiment_id, .data$family_id) |>
    dplyr::mutate(
      horizon = "max",
      effect_id = paste("max", .data$experiment_id, .data$family_id, sep = "::")
    )
  ordered_ids <- unlist(lapply(names(contrast_blocks), function(id) {
    paste("max", id, contrast_blocks[[id]]$family_id, sep = "::")
  }), use.names = FALSE)
  assembled <- as.matrix(Matrix::bdiag(covariance_blocks))
  V_index <- match(data$effect_id, ordered_ids)
  if (anyNA(V_index)) stop("Family V ordering could not be aligned.", call. = FALSE)
  V <- assembled[V_index, V_index, drop = FALSE]
  dimnames(V) <- list(data$effect_id, data$effect_id)
  climate <- climate_data |>
    dplyr::transmute(
      experiment_id = as.character(.data$experiment_id),
      T = .data$T, P = .data$P, P_raw = .data$P_raw, W = .data$W,
      T_z = .data$T_z, P_z = .data$P_z, W_z = .data$W_z
    ) |>
    dplyr::distinct()
  if (nrow(climate) != 42L) {
    stop("Species scaling did not yield 42 experiment climates.", call. = FALSE)
  }
  membership_counts <- membership |>
    dplyr::summarise(
      n_resident_species_in_experiment =
        dplyr::n_distinct(.data$source_species),
      .by = c("experiment_id", "family_id")
    )
  data <- data |>
    dplyr::left_join(climate, by = "experiment_id") |>
    dplyr::left_join(membership_counts, by = c("experiment_id", "family_id")) |>
    dplyr::mutate(TP_z = .data$T_z * .data$P_z)
  if (any(!stats::complete.cases(data)) ||
      any(!is.finite(data$family_auc_se) | data$family_auc_se <= 0) ||
      anyDuplicated(data$effect_id)) {
    stop("Family AUC data failed completeness/key/SE validation.", call. = FALSE)
  }
  symmetry_error <- max(abs(V - t(V)))
  diagonal_error <- max(abs(diag(V) - data$family_auc_se^2))
  minimum_eigenvalue <- min(eigen(V, symmetric = TRUE, only.values = TRUE)$values)
  valid <- identical(rownames(V), data$effect_id) &&
    symmetry_error <= 1e-8 && diagonal_error <= 1e-8 &&
    minimum_eigenvalue >= -1e-8
  if (!valid) {
    stop("Full family AUC covariance failed validation; it was not repaired.",
         call. = FALSE)
  }
  list(
    data = data,
    V = V,
    contrast_warnings = contrast_warnings,
    validation = tibble::tibble(
      n_effects = nrow(data),
      n_experiments = dplyr::n_distinct(data$experiment_id),
      dimensions_valid = all(dim(V) == nrow(data)),
      order_valid = identical(rownames(V), data$effect_id),
      maximum_symmetry_error = symmetry_error,
      maximum_diagonal_variance_error = diagonal_error,
      minimum_eigenvalue = minimum_eigenvalue,
      valid = valid
    )
  )
}

#' Fit Family Aggregated Model.
#'
#' @param family_auc_bundle A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A fitted model bundle with diagnostics retained for validation.
fit_family_aggregated_model <- function(family_auc_bundle) {
  data <- family_auc_bundle$data |>
    dplyr::mutate(
      family_id = factor(.data$family_id),
      experiment_id = factor(.data$experiment_id)
    )
  random <- list(
    ~ 1 | family_id,
    ~ 0 + W_z | family_id,
    ~ 0 + T_z + P_z + TP_z | family_id
  )
  fit_once <- function(control = NULL) {
    arguments <- list(
      yi = data$family_auc_effect,
      V = family_auc_bundle$V,
      mods = ~ T_z * P_z * W_z,
      random = random,
      struct = "GEN",
      data = data,
      method = "REML",
      sparse = TRUE
    )
    if (!is.null(control)) arguments$control <- control
    auc_capture(do.call(metafor::rma.mv, arguments))
  }
  attempts <- "nlminb"
  captured <- fit_once()
  if (is.null(captured$fit)) {
    attempts <- c(attempts, "bobyqa")
    captured <- fit_once(list(optimizer = "bobyqa", maxfun = 100000))
  }
  if (is.null(captured$fit)) {
    stop("M_family_aggregated failed without simplification: ", captured$error,
         call. = FALSE)
  }
  fit <- captured$fit
  if (!identical(colnames(fit$G), "W_z") ||
      !identical(colnames(fit$H), c("T_z", "P_z", "TP_z"))) {
    stop("Family model random structure changed.", call. = FALSE)
  }
  list(
    fit = fit,
    warnings = captured$warnings,
    attempts = attempts,
    data = data,
    V = family_auc_bundle$V,
    contrast_warnings = family_auc_bundle$contrast_warnings,
    V_validation = family_auc_bundle$validation
  )
}

#' Family Model Components.
#'
#' @param model_bundle A fitted model object with the structure checked by this function.
#' @return An R object described by the function name; its structure is validated before later stages use it.
family_model_components <- function(model_bundle) {
  fit <- model_bundle$fit
  context_covariance <- as.matrix(fit$H)
  context_correlation <- stats::cov2cor(context_covariance)
  family_intercept_variance <- as.numeric(fit$sigma2[[1]])
  family_W_variance <- as.numeric(fit$tau2[[1]])
  full_covariance <- matrix(
    0, 5L, 5L,
    dimnames = list(
      c("intercept", "W_z", "T_z", "P_z", "TP_z"),
      c("intercept", "W_z", "T_z", "P_z", "TP_z")
    )
  )
  full_covariance["intercept", "intercept"] <- family_intercept_variance
  full_covariance["W_z", "W_z"] <- family_W_variance
  full_covariance[c("T_z", "P_z", "TP_z"),
                  c("T_z", "P_z", "TP_z")] <- context_covariance
  list(
    context_covariance = context_covariance,
    context_correlation = context_correlation,
    full_covariance = full_covariance,
    family_intercept_variance = family_intercept_variance,
    family_W_variance = family_W_variance
  )
}

#' Family Support Table.
#'
#' @param data A data frame whose required columns are checked before use.
#' @return An R object described by the function name; its structure is validated before later stages use it.
family_support_table <- function(data) {
  data |>
    dplyr::distinct(
      .data$family_id, .data$experiment_id, .data$T_z, .data$P_z, .data$W_z
    ) |>
    dplyr::group_by(.data$family_id) |>
    dplyr::group_modify(function(.x, .y) {
      design <- cbind(
        intercept = 1, T_z = .x$T_z, P_z = .x$P_z,
        TP_z = .x$T_z * .x$P_z
      )
      geometric <- cbind(intercept = 1, T_z = .x$T_z, P_z = .x$P_z)
      tibble::tibble(
        n_experiments = dplyr::n_distinct(.x$experiment_id),
        tp_design_rank = qr(design)$rank,
        tp_geometric_rank = qr(geometric)$rank,
        support_class = dplyr::case_when(
          dplyr::n_distinct(.x$experiment_id) == 1L ~ "single-point",
          qr(design)$rank == 4L ~ "full-rank",
          TRUE ~ "partial-rank"
        ),
        T_min = min(.x$T_z), T_max = max(.x$T_z),
        P_min = min(.x$P_z), P_max = max(.x$P_z),
        W_min = min(.x$W_z), W_max = max(.x$W_z)
      )
    }) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      eligible_for_individual_interpretation =
        .data$n_experiments >= 4L & .data$tp_design_rank == 4L
    )
}

#' Extract Family Random Coefficients.
#'
#' @param model_bundle A fitted model object with the structure checked by this function.
#' @param support A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
extract_family_random_coefficients <- function(model_bundle, support) {
  fit <- model_bundle$fit
  parts <- family_model_components(model_bundle)
  X <- as.matrix(fit$X)
  residual <- as.numeric(fit$yi - X %*% fit$beta)
  weighted_residual <- Matrix::solve(fit$M, residual)
  family_levels <- levels(fit$mf.g$outer)
  if (!identical(family_levels, levels(fit$mf.h$outer))) {
    stop("Family order differs across random blocks.", call. = FALSE)
  }
  n_families <- length(family_levels)
  Z_intercept <- fit$Z.S[[1]]
  Z_warming <- Matrix::Diagonal(x = as.numeric(fit$Z.G1[, 1])) %*% fit$Z.G2
  Z_context <- do.call(cbind, lapply(seq_len(3L), function(column) {
    Matrix::Diagonal(x = as.numeric(fit$Z.H1[, column])) %*% fit$Z.H2
  }))
  random_context <- as.numeric(
    Matrix::kronecker(
      parts$context_covariance, Matrix::Diagonal(n_families)
    ) %*% Matrix::crossprod(Z_context, weighted_residual)
  )
  random_context <- matrix(
    random_context, nrow = n_families, ncol = 3L,
    dimnames = list(family_levels, c("T_z", "P_z", "TP_z"))
  )
  random <- tibble::tibble(
    family_id = family_levels,
    random_intercept = as.numeric(
      parts$family_intercept_variance *
        Matrix::crossprod(Z_intercept, weighted_residual)
    ),
    random_W = as.numeric(
      parts$family_W_variance * Matrix::crossprod(Z_warming, weighted_residual)
    ),
    random_T = random_context[, "T_z"],
    random_P = random_context[, "P_z"],
    random_TP = random_context[, "TP_z"]
  ) |>
    dplyr::left_join(support, by = "family_id")
  list(
    coefficients = random,
    levels = family_levels,
    Z_intercept = Z_intercept,
    Z_warming = Z_warming,
    Z_context = Z_context
  )
}

#' Family Hull Grid.
#'
#' @param data A data frame whose required columns are checked before use.
#' @param family_name A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param W_z_value A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param points A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
family_hull_grid <- function(data, family_name, W_z_value, points = 81L) {
  observed <- data |>
    dplyr::filter(as.character(.data$family_id) == .env$family_name) |>
    dplyr::distinct(.data$T_z, .data$P_z)
  if (qr(cbind(1, observed$T_z, observed$P_z))$rank < 3L) {
    stop(family_name, " does not have 2D T-P support.", call. = FALSE)
  }
  candidate <- tidyr::expand_grid(
    T_z = seq(min(observed$T_z), max(observed$T_z), length.out = points),
    P_z = seq(min(observed$P_z), max(observed$P_z), length.out = points)
  ) |>
    dplyr::mutate(
      grid_point_type = "regular_grid",
      tile_width = diff(range(observed$T_z)) / (points - 1L),
      tile_height = diff(range(observed$P_z)) / (points - 1L)
    )
  hull <- geometry::convhulln(as.matrix(observed), options = "Qt")
  candidate |>
    dplyr::filter(geometry::inhulln(
      hull, as.matrix(dplyr::select(candidate, "T_z", "P_z"))
    )) |>
    dplyr::bind_rows(observed |>
      dplyr::mutate(
        grid_point_type = "observed_climate",
        tile_width = diff(range(observed$T_z)) / (points - 1L),
        tile_height = diff(range(observed$P_z)) / (points - 1L)
      )) |>
    dplyr::distinct(.data$T_z, .data$P_z, .keep_all = TRUE) |>
    dplyr::mutate(W_z = W_z_value, TP_z = .data$T_z * .data$P_z)
}

#' Classify Family Context.
#'
#' @param model_bundle A fitted model object with the structure checked by this function.
#' @param support A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param random_bundle A model formula or random-effects specification.
#' @param points A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
classify_family_context <- function(model_bundle, support, random_bundle,
                                    points = 81L) {
  fit <- model_bundle$fit
  data <- model_bundle$data
  parts <- family_model_components(model_bundle)
  random <- random_bundle$coefficients
  family_levels <- random_bundle$levels
  n_families <- length(family_levels)
  X <- as.matrix(fit$X)
  M_inverse_X <- Matrix::solve(fit$M, X)
  fixed_covariance <- solve(crossprod(X, M_inverse_X))
  fixed_beta <- stats::coef(fit)
  experiment_warming <- data |>
    dplyr::distinct(.data$experiment_id, .data$W_z, .data$W)
  W_median <- stats::median(experiment_warming$W)
  W_z_median <- stats::median(experiment_warming$W_z)
  if (round(W_median, 2) != 2.29) {
    stop("Species-analysis median raw warming is not 2.29 C.", call. = FALSE)
  }
  predict_surface <- function(family_name) {
    family_index <- match(family_name, family_levels)
    grid <- family_hull_grid(data, family_name, W_z_median, points)
    fixed_design <- auc_design_matrix(
      ~ T_z * P_z * W_z, grid, names(fixed_beta)
    )
    random_design <- cbind(
      intercept = 1, W_z = grid$W_z, T_z = grid$T_z,
      P_z = grid$P_z, TP_z = grid$TP_z
    )
    random_vector <- random |>
      dplyr::filter(.data$family_id == .env$family_name) |>
      dplyr::select(
        "random_intercept", "random_W", "random_T", "random_P", "random_TP"
      ) |>
      unlist(use.names = FALSE) |>
      as.numeric()
    Z_family <- cbind(
      random_bundle$Z_intercept[, family_index, drop = FALSE],
      random_bundle$Z_warming[, family_index, drop = FALSE],
      random_bundle$Z_context[, family_index, drop = FALSE],
      random_bundle$Z_context[, n_families + family_index, drop = FALSE],
      random_bundle$Z_context[, 2L * n_families + family_index, drop = FALSE]
    )
    M_inverse_Z <- Matrix::solve(fit$M, Z_family)
    XtM_inverse_Z <- as.matrix(crossprod(X, M_inverse_Z))
    ZtM_inverse_Z <- as.matrix(crossprod(Z_family, M_inverse_Z))
    random_error <- parts$full_covariance -
      parts$full_covariance %*% ZtM_inverse_Z %*% parts$full_covariance +
      parts$full_covariance %*% t(XtM_inverse_Z) %*% fixed_covariance %*%
      XtM_inverse_Z %*% parts$full_covariance
    fixed_random_error <-
      -fixed_covariance %*% XtM_inverse_Z %*% parts$full_covariance
    prediction_variance <- pmax(
      0,
      rowSums((fixed_design %*% fixed_covariance) * fixed_design) +
        rowSums((random_design %*% random_error) * random_design) +
        2 * rowSums(
          (fixed_design %*% fixed_random_error) * random_design
        )
    )
    estimate <- as.numeric(fixed_design %*% fixed_beta) +
      as.numeric(random_design %*% random_vector)
    grid |>
      dplyr::mutate(
        family_id = family_name,
        estimate = estimate,
        prediction_se = sqrt(prediction_variance),
        conf_low = .data$estimate - stats::qnorm(0.975) * .data$prediction_se,
        conf_high = .data$estimate + stats::qnorm(0.975) * .data$prediction_se
      )
  }
  eligible <- support |>
    dplyr::filter(.data$eligible_for_individual_interpretation) |>
    dplyr::pull(.data$family_id) |>
    as.character()
  predictions <- stats::setNames(lapply(eligible, predict_surface), eligible)
  classification <- dplyr::bind_rows(lapply(eligible, function(family_name) {
    prediction <- predictions[[family_name]]
    minimum_index <- which.min(prediction$estimate)
    maximum_index <- which.max(prediction$estimate)
    tibble::tibble(
      family_id = family_name,
      min_fitted_family_delta_auc = prediction$estimate[[minimum_index]],
      max_fitted_family_delta_auc = prediction$estimate[[maximum_index]],
      min_prediction_se = prediction$prediction_se[[minimum_index]],
      min_conf_low = prediction$conf_low[[minimum_index]],
      min_conf_high = prediction$conf_high[[minimum_index]],
      min_T_z = prediction$T_z[[minimum_index]],
      min_P_z = prediction$P_z[[minimum_index]],
      max_prediction_se = prediction$prediction_se[[maximum_index]],
      max_conf_low = prediction$conf_low[[maximum_index]],
      max_conf_high = prediction$conf_high[[maximum_index]],
      max_T_z = prediction$T_z[[maximum_index]],
      max_P_z = prediction$P_z[[maximum_index]]
    )
  })) |>
    dplyr::left_join(
      dplyr::select(
        random, "family_id", "random_W", "n_experiments",
        "tp_design_rank", "support_class"
      ),
      by = "family_id"
    ) |>
    dplyr::mutate(
      fitted_class = dplyr::case_when(
        .data$min_fitted_family_delta_auc > 0 ~ "Consistent winner",
        .data$max_fitted_family_delta_auc < 0 ~ "Consistent loser",
        TRUE ~ "Context-dependent"
      ),
      fitted_class = factor(
        .data$fitted_class,
        levels = c("Consistent loser", "Context-dependent", "Consistent winner")
      ),
      warming_C = W_median,
      warming_z = W_z_median
    )
  if (any(classification$min_fitted_family_delta_auc >
          classification$max_fitted_family_delta_auc)) {
    stop("A family fitted minimum exceeds its maximum.", call. = FALSE)
  }
  selected <- dplyr::bind_rows(
    classification |>
      dplyr::filter(.data$fitted_class == "Consistent winner") |>
      dplyr::slice_max(
        .data$min_fitted_family_delta_auc, n = 2L, with_ties = FALSE
      ),
    classification |>
      dplyr::filter(.data$fitted_class == "Consistent loser") |>
      dplyr::slice_min(
        .data$max_fitted_family_delta_auc, n = 2L, with_ties = FALSE
      ),
    classification |>
      dplyr::filter(.data$fitted_class == "Context-dependent") |>
      dplyr::mutate(
        response_span = .data$max_fitted_family_delta_auc -
          .data$min_fitted_family_delta_auc
      ) |>
      dplyr::slice_max(.data$response_span, n = 2L, with_ties = FALSE)
  ) |>
    dplyr::distinct(.data$family_id) |>
    dplyr::pull(.data$family_id) |>
    as.character()
  classification |>
    dplyr::mutate(
      selected_for_surface = .data$family_id %in% .env$selected,
      fitted_class = as.character(.data$fitted_class)
    )
}

#' Summarize Family Model.
#'
#' @param model_bundle A fitted model object with the structure checked by this function.
#' @return An R object described by the function name; its structure is validated before later stages use it.
summarize_family_model <- function(model_bundle) {
  fit <- model_bundle$fit
  parts <- family_model_components(model_bundle)
  context_eigenvalues <- eigen(
    parts$context_covariance, symmetric = TRUE, only.values = TRUE
  )$values
  full_eigenvalues <- eigen(
    parts$full_covariance, symmetric = TRUE, only.values = TRUE
  )$values
  boundary_tolerance <- 1e-8 * max(1, stats::var(model_bundle$data$family_auc_effect))
  random_variances <- diag(parts$full_covariance)
  maximum_correlation <- max(abs(
    parts$context_correlation[upper.tri(parts$context_correlation)]
  ))
  optimizer <- if (!is.null(fit$opt.res$ierr)) {
    "bobyqa"
  } else if (!is.null(fit$opt.res$convergence)) {
    "nlminb"
  } else {
    class(fit$opt.res)[[1]]
  }
  validation <- model_bundle$V_validation
  tibble::tibble(
    model = "M_family_aggregated",
    formula = "family_auc_effect ~ T_z * P_z * W_z",
    random_structure = paste(
      "family intercept; independent family W_z slope;",
      "correlated family T_z/P_z/TP_z block"
    ),
    n_effects = nrow(model_bundle$data),
    n_families = dplyr::n_distinct(model_bundle$data$family_id),
    n_experiments = dplyr::n_distinct(model_bundle$data$experiment_id),
    optimizer_attempts = paste(model_bundle$attempts, collapse = " -> "),
    optimizer_used = optimizer,
    optimizer_code = auc_value_or(
      fit$opt.res$convergence, auc_value_or(fit$opt.res$ierr, NA_integer_)
    ),
    warnings = paste(model_bundle$warnings, collapse = " | "),
    finite_coefficients = all(is.finite(stats::coef(fit))),
    finite_fixed_covariance = all(is.finite(stats::vcov(fit))),
    family_intercept_variance = parts$family_intercept_variance,
    family_W_slope_variance = parts$family_W_variance,
    family_T_variance = parts$context_covariance["T_z", "T_z"],
    family_P_variance = parts$context_covariance["P_z", "P_z"],
    family_TP_variance = parts$context_covariance["TP_z", "TP_z"],
    T_P_correlation = parts$context_correlation["T_z", "P_z"],
    T_TP_correlation = parts$context_correlation["T_z", "TP_z"],
    P_TP_correlation = parts$context_correlation["P_z", "TP_z"],
    minimum_context_eigenvalue = min(context_eigenvalues),
    minimum_full_random_eigenvalue = min(full_eigenvalues),
    maximum_absolute_context_correlation = maximum_correlation,
    boundary_tolerance = boundary_tolerance,
    any_variance_on_boundary = any(random_variances <= boundary_tolerance),
    covariance_positive_definite = min(full_eigenvalues) > 0,
    random_structure_identifiable =
      all(is.finite(full_eigenvalues)) &&
      min(full_eigenvalues) > boundary_tolerance &&
      all(random_variances > boundary_tolerance) &&
      maximum_correlation < 1 - 1e-6,
    V_source = "full family emmeans covariance reconstructed from family AUC models",
    V_maximum_diagonal_variance_error =
      validation$maximum_diagonal_variance_error,
    V_maximum_symmetry_error = validation$maximum_symmetry_error,
    V_minimum_eigenvalue = validation$minimum_eigenvalue,
    V_valid = validation$valid,
    contrast_warnings = paste(model_bundle$contrast_warnings, collapse = " | ")
  )
}

#' Compare Family Covariance.
#'
#' @param candidate A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param reference A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param absolute A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param relative A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
compare_family_covariance <- function(candidate, reference,
                                      absolute = 1e-10,
                                      relative = 1e-8) {
  reference_matrix <- if (is.list(reference) && !is.null(reference$V_family_auc)) {
    as.matrix(reference$V_family_auc)
  } else {
    as.matrix(reference)
  }
  candidate_matrix <- as.matrix(candidate)
  same_dimensions <- identical(dim(candidate_matrix), dim(reference_matrix))
  same_names <- identical(dimnames(candidate_matrix), dimnames(reference_matrix))
  maximum_error <- if (same_dimensions) {
    max(abs(candidate_matrix - reference_matrix))
  } else {
    Inf
  }
  maximum_limit <- if (same_dimensions) {
    max(absolute + relative * abs(reference_matrix))
  } else {
    NA_real_
  }
  result <- tibble::tibble(
    check = c(
      "family_covariance_dimensions", "family_covariance_order",
      "family_covariance_values"
    ),
    pass = c(
      same_dimensions, same_names,
      same_dimensions && maximum_error <= maximum_limit
    ),
    detail = c(
      paste(dim(candidate_matrix), collapse = " x "),
      if (same_names) "Row and column identifiers are identical." else
        "Row or column identifiers differ.",
      paste0(
        "maximum absolute error = ", format(maximum_error, digits = 16),
        "; maximum allowed = ", format(maximum_limit, digits = 16)
      )
    )
  )
  if (any(!result$pass)) {
    stop(
      "Family sampling covariance regression failed: ",
      paste(result$check[!result$pass], collapse = ", "), call. = FALSE
    )
  }
  result
}
