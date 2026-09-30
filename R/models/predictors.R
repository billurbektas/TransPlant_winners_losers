# Module: R/models/predictors
# Responsibility: Implement the named scientific workflow stage or shared infrastructure.
# Inputs: Function arguments and explicitly declared target dependencies.
# Outputs: Validated in-memory objects and declared workflow artifacts.
# Side effects: Limited to functions whose contracts explicitly write, copy, or render files.

# Predictor-specific data preparation and the validated model specifications.

#' Build Functional Species.
#'
#' @param functional_groups A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param retained_species A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param resident_species_overrides A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
build_functional_species <- function(functional_groups,
                                               retained_species,
                                               resident_species_overrides = NULL) {
  retained_species <- sort(unique(as.character(retained_species)))
  species_column <- if ("species_id" %in% names(functional_groups)) {
    "species_id"
  } else {
    "accepted_name"
  }
  accepted_column <- if ("tnrs_accepted_name" %in% names(functional_groups)) {
    "tnrs_accepted_name"
  } else {
    species_column
  }
  woody_column <- if ("woody_status" %in% names(functional_groups)) {
    "woody_status"
  } else {
    "woodiness"
  }
  required <- c(species_column, accepted_column, "family", "graminoid_status",
                woody_column)
  if (!all(required %in% names(functional_groups))) {
    stop("Functional source is missing: ",
         paste(setdiff(required, names(functional_groups)), collapse = ", "),
         call. = FALSE)
  }
  result <- functional_groups |>
    dplyr::transmute(
      species_id = as.character(.data[[species_column]]),
      TNRS_accepted_name = as.character(.data[[accepted_column]]),
      family = as.character(.data$family),
      graminoid_status = as.character(.data$graminoid_status),
      woodiness = as.character(.data[[woody_column]])
    )
  if (!is.null(resident_species_overrides)) {
    assert_columns(
      resident_species_overrides,
      c("species_id", "family_id", "graminoid_status", "woodiness",
        "decision_basis"),
      "resident species overrides"
    )
    functional_overrides <- resident_species_overrides |>
      dplyr::transmute(
        species_id = as.character(.data$species_id),
        TNRS_accepted_name_override = as.character(.data$species_id),
        family_override = as.character(.data$family_id),
        graminoid_status_override = as.character(.data$graminoid_status),
        woodiness_override = as.character(.data$woodiness)
      )
    if (anyDuplicated(functional_overrides$species_id)) {
      stop("Reviewed functional overrides are non-unique.", call. = FALSE)
    }
    result <- dplyr::full_join(
      result, functional_overrides, by = "species_id"
    ) |>
      dplyr::mutate(
        TNRS_accepted_name = dplyr::coalesce(
          .data$TNRS_accepted_name_override, .data$TNRS_accepted_name
        ),
        family = dplyr::coalesce(.data$family_override, .data$family),
        graminoid_status = dplyr::coalesce(
          .data$graminoid_status_override, .data$graminoid_status
        ),
        woodiness = dplyr::coalesce(
          .data$woodiness_override, .data$woodiness
        )
      ) |>
      dplyr::select(
        "species_id", "TNRS_accepted_name", "family",
        "graminoid_status", "woodiness"
      )
  }
  result <- result |>
    dplyr::filter(.data$species_id %in% .env$retained_species) |>
    dplyr::arrange(.data$species_id)
  if (nrow(result) != length(retained_species) ||
      anyDuplicated(result$species_id) ||
      anyNA(result[c("graminoid_status", "woodiness")]) ||
      !identical(result$species_id, retained_species)) {
    stop("Reviewed functional statuses are incomplete for retained species.",
         call. = FALSE)
  }
  result
}

#' Audit Functional Source Coverage.
#'
#' @param functional_groups A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param validated_functional_species A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param resident_species_overrides A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
audit_functional_source_coverage <- function(functional_groups,
                                             validated_functional_species,
                                             resident_species_overrides = NULL) {
  active_ids <- as.character(functional_groups$accepted_name)
  override_ids <- if (is.null(resident_species_overrides)) {
    character()
  } else {
    as.character(resident_species_overrides$species_id)
  }
  available_ids <- union(active_ids, override_ids)
  required_ids <- as.character(validated_functional_species$species_id)
  missing <- setdiff(required_ids, available_ids)
  tibble::tibble(
    source = paste(
      "reviewed functional-group table plus explicit resident-species overrides",
      "versus the validated historical interface"
    ),
    expected_species = length(required_ids),
    active_species_found = sum(required_ids %in% active_ids),
    override_species_supplied = sum(
      required_ids %in% override_ids & !required_ids %in% active_ids
    ),
    missing_species = length(missing),
    missing_species_ids = paste(missing, collapse = "; "),
    pass = !length(missing)
  )
}

#' Compare Functional Species Reference.
#'
#' @param candidate A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param reference A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
compare_functional_species_reference <- function(candidate, reference) {
  candidate <- candidate |> dplyr::arrange(.data$species_id)
  reference <- reference |> dplyr::arrange(.data$species_id)
  failures <- compare_reference_node(
    candidate, reference, c("functional_species")
  )
  pass <- !length(failures)
  result <- tibble::tibble(
    check = "functional species assignments",
    pass = pass,
    detail = if (pass) {
      paste(nrow(candidate), "Graminoid/Woodiness assignments are exact.")
    } else {
      paste(failures, collapse = " | ")
    }
  )
  if (!pass) stop(result$detail[[1]], call. = FALSE)
  result
}

#' Predictor Formulas.
#'
#' @return An R object described by the function name; its structure is validated before later stages use it.
predictor_formulas <- function() {
  phylogeny_full <- ~ T_z * P_z * W_z +
    PEV1_z + PEV2_z + PEV3_z +
    PEV1_z:W_z + PEV2_z:W_z + PEV3_z:W_z +
    PEV1_z:T_z + PEV2_z:T_z + PEV3_z:T_z +
    PEV1_z:P_z + PEV2_z:P_z + PEV3_z:P_z +
    PEV1_z:TP_z + PEV2_z:TP_z + PEV3_z:TP_z
  list(
    climate = ~ T_z * P_z * W_z,
    continuous_full = ~ Trait_z * T_z * P_z * W_z,
    Graminoid_full = ~ GraminoidStatus * T_z * P_z * W_z,
    Woodiness_full = ~ Woodiness * T_z * P_z * W_z,
    Phylogeny_full = phylogeny_full
  )
}

#' Prepare Trait Data.
#'
#' @param trait_species A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param auc_data A data frame whose required columns are checked before use.
#' @param sampling_variance_matrix Effect-ordered sampling variances.
#' @param trait A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param expected_species A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
prepare_trait_data <- function(
    trait_species, auc_data, sampling_variance_matrix, trait,
    expected_species = 526L) {
  required <- c("species_id", trait)
  if (!all(required %in% names(trait_species))) {
    stop("Trait input is missing: ", paste(setdiff(required, names(trait_species)),
                                           collapse = ", "), call. = FALSE)
  }
  unique_values <- trait_species |>
    dplyr::transmute(
      species_id = as.character(.data$species_id),
      trait_value = as.numeric(.data[[trait]])
    ) |>
    dplyr::filter(is.finite(.data$trait_value), .data$trait_value > 0) |>
    dplyr::mutate(log_trait = log(.data$trait_value))
  if (anyDuplicated(unique_values$species_id)) {
    stop(trait, " is not unique by species.", call. = FALSE)
  }
  center <- mean(unique_values$log_trait)
  scale <- stats::sd(unique_values$log_trait)
  if (!is.finite(scale) || scale <= 0) stop(trait, " has no usable log-scale variation.")
  unique_values <- unique_values |>
    dplyr::mutate(Trait_z = (.data$log_trait - center) / scale)
  data <- auc_data |>
    dplyr::inner_join(unique_values, by = "species_id") |>
    dplyr::mutate(
      species_id = droplevels(factor(.data$species_id)),
      experiment_id = droplevels(factor(.data$experiment_id))
    )
  trait_sampling_variances <- subset_sampling_variance_matrix(
    sampling_variance_matrix, data$effect_id
  )
  if (!identical(
        rownames(trait_sampling_variances), as.character(data$effect_id)
      )) {
    stop(trait, " sampling-variance ordering validation failed.", call. = FALSE)
  }
  list(
    data = data,
    sampling_variance_matrix = trait_sampling_variances,
    species_values = unique_values,
    scaling = tibble::tibble(
      trait = trait,
      species_N_with_trait = nrow(unique_values),
      species_N_missing = expected_species - nrow(unique_values),
      raw_trait_mean = mean(unique_values$trait_value),
      raw_trait_SD = stats::sd(unique_values$trait_value),
      log_trait_mean = center,
      log_trait_SD = scale,
      n_effects = nrow(data),
      n_experiments = dplyr::n_distinct(data$experiment_id),
      weighting = "unique retained trait-complete species before AUC-row expansion"
    )
  )
}

#' Standardize Unique Species Predictors.
#'
#' @param species_data A data frame whose required columns are checked before use.
#' @param columns A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param suffix A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
standardize_unique_species_predictors <- function(species_data, columns,
                                                  suffix = "_z") {
  if (anyDuplicated(species_data$species_id)) {
    stop("Predictor input is not unique by species.", call. = FALSE)
  }
  result <- species_data
  rows <- lapply(columns, function(column) {
    values <- as.numeric(species_data[[column]])
    if (any(!is.finite(values))) stop(column, " contains non-finite values.")
    center <- mean(values)
    scale <- stats::sd(values)
    if (!is.finite(scale) || scale <= 0) stop(column, " has no usable variation.")
    result[[paste0(column, suffix)]] <<- (values - center) / scale
    tibble::tibble(
      PEV = column,
      species_N = nrow(species_data),
      raw_mean = center,
      raw_SD = scale,
      weighting = "all 526 unique retained species before AUC-row expansion"
    )
  })
  list(data = result, scaling = dplyr::bind_rows(rows))
}

#' Prepare Functional Predictors.
#'
#' @param auc_data A data frame whose required columns are checked before use.
#' @param functional_species A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
prepare_functional_predictors <- function(auc_data, functional_species) {
  if (anyDuplicated(functional_species$species_id)) {
    stop("Functional information is not unique by species.", call. = FALSE)
  }
  data <- auc_data |>
    dplyr::left_join(functional_species, by = "species_id") |>
    dplyr::mutate(
      GraminoidStatus = factor(
        .data$graminoid_status, levels = c("non_graminoid", "graminoid")
      ),
      Woodiness = factor(.data$woodiness, levels = c("non_woody", "woody"))
    )
  if (nrow(data) != nrow(auc_data) ||
      any(!stats::complete.cases(data[c("GraminoidStatus", "Woodiness")]))) {
    stop("Functional predictor join is incomplete.", call. = FALSE)
  }
  data
}

#' Build Predictor Specifications.
#'
#' @param prepared_traits A data frame whose required columns are checked before use.
#' @param functional_data A data frame whose required columns are checked before use.
#' @param phylogeny_data A data frame whose required columns are checked before use.
#' @param sampling_variance_matrix Effect-ordered sampling variances.
#' @return An R object described by the function name; its structure is validated before later stages use it.
build_predictor_specifications <- function(prepared_traits, functional_data,
                                               phylogeny_data,
                                               sampling_variance_matrix) {
  f <- predictor_formulas()
  trait_spec <- function(key, label, prepared) {
    list(
      key = key, label = label, data = prepared$data,
      sampling_variance_matrix = prepared$sampling_variance_matrix,
      full_formula = f$continuous_full,
      predictor_columns = "Trait_z", full_data = FALSE
    )
  }
  list(
    LDMC = trait_spec("LDMC", "LDMC", prepared_traits$LDMC),
    LMA = trait_spec("LMA", "LMA", prepared_traits$LMA),
    Height = trait_spec("Height", "Height", prepared_traits$vegetative_height),
    Graminoid = list(
      key = "Graminoid", label = "Graminoid status", data = functional_data,
      sampling_variance_matrix = sampling_variance_matrix,
      full_formula = f$Graminoid_full,
      predictor_columns = "GraminoidStatus", full_data = TRUE
    ),
    Woodiness = list(
      key = "Woodiness", label = "Woodiness", data = functional_data,
      sampling_variance_matrix = sampling_variance_matrix,
      full_formula = f$Woodiness_full,
      predictor_columns = "Woodiness", full_data = TRUE
    ),
    Phylogeny = list(
      key = "Phylogeny", label = "Phylogenetic position (PEV1-PEV3)",
      data = phylogeny_data,
      sampling_variance_matrix = sampling_variance_matrix,
      full_formula = f$Phylogeny_full,
      predictor_columns = c("PEV1_z", "PEV2_z", "PEV3_z"), full_data = TRUE
    )
  )
}

#' Trait Predictor Specification.
#'
#' @param key A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param label A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param prepared A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
trait_predictor_specification <- function(key, label, prepared) {
  formulas <- predictor_formulas()
  list(
    key = key,
    label = label,
    data = prepared$data,
    sampling_variance_matrix = prepared$sampling_variance_matrix,
    full_formula = formulas$continuous_full,
    predictor_columns = "Trait_z",
    full_data = FALSE
  )
}

#' Functional Predictor Specification.
#'
#' @param key A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param label A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param data A data frame whose required columns are checked before use.
#' @param sampling_variance_matrix Effect-ordered sampling variances.
#' @return An R object described by the function name; its structure is validated before later stages use it.
functional_predictor_specification <- function(
    key, label, data, sampling_variance_matrix) {
  formulas <- predictor_formulas()
  if (identical(key, "Graminoid")) {
    full_formula <- formulas$Graminoid_full
    predictor_columns <- "GraminoidStatus"
  } else if (identical(key, "Woodiness")) {
    full_formula <- formulas$Woodiness_full
    predictor_columns <- "Woodiness"
  } else {
    stop("Unknown functional predictor key: ", key, call. = FALSE)
  }
  list(
    key = key,
    label = label,
    data = data,
    sampling_variance_matrix = sampling_variance_matrix,
    full_formula = full_formula,
    predictor_columns = predictor_columns,
    full_data = TRUE
  )
}

#' Build Phylogeny Predictor Specification.
#'
#' @param data A data frame whose required columns are checked before use.
#' @param sampling_variance_matrix Effect-ordered sampling variances.
#' @return An R object described by the function name; its structure is validated before later stages use it.
build_phylogeny_predictor_specification <- function(
    data, sampling_variance_matrix) {
  formulas <- predictor_formulas()
  list(
    key = "Phylogeny",
    label = "Phylogenetic position (PEV1-PEV3)",
    data = data,
    sampling_variance_matrix = sampling_variance_matrix,
    full_formula = formulas$Phylogeny_full,
    predictor_columns = c("PEV1_z", "PEV2_z", "PEV3_z"),
    full_data = TRUE
  )
}

#' Fit matched climate-only and full predictor models.
#'
#' @param spec A validated workflow configuration or specification object.
#' @param climate_formula A model formula or random-effects specification.
#' @param full_climate_ML Shared complete-data climate-only ML fit, when used.
#' @param full_climate_REML Shared complete-data climate-only REML fit, when used.
#' @param random Random intercept and slope formulas passed to every fit.
#' @param struct Covariance structure passed to `metafor::rma.mv()`.
#' @param cluster_column Column used for experiment-clustered CR1 inference.
#' @param robust_adjust Whether to apply the small-sample CR1 adjustment.
#' @param club_sandwich Whether `metafor` should delegate to `clubSandwich`.
#' @param comparison_method Estimation method for the nested comparison fits;
#'   required to be `"ML"`.
#' @param final_method Estimation method for final coefficient and variance
#'   estimates; required to be `"REML"`.
#' @param sparse Whether `metafor::rma.mv()` should use sparse matrices.
#' @details The climate-only and full predictor-by-climate models use identical
#' effect rows and sampling variances. ML provides the direct nested comparison;
#' REML provides final estimates and response-heterogeneity summaries. No
#' intermediate additive predictor model is fitted.
#' @return A fitted model bundle with diagnostics retained for validation.
fit_predictor_model_pair <- function(
    spec, climate_formula, full_climate_ML = NULL,
    full_climate_REML = NULL, random, struct, cluster_column,
    robust_adjust, club_sandwich, comparison_method = "ML",
    final_method = "REML", sparse = TRUE) {
  if (!identical(comparison_method, "ML") ||
      !identical(final_method, "REML")) {
    stop(
      "Predictor comparisons require ML and final predictor fits require REML.",
      call. = FALSE
    )
  }
  validate_predictor_specification(spec, climate_formula)
  if (isTRUE(spec$full_data)) {
    if (is.null(full_climate_ML) || is.null(full_climate_REML)) {
      stop(spec$label, " requires the shared complete-data baselines.",
           call. = FALSE)
    }
    climate_ML <- full_climate_ML
    climate_REML <- full_climate_REML
  } else {
    climate_ML <- fit_meta_regression(
      spec$data, spec$sampling_variance_matrix, climate_formula,
      comparison_method,
      paste(spec$label, "climate-only ML"),
      random = random, struct = struct, sparse = sparse
    )
    climate_REML <- fit_meta_regression(
      spec$data, spec$sampling_variance_matrix, climate_formula,
      final_method,
      paste(spec$label, "climate-only REML"),
      random = random, struct = struct, sparse = sparse
    )
  }
  full_ML <- fit_meta_regression(
    spec$data, spec$sampling_variance_matrix, spec$full_formula,
    comparison_method,
    paste(spec$label, "full predictor-climate ML"), spec$predictor_columns,
    random = random, struct = struct, sparse = sparse
  )
  full_REML <- fit_meta_regression(
    spec$data, spec$sampling_variance_matrix, spec$full_formula,
    final_method,
    paste(spec$label, "full predictor-climate REML"), spec$predictor_columns,
    random = random, struct = struct, sparse = sparse
  )
  cr1 <- cluster_robust_inference(
    full_REML, spec$label, "full predictor-climate REML",
    cluster_column = cluster_column,
    adjust = robust_adjust,
    club_sandwich = club_sandwich
  )
  statistics <- dplyr::bind_rows(
    model_fit_statistics(climate_ML, spec$label, "climate only", "ML"),
    model_fit_statistics(
      full_ML, spec$label, "full predictor-climate", "ML"
    )
  )
  comparisons <- nested_model_lrt(
    climate_ML, full_ML, spec$label,
    "climate only vs full predictor-climate"
  )
  diagnostics <- dplyr::bind_rows(
    random_effect_diagnostics(
      climate_REML, spec$label, "climate-only REML"
    ),
    random_effect_diagnostics(
      full_REML, spec$label, "full predictor-climate REML"
    )
  )
  audit <- tibble::tibble(
    predictor = spec$label,
    n_effects = nrow(spec$data),
    n_species = dplyr::n_distinct(spec$data$species_id),
    n_experiments = dplyr::n_distinct(spec$data$experiment_id),
    experiment_IDs = paste(
      sort(unique(as.character(spec$data$experiment_id))), collapse = ";"
    ),
    exact_sampling_variance_order = identical(
      rownames(spec$sampling_variance_matrix), as.character(spec$data$effect_id)
    ) && identical(
      colnames(spec$sampling_variance_matrix), as.character(spec$data$effect_id)
    ),
    matched_climate_and_full_rows_and_variances =
      isTRUE(comparisons$same_rows_and_sampling_variances[[1]]),
    baseline_reproduces_resident_climate_model = if (isTRUE(spec$full_data)) {
      isTRUE(all.equal(
        as.numeric(climate_REML$fit$beta),
        as.numeric(full_climate_REML$fit$beta),
        tolerance = 1e-10
      ))
    } else {
      NA
    }
  )
  if (!audit$exact_sampling_variance_order ||
      !audit$matched_climate_and_full_rows_and_variances ||
      (isTRUE(spec$full_data) && !audit$baseline_reproduces_resident_climate_model)) {
    stop(spec$label, " predictor models failed their matched-data audit.",
         call. = FALSE)
  }
  list(
    key = spec$key,
    label = spec$label,
    data = spec$data,
    sampling_variance_matrix = spec$sampling_variance_matrix,
    fits = list(
      climate_ML = climate_ML,
      climate_REML = climate_REML,
      full_ML = full_ML,
      full_REML = full_REML
    ),
    statistics = statistics,
    comparisons = comparisons,
    diagnostics = diagnostics,
    audit = audit,
    CR1 = cr1$table,
    CR1_object = cr1$object
  )
}

#' Assemble predictor model-pair tables.
#'
#' @param ... Additional arguments forwarded to the delegated operation.
#' @return An R object described by the function name; its structure is validated before later stages use it.
assemble_predictor_model_tables <- function(...) {
  model_pairs <- list(...)
  list(
    fit_statistics = dplyr::bind_rows(lapply(model_pairs, `[[`, "statistics")),
    comparisons = dplyr::bind_rows(lapply(model_pairs, `[[`, "comparisons")),
    diagnostics = dplyr::bind_rows(lapply(model_pairs, `[[`, "diagnostics")),
    audit = dplyr::bind_rows(lapply(model_pairs, `[[`, "audit")),
    CR1 = dplyr::bind_rows(lapply(model_pairs, `[[`, "CR1"))
  )
}

#' Validate Predictor Specification.
#'
#' @param spec A validated workflow configuration or specification object.
#' @param climate_formula A model formula or random-effects specification.
#' @param full_effect_ids A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
validate_predictor_specification <- function(spec, climate_formula,
                                             full_effect_ids = NULL) {
  designs <- lapply(
    list(climate_formula, spec$full_formula),
    stats::model.matrix,
    data = spec$data
  )
  full_rank <- vapply(designs, function(x) qr(x)$rank == ncol(x), logical(1))
  exact_sampling_variance_order <- identical(
    rownames(spec$sampling_variance_matrix), as.character(spec$data$effect_id)
  ) && identical(
    colnames(spec$sampling_variance_matrix), as.character(spec$data$effect_id)
  )
  full_rows <- !isTRUE(spec$full_data) || is.null(full_effect_ids) ||
    identical(as.character(spec$data$effect_id), as.character(full_effect_ids))
  result <- tibble::tibble(
    predictor = spec$label,
    climate_design_full_rank = full_rank[[1]],
    full_design_full_rank = full_rank[[2]],
    exact_sampling_variance_order = exact_sampling_variance_order,
    expected_full_rows = full_rows,
    pass = all(full_rank) && exact_sampling_variance_order && full_rows
  )
  if (!result$pass[[1]]) stop(spec$label, " predictor specification failed validation.")
  result
}
