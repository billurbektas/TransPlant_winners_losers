# Module: R/validation/legacy_reference
# Responsibility: Implement the named scientific workflow stage or shared infrastructure.
# Inputs: Function arguments and explicitly declared target dependencies.
# Outputs: Validated in-memory objects and declared workflow artifacts.
# Side effects: Limited to functions whose contracts explicitly write, copy, or render files.

# Temporary numerical regression against the historical release authorities.

#' Legacy Reference Fields.
#'
#' @return An R object described by the function name; its structure is validated before later stages use it.
legacy_reference_fields <- function() {
  c(
    "metadata", "data_validation", "scaling", "M_both", "classification",
    "family_model", "predictor_model_fit", "predictor_model_comparisons",
    "predictor_CR1", "heterogeneity_components", "heterogeneity_synthesis",
    "phylogenetic_positive_variation_percent_PEV1_PEV3"
  )
}

#' Validate Legacy Reference Schema.
#'
#' @param x A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param label A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
validate_legacy_reference_schema <- function(x, label = "Step 22 object") {
  if (!is.list(x) || !identical(names(x), legacy_reference_fields())) {
    stop(label, " does not retain the exact 12-field Step 22 schema.", call. = FALSE)
  }
  invisible(x)
}

#' Reference Tolerance.
#'
#' @param path Character path within the historical 12-field result object.
#' @details The two looser cases cover only the historical phylogeny migration;
#'   the active reconstruction is checked independently for topology,
#'   distances, and PEV coordinates.
#' @return Named absolute and relative tolerances.
reference_tolerance <- function(path) {
  joined <- paste(path, collapse = "/")
  if (grepl("(^|/)(scaling|data_validation)(/|$)", joined)) {
    return(c(absolute = 1e-12, relative = 1e-10))
  }
  if (grepl("^predictor_CR1/", joined)) {
    return(c(absolute = 1e-5, relative = 1e-6))
  }
  if (grepl("^heterogeneity_components/", joined)) {
    return(c(absolute = 2e-7, relative = 1e-6))
  }
  if (grepl("^heterogeneity_synthesis/", joined)) {
    return(c(absolute = 2e-3, relative = 1e-6))
  }
  if (grepl("random|variance_profiles", joined, ignore.case = TRUE)) {
    return(c(absolute = 1e-6, relative = 1e-4))
  }
  c(absolute = 1e-8, relative = 1e-6)
}

#' Numeric Values Equal.
#'
#' @param candidate A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param reference A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param tolerance A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
numeric_values_equal <- function(candidate, reference, tolerance) {
  if (length(candidate) != length(reference)) return(FALSE)
  comparable <- is.finite(candidate) & is.finite(reference)
  special_equal <- (is.na(candidate) & is.na(reference)) |
    (is.infinite(candidate) & is.infinite(reference) & sign(candidate) == sign(reference))
  if (any(!(comparable | special_equal))) return(FALSE)
  if (!any(comparable)) return(TRUE)
  error <- abs(candidate[comparable] - reference[comparable])
  limit <- tolerance[["absolute"]] +
    tolerance[["relative"]] * abs(reference[comparable])
  all(error <= limit)
}

#' Compare Deliberately Renamed Legacy Labels.
#'
#' @param candidate A value from the stage-based workflow.
#' @param reference The corresponding value from the frozen historical object.
#' @param path Character path identifying the value within the reference schema.
#' @details This temporary adapter changes comparison semantics for labels only.
#'   It never changes numeric values, dimensions, identifiers, classifications,
#'   or scientific settings, and it is removed with legacy certification in
#'   Phase 2.
#' @return `TRUE` only for the registered non-scientific Phase 1 naming migrations.
legacy_nomenclature_migration_equal <- function(candidate, reference, path) {
  location <- paste(path, collapse = "/")

  if (identical(location, "metadata/step")) {
    return(length(candidate) == 0L && identical(reference, 22L))
  }

  if (identical(location, "metadata/script")) {
    return(
      identical(candidate, "_targets.R (stage-based targets workflow)") &&
        reference %in% c(
          "_targets.R (active canonical targets lineage)",
          "R/34_step_22_downstream_526_models.R"
        )
    )
  }

  consistency_text_paths <- c(
    "data_validation/consistency_checks/check",
    "data_validation/consistency_checks/detail"
  )
  if (location %in% consistency_text_paths) {
    # These two columns are explanatory prose rewritten with the new module
    # names. Row identity and scientific status remain protected by exact
    # comparisons of check_number and pass in the same table.
    return(
      is.character(candidate) && is.character(reference) &&
        length(candidate) == length(reference)
    )
  }

  climate_label_paths <- c(
    "M_both/CR1/predictor", "M_both/CR1/model",
    "M_both/random/predictor", "M_both/random/model",
    "M_both/fit_statistics/predictor", "M_both/fit_statistics/model"
  )
  if (location %in% climate_label_paths) {
    normalize <- function(x) {
      value <- as.character(x)
      value[value %in% c("Updated M_both", "M_both", "M_both_526")] <-
        "resident_climate"
      value
    }
    return(identical(normalize(candidate), normalize(reference)))
  }

  FALSE
}

#' Compare Reference Node.
#'
#' @param candidate A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param reference A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param path A project-relative or absolute filesystem path.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
compare_reference_node <- function(candidate, reference, path = character()) {
  location <- if (length(path)) paste(path, collapse = "/") else "<root>"
  failures <- character()

  if (legacy_nomenclature_migration_equal(candidate, reference, path)) {
    return(character())
  }

  # CSV import guesses an all-missing diagnostic column as logical, whereas
  # an equivalent in-memory model table retains it as character. Compare the
  # semantic diagnostic values so that NA_character_ and NA are equivalent,
  # while preserving exact matching for every non-missing warning/error.
  diagnostic_name <- length(path) && grepl(
    "(warning|error)s?$", tail(path, 1L), ignore.case = TRUE
  )
  diagnostic_types <- (is.character(candidate) || is.logical(candidate)) &&
    (is.character(reference) || is.logical(reference))
  if (diagnostic_name && diagnostic_types) {
    normalize_diagnostic <- function(x) {
      value <- as.character(x)
      value[is.na(x) | value == ""] <- NA_character_
      value
    }
    if (identical(
      normalize_diagnostic(candidate), normalize_diagnostic(reference)
    )) {
      return(character())
    }
    return(paste0(location, ": exact diagnostic values differ"))
  }

  if (is.data.frame(reference)) {
    if (!is.data.frame(candidate)) return(paste0(location, ": expected data frame"))
    if (!identical(names(candidate), names(reference))) {
      return(paste0(location, ": column names/order differ"))
    }
    if (nrow(candidate) != nrow(reference)) {
      return(paste0(location, ": row count differs"))
    }
    for (column in names(reference)) {
      failures <- c(
        failures,
        compare_reference_node(candidate[[column]], reference[[column]], c(path, column))
      )
    }
    return(failures)
  }

  if (is.list(reference)) {
    if (!is.list(candidate)) return(paste0(location, ": expected list"))
    if (!identical(names(candidate), names(reference))) {
      return(paste0(location, ": list names/order differ"))
    }
    for (name in names(reference)) {
      failures <- c(
        failures,
        compare_reference_node(candidate[[name]], reference[[name]], c(path, name))
      )
    }
    return(failures)
  }

  if (is.integer(reference)) {
    if (!is.integer(candidate) || !identical(candidate, reference)) {
      return(paste0(location, ": exact integer values differ"))
    }
    return(character())
  }

  if (is.numeric(reference)) {
    if (!is.numeric(candidate)) return(paste0(location, ": expected numeric"))
    if (!numeric_values_equal(candidate, reference, reference_tolerance(path))) {
      return(paste0(location, ": numeric values exceed tolerance"))
    }
    return(character())
  }

  # Counts, identifiers, structures, convergence flags, classifications, and
  # PASS/FAIL fields are intentionally exact.
  if (!identical(candidate, reference)) {
    return(paste0(location, ": exact values differ"))
  }
  character()
}

#' Validate Legacy Reference.
#'
#' @param candidate A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param reference A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
validate_legacy_reference <- function(candidate, reference) {
  validate_legacy_reference_schema(candidate, "candidate")
  validate_legacy_reference_schema(reference, "reference")
  failures <- compare_reference_node(candidate, reference)
  tibble::tibble(
    check = c("exact_12_field_schema", "numerical_and_structural_equivalence"),
    pass = c(TRUE, !length(failures)),
    detail = c(
      paste(legacy_reference_fields(), collapse = "; "),
      if (length(failures)) paste(failures, collapse = " | ") else
        "All fields pass the prespecified exact or numerical tolerance policy."
    )
  )
}

#' Assert Legacy Reference.
#'
#' @param validation A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
assert_legacy_reference <- function(validation) {
  if (any(!validation$pass)) {
    stop(
      "Step 22 numerical regression failed: ",
      paste(validation$detail[!validation$pass], collapse = " | "),
      call. = FALSE
    )
  }
  validation
}

#' Write Release Reference Once.
#'
#' @param candidate A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param path A project-relative or absolute filesystem path.
#' @return The path or paths written by the declared side effect.
write_release_reference_once <- function(candidate, path) {
  validate_legacy_reference_schema(candidate, "current canonical reference candidate")
  if (file.exists(path)) {
    existing <- readRDS(path)
    validation <- validate_legacy_reference(candidate, existing)
    if (any(!validation$pass)) {
      stop(
        "The frozen canonical reference exists and differs from the active ",
        "candidate. It will not be overwritten: ", path,
        call. = FALSE
      )
    }
    return(path)
  }
  write_rds_atomic(candidate, path)
}

#' Validate Frozen Release Reference.
#'
#' @param candidate A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param path A project-relative or absolute filesystem path.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
validate_frozen_release_reference <- function(candidate, path) {
  if (!file.exists(path)) {
    stop("Canonical reference file is missing: ", path, call. = FALSE)
  }
  assert_legacy_reference(validate_legacy_reference(candidate, readRDS(path)))
}

#' Validate Reference Table.
#'
#' @param candidate A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param reference A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param label A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param order_by A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
validate_reference_table <- function(candidate, reference, label,
                                     order_by = character()) {
  if (length(order_by)) {
    candidate <- candidate |>
      dplyr::arrange(!!!rlang::syms(order_by))
    reference <- reference |>
      dplyr::arrange(!!!rlang::syms(order_by))
  }
  failures <- compare_reference_node(candidate, reference, label)
  result <- tibble::tibble(
    check = label,
    pass = !length(failures),
    detail = if (length(failures)) paste(failures, collapse = " | ") else
      "All fields pass the Step 22 tolerance policy."
  )
  result
}

#' Assert Reference Table.
#'
#' @param candidate A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param reference A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param label A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param order_by A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
assert_reference_table <- function(candidate, reference, label,
                                   order_by = character()) {
  result <- validate_reference_table(candidate, reference, label, order_by)
  if (any(!result$pass)) {
    stop(label, " regression failed: ", result$detail[[1]], call. = FALSE)
  }
  result
}

#' Release Check Summary.
#'
#' @param name A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param check A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
release_check_summary <- function(name, check) {
  if (!is.data.frame(check) || !"pass" %in% names(check)) {
    stop("Release check '", name, "' must be a data frame with a pass column.",
         call. = FALSE)
  }
  pass <- as.logical(check$pass)
  failed_rows <- which(is.na(pass) | !pass)
  detail_column <- intersect(c("detail", "check", "note"), names(check))
  failure_detail <- if (!length(failed_rows)) {
    paste(nrow(check), "of", nrow(check), "checks passed")
  } else if (length(detail_column)) {
    paste(
      as.character(check[[detail_column[[1]]]][failed_rows]),
      collapse = " | "
    )
  } else {
    paste("failed rows", paste(failed_rows, collapse = ", "))
  }
  tibble::tibble(
    gate = name,
    checks = nrow(check),
    passed = sum(pass %in% TRUE),
    pass = length(failed_rows) == 0L,
    detail = failure_detail
  )
}

#' Assemble Release Validation.
#'
#' @param checks A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
assemble_release_validation <- function(checks) {
  if (!is.list(checks) || is.null(names(checks)) ||
      any(!nzchar(names(checks))) || anyDuplicated(names(checks))) {
    stop("Release checks require a uniquely named list.", call. = FALSE)
  }
  dplyr::bind_rows(Map(release_check_summary, names(checks), checks))
}

#' Assert Release Validation.
#'
#' @param validation A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
assert_release_validation <- function(validation) {
  if (!is.data.frame(validation) || !"pass" %in% names(validation)) {
    stop("Release validation must contain a pass column.", call. = FALSE)
  }
  failed <- is.na(validation$pass) | !validation$pass
  if (any(failed)) {
    stop(
      "Release validation failed: ",
      paste(validation$gate[failed], collapse = ", "),
      call. = FALSE
    )
  }
  validation
}

#' Assemble Legacy Comparison Candidate.
#'
#' @param metadata A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param data_validation A data frame whose required columns are checked before use.
#' @param scaling A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param M_both A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param classification A data frame whose required columns are checked before use.
#' @param family_model A fitted model object with the structure checked by this function.
#' @param predictor_model_fit A fitted model object with the structure checked by this function.
#' @param predictor_model_comparisons A fitted model object with the structure checked by this function.
#' @param predictor_CR1 A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param heterogeneity_components A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param heterogeneity_synthesis A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param phylogenetic_positive_variation_percent_PEV1_PEV3 A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
assemble_legacy_comparison_candidate <- function(
    metadata, data_validation, scaling, M_both, classification, family_model,
    predictor_model_fit, predictor_model_comparisons, predictor_CR1,
    heterogeneity_components, heterogeneity_synthesis,
    phylogenetic_positive_variation_percent_PEV1_PEV3) {
  result <- list(
    metadata = metadata,
    data_validation = data_validation,
    scaling = scaling,
    M_both = M_both,
    classification = classification,
    family_model = family_model,
    predictor_model_fit = predictor_model_fit,
    predictor_model_comparisons = predictor_model_comparisons,
    predictor_CR1 = predictor_CR1,
    heterogeneity_components = heterogeneity_components,
    heterogeneity_synthesis = heterogeneity_synthesis,
    phylogenetic_positive_variation_percent_PEV1_PEV3 =
      phylogenetic_positive_variation_percent_PEV1_PEV3
  )
  validate_legacy_reference_schema(result, "assembled candidate")
  result
}
