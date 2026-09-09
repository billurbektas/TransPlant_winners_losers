# Numerical regression against the immutable Step 22 authority.

step22_reference_fields <- function() {
  c(
    "metadata", "data_validation", "scaling", "M_both", "classification",
    "family_model", "predictor_model_fit", "predictor_model_comparisons",
    "predictor_CR1", "heterogeneity_components", "heterogeneity_synthesis",
    "phylogenetic_positive_variation_percent_PEV1_PEV3"
  )
}

validate_step22_schema <- function(x, label = "Step 22 object") {
  if (!is.list(x) || !identical(names(x), step22_reference_fields())) {
    stop(label, " does not retain the exact 12-field Step 22 schema.", call. = FALSE)
  }
  invisible(x)
}

reference_tolerance <- function(path) {
  joined <- paste(path, collapse = "/")
  if (grepl("(^|/)(scaling|data_validation)(/|$)", joined)) {
    return(c(absolute = 1e-12, relative = 1e-10))
  }
  if (grepl("random|variance_profiles", joined, ignore.case = TRUE)) {
    return(c(absolute = 1e-6, relative = 1e-4))
  }
  c(absolute = 1e-8, relative = 1e-6)
}

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

compare_reference_node <- function(candidate, reference, path = character()) {
  location <- if (length(path)) paste(path, collapse = "/") else "<root>"
  failures <- character()

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

validate_step22_reference <- function(candidate, reference) {
  validate_step22_schema(candidate, "candidate")
  validate_step22_schema(reference, "reference")
  failures <- compare_reference_node(candidate, reference)
  tibble::tibble(
    check = c("exact_12_field_schema", "numerical_and_structural_equivalence"),
    pass = c(TRUE, !length(failures)),
    detail = c(
      paste(step22_reference_fields(), collapse = "; "),
      if (length(failures)) paste(failures, collapse = " | ") else
        "All fields pass the prespecified exact or numerical tolerance policy."
    )
  )
}

assert_step22_reference <- function(validation) {
  if (any(!validation$pass)) {
    stop(
      "Step 22 numerical regression failed: ",
      paste(validation$detail[!validation$pass], collapse = " | "),
      call. = FALSE
    )
  }
  validation
}

assert_reference_table <- function(candidate, reference, label,
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
  if (any(!result$pass)) {
    stop(label, " regression failed: ", result$detail[[1]], call. = FALSE)
  }
  result
}

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

assemble_release_validation <- function(checks) {
  if (!is.list(checks) || is.null(names(checks)) ||
      any(!nzchar(names(checks))) || anyDuplicated(names(checks))) {
    stop("Release checks require a uniquely named list.", call. = FALSE)
  }
  dplyr::bind_rows(Map(release_check_summary, names(checks), checks))
}

assert_release_validation <- function(validation) {
  if (!is.data.frame(validation) || !"pass" %in% names(validation)) {
    stop("Release validation must contain a pass column.", call. = FALSE)
  }
  failed <- is.na(validation$pass) | !validation$pass
  if (any(failed)) {
    stop(
      "Canonical release validation failed: ",
      paste(validation$gate[failed], collapse = ", "),
      call. = FALSE
    )
  }
  validation
}

assemble_step22_candidate <- function(
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
  validate_step22_schema(result, "assembled candidate")
  result
}
