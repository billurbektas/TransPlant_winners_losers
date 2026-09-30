# Module: R/validation/legacy_phylogeny
# Responsibility: Implement the named scientific workflow stage or shared infrastructure.
# Inputs: Function arguments and explicitly declared target dependencies.
# Outputs: Validated in-memory objects and declared workflow artifacts.
# Side effects: Limited to functions whose contracts explicitly write, copy, or render files.

# Validation-only comparisons between the independently reconstructed
# phylogenetic branch and the frozen transition-era PEV interface.

#' Compare PEV Coordinates with a Historical Reference.
#'
#' @param candidate Reconstructed species-level PEV coordinates.
#' @param reference Frozen species-level PEV coordinates.
#' @param tolerance Maximum permitted absolute coordinate difference.
#' @return A two-row validation table for species order and coordinate values.
compare_pev_coordinates <- function(candidate, reference, tolerance = 1e-8) {
  axes <- c("PEV1", "PEV2", "PEV3")
  candidate <- orient_pev_axes(candidate[, c("species_id", axes)])
  reference <- orient_pev_axes(reference[, c("species_id", axes)])
  ids_equal <- identical(candidate$species_id, reference$species_id)
  values_equal <- ids_equal && all(
    abs(as.matrix(candidate[axes]) - as.matrix(reference[axes])) <= tolerance
  )
  tibble::tibble(
    check = c("PEV species order", "oriented PEV1-PEV3 coordinates"),
    pass = c(ids_equal, values_equal),
    detail = c(
      paste(nrow(candidate), "candidate and", nrow(reference), "reference rows"),
      paste("absolute tolerance", format(tolerance, scientific = TRUE))
    )
  )
}

#' Build a Historical PEV Validation Interface.
#'
#' @param recomputed Independently reconstructed species-level PEV coordinates.
#' @param reference Frozen historical PEV coordinates used only for comparison.
#' @param tolerance Maximum permitted absolute coordinate difference.
#' @return The reconstructed data and a one-row audit proving no frozen values were used for modeling.
build_legacy_pev_validation_interface <- function(recomputed, reference,
                                                   tolerance = 1e-8) {
  checks <- compare_pev_coordinates(recomputed, reference, tolerance)
  if (any(!checks$pass)) {
    stop(
      "Recomputed PEV coordinates exceed the historical-interface tolerance: ",
      paste(checks$check[!checks$pass], collapse = "; "), call. = FALSE
    )
  }
  axes <- c("PEV1", "PEV2", "PEV3")
  maximum_error <- max(
    abs(as.matrix(recomputed[axes]) - as.matrix(reference[axes]))
  )
  list(
    data = recomputed,
    audit = tibble::tibble(
      interface = "independently reconstructed PEV1-PEV3 modeling coordinates",
      recomputed_species = nrow(recomputed),
      maximum_recomputation_error = maximum_error,
      validation_tolerance = tolerance,
      historical_values_used_for_modeling = FALSE,
      reason = paste(
        "The frozen coordinates are comparator-only. The independently",
        "recomputed, response-independently oriented coordinates are supplied",
        "to the optimizer after the migration tolerance audit."
      ),
      pass = TRUE
    )
  )
}

#' Migration Key.
#'
#' @param data A data frame whose required columns are checked before use.
#' @param columns A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
migration_key <- function(data, columns) {
  if (!length(columns)) return(rep("", nrow(data)))
  do.call(
    paste,
    c(lapply(data[columns], function(x) {
      value <- as.character(x)
      value[is.na(value)] <- "<NA>"
      value
    }), sep = "\r")
  )
}

#' Migration Numeric Check.
#'
#' @param domain A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param candidate A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param reference A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param keys A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param columns A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param absolute_tolerance A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param relative_tolerance A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param rounding_digits A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param sign_columns A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param p_columns A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param alpha A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
migration_numeric_check <- function(
    domain, candidate, reference, keys, columns, absolute_tolerance,
    relative_tolerance, rounding_digits = 4L, sign_columns = character(),
    p_columns = character(), alpha = 0.05) {
  required <- unique(c(keys, columns, sign_columns, p_columns))
  assert_columns(candidate, required, paste(domain, "candidate"))
  assert_columns(reference, required, paste(domain, "reference"))
  candidate <- candidate[, required, drop = FALSE]
  reference <- reference[, required, drop = FALSE]
  candidate <- candidate[order(migration_key(candidate, keys)), , drop = FALSE]
  reference <- reference[order(migration_key(reference, keys)), , drop = FALSE]
  keys_equal <- identical(
    migration_key(candidate, keys), migration_key(reference, keys)
  )
  if (!keys_equal || nrow(candidate) != nrow(reference)) {
    return(tibble::tibble(
      domain = domain,
      values_compared = 0L,
      maximum_absolute_difference = Inf,
      absolute_tolerance = absolute_tolerance,
      relative_tolerance = relative_tolerance,
      sign_changed = NA,
      p_threshold_changed = NA,
      rounded_value_changed = NA,
      pass = FALSE,
      detail = "Row identities or ordering differ."
    ))
  }
  candidate_values <- unlist(candidate[columns], use.names = FALSE)
  reference_values <- unlist(reference[columns], use.names = FALSE)
  comparable <- is.finite(candidate_values) & is.finite(reference_values)
  special_equal <- (is.na(candidate_values) & is.na(reference_values)) |
    (is.infinite(candidate_values) & is.infinite(reference_values) &
       sign(candidate_values) == sign(reference_values))
  structure_equal <- all(comparable | special_equal)
  differences <- abs(candidate_values[comparable] - reference_values[comparable])
  limits <- absolute_tolerance +
    relative_tolerance * abs(reference_values[comparable])
  numerical_equal <- structure_equal && all(differences <= limits)
  maximum_error <- if (length(differences)) max(differences) else 0

  sign_changed <- FALSE
  if (length(sign_columns)) {
    candidate_sign <- unlist(candidate[sign_columns], use.names = FALSE)
    reference_sign <- unlist(reference[sign_columns], use.names = FALSE)
    finite_sign <- is.finite(candidate_sign) & is.finite(reference_sign)
    sign_changed <- any(
      sign(candidate_sign[finite_sign]) != sign(reference_sign[finite_sign])
    )
  }
  p_threshold_changed <- FALSE
  if (length(p_columns)) {
    candidate_p <- unlist(candidate[p_columns], use.names = FALSE)
    reference_p <- unlist(reference[p_columns], use.names = FALSE)
    finite_p <- is.finite(candidate_p) & is.finite(reference_p)
    p_threshold_changed <- any(
      (candidate_p[finite_p] < alpha) != (reference_p[finite_p] < alpha)
    )
  }
  rounded_value_changed <- any(
    round(candidate_values[comparable], rounding_digits) !=
      round(reference_values[comparable], rounding_digits)
  )
  pass <- numerical_equal && !sign_changed && !p_threshold_changed
  tibble::tibble(
    domain = domain,
    values_compared = as.integer(length(candidate_values)),
    maximum_absolute_difference = maximum_error,
    absolute_tolerance = absolute_tolerance,
    relative_tolerance = relative_tolerance,
    sign_changed = sign_changed,
    p_threshold_changed = p_threshold_changed,
    rounded_value_changed = rounded_value_changed,
    pass = pass,
    detail = paste(
      "keys exact:", keys_equal,
      "numerically equivalent:", numerical_equal,
      "rounding digits:", rounding_digits
    )
  )
}

#' Build Phylogeny Migration Equivalence.
#'
#' @param active_pev A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param historical_pev A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param active_sequence A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param historical_sequence A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param active_heterogeneity A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param historical_heterogeneity A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param rounding_digits A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param alpha A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
build_phylogeny_migration_equivalence <- function(
    active_pev, historical_pev, active_sequence, historical_sequence,
    active_heterogeneity, historical_heterogeneity,
    rounding_digits = 4L, alpha = 0.05) {
  axes <- c("PEV1", "PEV2", "PEV3")
  random_variance_columns <- c(
    "species_intercept_variance", "species_W_variance",
    "species_T_variance", "species_P_variance", "species_TP_variance"
  )
  random_correlation_columns <- c(
    "T_P_correlation", "T_TP_correlation", "P_TP_correlation",
    "maximum_absolute_context_correlation"
  )
  covariance_diagnostic_columns <- c(
    "minimum_full_covariance_eigenvalue"
  )
  heterogeneity_columns <- c("H_overall", "H_W", "H_context", "H_total")
  synthesis_columns <- c(
    "H_overall baseline", "H_overall predictor",
    "Overall heterogeneity resolved (%)", "H_W baseline", "H_W predictor",
    "Warming heterogeneity resolved (%)", "H_context baseline",
    "H_context predictor", "Context heterogeneity resolved (%)",
    "H_total baseline", "H_total predictor",
    "Total heterogeneity resolved (%)", "Delta AIC", "Delta AICc"
  )
  checks <- dplyr::bind_rows(
    migration_numeric_check(
      "oriented PEV1-PEV3 coordinates", active_pev, historical_pev,
      "species_id", axes, 1e-8, 1e-6, rounding_digits, axes
    ),
    migration_numeric_check(
      "phylogeny fixed-effect coefficients",
      active_sequence$CR1, historical_sequence$CR1,
      c("predictor", "model", "term"),
      "estimate",
      1e-5, 1e-4, rounding_digits, "estimate"
    ),
    migration_numeric_check(
      "phylogeny model-based standard errors, intervals, and p-values",
      active_sequence$CR1, historical_sequence$CR1,
      c("predictor", "model", "term"),
      c(
        "model_SE", "model_CI_low", "model_CI_high",
        "model_p_value"
      ),
      1e-5, 1e-4, rounding_digits, "estimate",
      "model_p_value", alpha
    ),
    migration_numeric_check(
      "phylogeny CR1 standard errors, intervals, p-values, and denominator df",
      active_sequence$CR1, historical_sequence$CR1,
      c("predictor", "model", "term"),
      c(
        "CR1_SE", "CR1_CI_low", "CR1_CI_high",
        "CR1_p_value", "CR1_denominator_df"
      ),
      1e-5, 1e-4, rounding_digits, "estimate", "CR1_p_value", alpha
    ),
    migration_numeric_check(
      "phylogeny random-effect variances",
      active_sequence$diagnostics, historical_sequence$diagnostics,
      c("predictor", "model"), random_variance_columns,
      1e-5, 1e-3, rounding_digits
    ),
    migration_numeric_check(
      "phylogeny random-effect correlations",
      active_sequence$diagnostics, historical_sequence$diagnostics,
      c("predictor", "model"), random_correlation_columns,
      1e-5, 1e-3, rounding_digits
    ),
    migration_numeric_check(
      "phylogeny covariance diagnostics",
      active_sequence$diagnostics, historical_sequence$diagnostics,
      c("predictor", "model"), covariance_diagnostic_columns,
      1e-5, 1e-3, rounding_digits
    ),
    migration_numeric_check(
      "phylogeny ML log-likelihood",
      active_sequence$statistics, historical_sequence$statistics,
      c("predictor", "model", "method"), "logLik",
      1e-7, 1e-6, rounding_digits
    ),
    migration_numeric_check(
      "phylogeny ML AIC",
      active_sequence$statistics, historical_sequence$statistics,
      c("predictor", "model", "method"), "AIC",
      1e-7, 1e-6, rounding_digits
    ),
    migration_numeric_check(
      "phylogeny ML AICc",
      active_sequence$statistics, historical_sequence$statistics,
      c("predictor", "model", "method"), "AICc",
      1e-7, 1e-6, rounding_digits
    ),
    migration_numeric_check(
      "phylogeny nested-model comparisons",
      active_sequence$comparisons, historical_sequence$comparisons,
      c("predictor", "comparison"),
      c("likelihood_ratio", "df", "p_value"),
      1e-7, 1e-6, rounding_digits, p_columns = "p_value", alpha = alpha
    ),
    migration_numeric_check(
      "phylogeny heterogeneity components",
      dplyr::filter(
        active_heterogeneity$components,
        .data$Predictor == active_sequence$label
      ),
      dplyr::filter(
        historical_heterogeneity$components,
        .data$Predictor == historical_sequence$label
      ),
      c("Predictor", "model_role"), heterogeneity_columns,
      1e-6, 1e-4, rounding_digits
    ),
    migration_numeric_check(
      "phylogeny heterogeneity-resolution percentages",
      dplyr::filter(
        active_heterogeneity$synthesis,
        .data$Predictor == active_sequence$label
      ),
      dplyr::filter(
        historical_heterogeneity$synthesis,
        .data$Predictor == historical_sequence$label
      ),
      "Predictor",
      grep("resolved", synthesis_columns, value = TRUE),
      2e-3, 1e-4, rounding_digits
    ),
    migration_numeric_check(
      "phylogeny heterogeneity synthesis levels and fit changes",
      dplyr::filter(
        active_heterogeneity$synthesis,
        .data$Predictor == active_sequence$label
      ),
      dplyr::filter(
        historical_heterogeneity$synthesis,
        .data$Predictor == historical_sequence$label
      ),
      "Predictor",
      setdiff(synthesis_columns, grep("resolved", synthesis_columns, value = TRUE)),
      2e-3, 1e-4, rounding_digits
    )
  )
  interpretation_changed <- any(
    !checks$pass | checks$sign_changed %in% TRUE |
      checks$p_threshold_changed %in% TRUE
  )
  dplyr::bind_rows(
    checks,
    tibble::tibble(
      domain = "scientific interpretation",
      values_compared = sum(checks$values_compared),
      maximum_absolute_difference = max(checks$maximum_absolute_difference),
      absolute_tolerance = NA_real_,
      relative_tolerance = NA_real_,
      sign_changed = any(checks$sign_changed %in% TRUE),
      p_threshold_changed = any(checks$p_threshold_changed %in% TRUE),
      rounded_value_changed = any(checks$rounded_value_changed %in% TRUE),
      pass = !interpretation_changed,
      detail = if (interpretation_changed) {
        "At least one prespecified numerical or interpretive migration criterion changed."
      } else {
        paste0(
          "No coefficient sign, p<", alpha,
          " decision, or scientific interpretation changed. Four-decimal",
          " display changes, if any, are reported separately."
        )
      }
    )
  )
}

#' Phylogeny Migration Release Check.
#'
#' @param audit A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
phylogeny_migration_release_check <- function(audit) {
  tibble::tibble(
    check = "active reconstructed phylogeny branch versus historical fixture refit",
    pass = all(audit$pass),
    detail = paste(
      sum(audit$pass), "of", nrow(audit),
      "migration equivalence criteria passed; fixtures were comparators only"
    )
  )
}

#' Format Phylogeny Migration Report.
#'
#' @param audit A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
format_phylogeny_migration_report <- function(audit) {
  format_value <- function(x) {
    ifelse(is.finite(x), format(x, digits = 10, scientific = TRUE), as.character(x))
  }
  rows <- vapply(seq_len(nrow(audit)), function(index) {
    paste0(
      "- ", audit$domain[[index]], ": ",
      if (isTRUE(audit$pass[[index]])) "PASS" else "FAIL",
      "; n=", audit$values_compared[[index]],
      "; max abs difference=", format_value(audit$maximum_absolute_difference[[index]]),
      "; sign changed=", audit$sign_changed[[index]],
      "; p<0.05 changed=", audit$p_threshold_changed[[index]],
      "; rounded value changed=", audit$rounded_value_changed[[index]],
      "; ", audit$detail[[index]]
    )
  }, character(1))
  c(
    "PHYLOGENY MIGRATION EQUIVALENCE REPORT",
    "======================================",
    "",
    paste(
      "The active branch rebuilds the 526-species tree, PCoA, orientation,",
      "standardization, phylogeny models, and heterogeneity from canonical inputs."
    ),
    paste(
      "Frozen transition PEV values are used only in a separate historical",
      "comparison refit and never as inputs to active scientific targets."
    ),
    "",
    rows,
    "",
    paste0(
      "OVERALL: ", if (all(audit$pass)) "PASS" else "FAIL",
      ". Scientific interpretation changed: ", !all(audit$pass), "."
    )
  )
}
