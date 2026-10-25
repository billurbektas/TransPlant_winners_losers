# Module: R/validation/release_summary
# Responsibility: Assemble current scientific release checks and readable model report.
# Inputs: Validated models, classifications, comparisons, and figures.
# Outputs: Release-validation result and current model text.
# Side effects: None; the publication target writes the returned text.

# Human-readable current model report assembled only from active scientific
# objects. Frozen reports and transition fixtures are not inputs.

#' Summarize one named release check
#'
#' @param name Name of the check.
#' @param check Data frame containing a logical `pass` column.
#' @return One row with the number of passing checks and failure details.
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

#' Collect named current-release checks
#'
#' @param checks Uniquely named list of validation data frames.
#' @return A table with one summary row per check.
assemble_release_validation <- function(checks) {
  if (!is.list(checks) || is.null(names(checks)) ||
      any(!nzchar(names(checks))) || anyDuplicated(names(checks))) {
    stop("Release checks require a uniquely named list.", call. = FALSE)
  }
  dplyr::bind_rows(Map(release_check_summary, names(checks), checks))
}

#' Require every current-release check to pass
#'
#' @param validation Release-check table with `gate` and logical `pass`.
#' @return The unchanged table; stops with failed gate names if any fail.
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

#' Report Number.
#'
#' @param x A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param digits A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
report_number <- function(x, digits = 7L) {
  ifelse(
    is.na(x), "NA",
    format(as.numeric(x), digits = digits, scientific = FALSE, trim = TRUE)
  )
}

#' Report Section.
#'
#' @param title A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param lines A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
report_section <- function(title, lines) {
  c("", title, strrep("-", nchar(title)), lines)
}

#' Current Release Checks.
#'
#' @param scope A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param predictor_validation A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param predictor_tables A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param family_model_summary A fitted model object with the structure checked by this function.
#' @param phylogeny_rebuild_audit A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param zero_distance_audit A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param heterogeneity_components A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param heterogeneity_synthesis A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
current_release_checks <- function(
    scope, predictor_validation, predictor_tables,
    family_model_summary, phylogeny_rebuild_audit, zero_distance_audit,
    heterogeneity_components, heterogeneity_synthesis) {
  tibble::tibble(
    check = c(
      "resident-model scope", "diagonal sampling variances",
      "predictor specifications", "predictor row/variance audits",
      "predictor convergence and covariance diagnostics",
      "family model diagnostics", "phylogeny reconstruction",
      "documented zero-distance placements",
      "heterogeneity components and synthesis"
    ),
    pass = c(
      scope$effects == 1604L && scope$species == 526L &&
        scope$experiments == 42L,
      isTRUE(scope$sampling_variance_validation$valid[[1]]),
      all(predictor_validation$pass),
      all(predictor_tables$audit$exact_sampling_variance_order &
            predictor_tables$audit$matched_climate_and_full_rows_and_variances &
            (is.na(predictor_tables$audit$baseline_reproduces_resident_climate_model) |
               predictor_tables$audit$baseline_reproduces_resident_climate_model)),
      all(predictor_tables$diagnostics$covariance_positive_definite &
            predictor_tables$diagnostics$covariance_identifiable &
            predictor_tables$diagnostics$finite_coefficients &
            predictor_tables$diagnostics$finite_fixed_covariance),
      isTRUE(family_model_summary$sampling_variance_valid[[1]]) &&
        isTRUE(family_model_summary$covariance_positive_definite[[1]]) &&
        isTRUE(family_model_summary$random_structure_identifiable[[1]]) &&
        isTRUE(family_model_summary$finite_coefficients[[1]]) &&
        isTRUE(family_model_summary$finite_fixed_covariance[[1]]),
      all(phylogeny_rebuild_audit$pass),
      all(zero_distance_audit$classification == "intentional shared placement"),
      nrow(heterogeneity_synthesis) == 6L &&
        nrow(heterogeneity_components) == 12L &&
        all(vapply(
          heterogeneity_components[vapply(
            heterogeneity_components, is.numeric, logical(1)
          )],
          function(x) all(is.finite(x)), logical(1)
        )) &&
        all(vapply(
          heterogeneity_synthesis[vapply(
            heterogeneity_synthesis, is.numeric, logical(1)
          )],
          function(x) all(is.finite(x)), logical(1)
        ))
    ),
    detail = c(
      paste(scope$effects, "effects;", scope$species, "species;",
            scope$experiments, "experiments"),
      paste(
        "minimum sampling variance",
        report_number(
          scope$sampling_variance_validation$minimum_sampling_variance
        )
      ),
      paste(sum(predictor_validation$pass), "of", nrow(predictor_validation), "checks passed"),
      paste(nrow(predictor_tables$audit), "predictor model pairs"),
      paste(nrow(predictor_tables$diagnostics), "REML diagnostic rows"),
      paste(family_model_summary$n_effects, "effects;",
            family_model_summary$n_families, "families"),
      paste(sum(phylogeny_rebuild_audit$pass), "of",
            nrow(phylogeny_rebuild_audit), "checks passed"),
      paste(nrow(zero_distance_audit), "zero-distance pairs"),
      paste(nrow(heterogeneity_components), "component rows;",
            nrow(heterogeneity_synthesis), "synthesis rows")
    )
  )
}

#' Format Classification Summary.
#'
#' @param summary A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
format_classification_summary <- function(summary) {
  paste0(
    summary$Resolution, ": eligible=", summary$`N eligible units`,
    "; winners=", summary$`Winner N`, " (", report_number(summary$`Winner %`), "%)",
    "; losers=", summary$`Loser N`, " (", report_number(summary$`Loser %`), "%)",
    "; context-dependent=", summary$`Context-dependent N`, " (",
    report_number(summary$`Context-dependent %`), "%)"
  )
}

#' Format Model Comparisons.
#'
#' @param comparisons A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
format_model_comparisons <- function(comparisons) {
  vapply(seq_len(nrow(comparisons)), function(index) {
    row <- comparisons[index, ]
    paste0(
      row$predictor, " | ", row$comparison,
      " | LR=", report_number(row$likelihood_ratio),
      " | df=", report_number(row$df),
      " | p=", report_number(row$p_value),
      " | valid=", row$valid_nested_comparison
    )
  }, character(1))
}

#' Format Heterogeneity Synthesis.
#'
#' @param synthesis A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
format_heterogeneity_synthesis <- function(synthesis) {
  vapply(seq_len(nrow(synthesis)), function(index) {
    row <- synthesis[index, ]
    paste0(
      row$Predictor,
      " | n=", row$`N effects`, "/", row$`N species`, "/", row$`N experiments`,
      " | H_overall=", report_number(row$`H_overall baseline`), " -> ",
      report_number(row$`H_overall predictor`), " (",
      report_number(row$`Overall heterogeneity resolved (%)`), "%)",
      " | H_W=", report_number(row$`H_W baseline`), " -> ",
      report_number(row$`H_W predictor`), " (",
      report_number(row$`Warming heterogeneity resolved (%)`), "%)",
      " | H_context=", report_number(row$`H_context baseline`), " -> ",
      report_number(row$`H_context predictor`), " (",
      report_number(row$`Context heterogeneity resolved (%)`), "%)",
      " | H_total=", report_number(row$`H_total baseline`), " -> ",
      report_number(row$`H_total predictor`), " (",
      report_number(row$`Total heterogeneity resolved (%)`), "%)",
      " | ", row$`ML climate-only vs full result`,
      " | Delta AIC (full - climate only)=",
      report_number(row$`Delta AIC full minus climate only`),
      " | Delta AICc (full - climate only)=",
      report_number(row$`Delta AICc full minus climate only`)
    )
  }, character(1))
}

#' Build Current Model Report.
#'
#' @param config A validated workflow configuration or specification object.
#' @param scope A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param climate_scaling A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param climate_model_results A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param species_summary A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param family_summary A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param trait_scaling A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param predictor_comparisons A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param heterogeneity_synthesis A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param release_checks A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
build_current_model_report <- function(
    config, scope, climate_scaling, climate_model_results, species_summary,
    family_summary, trait_scaling, predictor_comparisons,
    heterogeneity_synthesis, release_checks) {
  fixed_lines <- vapply(seq_len(nrow(climate_model_results$fixed)), function(index) {
    row <- climate_model_results$fixed[index, ]
    paste0(
      row$term, ": estimate=", report_number(row$estimate),
      "; SE=", report_number(row$model_SE),
      "; 95% CI=[", report_number(row$model_CI_low), ", ",
      report_number(row$model_CI_high), "]; p=", report_number(row$model_p_value)
    )
  }, character(1))
  cr1_lines <- vapply(seq_len(nrow(climate_model_results$CR1)), function(index) {
    row <- climate_model_results$CR1[index, ]
    paste0(
      row$term, ": estimate=", report_number(row$estimate),
      "; CR1 SE=", report_number(row$CR1_SE),
      "; CR1 95% CI=[", report_number(row$CR1_CI_low), ", ",
      report_number(row$CR1_CI_high), "]; CR1 p=",
      report_number(row$CR1_p_value), "; df=",
      report_number(row$CR1_denominator_df)
    )
  }, character(1))
  random_row <- climate_model_results$random[1, ]
  random_lines <- c(
    paste0("species intercept variance=", report_number(random_row$species_intercept_variance)),
    paste0("species warming-slope variance=", report_number(random_row$species_W_variance)),
    paste0("species T/P/TP variances=", paste(report_number(c(
      random_row$species_T_variance, random_row$species_P_variance,
      random_row$species_TP_variance
    )), collapse = ", ")),
    paste0("T-P/T-TP/P-TP correlations=", paste(report_number(c(
      random_row$T_P_correlation, random_row$T_TP_correlation,
      random_row$P_TP_correlation
    )), collapse = ", ")),
    paste0(
      "covariance positive definite=", random_row$covariance_positive_definite,
      "; identifiable=", random_row$covariance_identifiable,
      "; any boundary variance=", random_row$any_variance_at_boundary
    )
  )
  profile_lines <- vapply(seq_len(nrow(climate_model_results$variance_profiles)), function(index) {
    row <- climate_model_results$variance_profiles[index, ]
    paste0(
      row$component, ": estimate=", report_number(row$estimate),
      "; profile 95% CI=[", report_number(row$CI_low), ", ",
      report_number(row$CI_high), "]"
    )
  }, character(1))
  scaling_lines <- vapply(seq_len(nrow(climate_scaling)), function(index) {
    row <- climate_scaling[index, ]
    paste0(
      row$scaled_field, " <- ", row$source_field,
      ": center=", report_number(row$center),
      "; scale=", report_number(row$scale),
      "; weighting=", row$weighting, "; n=", row$n
    )
  }, character(1))
  trait_lines <- vapply(seq_len(nrow(trait_scaling)), function(index) {
    row <- trait_scaling[index, ]
    paste0(
      row$trait, ": species with trait=", row$species_N_with_trait,
      "; missing=", row$species_N_missing,
      "; log center=", report_number(row$log_trait_mean),
      "; log SD=", report_number(row$log_trait_SD),
      "; effects=", row$n_effects, "; weighting=", row$weighting
    )
  }, character(1))
  comparison_groups <- list(
    "Trait model comparisons" = c("LDMC", "LMA", "Height"),
    "Graminoid and Woodiness model comparisons" = c("Graminoid status", "Woodiness"),
    "PEV1-PEV3 model comparisons" = "Phylogenetic position (PEV1-PEV3)"
  )
  report <- c(
    "CURRENT RESIDENT AND PREDICTOR MODEL RESULTS",
    "============================================",
    "",
    "Generated from the current stage-based workflow; historical reports are not inputs.",
    report_section("Dataset and exclusions", c(
      paste("Effects:", scope$effects), paste("Species:", scope$species),
      paste("Experiments:", scope$experiments),
      paste("Excluded resident-model taxa:", paste(config$exclusions, collapse = "; ")),
      "The excluded taxa remain in the community winner/loser calculations."
    )),
    report_section("Climate scaling", scaling_lines),
    report_section("Resident climate fixed effects", fixed_lines),
    report_section("Resident climate experiment-clustered CR1 fixed effects", cr1_lines),
    report_section("Resident climate random-effect variances and correlations", random_lines),
    report_section("Resident climate variance-component profile intervals", profile_lines),
    report_section("Species classification", format_classification_summary(species_summary)),
    report_section("Family classification", format_classification_summary(family_summary)),
    report_section("Trait scaling", trait_lines)
  )
  for (title in names(comparison_groups)) {
    report <- c(
      report,
      report_section(
        title,
        format_model_comparisons(dplyr::filter(
          predictor_comparisons,
          .data$predictor %in% comparison_groups[[title]]
        ))
      )
    )
  }
  release_lines <- vapply(seq_len(nrow(release_checks)), function(index) {
    paste0(
      release_checks$check[[index]], ": ",
      if (isTRUE(release_checks$pass[[index]])) "PASS" else "FAIL",
      " — ", release_checks$detail[[index]]
    )
  }, character(1))
  c(
    report,
    report_section(
      "Response-scale heterogeneity-resolution synthesis",
      format_heterogeneity_synthesis(heterogeneity_synthesis)
    ),
    report_section("Scientific-data checks (not the final release gate)", release_lines),
    "",
    paste(
      "SCIENTIFIC-DATA CHECK STATUS:",
      if (all(release_checks$pass)) "PASS" else "FAIL"
    ),
    "Final release status is recorded separately in output/audits/release_validation.csv."
  )
}
