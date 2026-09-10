# Human-readable current downstream report assembled only from active canonical
# scientific objects. Frozen reports and transition fixtures are not inputs.

report_number <- function(x, digits = 7L) {
  ifelse(
    is.na(x), "NA",
    format(as.numeric(x), digits = digits, scientific = FALSE, trim = TRUE)
  )
}

report_section <- function(title, lines) {
  c("", title, strrep("-", nchar(title)), lines)
}

active_downstream_release_checks <- function(
    scope, consistency_checks, predictor_validation, predictor_tables,
    family_model_summary, phylogeny_rebuild_audit, zero_distance_audit,
    heterogeneity_components, heterogeneity_synthesis) {
  tibble::tibble(
    check = c(
      "downstream scope", "sampling covariance", "downstream consistency",
      "predictor specifications", "predictor row/covariance audits",
      "predictor convergence and covariance diagnostics",
      "family model diagnostics", "phylogeny reconstruction",
      "documented zero-distance placements",
      "heterogeneity components and synthesis"
    ),
    pass = c(
      scope$effects == 1604L && scope$species == 526L &&
        scope$experiments == 42L,
      isTRUE(scope$V_validation$valid[[1]]),
      all(consistency_checks$pass),
      all(predictor_validation$pass),
      all(predictor_tables$audit$exact_V_order &
            predictor_tables$audit$matched_M0_M1_full_rows &
            (is.na(predictor_tables$audit$baseline_reproduces_M_both) |
               predictor_tables$audit$baseline_reproduces_M_both)),
      all(predictor_tables$diagnostics$covariance_positive_definite &
            predictor_tables$diagnostics$covariance_identifiable &
            predictor_tables$diagnostics$finite_coefficients &
            predictor_tables$diagnostics$finite_fixed_covariance),
      isTRUE(family_model_summary$V_valid[[1]]) &&
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
      paste("minimum eigenvalue", report_number(scope$V_validation$minimum_eigenvalue)),
      paste(sum(consistency_checks$pass), "of", nrow(consistency_checks), "checks passed"),
      paste(sum(predictor_validation$pass), "of", nrow(predictor_validation), "checks passed"),
      paste(nrow(predictor_tables$audit), "predictor sequences"),
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

format_classification_summary <- function(summary) {
  paste0(
    summary$Resolution, ": eligible=", summary$`N eligible units`,
    "; winners=", summary$`Winner N`, " (", report_number(summary$`Winner %`), "%)",
    "; losers=", summary$`Loser N`, " (", report_number(summary$`Loser %`), "%)",
    "; context-dependent=", summary$`Context-dependent N`, " (",
    report_number(summary$`Context-dependent %`), "%)"
  )
}

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
      " | ", row$`ML M0 vs M1 result`,
      " | ", row$`ML M1 vs full/context result`,
      " | Delta AIC=", report_number(row$`Delta AIC`),
      " | Delta AICc=", report_number(row$`Delta AICc`)
    )
  }, character(1))
}

build_current_downstream_model_report <- function(
    config, scope, climate_scaling, mboth_bundle, species_summary,
    family_summary, trait_scaling, predictor_comparisons,
    heterogeneity_synthesis, consistency_checks, release_checks) {
  fixed_lines <- vapply(seq_len(nrow(mboth_bundle$fixed)), function(index) {
    row <- mboth_bundle$fixed[index, ]
    paste0(
      row$term, ": estimate=", report_number(row$estimate),
      "; SE=", report_number(row$model_SE),
      "; 95% CI=[", report_number(row$model_CI_low), ", ",
      report_number(row$model_CI_high), "]; p=", report_number(row$model_p_value)
    )
  }, character(1))
  cr1_lines <- vapply(seq_len(nrow(mboth_bundle$CR1)), function(index) {
    row <- mboth_bundle$CR1[index, ]
    paste0(
      row$term, ": estimate=", report_number(row$estimate),
      "; CR1 SE=", report_number(row$CR1_SE),
      "; CR1 95% CI=[", report_number(row$CR1_CI_low), ", ",
      report_number(row$CR1_CI_high), "]; CR1 p=",
      report_number(row$CR1_p_value), "; df=",
      report_number(row$CR1_denominator_df)
    )
  }, character(1))
  random_row <- mboth_bundle$random[1, ]
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
  profile_lines <- vapply(seq_len(nrow(mboth_bundle$variance_profiles)), function(index) {
    row <- mboth_bundle$variance_profiles[index, ]
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
    "CURRENT CANONICAL DOWNSTREAM MODEL RESULTS",
    "==========================================",
    "",
    "Generated from current canonical targets; frozen Step 22/transition reports are not inputs.",
    report_section("Dataset and exclusions", c(
      paste("Effects:", scope$effects), paste("Species:", scope$species),
      paste("Experiments:", scope$experiments),
      paste("Excluded downstream taxa:", paste(config$exclusions, collapse = "; ")),
      "The excluded taxa remain in upstream winner/loser calculations."
    )),
    report_section("Climate scaling", scaling_lines),
    report_section("M_both fixed effects", fixed_lines),
    report_section("M_both experiment-clustered CR1 fixed effects", cr1_lines),
    report_section("M_both random-effect variances and correlations", random_lines),
    report_section("M_both variance-component profile intervals", profile_lines),
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
  consistency_lines <- vapply(seq_len(nrow(consistency_checks)), function(index) {
    paste0(
      consistency_checks$check_number[[index]], ". ",
      consistency_checks$check[[index]], ": ",
      if (isTRUE(consistency_checks$pass[[index]])) "PASS" else "FAIL",
      " — ", consistency_checks$detail[[index]]
    )
  }, character(1))
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
    report_section("Canonical consistency checks", consistency_lines),
    report_section("Active release checks", release_lines),
    "",
    paste("FINAL STATUS:", if (all(release_checks$pass)) "PASS" else "FAIL")
  )
}
