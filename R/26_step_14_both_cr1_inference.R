# Step 14: experiment-clustered CR1 inference for the saved M_both ----------
#
# This script does not fit or update M_both. It reads the already fitted model
# and replaces only the model-based coefficient covariance used for inference
# with metafor's experiment-clustered CR1 sandwich covariance (adjust = TRUE,
# clubSandwich = FALSE). Point estimates and the fitted response surface are
# therefore unchanged.

required_packages <- c("dplyr", "metafor", "readr", "tibble")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "),
       call. = FALSE)
}

model_file <- file.path(
  "output", "step_by_step_modeling", "species_both",
  "auc_change_species_both_model.rds"
)
output_dir <- dirname(model_file)
if (!file.exists(model_file)) {
  stop("Saved M_both is required: ", model_file, call. = FALSE)
}

M_both <- readRDS(model_file)
if (!inherits(M_both, "rma.mv")) {
  stop("The saved object is not an rma.mv model.", call. = FALSE)
}
if (is.null(M_both$data) || !"experiment_id" %in% names(M_both$data)) {
  stop("The fitted model does not retain experiment_id in fitted-row order.",
       call. = FALSE)
}
if (length(M_both$data$experiment_id) != M_both$k) {
  stop("experiment_id length does not match the fitted model rows.",
       call. = FALSE)
}
if (anyNA(M_both$data$experiment_id)) {
  stop("experiment_id contains missing values.", call. = FALSE)
}

# CR1 inference only: robust.rma does not refit M_both or alter its estimates.
captured_warnings <- character()
M_both_cr1 <- withCallingHandlers(
  metafor::robust(
    M_both,
    cluster = M_both$data$experiment_id,
    adjust = TRUE,
    clubSandwich = FALSE
  ),
  warning = function(w) {
    captured_warnings <<- c(captured_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

term <- names(stats::coef(M_both))
ordinary_estimate <- as.numeric(stats::coef(M_both))
robust_estimate <- as.numeric(M_both_cr1$beta)
if (!isTRUE(all.equal(ordinary_estimate, robust_estimate,
                      tolerance = sqrt(.Machine$double.eps)))) {
  stop("CR1 unexpectedly changed the fitted fixed-effect estimates.",
       call. = FALSE)
}

ordinary_se <- as.numeric(M_both$se)
ordinary_ci_low <- as.numeric(M_both$ci.lb)
ordinary_ci_high <- as.numeric(M_both$ci.ub)
ordinary_p <- as.numeric(M_both$pval)
cr1_se <- as.numeric(M_both_cr1$se)
cr1_ci_low <- as.numeric(M_both_cr1$ci.lb)
cr1_ci_high <- as.numeric(M_both_cr1$ci.ub)
cr1_p <- as.numeric(M_both_cr1$pval)
cr1_df <- rep(as.numeric(M_both_cr1$ddf), length.out = length(term))

coefficient_table <- tibble::tibble(
  term = term,
  highlight_T_by_P_by_W = term == "T_z:P_z:W_z",
  estimate = ordinary_estimate,
  model_standard_error = ordinary_se,
  model_conf_low = ordinary_ci_low,
  model_conf_high = ordinary_ci_high,
  model_p_value = ordinary_p,
  CR1_standard_error = cr1_se,
  CR1_conf_low = cr1_ci_low,
  CR1_conf_high = cr1_ci_high,
  CR1_p_value = cr1_p,
  CR1_denominator_df = cr1_df
)

three_way <- coefficient_table |>
  dplyr::filter(.data$highlight_T_by_P_by_W)
if (nrow(three_way) != 1L) {
  stop("Could not identify exactly one T_z:P_z:W_z coefficient.",
       call. = FALSE)
}

summary_table <- tibble::tibble(
  model = "M_both",
  response = "Delta AUC (warmed minus origin)",
  model_refitted = FALSE,
  robust_method = "experiment-clustered CR1",
  robust_call = paste0(
    "metafor::robust(M_both, cluster = M_both$data$experiment_id, ",
    "adjust = TRUE, clubSandwich = FALSE)"
  ),
  n_effects = M_both$k,
  n_fixed_coefficients = M_both$p,
  n_experiment_clusters = dplyr::n_distinct(M_both$data$experiment_id),
  CR1_denominator_df = unique(cr1_df)[[1]],
  warnings = paste(unique(captured_warnings), collapse = " | "),
  fixed_estimates_identical = isTRUE(all.equal(
    ordinary_estimate, robust_estimate,
    tolerance = sqrt(.Machine$double.eps)
  )),
  three_way_term = three_way$term[[1]],
  three_way_estimate = three_way$estimate[[1]],
  three_way_model_standard_error = three_way$model_standard_error[[1]],
  three_way_model_conf_low = three_way$model_conf_low[[1]],
  three_way_model_conf_high = three_way$model_conf_high[[1]],
  three_way_model_p_value = three_way$model_p_value[[1]],
  three_way_CR1_standard_error = three_way$CR1_standard_error[[1]],
  three_way_CR1_conf_low = three_way$CR1_conf_low[[1]],
  three_way_CR1_conf_high = three_way$CR1_conf_high[[1]],
  three_way_CR1_p_value = three_way$CR1_p_value[[1]],
  point_surface_changed = FALSE,
  inference_changed_at_alpha_0_05 =
    three_way$model_conf_low[[1]] > 0 &&
      three_way$CR1_conf_low[[1]] <= 0 &&
      three_way$CR1_conf_high[[1]] >= 0,
  interpretation = paste(
    "CR1 changes uncertainty, not fitted values: the population-level fixed",
    "T-by-P-by-W point surface is exactly unchanged. The positive three-way",
    "coefficient has a model-based 95% CI excluding zero but an",
    "experiment-clustered CR1 95% CI spanning zero, so evidence for this",
    "individual interaction is weakened at alpha = 0.05."
  )
)

readr::write_csv(
  coefficient_table,
  file.path(output_dir, "both_fixed_coefficients_model_vs_experiment_CR1.csv")
)
readr::write_csv(
  summary_table,
  file.path(output_dir, "both_experiment_CR1_summary.csv")
)
saveRDS(
  M_both_cr1,
  file.path(output_dir, "auc_change_species_both_experiment_CR1.rds")
)

cat("\nSaved-M_both experiment-clustered CR1 inference\n")
cat("Model refitted: FALSE\n")
cat("Experiment clusters: ", summary_table$n_experiment_clusters, "\n",
    sep = "")
cat("CR1 denominator df: ", summary_table$CR1_denominator_df, "\n", sep = "")
print(coefficient_table, n = Inf, width = Inf)
cat("\nT_z:P_z:W_z: estimate = ",
    sprintf("%.8f", three_way$estimate),
    "; model SE = ", sprintf("%.8f", three_way$model_standard_error),
    ", 95% CI [", sprintf("%.8f", three_way$model_conf_low), ", ",
    sprintf("%.8f", three_way$model_conf_high), "]",
    "; CR1 SE = ", sprintf("%.8f", three_way$CR1_standard_error),
    ", 95% CI [", sprintf("%.8f", three_way$CR1_conf_low), ", ",
    sprintf("%.8f", three_way$CR1_conf_high), "]\n", sep = "")
cat(summary_table$interpretation, "\n")
cat("Outputs: ", normalizePath(output_dir), "\n", sep = "")
