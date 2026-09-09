# Step 9: M_both -------------------------------------------------------------
#
# Response: Delta AUC = warmed AUC minus origin AUC, maximum-duration horizon.
# Combines the previously fitted species warming-slope and background-context
# blocks. No experiment effects, traits, duration, colonizers, or other terms.

source("config/analysis.R")
source("R/pipeline/utils.R")
source("R/pipeline/frequentist_hierarchical_auc.R")

required_packages <- c("dplyr", "Matrix", "metafor", "readr", "tibble")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "),
       call. = FALSE)
}

horizon <- "max"
output_dir <- file.path("output", "step_by_step_modeling", "species_both")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

warming_model_file <- file.path(
  "output", "step_by_step_modeling", "auc_intercept_species",
  "auc_change_climate_species_slope_model.rds"
)
context_model_file <- file.path(
  "output", "step_by_step_modeling", "species_context",
  "auc_change_species_context_model.rds"
)
if (!file.exists(warming_model_file) || !file.exists(context_model_file)) {
  stop("Saved M_warming and M_context models are required.", call. = FALSE)
}
M_warming <- readRDS(warming_model_file)
M_context <- readRDS(context_model_file)

scaled_auc <- scale_resident_auc_data(
  read_resident_auc_data(here::here(), horizons = horizon)
)
auc_data <- scaled_auc$data |>
  dplyr::mutate(
    species_id = as.character(.data$species_id),
    TP_z = .data$T_z * .data$P_z
  )
auc_by_plot <- readRDS(file.path("output", "intermediate", "auc_by_plot.rds"))
sampling <- reconstruct_auc_covariance(auc_data, auc_by_plot)
if (!isTRUE(sampling$audit$full_covariance_valid[[1]])) {
  stop("The full known Delta-AUC sampling covariance did not validate.",
       call. = FALSE)
}
V_auc <- sampling$V

# Exact requested model. The three entries are separate species-level blocks:
# independent intercept, independent warming slope, and correlated T/P/T*P.
both_model_file <- file.path(output_dir, "auc_change_species_both_model.rds")
fit_capture_file <- file.path(output_dir, "fit_capture.rds")
if (file.exists(both_model_file) && file.exists(fit_capture_file)) {
  M_both <- readRDS(both_model_file)
  saved_capture <- readRDS(fit_capture_file)
  fit_attempt <- list(
    fit = M_both, warnings = saved_capture$warnings, error = saved_capture$error
  )
} else {
  fit_attempt <- auc_capture(
    metafor::rma.mv(
      yi = auc_effect,
      V = V_auc,
      mods = ~ T_z * P_z * W_z,
      random = list(
        ~ 1 | species_id,
        ~ 0 + W_z | species_id,
        ~ 0 + T_z + P_z + TP_z | species_id
      ),
      struct = "GEN",
      data = auc_data,
      method = "REML",
      sparse = TRUE
    )
  )
  if (is.null(fit_attempt$fit)) {
    stop("M_both failed: ", fit_attempt$error, call. = FALSE)
  }
  M_both <- fit_attempt$fit
  invisible(write_rds_atomic(M_both, both_model_file))
  invisible(write_rds_atomic(
    list(warnings = fit_attempt$warnings, error = fit_attempt$error),
    fit_capture_file
  ))
}

# Verify that metafor instantiated the intended separate blocks.
if (!identical(colnames(M_both$G), "W_z") ||
    !identical(colnames(M_both$H), c("T_z", "P_z", "TP_z"))) {
  stop(
    "Unexpected random blocks: G=", paste(colnames(M_both$G), collapse = ","),
    "; H=", paste(colnames(M_both$H), collapse = ","), call. = FALSE
  )
}
optimizer_code <- auc_value_or(M_both$opt.res$convergence, NA_integer_)
fixed_covariance <- stats::vcov(M_both)
converged <- !nzchar(fit_attempt$error) && !length(fit_attempt$warnings) &&
  isTRUE(optimizer_code == 0L) &&
  all(is.finite(stats::coef(M_both))) &&
  all(is.finite(fixed_covariance))

variance_values <- c(
  species_intercept = as.numeric(M_both$sigma2[[1]]),
  species_W = as.numeric(M_both$tau2[[1]]),
  species_T = as.numeric(M_both$gamma2[[1]]),
  species_P = as.numeric(M_both$gamma2[[2]]),
  species_TP = as.numeric(M_both$gamma2[[3]])
)
boundary_tolerance <- 1e-8 * max(1, stats::var(auc_data$auc_effect))

# Profile the four requested slope variances. Cache only an exact match to the
# fitted variance/correlation signature.
profile_specifications <- tibble::tribble(
  ~component, ~argument, ~index,
  "species_W", "tau2", 1L,
  "species_T", "gamma2", 1L,
  "species_P", "gamma2", 2L,
  "species_TP", "gamma2", 3L
)
profile_signature <- paste(
  format(
    c(variance_values, M_both$rho, M_both$phi),
    digits = 16, scientific = TRUE
  ),
  collapse = "|"
)
profile_cache_file <- file.path(output_dir, "variance_profile_intervals.rds")
profile_cache <- if (file.exists(profile_cache_file)) {
  readRDS(profile_cache_file)
} else {
  NULL
}
profile_cache_valid <- !is.null(profile_cache) &&
  identical(profile_cache$signature, profile_signature) &&
  nrow(profile_cache$table) == nrow(profile_specifications)

if (profile_cache_valid) {
  variance_profiles <- profile_cache$table
} else {
  profile_one <- function(i) {
    specification <- profile_specifications[i, ]
    arguments <- list(object = M_both)
    arguments[[specification$argument[[1]]]] <- specification$index[[1]]
    captured <- auc_capture(do.call(stats::confint, arguments))
    estimate <- variance_values[[specification$component[[1]]]]
    if (is.null(captured$fit)) {
      return(tibble::tibble(
        component = specification$component[[1]], estimate = estimate,
        conf_low = NA_real_, conf_high = NA_real_,
        warnings = paste(captured$warnings, collapse = " | "),
        error = captured$error
      ))
    }
    interval_name <- if (specification$argument[[1]] == "tau2") {
      "tau^2"
    } else {
      paste0("gamma^2.", specification$index[[1]])
    }
    interval <- captured$fit$random[interval_name, ]
    tibble::tibble(
      component = specification$component[[1]],
      estimate = as.numeric(interval[["estimate"]]),
      conf_low = as.numeric(interval[["ci.lb"]]),
      conf_high = as.numeric(interval[["ci.ub"]]),
      warnings = paste(captured$warnings, collapse = " | "), error = ""
    )
  }
  profile_cores <- suppressWarnings(as.integer(
    Sys.getenv("TRANSPLANT_STEP9_PROFILE_CORES", "4")
  ))
  detected_cores <- suppressWarnings(parallel::detectCores())
  if (!length(detected_cores) || !is.finite(detected_cores)) {
    detected_cores <- profile_cores
  }
  profile_cores <- max(1L, min(profile_cores, detected_cores))
  profile_results <- if (.Platform$OS.type != "windows" && profile_cores > 1L) {
    parallel::mclapply(
      seq_len(nrow(profile_specifications)), profile_one,
      mc.cores = profile_cores
    )
  } else {
    lapply(seq_len(nrow(profile_specifications)), profile_one)
  }
  variance_profiles <- dplyr::bind_rows(profile_results)
  invisible(write_rds_atomic(
    list(signature = profile_signature, table = variance_profiles),
    profile_cache_file
  ))
}
variance_profiles <- variance_profiles |>
  dplyr::mutate(
    warnings = dplyr::if_else(
      .data$warnings == "Extra argument ('progbar') disregarded.",
      "", .data$warnings, missing = ""
    ),
    boundary_tolerance = boundary_tolerance,
    finite_profile_interval = is.finite(.data$conf_low) &
      is.finite(.data$conf_high),
    effectively_on_boundary = .data$estimate <= .data$boundary_tolerance |
      !.data$finite_profile_interval |
      .data$conf_low <= .data$boundary_tolerance
  )

variance_table <- tibble::tibble(
  component = names(variance_values),
  estimate = unname(variance_values)
) |>
  dplyr::left_join(
    dplyr::select(variance_profiles, -"estimate"), by = "component"
  ) |>
  dplyr::mutate(
    effectively_on_boundary = dplyr::if_else(
      .data$component == "species_intercept",
      .data$estimate <= .env$boundary_tolerance,
      .data$effectively_on_boundary
    )
  )

# Full block-diagonal covariance: intercept, W, and correlated T/P/T*P.
context_covariance <- as.matrix(M_both$H)
context_correlation <- stats::cov2cor(context_covariance)
full_covariance <- matrix(
  0, nrow = 5L, ncol = 5L,
  dimnames = list(
    c("intercept", "W_z", "T_z", "P_z", "TP_z"),
    c("intercept", "W_z", "T_z", "P_z", "TP_z")
  )
)
full_covariance["intercept", "intercept"] <- variance_values[["species_intercept"]]
full_covariance["W_z", "W_z"] <- variance_values[["species_W"]]
full_covariance[c("T_z", "P_z", "TP_z"),
                c("T_z", "P_z", "TP_z")] <- context_covariance
covariance_eigenvalues <- eigen(
  full_covariance, symmetric = TRUE, only.values = TRUE
)$values
eigenvalue_table <- tibble::tibble(
  eigenvalue_index = seq_along(covariance_eigenvalues),
  eigenvalue = covariance_eigenvalues
)
correlation_table <- as.data.frame(
  as.table(context_correlation), responseName = "correlation"
) |>
  tibble::as_tibble() |>
  dplyr::rename(coefficient_1 = "Var1", coefficient_2 = "Var2") |>
  dplyr::filter(as.integer(.data$coefficient_1) <
                  as.integer(.data$coefficient_2))
max_absolute_context_correlation <- max(
  abs(context_correlation[upper.tri(context_correlation)])
)
covariance_positive_definite <- all(is.finite(covariance_eigenvalues)) &&
  min(covariance_eigenvalues) > boundary_tolerance
covariance_identifiable <- covariance_positive_definite &&
  all(variance_profiles$finite_profile_interval) &&
  !any(variance_profiles$effectively_on_boundary) &&
  max_absolute_context_correlation < 1 - 1e-6

fit_statistics_one <- function(model, model_name) {
  statistics <- metafor::fitstats(model)
  likelihood <- stats::logLik(model)
  estimation_method <- model$method
  tibble::tibble(
    model = model_name,
    logLik = as.numeric(likelihood),
    df = attr(likelihood, "df"),
    AIC = as.numeric(statistics["AIC", 1]),
    AICc = as.numeric(statistics["AICc", 1]),
    method = estimation_method
  )
}
model_comparison <- dplyr::bind_rows(
  fit_statistics_one(M_both, "M_both"),
  fit_statistics_one(M_context, "M_context"),
  fit_statistics_one(M_warming, "M_warming")
) |>
  dplyr::arrange(.data$AICc)

# Compare fixed T*P*W predictions at all 42 observed experiment climates. The
# response-scale heuristic calls a change material if it exceeds 10% of the
# observed Delta-AUC IQR or changes sign at more than 5% of experiment climates.
observed_climates <- auc_data |>
  dplyr::distinct(.data$experiment_id, .data$T_z, .data$P_z, .data$W_z) |>
  dplyr::mutate(TP_z = .data$T_z * .data$P_z)
fixed_prediction <- function(model) {
  design <- auc_design_matrix(
    ~ T_z * P_z * W_z, observed_climates, names(stats::coef(model))
  )
  as.numeric(design %*% stats::coef(model))
}
fixed_surface <- observed_climates |>
  dplyr::mutate(
    M_both = fixed_prediction(M_both),
    M_context = fixed_prediction(M_context),
    M_warming = fixed_prediction(M_warming)
  )
material_threshold <- 0.10 * stats::IQR(auc_data$auc_effect)
surface_comparison_one <- function(reference_name) {
  reference <- fixed_surface[[reference_name]]
  difference <- fixed_surface$M_both - reference
  sign_disagreement <- sign(fixed_surface$M_both) != sign(reference)
  tibble::tibble(
    comparison = paste0("M_both vs ", reference_name),
    n_observed_experiment_climates = nrow(fixed_surface),
    pearson_correlation = stats::cor(fixed_surface$M_both, reference),
    median_absolute_difference = stats::median(abs(difference)),
    root_mean_square_difference = sqrt(mean(difference^2)),
    maximum_absolute_difference = max(abs(difference)),
    n_sign_disagreements = sum(sign_disagreement),
    percent_sign_disagreement = 100 * mean(sign_disagreement),
    material_difference_threshold = material_threshold,
    materially_changed_descriptive =
      max(abs(difference)) > material_threshold |
      mean(sign_disagreement) > 0.05
  )
}
fixed_surface_comparison <- dplyr::bind_rows(
  surface_comparison_one("M_context"),
  surface_comparison_one("M_warming")
)
fixed_coefficient_comparison <- tibble::tibble(
  term = names(stats::coef(M_both)),
  M_both = as.numeric(stats::coef(M_both)),
  M_context = as.numeric(stats::coef(M_context)),
  M_warming = as.numeric(stats::coef(M_warming)),
  both_minus_context = .data$M_both - .data$M_context,
  both_minus_warming = .data$M_both - .data$M_warming
)

model_summary <- tibble::tibble(
  converged_without_warning = converged,
  optimizer_code = optimizer_code,
  warnings = paste(fit_attempt$warnings, collapse = " | "),
  species_intercept_variance = variance_values[["species_intercept"]],
  species_W_variance = variance_values[["species_W"]],
  species_T_variance = variance_values[["species_T"]],
  species_P_variance = variance_values[["species_P"]],
  species_TP_variance = variance_values[["species_TP"]],
  max_absolute_context_correlation = max_absolute_context_correlation,
  minimum_covariance_eigenvalue = min(covariance_eigenvalues),
  covariance_positive_definite = covariance_positive_definite,
  covariance_identifiable = covariance_identifiable,
  any_profiled_variance_on_boundary =
    any(variance_profiles$effectively_on_boundary),
  sampling_covariance_source = sampling$selected_v_type
)

invisible(write_csv_atomic(model_summary, file.path(output_dir, "model_summary.csv")))
invisible(write_csv_atomic(variance_table, file.path(output_dir, "variance_components.csv")))
invisible(write_csv_atomic(variance_profiles, file.path(output_dir, "variance_profile_intervals.csv")))
invisible(write_csv_atomic(correlation_table, file.path(output_dir, "context_correlations.csv")))
invisible(write_csv_atomic(eigenvalue_table, file.path(output_dir, "covariance_eigenvalues.csv")))
invisible(write_csv_atomic(model_comparison, file.path(output_dir, "model_fit_comparison.csv")))
invisible(write_csv_atomic(fixed_coefficient_comparison,
                           file.path(output_dir, "fixed_coefficient_comparison.csv")))
invisible(write_csv_atomic(fixed_surface,
                           file.path(output_dir, "fixed_surface_observed_climates.csv")))
invisible(write_csv_atomic(fixed_surface_comparison,
                           file.path(output_dir, "fixed_surface_comparison.csv")))

cat("M_both completed.\n")
cat("Converged without warnings: ", converged, "\n", sep = "")
cat("Optimizer code: ", optimizer_code, "\n", sep = "")
cat("Warnings: ", ifelse(length(fit_attempt$warnings),
                           paste(fit_attempt$warnings, collapse = " | "), "none"),
    "\n", sep = "")
cat("\nVariance components:\n")
print(variance_table)
cat("\nT/P/T*P correlations:\n")
print(correlation_table)
cat("\nCovariance eigenvalues:\n")
print(eigenvalue_table)
cat("Positive definite: ", covariance_positive_definite, "\n", sep = "")
cat("Identifiable by profiles/eigenvalues/correlations: ",
    covariance_identifiable, "\n", sep = "")
cat("\nREML fit statistics:\n")
print(model_comparison)
cat("\nFixed-surface comparisons:\n")
print(fixed_surface_comparison)
cat("Outputs: ", normalizePath(output_dir), "\n", sep = "")
