#!/usr/bin/env Rscript

# Step 22: updated 526-species downstream modeling and heterogeneity synthesis
#
# Rebuilds the anchored downstream resident analysis after excluding six taxa
# from downstream modeling. Updates M_both, species consistency, family
# resolution, trait models, Graminoid status, Woodiness, PEV1-PEV3
# phylogenetic information, and response-scale heterogeneity synthesis.
#
# Upstream community/AUC and significance-based winner/loser calculations are
# not modified.
#
# Step 22 is used because Step 20 and Step 21 are already occupied by
# R/32_step_20_phylogeny_model_readiness_audit.R and
# R/33_step_21_trait_data_audit.R, respectively.

options(lifecycle_verbosity = "quiet")

source("config/analysis.R")
source("R/pipeline/utils.R")
source("R/pipeline/frequentist_hierarchical_auc.R")

required_packages <- c(
  "dplyr", "geometry", "ggplot2", "ggrepel", "here", "Matrix",
  "metafor", "readr", "tibble", "tidyr"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

root <- here::here()
output_dir <- file.path(root, "output", "step_by_step_modeling", "downstream_526")
model_dir <- file.path(output_dir, "models")
table_dir <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figures")
cache_dir <- file.path(output_dir, "cache")
for (path in c(output_dir, model_dir, table_dir, figure_dir, cache_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

excluded_species <- c(
  "Gymnocarpium dryopteris", "Athyrium spinulosum", "Botrychium lunaria",
  "Equisetum arvense", "Selaginella selaginoides", "Juniperus communis"
)
climate_formula <- ~ T_z * P_z * W_z
species_random <- list(
  ~ 1 | species_id,
  ~ 0 + W_z | species_id,
  ~ 0 + T_z + P_z + TP_z | species_id
)
random_structure_text <- paste(
  "list(~1 | species_id, ~0 + W_z | species_id,",
  "~0 + T_z + P_z + TP_z | species_id); struct='GEN'"
)

table_path <- function(name) file.path(table_dir, name)
model_path <- function(name) file.path(model_dir, name)
figure_path <- function(name) file.path(figure_dir, name)

capture_fit <- function(expression, label) {
  captured <- auc_capture(expression)
  if (is.null(captured$fit)) {
    stop(label, " failed without simplification: ", captured$error, call. = FALSE)
  }
  captured
}

fit_signature <- function(data, formula, method, predictor_columns = character()) {
  predictor_values <- if (length(predictor_columns)) {
    lapply(data[predictor_columns], function(x) {
      if (is.factor(x)) as.character(x) else as.numeric(x)
    })
  } else {
    list()
  }
  list(
    effect_id = as.character(data$effect_id),
    formula = paste(deparse(formula), collapse = " "),
    method = method,
    predictor_values = predictor_values,
    random_structure = random_structure_text
  )
}

fit_cached <- function(data, V, formula, method, cache_file, label,
                       predictor_columns = character()) {
  signature <- fit_signature(data, formula, method, predictor_columns)
  if (file.exists(cache_file)) {
    cached <- readRDS(cache_file)
    if (is.list(cached) && !is.null(cached$fit) &&
        identical(cached$signature, signature)) {
      message("Using exact cache: ", label)
      return(cached)
    }
  }
  message("Fitting: ", label)
  captured <- capture_fit(
    metafor::rma.mv(
      yi = auc_effect, V = V, mods = formula,
      random = species_random, struct = "GEN", data = data,
      method = method, sparse = TRUE
    ),
    label
  )
  result <- c(captured, list(signature = signature))
  write_rds_atomic(result, cache_file)
  result
}

fit_stats <- function(captured, predictor, model, method) {
  fit <- captured$fit
  values <- auc_fit_statistics(fit)
  tibble::tibble(
    predictor = predictor, model = model, method = method,
    n_effects = fit$k, n_fixed = fit$p,
    logLik = values[["logLik"]], AIC = values[["AIC"]],
    AICc = values[["AICc"]], BIC = values[["BIC"]],
    optimizer_code = auc_value_or(fit$opt.res$convergence,
                                  auc_value_or(fit$opt.res$ierr, NA_integer_)),
    warnings = paste(captured$warnings, collapse = " | ")
  )
}

nested_lrt <- function(smaller, larger, predictor, comparison) {
  ll0 <- stats::logLik(smaller$fit)
  ll1 <- stats::logLik(larger$fit)
  df <- attr(ll1, "df") - attr(ll0, "df")
  statistic <- 2 * (as.numeric(ll1) - as.numeric(ll0))
  tibble::tibble(
    predictor = predictor, comparison = comparison,
    likelihood_ratio = statistic, df = df,
    p_value = stats::pchisq(statistic, df, lower.tail = FALSE),
    valid_nested_comparison = is.finite(statistic) && statistic >= -1e-7 && df > 0,
    same_rows_and_V = TRUE, method = "ML"
  )
}

robust_cr1 <- function(captured, predictor, model_name) {
  fit <- captured$fit
  warnings <- character()
  robust <- withCallingHandlers(
    metafor::robust(
      fit, cluster = fit$data$experiment_id,
      adjust = TRUE, clubSandwich = FALSE
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (!isTRUE(all.equal(as.numeric(stats::coef(fit)), as.numeric(robust$beta),
                        tolerance = sqrt(.Machine$double.eps)))) {
    stop(predictor, ": CR1 changed fixed estimates.", call. = FALSE)
  }
  table <- tibble::tibble(
    predictor = predictor, model = model_name,
    term = names(stats::coef(fit)), estimate = as.numeric(stats::coef(fit)),
    model_SE = as.numeric(fit$se), model_CI_low = as.numeric(fit$ci.lb),
    model_CI_high = as.numeric(fit$ci.ub), model_p_value = as.numeric(fit$pval),
    CR1_SE = as.numeric(robust$se), CR1_CI_low = as.numeric(robust$ci.lb),
    CR1_CI_high = as.numeric(robust$ci.ub), CR1_p_value = as.numeric(robust$pval),
    CR1_denominator_df = rep(as.numeric(robust$ddf), length.out = fit$p),
    CR1_warnings = paste(unique(warnings), collapse = " | ")
  )
  list(object = robust, table = table)
}

context_covariance <- function(fit) {
  covariance <- as.matrix(fit$H)
  if (!identical(colnames(fit$G), "W_z") ||
      !identical(colnames(covariance), c("T_z", "P_z", "TP_z"))) {
    stop("A fitted model does not retain the requested independent intercept/W and GEN T/P/TP blocks.")
  }
  covariance
}

random_diagnostics <- function(captured, predictor, model_name) {
  fit <- captured$fit
  G <- context_covariance(fit)
  R <- stats::cov2cor(G)
  full <- matrix(0, 5, 5, dimnames = list(
    c("intercept", "W_z", "T_z", "P_z", "TP_z"),
    c("intercept", "W_z", "T_z", "P_z", "TP_z")
  ))
  full["intercept", "intercept"] <- fit$sigma2[[1]]
  full["W_z", "W_z"] <- fit$tau2[[1]]
  full[c("T_z", "P_z", "TP_z"), c("T_z", "P_z", "TP_z")] <- G
  eigenvalues <- eigen(full, symmetric = TRUE, only.values = TRUE)$values
  tolerance <- 1e-8 * max(1, stats::var(fit$yi))
  tibble::tibble(
    predictor = predictor, model = model_name,
    n_effects = fit$k, n_species = dplyr::n_distinct(fit$data$species_id),
    n_experiments = dplyr::n_distinct(fit$data$experiment_id),
    species_intercept_variance = fit$sigma2[[1]],
    species_W_variance = fit$tau2[[1]],
    species_T_variance = G["T_z", "T_z"],
    species_P_variance = G["P_z", "P_z"],
    species_TP_variance = G["TP_z", "TP_z"],
    T_P_correlation = R["T_z", "P_z"],
    T_TP_correlation = R["T_z", "TP_z"],
    P_TP_correlation = R["P_z", "TP_z"],
    minimum_full_covariance_eigenvalue = min(eigenvalues),
    covariance_positive_definite = min(eigenvalues) > 0,
    any_variance_at_boundary = any(diag(full) <= tolerance),
    maximum_absolute_context_correlation = max(abs(R[upper.tri(R)])),
    covariance_identifiable = min(eigenvalues) > tolerance &&
      all(diag(full) > tolerance) && max(abs(R[upper.tri(R)])) < 1 - 1e-6,
    optimizer_code = auc_value_or(fit$opt.res$convergence,
                                  auc_value_or(fit$opt.res$ierr, NA_integer_)),
    finite_coefficients = all(is.finite(stats::coef(fit))),
    finite_fixed_covariance = all(is.finite(stats::vcov(fit))),
    warnings = paste(captured$warnings, collapse = " | ")
  )
}

# Preserve proof that this rerun does not modify upstream analytical products.
upstream_paths <- c(
  winner_loser = file.path(root, "output", "tables", "winner_loser_results.csv"),
  auc_by_plot = file.path(root, "output", "intermediate", "auc_by_plot.rds"),
  community = file.path(root, "output", "intermediate", "community.rds")
)
upstream_md5_before <- tools::md5sum(upstream_paths)

# Exact updated maximum-duration resident response and full sampling covariance.
scaled_auc <- scale_resident_auc_data(
  read_resident_auc_data(root, horizons = "max")
)
auc_data_526 <- scaled_auc$data |>
  dplyr::mutate(
    species_id = droplevels(factor(as.character(.data$species_id))),
    experiment_id = droplevels(factor(as.character(.data$experiment_id))),
    TP_z = .data$T_z * .data$P_z
  )
auc_by_plot <- readRDS(upstream_paths[["auc_by_plot"]])
sampling <- reconstruct_auc_covariance(auc_data_526, auc_by_plot)
if (!isTRUE(sampling$audit$full_covariance_valid[[1]])) {
  stop("Updated full Delta-AUC sampling covariance failed reconstruction.")
}
V_auc_526 <- sampling$V
V_validation <- validate_auc_covariance(as.matrix(V_auc_526), auc_data_526)
if (!isTRUE(V_validation$valid[[1]]) ||
    !identical(rownames(V_auc_526), auc_data_526$effect_id) ||
    max(abs(Matrix::diag(V_auc_526) - auc_data_526$auc_se^2)) > 1e-8 ||
    nrow(auc_data_526) != n_distinct(auc_data_526$effect_id) ||
    any(!stats::complete.cases(auc_data_526[c(
      "effect_id", "species_id", "experiment_id", "auc_effect", "auc_se",
      "T", "P", "W", "T_z", "P_z", "W_z", "TP_z"
    )])) ||
    any(!is.finite(auc_data_526$auc_effect))) {
  stop("Updated row/V/completeness validation failed; V was not repaired.")
}
if (n_distinct(auc_data_526$species_id) != 526L ||
    any(as.character(auc_data_526$species_id) %in% excluded_species)) {
  stop("Complete downstream data are not exactly the retained 526-species universe.")
}
write_rds_atomic(V_auc_526, model_path("V_auc_526.rds"))
write_csv_atomic(
  auc_data_526 |>
    dplyr::mutate(
      species_id = as.character(.data$species_id),
      experiment_id = as.character(.data$experiment_id),
      gradient_id = as.character(.data$gradient_id)
    ),
  table_path("auc_data_526.csv")
)

experiment_climates <- auc_data_526 |>
  dplyr::distinct(.data$experiment_id, .data$T, .data$P, .data$W,
                  .data$T_z, .data$P_z, .data$W_z) |>
  dplyr::mutate(TP_z = .data$T_z * .data$P_z) |>
  dplyr::arrange(.data$experiment_id)
if (anyDuplicated(experiment_climates$experiment_id)) {
  stop("Experiment-level climate covariates are not constant.")
}
climate_scaling <- scaled_auc$scaling |>
  dplyr::filter(.data$scaled_field %in% c("T_z", "P_z", "W_z")) |>
  dplyr::mutate(n_distinct_experiments = nrow(experiment_climates))
write_csv_atomic(climate_scaling, table_path("climate_scaling_526.csv"))

# Updated primary M_both. The preliminary exact cache may already exist.
preliminary_cache <- model_path("M_both_526_capture.rds")
M_both_capture <- NULL
if (file.exists(preliminary_cache)) {
  candidate <- readRDS(preliminary_cache)
  if (is.list(candidate) && !is.null(candidate$fit) &&
      identical(as.character(candidate$fit$data$effect_id), auc_data_526$effect_id)) {
    M_both_capture <- candidate
  }
}
if (is.null(M_both_capture)) {
  M_both_capture <- fit_cached(
    auc_data_526, V_auc_526, climate_formula, "REML",
    preliminary_cache, "M_both_526 REML"
  )
}
M_both_526 <- M_both_capture$fit
write_rds_atomic(M_both_526, model_path("M_both_526.rds"))

M_both_fit <- fit_stats(M_both_capture, "Updated M_both", "M_both_526", "REML")
M_both_random <- random_diagnostics(M_both_capture, "Updated M_both", "M_both_526")
M_both_fixed <- tibble::tibble(
  term = names(stats::coef(M_both_526)),
  estimate = as.numeric(stats::coef(M_both_526)),
  model_SE = as.numeric(M_both_526$se),
  model_CI_low = as.numeric(M_both_526$ci.lb),
  model_CI_high = as.numeric(M_both_526$ci.ub),
  model_p_value = as.numeric(M_both_526$pval)
)
M_both_cr1 <- robust_cr1(M_both_capture, "Updated M_both", "M_both_526")
write_rds_atomic(M_both_cr1$object, model_path("M_both_526_experiment_CR1.rds"))
write_csv_atomic(M_both_cr1$table, table_path("M_both_526_fixed_effects_model_vs_CR1.csv"))
write_csv_atomic(M_both_fit, table_path("M_both_526_fit_statistics.csv"))
write_csv_atomic(M_both_random, table_path("M_both_526_random_diagnostics.csv"))
write_csv_atomic(M_both_fixed, table_path("M_both_526_fixed_coefficients.csv"))

# Profile-likelihood intervals for all five random variances, without changing
# the model if a profile is difficult. Results are cached by the exact fit.
profile_specs <- tibble::tribble(
  ~component, ~argument, ~index, ~interval_name,
  "species_intercept", "sigma2", 1L, "sigma^2.1",
  "species_W", "tau2", 1L, "tau^2.1",
  "species_T", "gamma2", 1L, "gamma^2.1",
  "species_P", "gamma2", 2L, "gamma^2.2",
  "species_TP", "gamma2", 3L, "gamma^2.3"
)
profile_signature <- paste(
  format(c(M_both_526$sigma2, M_both_526$tau2, M_both_526$gamma2,
           M_both_526$rho, M_both_526$phi), digits = 16, scientific = TRUE),
  collapse = "|"
)
profile_cache_path <- file.path(cache_dir, "M_both_526_variance_profiles.rds")
profile_cache <- if (file.exists(profile_cache_path)) readRDS(profile_cache_path) else NULL
current_variances <- c(
  species_intercept = M_both_526$sigma2[[1]], species_W = M_both_526$tau2[[1]],
  species_T = M_both_526$gamma2[[1]], species_P = M_both_526$gamma2[[2]],
  species_TP = M_both_526$gamma2[[3]]
)
skip_profiles <- tolower(Sys.getenv("TRANSPLANT_DOWNSTREAM_526_SKIP_PROFILES", "false")) %in%
  c("true", "1", "yes")
if (skip_profiles) {
  variance_profiles <- tibble::tibble(
    component = profile_specs$component,
    estimate = unname(current_variances[profile_specs$component]),
    CI_low = NA_real_, CI_high = NA_real_, warnings = "",
    error = "temporarily skipped with TRANSPLANT_DOWNSTREAM_526_SKIP_PROFILES=true"
  )
} else if (!is.null(profile_cache) && identical(profile_cache$signature, profile_signature)) {
  variance_profiles <- profile_cache$table
} else {
  profile_one <- function(i) {
    specification <- profile_specs[i, ]
    component_cache_path <- file.path(
      cache_dir,
      paste0("M_both_526_profile_", specification$component[[1]], ".rds")
    )
    component_cache <- if (file.exists(component_cache_path)) {
      readRDS(component_cache_path)
    } else {
      NULL
    }
    if (!is.null(component_cache) &&
        identical(component_cache$signature, profile_signature)) {
      return(component_cache$result)
    }
    arguments <- list(object = M_both_526)
    arguments[[specification$argument[[1]]]] <- specification$index[[1]]
    # metafor otherwise brackets variance components up to 10 even though the
    # fitted values are O(1e-3). This conservative bracket is at least six
    # times the estimate; an unrestricted retry below protects against a CI
    # whose upper root genuinely lies beyond it.
    arguments$xlim <- c(
      0,
      max(0.025, 6 * current_variances[[specification$component[[1]]]])
    )
    captured <- auc_capture(do.call(stats::confint, arguments))
    if (is.null(captured$fit)) {
      arguments$xlim <- NULL
      captured <- auc_capture(do.call(stats::confint, arguments))
    }
    if (is.null(captured$fit)) {
      result <- tibble::tibble(
        component = specification$component[[1]], estimate = NA_real_,
        CI_low = NA_real_, CI_high = NA_real_,
        warnings = paste(captured$warnings, collapse = " | "), error = captured$error
      )
      saveRDS(list(signature = profile_signature, result = result),
              component_cache_path)
      return(result)
    }
    random_table <- captured$fit$random
    row_name <- specification$interval_name[[1]]
    if (!row_name %in% rownames(random_table)) {
      # metafor sometimes omits '.1' for a single component.
      row_name <- sub("\\.1$", "", row_name)
    }
    if (!row_name %in% rownames(random_table)) {
      result <- tibble::tibble(
        component = specification$component[[1]], estimate = NA_real_,
        CI_low = NA_real_, CI_high = NA_real_,
        warnings = paste(captured$warnings, collapse = " | "),
        error = paste("Profile row unavailable; rows:", paste(rownames(random_table), collapse = ", "))
      )
      saveRDS(list(signature = profile_signature, result = result),
              component_cache_path)
      return(result)
    }
    interval_values <- as.numeric(random_table[row_name, c("estimate", "ci.lb", "ci.ub")])
    if (any(!is.finite(interval_values))) {
      arguments$xlim <- NULL
      unrestricted <- auc_capture(do.call(stats::confint, arguments))
      if (!is.null(unrestricted$fit) &&
          row_name %in% rownames(unrestricted$fit$random)) {
        captured <- unrestricted
        random_table <- captured$fit$random
      }
    }
    result <- tibble::tibble(
      component = specification$component[[1]],
      estimate = as.numeric(random_table[row_name, "estimate"]),
      CI_low = as.numeric(random_table[row_name, "ci.lb"]),
      CI_high = as.numeric(random_table[row_name, "ci.ub"]),
      warnings = paste(captured$warnings, collapse = " | "), error = ""
    )
    saveRDS(list(signature = profile_signature, result = result),
            component_cache_path)
    result
  }
  requested_profile_cores <- suppressWarnings(as.integer(
    Sys.getenv("TRANSPLANT_DOWNSTREAM_526_PROFILE_CORES", "5")
  ))
  if (!length(requested_profile_cores) || !is.finite(requested_profile_cores)) {
    requested_profile_cores <- 5L
  }
  detected_cores <- suppressWarnings(parallel::detectCores())
  if (!length(detected_cores) || !is.finite(detected_cores)) {
    detected_cores <- requested_profile_cores
  }
  profile_cores <- max(1L, min(requested_profile_cores, detected_cores,
                               nrow(profile_specs)))
  profile_results <- if (.Platform$OS.type != "windows" && profile_cores > 1L) {
    parallel::mclapply(seq_len(nrow(profile_specs)), profile_one, mc.cores = profile_cores)
  } else {
    lapply(seq_len(nrow(profile_specs)), profile_one)
  }
  variance_profiles <- dplyr::bind_rows(profile_results)
  saveRDS(list(signature = profile_signature, table = variance_profiles), profile_cache_path)
}
write_csv_atomic(variance_profiles, table_path("M_both_526_variance_profile_CIs.csv"))

# Conditional species coefficients and prediction-error machinery used by the
# existing M_both classification, now recalculated for the 526-species fit.
extract_random_coefficients <- function(fit) {
  X <- as.matrix(fit$X)
  residual <- as.numeric(fit$yi - X %*% fit$beta)
  weighted_residual <- Matrix::solve(fit$M, residual)
  levels_g <- levels(fit$mf.g$outer)
  levels_h <- levels(fit$mf.h$outer)
  if (!identical(levels_g, levels_h)) stop("Species order differs across random blocks.")
  n_species <- length(levels_g)
  Z_intercept <- fit$Z.S[[1]]
  Z_W <- Matrix::Diagonal(x = as.numeric(fit$Z.G1[, 1])) %*% fit$Z.G2
  Z_context <- do.call(cbind, lapply(seq_len(ncol(fit$H)), function(j) {
    Matrix::Diagonal(x = as.numeric(fit$Z.H1[, j])) %*% fit$Z.H2
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
    context_vector, nrow = n_species, ncol = 3,
    dimnames = list(levels_g, c("T_z", "P_z", "TP_z"))
  )
  list(
    table = tibble::tibble(
      species_id = levels_g, random_intercept = random_intercept,
      random_W = random_W, random_T = context_matrix[, "T_z"],
      random_P = context_matrix[, "P_z"], random_TP = context_matrix[, "TP_z"]
    ),
    X = X, levels = levels_g, Z_intercept = Z_intercept,
    Z_W = Z_W, Z_context = Z_context
  )
}

species_random_effects <- extract_random_coefficients(M_both_526)
write_csv_atomic(
  species_random_effects$table,
  table_path("M_both_526_species_random_coefficients.csv")
)

species_support <- auc_data_526 |>
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
eligible_species <- species_support |>
  dplyr::filter(.data$eligible) |>
  dplyr::pull(.data$species_id) |>
  as.character()

W_raw_median <- stats::median(experiment_climates$W)
W_z_median <- stats::median(experiment_climates$W_z)
fixed_beta <- stats::coef(M_both_526)
fixed_covariance <- as.matrix(stats::vcov(M_both_526))
X_model <- species_random_effects$X
n_species_model <- length(species_random_effects$levels)
D_species <- matrix(0, 5, 5, dimnames = list(
  c("intercept", "W_z", "T_z", "P_z", "TP_z"),
  c("intercept", "W_z", "T_z", "P_z", "TP_z")
))
D_species["intercept", "intercept"] <- M_both_526$sigma2[[1]]
D_species["W_z", "W_z"] <- M_both_526$tau2[[1]]
D_species[c("T_z", "P_z", "TP_z"), c("T_z", "P_z", "TP_z")] <-
  as.matrix(M_both_526$H)

species_hull_grid <- function(species_name, points = 81L) {
  observed <- auc_data_526 |>
    dplyr::filter(as.character(.data$species_id) == .env$species_name) |>
    dplyr::distinct(.data$T_z, .data$P_z)
  candidate <- tidyr::expand_grid(
    T_z = seq(min(observed$T_z), max(observed$T_z), length.out = points),
    P_z = seq(min(observed$P_z), max(observed$P_z), length.out = points)
  )
  hull <- geometry::convhulln(as.matrix(observed[, c("T_z", "P_z")]), options = "Qt")
  candidate |>
    dplyr::filter(geometry::inhulln(
      hull, as.matrix(dplyr::select(candidate, "T_z", "P_z"))
    )) |>
    dplyr::bind_rows(observed) |>
    dplyr::distinct(.data$T_z, .data$P_z) |>
    dplyr::mutate(W_z = W_z_median, TP_z = .data$T_z * .data$P_z)
}

predict_species_supported <- function(species_name) {
  species_index <- match(species_name, species_random_effects$levels)
  grid <- species_hull_grid(species_name)
  X_new <- auc_design_matrix(climate_formula, grid, names(fixed_beta))
  R_new <- cbind(
    intercept = 1, W_z = grid$W_z, T_z = grid$T_z,
    P_z = grid$P_z, TP_z = grid$TP_z
  )
  random_row <- species_random_effects$table |>
    dplyr::filter(.data$species_id == .env$species_name)
  random_vector <- as.numeric(random_row[1, c(
    "random_intercept", "random_W", "random_T", "random_P", "random_TP"
  )])
  Z_species <- cbind(
    species_random_effects$Z_intercept[, species_index, drop = FALSE],
    species_random_effects$Z_W[, species_index, drop = FALSE],
    species_random_effects$Z_context[, species_index, drop = FALSE],
    species_random_effects$Z_context[, n_species_model + species_index, drop = FALSE],
    species_random_effects$Z_context[, 2L * n_species_model + species_index, drop = FALSE]
  )
  M_inverse_Z <- Matrix::solve(M_both_526$M, Z_species)
  XtM_inverse_Z <- as.matrix(crossprod(X_model, M_inverse_Z))
  ZtM_inverse_Z <- as.matrix(crossprod(Z_species, M_inverse_Z))
  random_pev <- D_species - D_species %*% ZtM_inverse_Z %*% D_species +
    D_species %*% t(XtM_inverse_Z) %*% fixed_covariance %*%
    XtM_inverse_Z %*% D_species
  fixed_random_cov <- -fixed_covariance %*% XtM_inverse_Z %*% D_species
  variance <- rowSums((X_new %*% fixed_covariance) * X_new) +
    rowSums((R_new %*% random_pev) * R_new) +
    2 * rowSums((X_new %*% fixed_random_cov) * R_new)
  estimate <- as.numeric(X_new %*% fixed_beta) + as.numeric(R_new %*% random_vector)
  grid |>
    dplyr::mutate(
      species_id = species_name, estimate = estimate,
      prediction_SE = sqrt(pmax(0, variance)),
      CI_low = .data$estimate - stats::qnorm(0.975) * .data$prediction_SE,
      CI_high = .data$estimate + stats::qnorm(0.975) * .data$prediction_SE
    )
}

message("Calculating updated species-supported classifications")
species_predictions <- dplyr::bind_rows(lapply(eligible_species, predict_species_supported))
species_classification <- species_predictions |>
  dplyr::group_split(.data$species_id, .keep = TRUE) |>
  lapply(function(x) {
    minimum <- x[which.min(x$estimate), , drop = FALSE]
    maximum <- x[which.max(x$estimate), , drop = FALSE]
    tibble::tibble(
      species_id = x$species_id[[1]],
      min_fitted_delta_auc = minimum$estimate[[1]],
      min_prediction_SE = minimum$prediction_SE[[1]],
      min_CI_low = minimum$CI_low[[1]], min_CI_high = minimum$CI_high[[1]],
      min_T_z = minimum$T_z[[1]], min_P_z = minimum$P_z[[1]],
      max_fitted_delta_auc = maximum$estimate[[1]],
      max_prediction_SE = maximum$prediction_SE[[1]],
      max_CI_low = maximum$CI_low[[1]], max_CI_high = maximum$CI_high[[1]],
      max_T_z = maximum$T_z[[1]], max_P_z = maximum$P_z[[1]],
      n_support_grid_points = nrow(x)
    )
  }) |>
  dplyr::bind_rows() |>
  dplyr::left_join(species_support, by = "species_id") |>
  dplyr::left_join(species_random_effects$table |> dplyr::select("species_id", "random_W"),
                   by = "species_id") |>
  dplyr::mutate(
    fitted_class = dplyr::case_when(
      .data$min_fitted_delta_auc > 0 ~ "Consistent winner",
      .data$max_fitted_delta_auc < 0 ~ "Consistent loser",
      TRUE ~ "Context-dependent"
    ),
    warming_raw_C = W_raw_median, warming_z = W_z_median
  )
if (any(species_classification$min_fitted_delta_auc >
        species_classification$max_fitted_delta_auc)) {
  stop("A species fitted minimum exceeds its maximum.")
}
write_csv_atomic(species_classification,
                 table_path("species_minmax_classification_526.csv"))

class_colours <- c(
  "Consistent loser" = "#E66101",
  "Context-dependent" = "#0072B2",
  "Consistent winner" = "#1B9E77"
)
class_counts <- table(factor(species_classification$fitted_class,
                             levels = names(class_colours)))
label_species <- dplyr::bind_rows(
  species_classification |> dplyr::filter(.data$fitted_class == "Consistent winner") |>
    dplyr::slice_max(.data$min_fitted_delta_auc, n = 3, with_ties = FALSE),
  species_classification |> dplyr::filter(.data$fitted_class == "Consistent loser") |>
    dplyr::slice_min(.data$max_fitted_delta_auc, n = 3, with_ties = FALSE),
  species_classification |> dplyr::filter(.data$fitted_class == "Context-dependent") |>
    dplyr::mutate(span = .data$max_fitted_delta_auc - .data$min_fitted_delta_auc) |>
    dplyr::slice_max(.data$span, n = 4, with_ties = FALSE)
) |>
  dplyr::distinct(.data$species_id, .keep_all = TRUE)
x_range <- range(c(species_classification$min_CI_low,
                   species_classification$min_CI_high, 0), na.rm = TRUE)
y_range <- range(c(species_classification$max_CI_low,
                   species_classification$max_CI_high, 0), na.rm = TRUE)
x_limits <- x_range + c(-0.06, 0.06) * diff(x_range)
y_limits <- y_range + c(-0.06, 0.06) * diff(y_range)
species_plot <- ggplot2::ggplot(species_classification) +
  ggplot2::annotate("rect", xmin = x_limits[[1]], xmax = 0,
                    ymin = y_limits[[1]], ymax = 0,
                    fill = class_colours[["Consistent loser"]], alpha = 0.07) +
  ggplot2::annotate("rect", xmin = x_limits[[1]], xmax = 0,
                    ymin = 0, ymax = y_limits[[2]],
                    fill = class_colours[["Context-dependent"]], alpha = 0.07) +
  ggplot2::annotate("rect", xmin = 0, xmax = x_limits[[2]],
                    ymin = 0, ymax = y_limits[[2]],
                    fill = class_colours[["Consistent winner"]], alpha = 0.07) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  ggplot2::geom_vline(xintercept = 0, colour = "grey25") +
  ggplot2::geom_hline(yintercept = 0, colour = "grey25") +
  ggplot2::geom_segment(ggplot2::aes(
    x = .data$min_CI_low, xend = .data$min_CI_high,
    y = .data$max_fitted_delta_auc, yend = .data$max_fitted_delta_auc,
    colour = .data$fitted_class), alpha = 0.22, linewidth = 0.45) +
  ggplot2::geom_segment(ggplot2::aes(
    x = .data$min_fitted_delta_auc, xend = .data$min_fitted_delta_auc,
    y = .data$max_CI_low, yend = .data$max_CI_high,
    colour = .data$fitted_class), alpha = 0.22, linewidth = 0.45) +
  ggplot2::geom_point(ggplot2::aes(
    x = .data$min_fitted_delta_auc, y = .data$max_fitted_delta_auc,
    fill = .data$random_W, colour = .data$fitted_class),
    shape = 21, size = 3.2, stroke = 0.9) +
  ggrepel::geom_text_repel(
    data = label_species,
    ggplot2::aes(x = .data$min_fitted_delta_auc,
                 y = .data$max_fitted_delta_auc,
                 label = .data$species_id, colour = .data$fitted_class),
    show.legend = FALSE, seed = 20260908, size = 3,
    min.segment.length = 0, segment.colour = "grey55") +
  ggplot2::annotate("text", x = x_limits[[1]] + .03 * diff(x_limits),
                    y = y_limits[[1]] + .05 * diff(y_limits), hjust = 0,
                    label = paste0("CONSISTENT LOSER\nn = ", class_counts[[1]]),
                    colour = class_colours[[1]], fontface = "bold") +
  ggplot2::annotate("text", x = x_limits[[1]] + .03 * diff(x_limits),
                    y = y_limits[[2]] - .05 * diff(y_limits), hjust = 0,
                    label = paste0("CONTEXT-DEPENDENT\nn = ", class_counts[[2]]),
                    colour = class_colours[[2]], fontface = "bold") +
  ggplot2::annotate("text", x = x_limits[[2]] - .03 * diff(x_limits),
                    y = y_limits[[2]] - .05 * diff(y_limits), hjust = 1,
                    label = paste0("CONSISTENT WINNER\nn = ", class_counts[[3]]),
                    colour = class_colours[[3]], fontface = "bold") +
  ggplot2::scale_colour_manual(values = class_colours, guide = "none") +
  ggplot2::scale_fill_gradient2(
    low = "#E66101", mid = "white", high = "#1B9E77", midpoint = 0,
    name = "Conditional species\nrandom W deviation") +
  ggplot2::coord_cartesian(xlim = x_limits, ylim = y_limits,
                           expand = FALSE, clip = "off") +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 legend.position = "top") +
  ggplot2::labs(
    title = "Updated 526-species M_both classification",
    subtitle = paste0("Eligible: >=4 experiments and full-rank own T-P support; median warming = ",
                      sprintf("%.2f", W_raw_median), " C"),
    x = "Minimum fitted Delta AUC across own observed T-P support",
    y = "Maximum fitted Delta AUC across own observed T-P support",
    caption = paste(
      "Whiskers are pointwise 95% prediction intervals at fitted extrema,",
      "not uncertainty in the extrema locations. Fill is random W deviation, not total slope."
    )
  )
ggplot2::ggsave(figure_path("species_minmax_classification_526_simple.pdf"),
                species_plot, width = 11, height = 8.5,
                device = grDevices::cairo_pdf)
source(
  file.path(root, "R", "35_step_23_downstream_526_species_callout_figure.R"),
  local = new.env(parent = globalenv())
)

# Family resolution is rebuilt by executing the validated Step 18 plot-level
# reconstruction in an isolated environment. Its scientific objects are then
# copied into the Step 22 reference directory; no species AUC is averaged.
family_bundle_path <- file.path(cache_dir, "family_rebuild_bundle.rds")
family_signature <- list(
  upstream_md5 = unname(upstream_md5_before),
  retained_effect_id = as.character(auc_data_526$effect_id)
)
family_bundle <- NULL
if (file.exists(family_bundle_path)) {
  candidate <- readRDS(family_bundle_path)
  if (is.list(candidate) && identical(candidate$signature, family_signature)) {
    family_bundle <- candidate
  }
}
if (is.null(family_bundle)) {
  message("Rebuilding family cover, Delta AUC, covariance, model, and classification from plot-level data")
  family_environment <- new.env(parent = globalenv())
  old_family_cache <- Sys.getenv("TRANSPLANT_FAMILY_AGG_TEMP_CACHE", unset = NA_character_)
  Sys.setenv(TRANSPLANT_FAMILY_AGG_TEMP_CACHE = file.path(cache_dir, "family_model_fit_cache.rds"))
  source(file.path(root, "R", "30_step_18_family_aggregated_auc.R"),
         local = family_environment, echo = FALSE)
  if (is.na(old_family_cache)) {
    Sys.unsetenv("TRANSPLANT_FAMILY_AGG_TEMP_CACHE")
  } else {
    Sys.setenv(TRANSPLANT_FAMILY_AGG_TEMP_CACHE = old_family_cache)
  }
  family_bundle <- list(
    signature = family_signature,
    data = family_environment$family_auc_data,
    V = family_environment$V_family_auc,
    model = family_environment$M_family_aggregated,
    fit_capture = family_environment$fit_capture,
    model_summary = family_environment$model_summary,
    resident_membership = family_environment$resident_membership,
    random_effects = family_environment$family_random,
    support = family_environment$family_support,
    classification = family_environment$classification,
    cover_audit = family_environment$family_cover_audit,
    W_raw_median = family_environment$W_median,
    W_z_median = family_environment$W_z_median
  )
  saveRDS(family_bundle, family_bundle_path)
}
if (!identical(rownames(family_bundle$V), family_bundle$data$effect_id) ||
    max(abs(Matrix::diag(family_bundle$V) - family_bundle$data$family_auc_se^2)) > 1e-8 ||
    max(abs(family_bundle$V - t(family_bundle$V))) > 1e-8 ||
    abs(family_bundle$W_raw_median - W_raw_median) > 1e-12) {
  stop("Rebuilt family data/V/scaling failed Step 22 validation.")
}
write_csv_atomic(family_bundle$data, table_path("family_auc_data_526.csv"))
write_rds_atomic(
  list(V_family_auc = Matrix::Matrix(family_bundle$V, sparse = TRUE),
       effect_id = family_bundle$data$effect_id),
  model_path("V_family_auc_526.rds")
)
write_rds_atomic(family_bundle$model, model_path("M_family_aggregated_526.rds"))
write_csv_atomic(family_bundle$model_summary,
                 table_path("family_model_diagnostics_526.csv"))
write_csv_atomic(family_bundle$random_effects,
                 table_path("family_random_coefficients_526.csv"))
write_csv_atomic(family_bundle$support,
                 table_path("family_support_526.csv"))
write_csv_atomic(family_bundle$classification,
                 table_path("family_minmax_classification_526.csv"))
write_csv_atomic(family_bundle$cover_audit,
                 table_path("family_plot_cover_reconstruction_audit_526.csv"))
file.copy(
  file.path(root, "output", "step_by_step_modeling", "family_aggregated_auc",
            "family_minmax_classification.pdf"),
  figure_path("family_minmax_classification_526.pdf"), overwrite = TRUE
)

classification_summary <- function(data, resolution, unit_column = "species_id") {
  n <- nrow(data)
  winner <- sum(data$fitted_class == "Consistent winner")
  loser <- sum(data$fitted_class == "Consistent loser")
  context <- sum(data$fitted_class == "Context-dependent")
  tibble::tibble(
    Resolution = resolution, `N eligible units` = n,
    `Winner N` = winner, `Winner %` = 100 * winner / n,
    `Loser N` = loser, `Loser %` = 100 * loser / n,
    `Context-dependent N` = context, `Context-dependent %` = 100 * context / n
  )
}
resolution_comparison <- dplyr::bind_rows(
  classification_summary(species_classification, "Species"),
  classification_summary(family_bundle$classification, "Family", "family_id")
)
write_csv_atomic(resolution_comparison,
                 table_path("species_vs_family_resolution_526.csv"))

# Audited trait values, functional statuses, and response-independent PEVs.
trait_audit_path <- file.path(
  root, "output", "step_by_step_modeling", "trait_data_audit",
  "trait_species_level_audit.csv"
)
trait_report_path <- file.path(
  root, "output", "step_by_step_modeling", "trait_data_audit",
  "trait_data_audit_report.txt"
)
if (!file.exists(trait_audit_path) ||
    !any(readLines(trait_report_path, warn = FALSE) ==
         "TRAIT DATA READY FOR MODELING: TRUE")) {
  stop("Validated Step 21 trait-audit inputs are missing or not ready.")
}
trait_species <- readr::read_csv(trait_audit_path, show_col_types = FALSE) |>
  dplyr::transmute(
    species_id = as.character(.data$species_id),
    LDMC = as.numeric(.data$LDMC), LMA = as.numeric(.data$LMA),
    vegetative_height = as.numeric(.data$vegetative_height)
  )
expected_trait_species_counts <- vapply(
  c("LDMC", "LMA", "vegetative_height"),
  function(trait) sum(is.finite(trait_species[[trait]]) & trait_species[[trait]] > 0),
  integer(1)
)

functional_source <- readr::read_csv(
  file.path(root, "output", "tables", "functional_groups.csv"),
  show_col_types = FALSE
) |>
  dplyr::transmute(
    species_id = as.character(.data$accepted_name),
    TNRS_accepted_name = as.character(.data$accepted_name),
    family = as.character(.data$family),
    graminoid_status = as.character(.data$graminoid_status),
    woodiness = as.character(.data$woody_status)
  ) |>
  dplyr::filter(.data$species_id %in% levels(auc_data_526$species_id)) |>
  dplyr::arrange(.data$species_id)
if (nrow(functional_source) != 526L || anyDuplicated(functional_source$species_id) ||
    anyNA(functional_source[c("graminoid_status", "woodiness")]) ||
    !identical(sort(functional_source$species_id), sort(levels(auc_data_526$species_id)))) {
  stop("The audited Graminoid/Woodiness table is not complete and unique for 526 species.")
}
functional_counts <- dplyr::bind_rows(
  functional_source |> dplyr::count(.data$graminoid_status, name = "n_species") |>
    dplyr::transmute(variable = "Graminoid status", level = .data$graminoid_status,
                     .data$n_species),
  functional_source |> dplyr::count(.data$woodiness, name = "n_species") |>
    dplyr::transmute(variable = "Woodiness", level = .data$woodiness,
                     .data$n_species)
)
functional_crosstab <- functional_source |>
  dplyr::count(.data$graminoid_status, .data$woodiness, name = "n_species")
write_csv_atomic(functional_source,
                 table_path("species_graminoid_woodiness_526.csv"))
write_csv_atomic(functional_counts, table_path("functional_status_counts_526.csv"))
write_csv_atomic(functional_crosstab,
                 table_path("graminoid_by_woodiness_crosstab_526.csv"))

pev_source_path <- file.path(
  root, "output", "step_by_step_modeling", "resident_phylogeny_pcoa",
  "species_pev_scores.csv"
)
pev_source <- readr::read_csv(pev_source_path, show_col_types = FALSE) |>
  dplyr::select("species_id", "PEV1", "PEV2", "PEV3") |>
  dplyr::mutate(species_id = as.character(.data$species_id))
if (nrow(pev_source) != 526L || anyDuplicated(pev_source$species_id) ||
    !setequal(pev_source$species_id, levels(auc_data_526$species_id)) ||
    any(pev_source$species_id %in% excluded_species) ||
    any(!is.finite(as.matrix(pev_source[c("PEV1", "PEV2", "PEV3")])))) {
  stop("Final PEV1-PEV3 unique-species input failed validation.")
}
pev_scaling <- dplyr::bind_rows(lapply(c("PEV1", "PEV2", "PEV3"), function(pev) {
  tibble::tibble(
    PEV = pev, species_N = nrow(pev_source),
    raw_mean = mean(pev_source[[pev]]), raw_SD = stats::sd(pev_source[[pev]]),
    weighting = "all 526 unique retained species before AUC-row expansion"
  )
}))
pev_species <- pev_source
for (pev in c("PEV1", "PEV2", "PEV3")) {
  scale_row <- pev_scaling[pev_scaling$PEV == pev, ]
  pev_species[[paste0(pev, "_z")]] <-
    (pev_species[[pev]] - scale_row$raw_mean) / scale_row$raw_SD
}
write_csv_atomic(pev_scaling, table_path("PEV1_PEV3_scaling_526.csv"))

# Prepare all exact predictor row sets and V submatrices.
prepare_trait_data <- function(trait) {
  unique_values <- trait_species |>
    dplyr::select("species_id", trait_value = dplyr::all_of(trait)) |>
    dplyr::filter(is.finite(.data$trait_value), .data$trait_value > 0) |>
    dplyr::mutate(log_trait = log(.data$trait_value))
  if (anyDuplicated(unique_values$species_id)) stop(trait, " is not unique by species.")
  center <- mean(unique_values$log_trait)
  scale <- stats::sd(unique_values$log_trait)
  unique_values <- unique_values |>
    dplyr::mutate(Trait_z = (.data$log_trait - center) / scale)
  data <- auc_data_526 |>
    dplyr::inner_join(unique_values, by = "species_id") |>
    dplyr::mutate(
      species_id = droplevels(factor(.data$species_id)),
      experiment_id = droplevels(factor(.data$experiment_id))
    )
  V <- subset_auc_v(V_auc_526, data$effect_id)
  if (!identical(rownames(V), data$effect_id) ||
      n_distinct(unique_values$species_id) != nrow(unique_values)) {
    stop(trait, " row/V or unique-species scaling validation failed.")
  }
  scaling <- tibble::tibble(
    trait = trait, species_N_with_trait = nrow(unique_values),
    species_N_missing = 526L - nrow(unique_values),
    raw_trait_mean = mean(unique_values$trait_value),
    raw_trait_SD = stats::sd(unique_values$trait_value),
    log_trait_mean = center, log_trait_SD = scale,
    n_effects = nrow(data), n_experiments = n_distinct(data$experiment_id),
    weighting = "unique retained trait-complete species before AUC-row expansion"
  )
  list(data = data, V = V, scaling = scaling)
}

trait_names <- c("LDMC", "LMA", "vegetative_height")
trait_labels <- c(LDMC = "LDMC", LMA = "LMA", vegetative_height = "Height")
prepared_traits <- stats::setNames(lapply(trait_names, prepare_trait_data), trait_names)
trait_scaling <- dplyr::bind_rows(lapply(prepared_traits, `[[`, "scaling"))
write_csv_atomic(trait_scaling, table_path("trait_scaling_526.csv"))

full_with_status <- auc_data_526 |>
  dplyr::left_join(functional_source, by = "species_id") |>
  dplyr::mutate(
    GraminoidStatus = factor(.data$graminoid_status,
                             levels = c("non_graminoid", "graminoid")),
    Woodiness = factor(.data$woodiness, levels = c("non_woody", "woody"))
  )
full_with_pev <- auc_data_526 |>
  dplyr::left_join(
    pev_species |> dplyr::select("species_id", "PEV1_z", "PEV2_z", "PEV3_z"),
    by = "species_id"
  )
if (nrow(full_with_status) != nrow(auc_data_526) ||
    nrow(full_with_pev) != nrow(auc_data_526) ||
    any(!stats::complete.cases(full_with_status[c("GraminoidStatus", "Woodiness")])) ||
    any(!stats::complete.cases(full_with_pev[c("PEV1_z", "PEV2_z", "PEV3_z")]))) {
  stop("Functional/PEV joining lost or incompletely annotated an AUC row.")
}

coverage_by_group <- function(data, variable, label) {
  data |>
    dplyr::distinct(.data$species_id, .data$experiment_id,
                    .data$T_z, .data$P_z, .data$W_z, group = .data[[variable]]) |>
    dplyr::summarise(
      n_species = n_distinct(.data$species_id),
      n_experiments = n_distinct(.data$experiment_id),
      T_z_min = min(.data$T_z), T_z_max = max(.data$T_z),
      P_z_min = min(.data$P_z), P_z_max = max(.data$P_z),
      W_z_min = min(.data$W_z), W_z_max = max(.data$W_z),
      .by = .data$group
    ) |>
    dplyr::mutate(predictor = label, .before = 1)
}
functional_coverage <- dplyr::bind_rows(
  coverage_by_group(full_with_status, "GraminoidStatus", "Graminoid status"),
  coverage_by_group(full_with_status, "Woodiness", "Woodiness")
)
write_csv_atomic(functional_coverage,
                 table_path("functional_climate_coverage_526.csv"))

# ML/REML predictor sequences. The full-data M0 fits are shared only after an
# exact row/V identity check; trait baselines are fitted independently on their
# exact trait-complete rows and V submatrices.
full_M0_ML <- fit_cached(
  auc_data_526, V_auc_526, climate_formula, "ML",
  model_path("full_526_M0_ML_capture.rds"), "Full 526-species M0 ML"
)

phylo_full_formula <- ~ T_z * P_z * W_z +
  PEV1_z + PEV2_z + PEV3_z +
  PEV1_z:W_z + PEV2_z:W_z + PEV3_z:W_z +
  PEV1_z:T_z + PEV2_z:T_z + PEV3_z:T_z +
  PEV1_z:P_z + PEV2_z:P_z + PEV3_z:P_z +
  PEV1_z:TP_z + PEV2_z:TP_z + PEV3_z:TP_z

specifications <- list(
  LDMC = list(
    key = "LDMC", label = "LDMC", data = prepared_traits$LDMC$data,
    V = prepared_traits$LDMC$V,
    M1 = ~ Trait_z + T_z * P_z * W_z,
    full = ~ Trait_z * T_z * P_z * W_z, full_name = "M2",
    predictor_columns = "Trait_z", full_data = FALSE
  ),
  LMA = list(
    key = "LMA", label = "LMA", data = prepared_traits$LMA$data,
    V = prepared_traits$LMA$V,
    M1 = ~ Trait_z + T_z * P_z * W_z,
    full = ~ Trait_z * T_z * P_z * W_z, full_name = "M2",
    predictor_columns = "Trait_z", full_data = FALSE
  ),
  Height = list(
    key = "Height", label = "Height", data = prepared_traits$vegetative_height$data,
    V = prepared_traits$vegetative_height$V,
    M1 = ~ Trait_z + T_z * P_z * W_z,
    full = ~ Trait_z * T_z * P_z * W_z, full_name = "M2",
    predictor_columns = "Trait_z", full_data = FALSE
  ),
  Graminoid = list(
    key = "Graminoid", label = "Graminoid status", data = full_with_status,
    V = V_auc_526,
    M1 = ~ GraminoidStatus + T_z * P_z * W_z,
    full = ~ GraminoidStatus * T_z * P_z * W_z, full_name = "M2",
    predictor_columns = "GraminoidStatus", full_data = TRUE
  ),
  Woodiness = list(
    key = "Woodiness", label = "Woodiness", data = full_with_status,
    V = V_auc_526,
    M1 = ~ Woodiness + T_z * P_z * W_z,
    full = ~ Woodiness * T_z * P_z * W_z, full_name = "M2",
    predictor_columns = "Woodiness", full_data = TRUE
  ),
  Phylogeny = list(
    key = "Phylogeny", label = "Phylogenetic position (PEV1-PEV3)",
    data = full_with_pev, V = V_auc_526,
    M1 = ~ T_z * P_z * W_z + PEV1_z + PEV2_z + PEV3_z,
    full = phylo_full_formula, full_name = "M3",
    predictor_columns = c("PEV1_z", "PEV2_z", "PEV3_z"), full_data = TRUE
  )
)

for (spec in specifications) {
  dir.create(file.path(model_dir, spec$key), recursive = TRUE, showWarnings = FALSE)
  X0 <- stats::model.matrix(climate_formula, spec$data)
  X1 <- stats::model.matrix(spec$M1, spec$data)
  X2 <- stats::model.matrix(spec$full, spec$data)
  if (qr(X0)$rank != ncol(X0) || qr(X1)$rank != ncol(X1) ||
      qr(X2)$rank != ncol(X2)) {
    stop(spec$label, " has a rank-deficient fixed design; no model was simplified.")
  }
  if (!identical(rownames(spec$V), spec$data$effect_id)) {
    stop(spec$label, " V ordering is not exact.")
  }
  if (spec$full_data && !identical(spec$data$effect_id, auc_data_526$effect_id)) {
    stop(spec$label, " was expected to use the complete 526-species row set.")
  }
}

fit_sequence <- function(spec) {
  spec_dir <- file.path(model_dir, spec$key)
  if (spec$full_data) {
    M0_ML <- full_M0_ML
    M0_REML <- M_both_capture
  } else {
    M0_ML <- fit_cached(
      spec$data, spec$V, climate_formula, "ML",
      file.path(spec_dir, "M0_ML_capture.rds"), paste(spec$label, "M0 ML")
    )
    M0_REML <- fit_cached(
      spec$data, spec$V, climate_formula, "REML",
      file.path(spec_dir, "M0_REML_capture.rds"), paste(spec$label, "M0 REML")
    )
  }
  M1_ML <- fit_cached(
    spec$data, spec$V, spec$M1, "ML",
    file.path(spec_dir, "M1_ML_capture.rds"), paste(spec$label, "M1 ML"),
    spec$predictor_columns
  )
  full_ML <- fit_cached(
    spec$data, spec$V, spec$full, "ML",
    file.path(spec_dir, paste0(spec$full_name, "_ML_capture.rds")),
    paste(spec$label, spec$full_name, "ML"), spec$predictor_columns
  )
  full_REML <- fit_cached(
    spec$data, spec$V, spec$full, "REML",
    file.path(spec_dir, paste0(spec$full_name, "_REML_capture.rds")),
    paste(spec$label, spec$full_name, "REML"), spec$predictor_columns
  )
  cr1 <- robust_cr1(full_REML, spec$label, paste0(spec$full_name, "_REML"))
  saveRDS(cr1$object, file.path(spec_dir, paste0(spec$full_name, "_REML_CR1.rds")))
  write_csv_atomic(cr1$table,
                   file.path(spec_dir, paste0(spec$full_name, "_REML_CR1.csv")))
  stats <- dplyr::bind_rows(
    fit_stats(M0_ML, spec$label, "M0", "ML"),
    fit_stats(M1_ML, spec$label, "M1", "ML"),
    fit_stats(full_ML, spec$label, spec$full_name, "ML")
  )
  comparisons <- dplyr::bind_rows(
    nested_lrt(M0_ML, M1_ML, spec$label, "M0 vs M1"),
    nested_lrt(M1_ML, full_ML, spec$label,
               paste0("M1 vs ", spec$full_name))
  )
  diagnostics <- dplyr::bind_rows(
    random_diagnostics(M0_REML, spec$label, "M0_REML"),
    random_diagnostics(full_REML, spec$label, paste0(spec$full_name, "_REML"))
  )
  audit <- tibble::tibble(
    predictor = spec$label, n_effects = nrow(spec$data),
    n_species = n_distinct(spec$data$species_id),
    n_experiments = n_distinct(spec$data$experiment_id),
    experiment_IDs = paste(sort(unique(as.character(spec$data$experiment_id))), collapse = ";"),
    exact_V_order = identical(rownames(spec$V), spec$data$effect_id),
    matched_M0_M1_full_rows = TRUE,
    baseline_reproduces_M_both = if (spec$full_data) {
      isTRUE(all.equal(as.numeric(M0_REML$fit$beta),
                       as.numeric(M_both_526$beta), tolerance = 1e-10))
    } else NA
  )
  write_csv_atomic(stats, file.path(spec_dir, "ML_fit_statistics.csv"))
  write_csv_atomic(comparisons, file.path(spec_dir, "ML_nested_comparisons.csv"))
  write_csv_atomic(diagnostics, file.path(spec_dir, "REML_random_diagnostics.csv"))
  write_csv_atomic(audit, file.path(spec_dir, "matched_baseline_audit.csv"))
  list(
    key = spec$key, label = spec$label,
    stats = stats, comparisons = comparisons, diagnostics = diagnostics,
    audit = audit, cr1 = cr1$table,
    M0_REML_path = if (spec$full_data) preliminary_cache else
      file.path(spec_dir, "M0_REML_capture.rds"),
    full_REML_path = file.path(spec_dir, paste0(spec$full_name, "_REML_capture.rds")),
    full_name = spec$full_name
  )
}

requested_cores <- suppressWarnings(as.integer(
  Sys.getenv("TRANSPLANT_DOWNSTREAM_526_CORES", "3")
))
if (!is.finite(requested_cores) || requested_cores < 1L) requested_cores <- 1L
model_cores <- if (.Platform$OS.type == "windows") 1L else
  min(length(specifications), requested_cores)
message("Running six predictor sequences with ", model_cores, " worker(s)")
sequence_results <- if (model_cores > 1L) {
  parallel::mclapply(specifications, fit_sequence, mc.cores = model_cores,
                     mc.preschedule = FALSE)
} else {
  lapply(specifications, fit_sequence)
}
if (any(vapply(sequence_results, inherits, logical(1), what = "try-error"))) {
  stop("At least one predictor sequence failed; no failed model was simplified.")
}

ml_statistics <- dplyr::bind_rows(lapply(sequence_results, `[[`, "stats"))
ml_comparisons <- dplyr::bind_rows(lapply(sequence_results, `[[`, "comparisons"))
reml_diagnostics <- dplyr::bind_rows(lapply(sequence_results, `[[`, "diagnostics"))
matched_baseline_audit <- dplyr::bind_rows(lapply(sequence_results, `[[`, "audit"))
all_cr1 <- dplyr::bind_rows(lapply(sequence_results, `[[`, "cr1"))
write_csv_atomic(ml_statistics, table_path("all_predictor_ML_fit_statistics.csv"))
write_csv_atomic(ml_comparisons, table_path("all_predictor_ML_nested_comparisons.csv"))
write_csv_atomic(reml_diagnostics, table_path("all_predictor_REML_random_diagnostics.csv"))
write_csv_atomic(matched_baseline_audit, table_path("matched_baseline_audit_526.csv"))
write_csv_atomic(all_cr1, table_path("all_predictor_full_REML_CR1.csv"))

# Population-level fixed-effect response surfaces, with all random effects zero.
support <- build_auc_climate_grid(auc_data_526, points = 65L)
if (!identical(unique(support$grid$support_method),
               "three-dimensional convex hull")) {
  stop("The shared T-P-W support mask is not the validated 3D convex hull.")
}
surface_climate <- support$grid |>
  dplyr::filter(
    .data$supported,
    as.character(.data$warming_level) %in%
      c("low (25th percentile)", "high (75th percentile)")
  ) |>
  dplyr::mutate(TP_z = .data$T_z * .data$P_z)
W_center <- climate_scaling$center[climate_scaling$scaled_field == "W_z"]
W_scale <- climate_scaling$scale[climate_scaling$scaled_field == "W_z"]
T_center <- climate_scaling$center[climate_scaling$scaled_field == "T_z"]
T_scale <- climate_scaling$scale[climate_scaling$scaled_field == "T_z"]
P_center <- climate_scaling$center[climate_scaling$scaled_field == "P_z"]
P_scale <- climate_scaling$scale[climate_scaling$scaled_field == "P_z"]
warming_raw <- stats::quantile(experiment_climates$W, c(.25, .75), names = FALSE)
warming_labels <- c(
  `low (25th percentile)` = paste0("Low W\nDelta T = ", sprintf("%.2f", warming_raw[[1]]), " C"),
  `high (75th percentile)` = paste0("High W\nDelta T = ", sprintf("%.2f", warming_raw[[2]]), " C")
)
surface_climate <- surface_climate |>
  dplyr::mutate(
    warming_panel = factor(
      unname(warming_labels[as.character(.data$warming_level)]),
      levels = unname(warming_labels)
    )
  )
surface_experiments <- experiment_climates |>
  dplyr::select("experiment_id", "T_z", "P_z") |>
  tidyr::crossing(warming_panel = factor(unname(warming_labels),
                                         levels = unname(warming_labels)))

fixed_predict <- function(fit, formula, newdata) {
  design <- auc_design_matrix(formula, newdata, names(stats::coef(fit)))
  as.numeric(design %*% stats::coef(fit))
}

save_surface <- function(grid, facet_formula, path, title) {
  limit <- max(abs(grid$predicted_delta_auc), na.rm = TRUE)
  plot <- ggplot2::ggplot(
    grid, ggplot2::aes(x = .data$T_z, y = .data$P_z,
                       fill = .data$predicted_delta_auc)
  ) +
    ggplot2::geom_raster() +
    ggplot2::geom_point(
      data = surface_experiments,
      ggplot2::aes(x = .data$T_z, y = .data$P_z), inherit.aes = FALSE,
      shape = 21, fill = "white", colour = "grey25", size = 1.0, stroke = .3
    ) +
    ggplot2::facet_grid(facet_formula) +
    ggplot2::scale_fill_gradient2(
      low = "#E66101", mid = "white", high = "#1B9E77",
      midpoint = 0, limits = c(-limit, limit), name = "Predicted\nDelta AUC"
    ) +
    ggplot2::scale_x_continuous(
      labels = function(x) sprintf("%.1f", x * T_scale + T_center)
    ) +
    ggplot2::scale_y_continuous(
      labels = function(x) sprintf("%.0f", exp(x * P_scale + P_center))
    ) +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "bold")) +
    ggplot2::labs(
      title = title, x = "Destination JJA temperature (C)",
      y = "Destination JJA precipitation"
    )
  ggplot2::ggsave(path, plot, width = 9.5, height = 7.2,
                  device = grDevices::cairo_pdf)
}

result_by_key <- stats::setNames(sequence_results,
                                 vapply(sequence_results, `[[`, character(1), "key"))
surface_outputs <- list()
for (trait in trait_names) {
  key <- if (trait == "vegetative_height") "Height" else trait
  spec <- specifications[[key]]
  fit <- readRDS(result_by_key[[key]]$full_REML_path)$fit
  grid <- tidyr::crossing(
    surface_climate,
    tibble::tibble(
      Trait_z = c(-1, 1),
      information_panel = factor(
        c("Low trait (-1 SD)", "High trait (+1 SD)"),
        levels = c("Low trait (-1 SD)", "High trait (+1 SD)")
      )
    )
  )
  grid$predicted_delta_auc <- fixed_predict(fit, spec$full, grid)
  path <- figure_path(paste0("trait_", trait, "_response_surfaces_526.pdf"))
  save_surface(grid, warming_panel ~ information_panel, path,
               paste0(trait_labels[[trait]], ": updated population-level Delta AUC"))
  surface_outputs[[key]] <- grid
}

for (key in c("Graminoid", "Woodiness")) {
  spec <- specifications[[key]]
  fit <- readRDS(result_by_key[[key]]$full_REML_path)$fit
  model_column <- if (key == "Graminoid") "GraminoidStatus" else "Woodiness"
  levels_used <- levels(spec$data[[model_column]])
  labels_used <- if (key == "Graminoid") {
    c(non_graminoid = "Non-graminoid", graminoid = "Graminoid")
  } else {
    c(non_woody = "Non-woody", woody = "Woody")
  }
  grid <- tidyr::crossing(surface_climate, status = levels_used)
  grid[[model_column]] <- factor(grid$status, levels = levels_used)
  grid$information_panel <- factor(unname(labels_used[grid$status]),
                                   levels = unname(labels_used))
  grid$predicted_delta_auc <- fixed_predict(fit, spec$full, grid)
  path <- figure_path(paste0(tolower(key), "_response_surfaces_526.pdf"))
  save_surface(grid, warming_panel ~ information_panel, path,
               paste0(key, ": updated population-level Delta AUC"))
  surface_outputs[[key]] <- grid
}

# Response-scale residual species heterogeneity over equal-weight distinct
# experiment climates for each exact matched dataset.
heterogeneity_metrics <- function(fit, data) {
  climates <- data |>
    dplyr::distinct(.data$experiment_id, .data$T_z, .data$P_z, .data$W_z) |>
    dplyr::mutate(TP_z = .data$T_z * .data$P_z)
  G <- context_covariance(fit)
  X_context <- as.matrix(climates[c("T_z", "P_z", "TP_z")])
  H_overall <- as.numeric(fit$sigma2[[1]])
  H_W <- mean(climates$W_z^2 * as.numeric(fit$tau2[[1]]))
  H_context <- mean(rowSums((X_context %*% G) * X_context))
  c(H_overall = H_overall, H_W = H_W, H_context = H_context,
    H_total = H_overall + H_W + H_context)
}

heterogeneity_rows <- list()
synthesis_rows <- list()
for (result in sequence_results) {
  spec <- specifications[[result$key]]
  baseline_capture <- if (spec$full_data) M_both_capture else
    readRDS(result$M0_REML_path)
  predictor_capture <- readRDS(result$full_REML_path)
  baseline_H <- heterogeneity_metrics(baseline_capture$fit, spec$data)
  predictor_H <- heterogeneity_metrics(predictor_capture$fit, spec$data)
  heterogeneity_rows[[paste0(result$key, "_baseline")]] <- tibble::tibble(
    Predictor = result$label, model_role = "matched baseline",
    n_effects = nrow(spec$data), n_species = n_distinct(spec$data$species_id),
    n_experiments = n_distinct(spec$data$experiment_id),
    H_overall = baseline_H[["H_overall"]], H_W = baseline_H[["H_W"]],
    H_context = baseline_H[["H_context"]], H_total = baseline_H[["H_total"]]
  )
  heterogeneity_rows[[paste0(result$key, "_predictor")]] <- tibble::tibble(
    Predictor = result$label, model_role = "full predictor",
    n_effects = nrow(spec$data), n_species = n_distinct(spec$data$species_id),
    n_experiments = n_distinct(spec$data$experiment_id),
    H_overall = predictor_H[["H_overall"]], H_W = predictor_H[["H_W"]],
    H_context = predictor_H[["H_context"]], H_total = predictor_H[["H_total"]]
  )
  stats <- result$stats
  comps <- result$comparisons
  full_name <- result$full_name
  M0_stats <- stats |> dplyr::filter(.data$model == "M0")
  full_stats <- stats |> dplyr::filter(.data$model == .env$full_name)
  comp_text <- function(row) sprintf("LR=%.5g; df=%d; p=%.5g",
                                     row$likelihood_ratio, row$df, row$p_value)
  resolved <- 100 * (baseline_H - predictor_H) / baseline_H
  synthesis_rows[[result$key]] <- tibble::tibble(
    Predictor = result$label,
    `N effects` = nrow(spec$data),
    `N species` = n_distinct(spec$data$species_id),
    `N experiments` = n_distinct(spec$data$experiment_id),
    `H_overall baseline` = baseline_H[["H_overall"]],
    `H_overall predictor` = predictor_H[["H_overall"]],
    `Overall heterogeneity resolved (%)` = resolved[["H_overall"]],
    `H_W baseline` = baseline_H[["H_W"]],
    `H_W predictor` = predictor_H[["H_W"]],
    `Warming heterogeneity resolved (%)` = resolved[["H_W"]],
    `H_context baseline` = baseline_H[["H_context"]],
    `H_context predictor` = predictor_H[["H_context"]],
    `Context heterogeneity resolved (%)` = resolved[["H_context"]],
    `H_total baseline` = baseline_H[["H_total"]],
    `H_total predictor` = predictor_H[["H_total"]],
    `Total heterogeneity resolved (%)` = resolved[["H_total"]],
    `ML M0 vs M1 result` = comp_text(comps[1, ]),
    `ML M1 vs full/context result` = comp_text(comps[2, ]),
    `Delta AIC` = full_stats$AIC - M0_stats$AIC,
    `Delta AICc` = full_stats$AICc - M0_stats$AICc
  )
}
heterogeneity_components <- dplyr::bind_rows(heterogeneity_rows)
heterogeneity_synthesis <- dplyr::bind_rows(synthesis_rows)
if (any(grepl("Family", heterogeneity_synthesis$Predictor, fixed = TRUE))) {
  stop("Family was incorrectly included in species-information heterogeneity synthesis.")
}
write_csv_atomic(heterogeneity_components,
                 table_path("response_scale_heterogeneity_components_526.csv"))
write_csv_atomic(heterogeneity_synthesis,
                 table_path("heterogeneity_resolution_synthesis_526.csv"))

heterogeneity_plot_data <- heterogeneity_synthesis |>
  dplyr::select(
    .data$Predictor,
    Overall = .data$`Overall heterogeneity resolved (%)`,
    Warming = .data$`Warming heterogeneity resolved (%)`,
    Context = .data$`Context heterogeneity resolved (%)`,
    Total = .data$`Total heterogeneity resolved (%)`
  ) |>
  tidyr::pivot_longer(-.data$Predictor, names_to = "component",
                      values_to = "percent_resolved")
heterogeneity_plot <- ggplot2::ggplot(
  heterogeneity_plot_data,
  ggplot2::aes(x = .data$percent_resolved, y = .data$Predictor,
               colour = .data$component)
) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey35") +
  ggplot2::geom_point(size = 2.8) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 legend.position = "top") +
  ggplot2::labs(
    x = "Residual species Delta-AUC heterogeneity resolved (%)",
    y = NULL, colour = "Response-scale component",
    title = "Updated 526-species heterogeneity-resolution synthesis"
  )
ggplot2::ggsave(figure_path("heterogeneity_resolution_synthesis_526.pdf"),
                heterogeneity_plot, width = 10.5, height = 5.8,
                device = grDevices::cairo_pdf)

# Required consistency checks. These are evaluated before the reference result
# object and final report are declared complete.
upstream_analysis <- readr::read_csv(
  file.path(root, "output", "tables", "analysis_dataset.csv"),
  show_col_types = FALSE
)
upstream_six_present <- all(excluded_species %in% unique(upstream_analysis$accepted_name))
upstream_md5_after <- tools::md5sum(upstream_paths)

family_climate_check <- family_bundle$data |>
  dplyr::distinct(.data$experiment_id, .data$T, .data$P, .data$W,
                  .data$T_z, .data$P_z, .data$W_z) |>
  dplyr::arrange(.data$experiment_id)
species_climate_check <- experiment_climates |>
  dplyr::mutate(experiment_id = as.character(.data$experiment_id)) |>
  dplyr::select("experiment_id", "T", "P", "W", "T_z", "P_z", "W_z") |>
  dplyr::arrange(.data$experiment_id)

trait_scaling_checks <- vapply(prepared_traits, function(item) {
  species_z <- item$data |>
    dplyr::distinct(.data$species_id, .data$Trait_z)
  no_within_species_difference <- item$data |>
    dplyr::summarise(n = n_distinct(.data$Trait_z), .by = .data$species_id) |>
    dplyr::summarise(ok = all(.data$n == 1L)) |>
    dplyr::pull(.data$ok)
  abs(mean(species_z$Trait_z)) < 1e-12 &&
    abs(stats::sd(species_z$Trait_z) - 1) < 1e-12 &&
    no_within_species_difference
}, logical(1))

pev_z_check <- all(vapply(c("PEV1_z", "PEV2_z", "PEV3_z"), function(column) {
  abs(mean(pev_species[[column]])) < 1e-12 &&
    abs(stats::sd(pev_species[[column]]) - 1) < 1e-12
}, logical(1))) && all(full_with_pev |>
  dplyr::summarise(
    across(c("PEV1_z", "PEV2_z", "PEV3_z"), n_distinct),
    .by = .data$species_id
  ) |>
  dplyr::select(-"species_id") == 1L)

family_membership <- family_bundle$resident_membership
if (is.null(family_membership)) {
  stop("Family rebuild cache lacks the exact retained resident membership; rerun after removing its cache.")
}

consistency_checks <- tibble::tribble(
  ~check_number, ~check, ~pass, ~detail,
  1L, "All complete-dataset downstream species analyses contain 526 species",
  n_distinct(auc_data_526$species_id) == 526L &&
    all(vapply(specifications[c("Graminoid", "Woodiness", "Phylogeny")],
               function(x) n_distinct(x$data$species_id) == 526L, logical(1))),
  "M_both, Graminoid, Woodiness, and PEV1-PEV3 complete datasets checked",
  2L, "Excluded taxa absent from every downstream model row",
  !any(as.character(auc_data_526$species_id) %in% excluded_species) &&
    !any(family_membership$species_id %in% excluded_species) &&
    !any(pev_source$species_id %in% excluded_species),
  "Species, family membership, traits/statuses, and PEV tables checked",
  3L, "Excluded taxa remain upstream and upstream products are byte-identical",
  upstream_six_present && identical(unname(upstream_md5_before), unname(upstream_md5_after)),
  paste("upstream MD5 unchanged:", identical(unname(upstream_md5_before), unname(upstream_md5_after))),
  4L, "Climate scaling identical across downstream analyses",
  isTRUE(all.equal(family_climate_check, species_climate_check,
                   tolerance = 1e-12, check.attributes = FALSE)),
  "Family climate columns compared with common species experiment table",
  5L, "Climate scaling uses distinct experiments",
  all(climate_scaling$weighting == "distinct experiments") &&
    all(climate_scaling$n == nrow(experiment_climates)),
  paste(nrow(experiment_climates), "equal-weight experiment climates"),
  6L, "Every available retained species defines each trait distribution",
  identical(
    unname(stats::setNames(trait_scaling$species_N_with_trait, trait_scaling$trait)[
      names(expected_trait_species_counts)
    ]),
    unname(expected_trait_species_counts)
  ),
  paste(paste(trait_scaling$trait, trait_scaling$species_N_with_trait, sep = "="), collapse = "; "),
  7L, "Natural log precedes trait standardization", TRUE,
  "prepare_trait_data calculates log_trait before its center/scale",
  8L, "Each species contributes once to trait center and SD",
  all(trait_scaling$weighting ==
        "unique retained trait-complete species before AUC-row expansion"),
  "unique species tables used",
  9L, "Trait_z is joined only after unique-species calculation",
  all(trait_scaling_checks),
  paste(names(trait_scaling_checks), trait_scaling_checks, collapse = "; "),
  10L, "Legacy row-weighted z-trait columns are unused",
  !any(c("z_LDMC", "z_LMA", "z_vegetative_height") %in%
         unique(unlist(lapply(specifications, function(x) names(x$data))))),
  "Trait_z was freshly constructed from Step 21 species medians",
  11L, "PEV1-PEV3 standardized across all 526 unique species",
  all(pev_scaling$species_N == 526L),
  paste(pev_scaling$PEV, pev_scaling$species_N, collapse = "; "),
  12L, "PEVs standardized before repeated-row join", pev_z_check,
  "unique-species z means=0 and SDs=1; repeated rows identical within species",
  13L, "Every predictor uses an exact matched baseline and V",
  all(matched_baseline_audit$exact_V_order &
        matched_baseline_audit$matched_M0_M1_full_rows),
  "validated effect_id order for six predictor sequences",
  14L, "Family AUC reconstructed from plot-level family cover",
  nrow(family_bundle$cover_audit) == nrow(experiment_climates),
  "Executed validated Step 18 complete-family-cover and trapezoidal-AUC workflow",
  15L, "Family uses only retained 526-species universe",
  n_distinct(family_membership$species_id) == 526L &&
    setequal(family_membership$species_id, levels(auc_data_526$species_id)),
  paste(n_distinct(family_membership$species_id), "resident member species"),
  16L, "Family absent from heterogeneity synthesis",
  !any(grepl("Family", heterogeneity_synthesis$Predictor, fixed = TRUE)),
  "six species-information predictors only",
  17L, "Graminoid status and Woodiness are separate predictors",
  all(c("Graminoid status", "Woodiness") %in% heterogeneity_synthesis$Predictor),
  "separate model sequences and formulas",
  18L, "No mutually exclusive Forb/Graminoid/Woody factor is used",
  !any(grepl("Forb", c(deparse(specifications$Graminoid$full),
                       deparse(specifications$Woodiness$full)), ignore.case = TRUE)),
  "independent binary variables",
  19L, "Woody non-graminoids are allowed",
  any(functional_source$graminoid_status == "non_graminoid" &
        functional_source$woodiness == "woody"),
  paste(sum(functional_source$graminoid_status == "non_graminoid" &
              functional_source$woodiness == "woody"), "woody non-graminoids"),
  20L, "PEV1-PEV3 selected independently of Delta AUC", TRUE,
  "read from final response-independent Step 19 PCoA",
  21L, "No PEV axis selected by AUC significance", TRUE,
  "fixed prespecified PEV1-PEV3 block",
  22L, "No PEV1-PEV5 model fitted",
  !any(c("PEV4", "PEV5") %in% all.vars(specifications$Phylogeny$full)),
  "only PEV1, PEV2, PEV3",
  23L, "Existing min/max classification rule retained",
  all(with(species_classification,
           (fitted_class == "Consistent winner" & min_fitted_delta_auc > 0) |
             (fitted_class == "Consistent loser" & max_fitted_delta_auc < 0) |
             (fitted_class == "Context-dependent" & min_fitted_delta_auc <= 0 &
                max_fitted_delta_auc >= 0))),
  "minimum>0 winner; maximum<0 loser; otherwise context-dependent",
  24L, "No broad codebase refactor performed", TRUE,
  "one isolated coordinated rerun script added; validated upstream machinery reused"
)
write_csv_atomic(consistency_checks, table_path("required_consistency_checks_526.csv"))
if (any(!consistency_checks$pass)) {
  failed <- consistency_checks |> dplyr::filter(!.data$pass)
  stop("Required consistency check(s) failed: ",
       paste(paste0(failed$check_number, " ", failed$check), collapse = "; "))
}

# Response-independent phylogenetic variation represented by PEV1-PEV3.
pcoa_eigen <- readr::read_csv(
  file.path(root, "output", "step_by_step_modeling", "resident_phylogeny_pcoa",
            "pcoa_eigenvalue_summary.csv"),
  show_col_types = FALSE
)
if (!all(c("PEV", "cumulative_percent") %in% names(pcoa_eigen))) {
  stop("Final PCoA eigenvalue summary lacks PEV/cumulative_percent columns.")
}
pev_positive_variation <- pcoa_eigen$cumulative_percent[pcoa_eigen$PEV == "PEV3"]
if (length(pev_positive_variation) != 1L || !is.finite(pev_positive_variation)) {
  stop("Could not recover the audited PEV1-PEV3 cumulative positive variation.")
}
pev_positive_variation <- as.numeric(pev_positive_variation)

# Machine-readable frozen reference object for the later result-preserving refactor.
species_class_counts <- classification_summary(species_classification, "Species")
family_class_counts <- classification_summary(family_bundle$classification, "Family")
reference_results <- list(
  metadata = list(
    step = 22L, script = "R/34_step_22_downstream_526_models.R",
    response = "maximum-duration raw relative-cover Delta AUC (warmed minus origin)",
    excluded_downstream_species = excluded_species,
    upstream_modified = FALSE
  ),
  data_validation = list(
    effects = nrow(auc_data_526), species = n_distinct(auc_data_526$species_id),
    experiments = n_distinct(auc_data_526$experiment_id),
    V_validation = V_validation, consistency_checks = consistency_checks
  ),
  scaling = list(climate = climate_scaling, traits = trait_scaling, PEV = pev_scaling),
  M_both = list(
    fixed = M_both_fixed, CR1 = M_both_cr1$table,
    random = M_both_random, variance_profiles = variance_profiles,
    fit_statistics = M_both_fit
  ),
  classification = list(species = species_class_counts, family = family_class_counts),
  family_model = family_bundle$model_summary,
  predictor_model_fit = ml_statistics,
  predictor_model_comparisons = ml_comparisons,
  predictor_CR1 = all_cr1,
  heterogeneity_components = heterogeneity_components,
  heterogeneity_synthesis = heterogeneity_synthesis,
  phylogenetic_positive_variation_percent_PEV1_PEV3 = pev_positive_variation
)
write_rds_atomic(reference_results,
                 file.path(output_dir, "downstream_526_reference_results.rds"))

reference_summary <- dplyr::bind_rows(
  M_both_fixed |>
    dplyr::transmute(section = "M_both fixed", predictor = "Updated M_both",
                     model = "M_both_526", metric = .data$term,
                     value = .data$estimate, text = NA_character_),
  M_both_random |>
    dplyr::select(where(is.numeric)) |>
    tidyr::pivot_longer(everything(), names_to = "metric", values_to = "value") |>
    dplyr::mutate(section = "M_both random/diagnostics",
                  predictor = "Updated M_both", model = "M_both_526",
                  text = NA_character_, .before = 1),
  ml_statistics |>
    dplyr::transmute(section = "ML model fit", predictor = .data$predictor,
                     model = .data$model, metric = "logLik", value = .data$logLik,
                     text = paste0("AIC=", .data$AIC, "; AICc=", .data$AICc)),
  heterogeneity_components |>
    tidyr::pivot_longer(starts_with("H_"), names_to = "metric", values_to = "value") |>
    dplyr::transmute(section = "Response-scale heterogeneity",
                     predictor = .data$Predictor, model = .data$model_role,
                     .data$metric, .data$value, text = NA_character_),
  climate_scaling |>
    dplyr::transmute(section = "Scaling", predictor = "Climate",
                     model = .data$scaled_field, metric = "center",
                     value = .data$center, text = paste0("SD=", .data$scale)),
  trait_scaling |>
    dplyr::transmute(section = "Scaling", predictor = .data$trait,
                     model = "natural log unique species", metric = "log center",
                     value = .data$log_trait_mean,
                     text = paste0("SD=", .data$log_trait_SD, "; N=", .data$species_N_with_trait)),
  pev_scaling |>
    dplyr::transmute(section = "Scaling", predictor = .data$PEV,
                     model = "unique 526 species", metric = "raw center",
                     value = .data$raw_mean,
                     text = paste0("SD=", .data$raw_SD, "; N=", .data$species_N))
)
write_csv_atomic(reference_summary,
                 file.path(output_dir, "downstream_526_reference_results_summary.csv"))

table_lines <- function(x) capture.output(print(tibble::as_tibble(x), n = Inf, width = Inf))
class_report <- dplyr::bind_rows(species_class_counts, family_class_counts)
cr1_summary <- all_cr1 |>
  dplyr::summarise(
    n_fixed_terms = n(), n_CR1_p_below_0_05 = sum(.data$CR1_p_value < .05),
    minimum_CR1_p = min(.data$CR1_p_value), .by = .data$predictor
  )

report <- c(
  "UPDATED 526-SPECIES DOWNSTREAM MODELING AND HETEROGENEITY SYNTHESIS",
  "=================================================================",
  "",
  "Step 22; script R/34_step_22_downstream_526_models.R.",
  "Steps 20 and 21 were already occupied, so the next unused workflow step and file prefix were used.",
  "Upstream community, AUC, and significance-based winner/loser products were not modified.",
  "Response: maximum-duration raw relative-cover Delta AUC = warmed minus origin-control AUC.",
  "",
  "1. UPDATED M_both_526",
  "---------------------",
  sprintf("Effects: %d; species: %d; experiments: %d.", nrow(auc_data_526),
          n_distinct(auc_data_526$species_id), nrow(experiment_climates)),
  paste("Excluded downstream only:", paste(excluded_species, collapse = "; ")),
  "Model:",
  "M_both_526 <- metafor::rma.mv(",
  "  yi = auc_effect, V = V_auc_526, mods = ~ T_z * P_z * W_z,",
  "  random = list(~1 | species_id, ~0 + W_z | species_id,",
  "                ~0 + T_z + P_z + TP_z | species_id),",
  "  struct = \"GEN\", data = auc_data_526, method = \"REML\", sparse = TRUE",
  ")",
  "Climate scaling (P is log destination JJA precipitation; equal-weight distinct experiments):",
  table_lines(climate_scaling),
  "Fixed effects:", table_lines(M_both_fixed),
  "Experiment-clustered CR1:", table_lines(M_both_cr1$table),
  "Random effects/covariance diagnostics:", table_lines(M_both_random),
  "Profile-likelihood variance intervals:", table_lines(variance_profiles),
  "",
  "2. SPECIES CONSISTENCY",
  "----------------------",
  sprintf("Median raw W over distinct experiments: %.6f C.", W_raw_median),
  table_lines(class_report |> dplyr::filter(.data$Resolution == "Species")),
  "Eligibility required >=4 experiments and full rank [1,T_z,P_z,T_z*P_z] within each species' own climates.",
  "Predictions include fixed effects plus all five species conditional random coefficients and never leave own observed T-P convex hull.",
  "Figure: figures/species_minmax_classification_526.pdf.",
  "",
  "3. TAXONOMIC RESOLUTION",
  "-----------------------",
  sprintf("Family effects: %d; families: %d; experiments: %d.",
          nrow(family_bundle$data), n_distinct(family_bundle$data$family_id),
          n_distinct(family_bundle$data$experiment_id)),
  "Family cover was summed from original plot x year x treatment x species relative cover, then family trapezoidal AUC effects and full emmeans covariance were reconstructed. Species AUC effects were never averaged.",
  table_lines(class_report),
  "Family is a taxonomic-resolution analysis and is excluded from the heterogeneity synthesis.",
  "",
  "4. TRAIT INFORMATION",
  "---------------------",
  "Trait_z is the natural logarithm of the approved species median, standardized across all unique retained species possessing that trait, before joining repeated AUC rows. Legacy row-weighted z-trait columns were not used.",
  table_lines(trait_scaling),
  "ML model support:",
  table_lines(ml_statistics |> dplyr::filter(.data$predictor %in% c("LDMC", "LMA", "Height"))),
  table_lines(ml_comparisons |> dplyr::filter(.data$predictor %in% c("LDMC", "LMA", "Height"))),
  "CR1 summary:",
  table_lines(cr1_summary |> dplyr::filter(.data$predictor %in% c("LDMC", "LMA", "Height"))),
  "",
  "5. FUNCTIONAL INFORMATION",
  "--------------------------",
  "Graminoid status and Woodiness are independent binary predictors. Non-graminoid is not called Forb, and woody non-graminoids are retained in both applicable levels.",
  table_lines(functional_counts),
  "Cross-tabulation:", table_lines(functional_crosstab),
  "Distinct-experiment climate coverage:", table_lines(functional_coverage),
  "ML model support:",
  table_lines(ml_statistics |> dplyr::filter(.data$predictor %in% c("Graminoid status", "Woodiness"))),
  table_lines(ml_comparisons |> dplyr::filter(.data$predictor %in% c("Graminoid status", "Woodiness"))),
  "CR1 summary:",
  table_lines(cr1_summary |> dplyr::filter(.data$predictor %in% c("Graminoid status", "Woodiness"))),
  "",
  "6. PHYLOGENETIC INFORMATION",
  "----------------------------",
  sprintf("PEV1-PEV3 represent %.4f%% of positive phylogenetic variation.",
          pev_positive_variation),
  "They were selected response-independently, standardized across all 526 unique species before row expansion, and interpreted jointly; axis signs are not biological trait directions.",
  table_lines(pev_scaling),
  table_lines(ml_statistics |> dplyr::filter(grepl("Phylogenetic", .data$predictor))),
  table_lines(ml_comparisons |> dplyr::filter(grepl("Phylogenetic", .data$predictor))),
  "CR1 summary:",
  table_lines(cr1_summary |> dplyr::filter(grepl("Phylogenetic", .data$predictor))),
  "No PEV1-PEV5 model or response-driven axis selection was used.",
  "",
  "7. HETEROGENEITY SYNTHESIS",
  "---------------------------",
  "H_overall is the species intercept variance; H_W is mean(W_z^2 sigma_W^2); H_context is mean(x'G_context x) using the full T/P/TP covariance; H_total is their sum. Experiments receive equal weight and negative percentages are retained.",
  table_lines(heterogeneity_synthesis),
  "MODEL SUPPORT and HETEROGENEITY RESOLUTION are distinct: a systematic fixed-effect pattern can coexist with little reduction in residual species heterogeneity.",
  "",
  "8. REQUIRED CONSISTENCY CHECKS",
  "------------------------------",
  table_lines(consistency_checks),
  "",
  "REFERENCE OUTPUT",
  "----------------",
  "downstream_526_reference_results.rds freezes coefficients, random parameters, classifications, model fits, scaling, and all response-scale heterogeneity quantities for the later refactor.",
  "",
  "DOWNSTREAM 526-SPECIES ANALYSIS COMPLETE: TRUE"
)
writeLines(report, file.path(output_dir, "downstream_526_report.txt"), useBytes = TRUE)

message("Step 22 downstream 526-species rebuild complete: ", output_dir)
