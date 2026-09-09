# Step 16: isolated trait models using M_both as the baseline ---------------
#
# Fits LDMC, LMA, and vegetative height separately. Within each trait, all
# models use exactly the same complete-case rows, known sampling covariance,
# climate scaling, trait scaling, and M_both species random structure.
#
# ML:   M0 = T * P * W
#       M1 = Trait + T * P * W
#       M2 = Trait * T * P * W
# REML: M0 and M2 for final variance-component comparison; M2 for estimates.

source("config/analysis.R")
source("R/pipeline/utils.R")
source("R/pipeline/frequentist_hierarchical_auc.R")

required_packages <- c(
  "dplyr", "geometry", "ggplot2", "here", "Matrix", "metafor",
  "parallel", "readr", "tibble", "tidyr"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "),
       call. = FALSE)
}

traits <- c("LDMC", "LMA", "vegetative_height")
trait_labels <- c(
  LDMC = "LDMC", LMA = "LMA", vegetative_height = "Vegetative height"
)
model_dir <- file.path("output", "step_by_step_modeling", "species_both_traits")
figure_dir <- file.path("output", "frequentist_hierarchical_auc", "figures")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

formulas <- list(
  M0 = ~ T_z * P_z * W_z,
  M1 = ~ Trait_z + T_z * P_z * W_z,
  M2 = ~ Trait_z * T_z * P_z * W_z
)
random_structure_text <- paste(
  "list(~1 | species_id, ~0 + W_z | species_id,",
  "~0 + T_z + P_z + TP_z | species_id); struct = GEN"
)

# Reconstruct the same maximum-duration resident data and known sampling V
# used to fit M_both. No category/significance field enters the response.
scaled_auc <- scale_resident_auc_data(
  read_resident_auc_data(here::here(), horizons = "max")
)
auc_data <- scaled_auc$data |>
  dplyr::mutate(
    species_id = factor(as.character(.data$species_id)),
    experiment_id = factor(as.character(.data$experiment_id)),
    TP_z = .data$T_z * .data$P_z
  )
auc_by_plot <- readRDS(file.path("output", "intermediate", "auc_by_plot.rds"))
sampling <- reconstruct_auc_covariance(auc_data, auc_by_plot)
if (!isTRUE(sampling$audit$full_covariance_valid[[1]])) {
  stop("The reconstructed full Delta-AUC sampling covariance is invalid.",
       call. = FALSE)
}
V_auc <- sampling$V

prepare_trait <- function(trait) {
  prepared <- trait_complete_data(auc_data, trait)
  data <- prepared$data |>
    dplyr::mutate(
      Trait_z = .data$trait_z,
      species_id = droplevels(factor(as.character(.data$species_id))),
      experiment_id = droplevels(factor(as.character(.data$experiment_id))),
      TP_z = .data$T_z * .data$P_z
    )
  V <- subset_auc_v(V_auc, data$effect_id)
  validation <- validate_auc_covariance(as.matrix(V), data)
  if (!isTRUE(validation$valid[[1]]) || nrow(data) != length(data$effect_id) ||
      anyDuplicated(data$effect_id)) {
    stop(trait, ": trait rows or covariance failed exact alignment.",
         call. = FALSE)
  }
  design_audit <- lapply(formulas, function(formula) {
    X <- stats::model.matrix(formula, data)
    c(rank = qr(X)$rank, columns = ncol(X), condition_number = kappa(X))
  })
  if (any(vapply(design_audit, function(x) x[["rank"]] != x[["columns"]],
                 logical(1)))) {
    stop(trait, ": at least one fixed-effect design is rank-deficient.",
         call. = FALSE)
  }
  list(
    data = data, V = V, scaling = prepared$scaling,
    design_audit = design_audit, covariance_validation = validation
  )
}

prepared_traits <- stats::setNames(lapply(traits, prepare_trait), traits)
trait_audit <- dplyr::bind_rows(lapply(traits, function(trait) {
  item <- prepared_traits[[trait]]
  tibble::tibble(
    trait = trait,
    n_effects = nrow(item$data),
    n_species = dplyr::n_distinct(item$data$species_id),
    n_experiments = dplyr::n_distinct(item$data$experiment_id),
    first_effect_id = item$data$effect_id[[1]],
    last_effect_id = item$data$effect_id[[nrow(item$data)]],
    trait_transform = item$scaling$transform,
    trait_center_log_scale = item$scaling$center,
    trait_scale_log_scale = item$scaling$scale,
    trait_scaling_n_species = item$scaling$n,
    trait_scaling_weighting = item$scaling$weighting,
    V_order_exact = identical(rownames(item$V), item$data$effect_id),
    V_symmetric = max(abs(item$V - Matrix::t(item$V))) <= 1e-8,
    M0_design_rank = item$design_audit$M0[["rank"]],
    M0_design_columns = item$design_audit$M0[["columns"]],
    M1_design_rank = item$design_audit$M1[["rank"]],
    M1_design_columns = item$design_audit$M1[["columns"]],
    M2_design_rank = item$design_audit$M2[["rank"]],
    M2_design_columns = item$design_audit$M2[["columns"]]
  )
}))
invisible(write_csv_atomic(
  trait_audit, file.path(model_dir, "trait_row_covariance_scaling_audit.csv")
))

fit_signature <- function(data, formula, method) {
  list(
    effect_id = data$effect_id,
    Trait_z = data$Trait_z,
    formula = paste(deparse(formula), collapse = " "),
    method = method,
    random_structure = random_structure_text
  )
}

fit_cached <- function(data, V, formula, method, cache_file) {
  signature <- fit_signature(data, formula, method)
  if (file.exists(cache_file)) {
    cached <- readRDS(cache_file)
    if (is.list(cached) && !is.null(cached$fit) &&
        identical(cached$signature, signature)) {
      return(cached)
    }
  }
  captured <- auc_capture(
    metafor::rma.mv(
      yi = auc_effect,
      V = V,
      mods = formula,
      random = list(
        ~ 1 | species_id,
        ~ 0 + W_z | species_id,
        ~ 0 + T_z + P_z + TP_z | species_id
      ),
      struct = "GEN",
      data = data,
      method = method,
      sparse = TRUE
    )
  )
  if (is.null(captured$fit)) {
    stop("Trait model failed [", basename(cache_file), "]: ",
         captured$error, call. = FALSE)
  }
  result <- list(
    fit = captured$fit, warnings = captured$warnings,
    error = captured$error, signature = signature
  )
  invisible(write_rds_atomic(result, cache_file))
  result
}

model_diagnostics <- function(fit_result, trait, model_id, method) {
  fit <- fit_result$fit
  intended_blocks <- identical(colnames(fit$G), "W_z") &&
    identical(colnames(fit$H), c("T_z", "P_z", "TP_z"))
  tibble::tibble(
    trait = trait, model = model_id, method = method,
    n_effects = fit$k, n_fixed = fit$p,
    optimizer_convergence = auc_value_or(fit$opt.res$convergence, NA_integer_),
    finite_coefficients = all(is.finite(stats::coef(fit))),
    finite_fixed_covariance = all(is.finite(stats::vcov(fit))),
    intended_random_blocks = intended_blocks,
    warnings = paste(fit_result$warnings, collapse = " | "),
    error = fit_result$error
  )
}

fit_stats <- function(fit, trait, model_id) {
  statistics <- auc_fit_statistics(fit)
  tibble::tibble(
    trait = trait, model = model_id, method = "ML",
    n_effects = fit$k, n_fixed = fit$p,
    logLik = statistics[["logLik"]], AIC = statistics[["AIC"]],
    AICc = statistics[["AICc"]], BIC = statistics[["BIC"]]
  )
}

nested_lrt <- function(smaller, larger, trait, comparison) {
  ll_small <- stats::logLik(smaller)
  ll_large <- stats::logLik(larger)
  df_small <- attr(ll_small, "df")
  df_large <- attr(ll_large, "df")
  df_difference <- df_large - df_small
  statistic <- 2 * (as.numeric(ll_large) - as.numeric(ll_small))
  tibble::tibble(
    trait = trait, comparison = comparison,
    smaller_model = sub(" vs.*", "", comparison),
    larger_model = sub(".*vs ", "", comparison),
    likelihood_ratio = statistic, df = df_difference,
    p_value = stats::pchisq(statistic, df_difference, lower.tail = FALSE),
    same_rows_and_V = TRUE, method = "ML",
    valid_nested_comparison = is.finite(statistic) && statistic >= -1e-7 &&
      df_difference > 0
  )
}

variance_components_both <- function(fit, trait, model_id) {
  values <- c(
    species_intercept = as.numeric(fit$sigma2[[1]]),
    species_W_slope = as.numeric(fit$tau2[[1]]),
    species_T = as.numeric(fit$gamma2[[1]]),
    species_P = as.numeric(fit$gamma2[[2]]),
    species_TP = as.numeric(fit$gamma2[[3]])
  )
  tibble::tibble(
    trait = trait, model = model_id,
    component = names(values), variance = as.numeric(values)
  )
}

fit_trait_models <- function(trait) {
  message("Fitting isolated M_both trait sequence: ", trait)
  item <- prepared_traits[[trait]]
  trait_dir <- file.path(model_dir, trait)
  dir.create(trait_dir, recursive = TRUE, showWarnings = FALSE)
  specifications <- tibble::tribble(
    ~cache_name, ~model_id, ~method,
    "M0_ML", "M0", "ML",
    "M1_ML", "M1", "ML",
    "M2_ML", "M2", "ML",
    "M0_REML", "M0", "REML",
    "M2_REML", "M2", "REML"
  )
  fitted <- list()
  diagnostics <- list()
  for (index in seq_len(nrow(specifications))) {
    specification <- specifications[index, ]
    cache_name <- specification$cache_name[[1]]
    model_id <- specification$model_id[[1]]
    method <- specification$method[[1]]
    message("  ", trait, " ", cache_name)
    fitted[[cache_name]] <- fit_cached(
      item$data, item$V, formulas[[model_id]], method,
      file.path(trait_dir, paste0(cache_name, ".rds"))
    )
    diagnostics[[cache_name]] <- model_diagnostics(
      fitted[[cache_name]], trait, model_id, method
    )
  }

  ml_statistics <- dplyr::bind_rows(
    fit_stats(fitted$M0_ML$fit, trait, "M0"),
    fit_stats(fitted$M1_ML$fit, trait, "M1"),
    fit_stats(fitted$M2_ML$fit, trait, "M2")
  )
  comparisons <- dplyr::bind_rows(
    nested_lrt(
      fitted$M0_ML$fit, fitted$M1_ML$fit, trait, "M0 vs M1"
    ),
    nested_lrt(
      fitted$M1_ML$fit, fitted$M2_ML$fit, trait, "M1 vs M2"
    )
  )

  variances <- dplyr::bind_rows(
    variance_components_both(fitted$M0_REML$fit, trait, "M0_REML"),
    variance_components_both(fitted$M2_REML$fit, trait, "M2_REML")
  )
  variance_reduction <- variances |>
    tidyr::pivot_wider(names_from = "model", values_from = "variance") |>
    dplyr::mutate(
      percent_reduction = 100 * (.data$M0_REML - .data$M2_REML) /
        .data$M0_REML,
      interpretation = "descriptive variance change; not a formal R-squared"
    )

  # CR1 is an inference calculation on the final REML M2, not a refit.
  cr1_warnings <- character()
  M2_cr1 <- withCallingHandlers(
    metafor::robust(
      fitted$M2_REML$fit,
      cluster = fitted$M2_REML$fit$data$experiment_id,
      adjust = TRUE, clubSandwich = FALSE
    ),
    warning = function(condition) {
      cr1_warnings <<- c(cr1_warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  invisible(write_rds_atomic(
    M2_cr1, file.path(trait_dir, "M2_REML_experiment_CR1.rds")
  ))
  ordinary <- fitted$M2_REML$fit
  coefficient_terms <- names(stats::coef(ordinary))
  coefficient_table <- tibble::tibble(
    trait = trait,
    term = coefficient_terms,
    trait_related = grepl("Trait_z", coefficient_terms, fixed = TRUE),
    estimate = as.numeric(stats::coef(ordinary)),
    model_standard_error = as.numeric(ordinary$se),
    model_conf_low = as.numeric(ordinary$ci.lb),
    model_conf_high = as.numeric(ordinary$ci.ub),
    model_p_value = as.numeric(ordinary$pval),
    CR1_standard_error = as.numeric(M2_cr1$se),
    CR1_conf_low = as.numeric(M2_cr1$ci.lb),
    CR1_conf_high = as.numeric(M2_cr1$ci.ub),
    CR1_p_value = as.numeric(M2_cr1$pval),
    CR1_denominator_df = rep(
      as.numeric(M2_cr1$ddf), length.out = length(stats::coef(ordinary))
    ),
    CR1_warnings = paste(unique(cr1_warnings), collapse = " | ")
  )
  if (!isTRUE(all.equal(
    as.numeric(stats::coef(ordinary)), as.numeric(M2_cr1$beta),
    tolerance = sqrt(.Machine$double.eps)
  ))) {
    stop(trait, ": CR1 unexpectedly changed M2 coefficients.", call. = FALSE)
  }

  list(
    audit = dplyr::bind_rows(diagnostics),
    ml_statistics = ml_statistics, comparisons = comparisons,
    variances = variances, variance_reduction = variance_reduction,
    coefficients = coefficient_table
  )
}

requested_cores <- suppressWarnings(as.integer(
  Sys.getenv("TRANSPLANT_TRAIT_CORES", "3")
))
if (!is.finite(requested_cores) || requested_cores < 1L) requested_cores <- 1L
trait_cores <- min(length(traits), requested_cores)
if (.Platform$OS.type == "windows") trait_cores <- 1L
results <- if (trait_cores > 1L) {
  parallel::mclapply(traits, fit_trait_models, mc.cores = trait_cores,
                     mc.preschedule = FALSE)
} else {
  lapply(traits, fit_trait_models)
}
names(results) <- traits
if (any(vapply(results, inherits, logical(1), what = "try-error"))) {
  stop("At least one trait sequence failed; inspect worker output.", call. = FALSE)
}

model_audit <- dplyr::bind_rows(lapply(results, `[[`, "audit"))
ml_statistics <- dplyr::bind_rows(lapply(results, `[[`, "ml_statistics"))
comparisons <- dplyr::bind_rows(lapply(results, `[[`, "comparisons"))
variance_components <- dplyr::bind_rows(lapply(results, `[[`, "variances"))
variance_reduction <- dplyr::bind_rows(lapply(results, `[[`, "variance_reduction"))
coefficient_table <- dplyr::bind_rows(lapply(results, `[[`, "coefficients"))

tables <- list(
  trait_model_diagnostics = model_audit,
  trait_ml_model_fit_statistics = ml_statistics,
  trait_ml_nested_comparisons = comparisons,
  trait_REML_variance_components = variance_components,
  trait_REML_variance_reduction = variance_reduction,
  trait_M2_REML_fixed_effects_model_vs_CR1 = coefficient_table
)
for (name in names(tables)) {
  invisible(write_csv_atomic(
    tables[[name]], file.path(model_dir, paste0(name, ".csv"))
  ))
}

# High-minus-low trait contrasts over the same 3D observed climate support.
support <- build_auc_climate_grid(auc_data, points = 80L)
if (!identical(unique(support$grid$support_method),
               "three-dimensional convex hull")) {
  stop("Trait surfaces require the three-dimensional climate hull.",
       call. = FALSE)
}
climate_grid <- support$grid |>
  dplyr::filter(.data$supported)

scaling_value <- function(field, statistic) {
  value <- scaled_auc$scaling[[statistic]][
    scaled_auc$scaling$scaled_field == field
  ]
  if (length(value) != 1L || !is.finite(value)) {
    stop("Missing scaling value for ", field, call. = FALSE)
  }
  value[[1]]
}
T_center <- scaling_value("T_z", "center")
T_scale <- scaling_value("T_z", "scale")
P_center <- scaling_value("P_z", "center")
P_scale <- scaling_value("P_z", "scale")
W_center <- scaling_value("W_z", "center")
W_scale <- scaling_value("W_z", "scale")

warming_experiments <- auc_data |>
  dplyr::distinct(.data$experiment_id, .data$W_z, .data$W)
warming_z_values <- as.numeric(stats::quantile(
  warming_experiments$W_z, c(0.25, 0.50, 0.75)
))
warming_C_values <- as.numeric(stats::quantile(
  warming_experiments$W, c(0.25, 0.50, 0.75)
))
warming_labels <- paste0(
  c("Low warming", "Median warming", "High warming"),
  "\nDelta T = ", sprintf("%.2f", warming_C_values), " °C"
)
climate_grid <- climate_grid |>
  dplyr::mutate(
    warming_panel = factor(
      warming_labels[match(.data$W_z, warming_z_values)],
      levels = warming_labels
    ),
    destination_T_JJA = .data$T_z * T_scale + T_center,
    destination_P_JJA = exp(.data$P_z * P_scale + P_center),
    warming_C = .data$W_z * W_scale + W_center,
    TP_z = .data$T_z * .data$P_z
  )

trait_surface <- function(trait) {
  trait_dir <- file.path(model_dir, trait)
  model <- readRDS(file.path(trait_dir, "M2_REML.rds"))$fit
  robust <- readRDS(file.path(trait_dir, "M2_REML_experiment_CR1.rds"))
  high <- climate_grid |>
    dplyr::mutate(Trait_z = 1)
  low <- climate_grid |>
    dplyr::mutate(Trait_z = -1)
  X_high <- auc_design_matrix(formulas$M2, high, names(stats::coef(model)))
  X_low <- auc_design_matrix(formulas$M2, low, names(stats::coef(model)))
  contrast <- X_high - X_low
  estimate <- as.numeric(contrast %*% stats::coef(model))
  model_vcov <- as.matrix(stats::vcov(model))
  robust_vcov <- as.matrix(robust$vb)
  model_se <- sqrt(pmax(
    0, rowSums((contrast %*% model_vcov) * contrast)
  ))
  CR1_se <- sqrt(pmax(
    0, rowSums((contrast %*% robust_vcov) * contrast)
  ))
  CR1_df <- unique(as.numeric(robust$ddf))[[1]]
  climate_grid |>
    dplyr::mutate(
      trait = trait,
      contrast = "Trait_z +1 SD minus Trait_z -1 SD",
      high_minus_low_delta_auc = estimate,
      model_standard_error = model_se,
      model_conf_low = estimate - stats::qnorm(0.975) * model_se,
      model_conf_high = estimate + stats::qnorm(0.975) * model_se,
      CR1_standard_error = CR1_se,
      CR1_conf_low = estimate - stats::qt(0.975, CR1_df) * CR1_se,
      CR1_conf_high = estimate + stats::qt(0.975, CR1_df) * CR1_se,
      CR1_denominator_df = CR1_df,
      interpretation = paste(
        "Positive means high trait is more winner-like; negative means high",
        "trait is more loser-like; species random effects are zero."
      )
    )
}

trait_grids <- dplyr::bind_rows(lapply(traits, trait_surface))
invisible(write_csv_atomic(
  trait_grids,
  file.path(model_dir, "trait_high_minus_low_prediction_grids.csv")
))
for (trait in traits) {
  invisible(write_csv_atomic(
    dplyr::filter(trait_grids, .data$trait == .env$trait),
    file.path(model_dir, trait,
              paste0(trait, "_high_minus_low_prediction_grid.csv"))
  ))
}

experiments <- auc_data |>
  dplyr::distinct(.data$experiment_id, .data$T_z, .data$P_z)
raw_T_labels <- function(x) sprintf("%.1f", x * T_scale + T_center)
raw_P_labels <- function(x) sprintf("%.0f", exp(x * P_scale + P_center))
global_limit <- max(abs(trait_grids$high_minus_low_delta_auc), na.rm = TRUE)

plot_trait_surface <- function(trait) {
  grid <- dplyr::filter(trait_grids, .data$trait == .env$trait)
  zero_contour <- grid |>
    dplyr::summarise(
      crosses_zero = min(.data$high_minus_low_delta_auc) <= 0 &&
        max(.data$high_minus_low_delta_auc) >= 0,
      .by = "warming_panel"
    ) |>
    dplyr::filter(.data$crosses_zero) |>
    dplyr::select("warming_panel") |>
    dplyr::inner_join(grid, by = "warming_panel")
  plot <- ggplot2::ggplot(
    grid,
    ggplot2::aes(
      x = .data$T_z, y = .data$P_z,
      fill = .data$high_minus_low_delta_auc
    )
  ) +
    ggplot2::geom_raster() +
    ggplot2::geom_contour(
      data = zero_contour,
      ggplot2::aes(
        x = .data$T_z, y = .data$P_z,
        z = .data$high_minus_low_delta_auc
      ),
      inherit.aes = FALSE, breaks = 0, colour = "grey15", linewidth = 0.4
    ) +
    ggplot2::geom_point(
      data = experiments,
      ggplot2::aes(x = .data$T_z, y = .data$P_z),
      inherit.aes = FALSE, shape = 21, fill = "white", colour = "grey25",
      size = 1.3, stroke = 0.35
    ) +
    ggplot2::facet_wrap(~ warming_panel, nrow = 1) +
    ggplot2::scale_fill_gradient2(
      low = "#E66101", mid = "white", high = "#1B9E77", midpoint = 0,
      limits = c(-global_limit, global_limit)
    ) +
    ggplot2::scale_x_continuous(labels = raw_T_labels) +
    ggplot2::scale_y_continuous(labels = raw_P_labels) +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    ) +
    ggplot2::labs(
      title = paste0(
        trait_labels[[trait]], ": high-trait minus low-trait Delta AUC"
      ),
      subtitle = paste0(
        "+1 SD minus -1 SD for the hypothetical average species; all species ",
        "random effects are zero; 3D observed climate support only."
      ),
      x = "Destination JJA temperature (°C)",
      y = "Destination JJA precipitation",
      fill = "High - low\nDelta AUC"
    )
  stem <- file.path(
    figure_dir,
    paste0("step_16_trait_", trait, "_high_minus_low_surface")
  )
  ggplot2::ggsave(
    paste0(stem, ".pdf"), plot, width = 11.5, height = 4.8,
    device = grDevices::cairo_pdf
  )
  ggplot2::ggsave(
    paste0(stem, ".png"), plot, width = 11.5, height = 4.8,
    dpi = 220, bg = "white"
  )
  stem
}
figure_stems <- stats::setNames(lapply(traits, plot_trait_surface), traits)

surface_summary <- trait_grids |>
  dplyr::summarise(
    n_supported_grid_cells = dplyr::n(),
    min_high_minus_low_delta_auc = min(.data$high_minus_low_delta_auc),
    max_high_minus_low_delta_auc = max(.data$high_minus_low_delta_auc),
    mean_high_minus_low_delta_auc = mean(.data$high_minus_low_delta_auc),
    crosses_zero = min(.data$high_minus_low_delta_auc) <= 0 &&
      max(.data$high_minus_low_delta_auc) >= 0,
    .by = c("trait", "warming_panel", "warming_C")
  )
invisible(write_csv_atomic(
  surface_summary, file.path(model_dir, "trait_surface_summary.csv")
))

cat("\nIsolated M_both trait step completed.\n")
cat("Random structure unchanged: ", random_structure_text, "\n", sep = "")
cat("Trait model output: ", normalizePath(model_dir), "\n", sep = "")
cat("Figures:\n")
for (trait in traits) {
  cat("  ", trait, ": ", normalizePath(paste0(figure_stems[[trait]], ".pdf")),
      "\n", sep = "")
}
