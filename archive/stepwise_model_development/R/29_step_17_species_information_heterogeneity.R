# Step 17: how much species information resolves response heterogeneity -------
#
# All headline heterogeneity quantities are among-species variances in
# predicted Delta AUC, averaged with equal weight over the 42 distinct observed
# experiment climates. Each predictor is compared with an M_both baseline fit
# to exactly the same effects.

source("config/analysis.R")
source("R/pipeline/utils.R")
source("R/pipeline/frequentist_hierarchical_auc.R")
source("R/pipeline/additional_analysis.R")

required_packages <- c(
  "dplyr", "geometry", "ggplot2", "here", "Matrix", "metafor",
  "readr", "tibble", "tidyr"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "),
       call. = FALSE)
}

output_dir <- file.path(
  "output", "step_by_step_modeling", "species_information_heterogeneity"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_files <- list(
  headline = file.path(output_dir, "heterogeneity_resolved_headline.csv"),
  detailed = file.path(output_dir, "random_effect_response_heterogeneity.csv"),
  functional_audit = file.path(output_dir, "functional_status_audit.csv"),
  functional_assignments = file.path(
    output_dir, "resident_functional_status_assignments.csv"
  ),
  family_audit = file.path(output_dir, "family_audit.csv"),
  synthesis_figure = file.path(output_dir, "heterogeneity_resolved_synthesis.pdf"),
  trait_LDMC = file.path(output_dir, "trait_LDMC_response_surfaces.pdf"),
  trait_LMA = file.path(output_dir, "trait_LMA_response_surfaces.pdf"),
  trait_height = file.path(
    output_dir, "trait_vegetative_height_response_surfaces.pdf"
  ),
  graminoid_figure = file.path(
    output_dir, "graminoid_status_response_surfaces.pdf"
  ),
  woody_figure = file.path(
    output_dir, "woody_status_response_surfaces.pdf"
  ),
  family_complete = file.path(
    output_dir, "family_response_surfaces_all_families.pdf"
  ),
  family_clean = file.path(
    output_dir, "family_response_surfaces_at_least_5_species.pdf"
  ),
  report = file.path(output_dir, "species_information_heterogeneity_report.txt")
)
temporary_cache_dir <- file.path(
  "/private", "tmp", "transplant_species_information_heterogeneity_cache"
)
dir.create(temporary_cache_dir, recursive = TRUE, showWarnings = FALSE)

load_or_fit_temporary <- function(cache_name, effect_ids, fitter) {
  cache_path <- file.path(temporary_cache_dir, paste0(cache_name, ".rds"))
  if (file.exists(cache_path)) {
    cached <- readRDS(cache_path)
    if (is.list(cached) && !is.null(cached$fit) &&
        identical(as.character(cached$fit$data$effect_id), effect_ids)) {
      return(cached)
    }
  }
  fitted <- fitter()
  saveRDS(fitted, cache_path)
  fitted
}

# Common maximum-duration resident response, scaling, effect ordering, and V.
scaled_auc <- scale_resident_auc_data(
  read_resident_auc_data(here::here(), horizons = "max")
)
auc_data <- scaled_auc$data |>
  dplyr::mutate(
    TP_z = .data$T_z * .data$P_z
  )
auc_by_plot <- readRDS(file.path("output", "intermediate", "auc_by_plot.rds"))
sampling <- reconstruct_auc_covariance(auc_data, auc_by_plot)
if (!isTRUE(sampling$audit$full_covariance_valid[[1]])) {
  stop("The reconstructed full Delta-AUC sampling covariance is invalid.",
       call. = FALSE)
}
V_auc <- sampling$V
if (!identical(rownames(V_auc), auc_data$effect_id)) {
  stop("V_auc is not in exact effect_id order.", call. = FALSE)
}
auc_data <- auc_data |>
  dplyr::mutate(
    species_id = as.character(.data$species_id),
    experiment_id = as.character(.data$experiment_id)
  )

experiment_climates <- auc_data |>
  dplyr::distinct(
    .data$experiment_id, .data$T_z, .data$P_z, .data$W_z
  ) |>
  dplyr::mutate(TP_z = .data$T_z * .data$P_z) |>
  dplyr::arrange(.data$experiment_id)
if (nrow(experiment_climates) != 42L) {
  stop("Expected exactly 42 distinct observed experiment climates; found ",
       nrow(experiment_climates), ".", call. = FALSE)
}

climate_formula <- ~ T_z * P_z * W_z
trait_formula <- ~ Trait_z * T_z * P_z * W_z
graminoid_formula <- ~ GraminoidStatus * T_z * P_z * W_z
woody_formula <- ~ WoodyStatus * T_z * P_z * W_z
species_random <- list(
  ~ 1 | species_id,
  ~ 0 + W_z | species_id,
  ~ 0 + T_z + P_z + TP_z | species_id
)

# A five-coefficient GEN block with structural zero correlations is exactly a
# separate intercept, separate W slope, and correlated T/P/TP block. Metafor
# allows at most two inner|outer formulas, so the family model uses one such
# block for family and one for species.
independent_block_correlations <- c(rep(0, 7), NA, NA, NA)
family_species_random <- list(
  ~ 1 + W_z + T_z + P_z + TP_z | family_id,
  ~ 1 + W_z + T_z + P_z + TP_z | species_id
)

capture_fit <- function(expression, label) {
  captured <- auc_capture(expression)
  if (is.null(captured$fit)) {
    stop(label, " failed: ", captured$error, call. = FALSE)
  }
  captured
}

fit_species_model <- function(data, V, formula, label) {
  data <- data |>
    dplyr::mutate(
      species_id = droplevels(factor(.data$species_id)),
      experiment_id = droplevels(factor(.data$experiment_id))
    )
  capture_fit(
    metafor::rma.mv(
      yi = auc_effect, V = V, mods = formula,
      random = species_random, struct = "GEN", data = data,
      method = "REML", sparse = TRUE
    ),
    label
  )
}

fit_family_model <- function(data, V) {
  data <- data |>
    dplyr::mutate(
      family_id = droplevels(factor(.data$family_id)),
      species_id = droplevels(factor(.data$species_id)),
      experiment_id = droplevels(factor(.data$experiment_id))
    )
  capture_fit(
    metafor::rma.mv(
      yi = auc_effect, V = V, mods = climate_formula,
      random = family_species_random, struct = c("GEN", "GEN"),
      rho = independent_block_correlations,
      phi = independent_block_correlations,
      data = data, method = "REML", sparse = TRUE,
      control = list(optimizer = "bobyqa", maxfun = 100000)
    ),
    "Family hierarchical model"
  )
}

safe_correlation <- function(covariance, first, second) {
  denominator <- sqrt(covariance[first, first] * covariance[second, second])
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  covariance[first, second] / denominator
}

heterogeneity_metrics <- function(intercept_variance, W_variance,
                                  context_covariance) {
  context_design <- as.matrix(
    experiment_climates[, c("T_z", "P_z", "TP_z")]
  )
  context_by_experiment <- rowSums(
    (context_design %*% context_covariance) * context_design
  )
  H_overall <- as.numeric(intercept_variance)
  H_W <- mean(experiment_climates$W_z^2 * as.numeric(W_variance))
  H_context <- mean(context_by_experiment)
  c(
    H_overall = H_overall,
    H_W = H_W,
    H_context = H_context,
    H_total = H_overall + H_W + H_context
  )
}

diagnostics <- function(fit_capture, intercept_variance, W_variance,
                        covariance, label) {
  fit <- fit_capture$fit
  values <- eigen(covariance, symmetric = TRUE, only.values = TRUE)$values
  boundary_tolerance <- 1e-8 * max(1, stats::var(fit$yi))
  optimizer_code <- auc_value_or(
    fit$opt.res$convergence, auc_value_or(fit$opt.res$ierr, NA_integer_)
  )
  optimizer <- if (!is.null(fit$opt.res$ierr)) {
    "bobyqa"
  } else if (!is.null(fit$opt.res$convergence)) {
    "nlminb"
  } else {
    class(fit$opt.res)[[1]]
  }
  correlations <- stats::cov2cor(covariance)
  max_absolute_correlation <- max(abs(correlations[upper.tri(correlations)]))
  any_boundary <- any(
    c(intercept_variance, W_variance, diag(covariance)) <= boundary_tolerance
  )
  positive_definite <- all(is.finite(values)) && min(values) > 0
  list(
    label = label,
    optimizer = optimizer,
    optimizer_code = optimizer_code,
    warnings = paste(fit_capture$warnings, collapse = " | "),
    finite = all(is.finite(stats::coef(fit))) &&
      all(is.finite(stats::vcov(fit))) && all(is.finite(covariance)),
    min_eigenvalue = min(values),
    positive_definite = positive_definite,
    max_absolute_correlation = max_absolute_correlation,
    boundary = any_boundary,
    identifiable = positive_definite && !any_boundary &&
      max_absolute_correlation < 1 - 1e-6
  )
}

detail_row <- function(predictor, role, level, fit_capture,
                       intercept_variance, W_variance, context_covariance,
                       n_effects, n_species) {
  fit <- fit_capture$fit
  metrics <- heterogeneity_metrics(
    intercept_variance, W_variance, context_covariance
  )
  diagnostic <- diagnostics(
    fit_capture, intercept_variance, W_variance, context_covariance,
                            paste(predictor, role, level))
  tibble::tibble(
    species_information = predictor,
    model_role = role,
    heterogeneity_level = level,
    n_effects = n_effects,
    n_species = n_species,
    n_experiments = nrow(experiment_climates),
    species_intercept_variance = as.numeric(intercept_variance),
    species_W_slope_variance = as.numeric(W_variance),
    species_T_variance = context_covariance["T_z", "T_z"],
    species_P_variance = context_covariance["P_z", "P_z"],
    species_TP_variance = context_covariance["TP_z", "TP_z"],
    T_P_correlation = safe_correlation(context_covariance, "T_z", "P_z"),
    T_TP_correlation = safe_correlation(context_covariance, "T_z", "TP_z"),
    P_TP_correlation = safe_correlation(context_covariance, "P_z", "TP_z"),
    H_overall = metrics[["H_overall"]],
    H_W = metrics[["H_W"]],
    H_context = metrics[["H_context"]],
    H_total = metrics[["H_total"]],
    optimizer = diagnostic$optimizer,
    optimizer_code = diagnostic$optimizer_code,
    warnings = diagnostic$warnings,
    finite_fit = diagnostic$finite,
    context_min_eigenvalue = diagnostic$min_eigenvalue,
    context_positive_definite = diagnostic$positive_definite,
    max_absolute_context_correlation = diagnostic$max_absolute_correlation,
    any_random_variance_at_boundary = diagnostic$boundary,
    random_structure_identifiable = diagnostic$identifiable
  )
}

extract_separate_species <- function(fit_capture, predictor, role, data) {
  fit <- fit_capture$fit
  if (!identical(colnames(fit$G), "W_z") ||
      !identical(colnames(fit$H), c("T_z", "P_z", "TP_z"))) {
    stop(predictor, " ", role, ": unexpected species random structure.",
         call. = FALSE)
  }
  detail_row(
    predictor, role, "species", fit_capture,
    fit$sigma2[[1]], fit$tau2[[1]], as.matrix(fit$H),
    nrow(data), dplyr::n_distinct(data$species_id)
  )
}

percent_resolved_row <- function(predictor, baseline_row, predictor_row) {
  metric_names <- c("H_overall", "H_W", "H_context", "H_total")
  baseline <- unlist(baseline_row[1, metric_names], use.names = TRUE)
  residual <- unlist(predictor_row[1, metric_names], use.names = TRUE)
  resolved <- 100 * (baseline - residual) / baseline
  tibble::tibble(
    `Species information` = predictor,
    `Overall response heterogeneity resolved (%)` = resolved[["H_overall"]],
    `Warming-related response heterogeneity resolved (%)` = resolved[["H_W"]],
    `Background-climate response heterogeneity resolved (%)` =
      resolved[["H_context"]],
    `Total response heterogeneity resolved (%)` = resolved[["H_total"]]
  )
}

# Taxonomic and functional-status mapping: one record per accepted species.
# The two binary statuses are derived independently and are both available for
# every angiosperm retained by the downstream analysis filter.
taxonomy_source <- readr::read_csv(
  file.path("output", "tables", "analysis_dataset.csv"),
  show_col_types = FALSE
) |>
  dplyr::transmute(
    accepted_name = as.character(.data$accepted_name),
    family = trimws(as.character(.data$family))
  ) |>
  dplyr::distinct()
functional_overrides <- readr::read_csv(
  file.path("config", "functional_group_overrides.csv"),
  show_col_types = FALSE,
  col_types = readr::cols(.default = readr::col_character())
)
taxonomy <- derive_functional_groups(taxonomy_source, functional_overrides) |>
  dplyr::transmute(
    species_id = as.character(.data$accepted_name),
    family_id = trimws(as.character(.data$family)),
    graminoid_status = as.character(.data$graminoid_status),
    woody_status = as.character(.data$woody_status),
    graminoid_review_basis = as.character(.data$graminoid_review_basis),
    woody_review_basis = as.character(.data$woody_review_basis)
  )
if (anyDuplicated(taxonomy$species_id)) {
  stop("Taxonomy/function mapping is not unique by accepted species.",
       call. = FALSE)
}
auc_taxonomy <- auc_data |>
  dplyr::left_join(taxonomy, by = "species_id")
functional_assignments <- auc_taxonomy |>
  dplyr::distinct(
    .data$species_id, .data$family_id, .data$graminoid_status,
    .data$woody_status, .data$graminoid_review_basis,
    .data$woody_review_basis
  ) |>
  dplyr::arrange(.data$species_id)
invisible(write_csv_atomic(
  functional_assignments, output_files$functional_assignments
))

status_audit_one <- function(data, variable, levels) {
  data |>
    dplyr::mutate(
      status = dplyr::coalesce(as.character(.data[[variable]]), "missing")
    ) |>
    dplyr::summarise(
      n_effects = dplyr::n(),
      n_species = dplyr::n_distinct(.data$species_id),
      n_experiments = dplyr::n_distinct(.data$experiment_id),
      .by = "status"
    ) |>
    tidyr::complete(
      status = c(levels, "missing"),
      fill = list(n_effects = 0L, n_species = 0L, n_experiments = 0L)
    ) |>
    dplyr::mutate(
      status_variable = variable,
      included_in_model = .data$status %in% .env$levels,
      handling = dplyr::if_else(
        .data$included_in_model,
        "included as its own binary level",
        "excluded because status is missing"
      ),
      .before = 1L
    )
}
functional_status_audit <- dplyr::bind_rows(
  status_audit_one(
    auc_taxonomy, "graminoid_status", c("non_graminoid", "graminoid")
  ),
  status_audit_one(
    auc_taxonomy, "woody_status", c("non_woody", "woody")
  )
)
invisible(write_csv_atomic(
  functional_status_audit, output_files$functional_audit
))

family_known <- !is.na(auc_taxonomy$family_id) & auc_taxonomy$family_id != ""
family_counts <- auc_taxonomy |>
  dplyr::filter(.env$family_known) |>
  dplyr::summarise(
    n_species = dplyr::n_distinct(.data$species_id),
    n_species_experiment_effects = dplyr::n(),
    n_experiments = dplyr::n_distinct(.data$experiment_id),
    .by = "family_id"
  ) |>
  dplyr::arrange(dplyr::desc(.data$n_species), .data$family_id) |>
  dplyr::mutate(
    included_in_model = TRUE,
    included_in_cleaner_figure = .data$n_species >= 5L,
    handling = dplyr::if_else(
      .data$n_species >= 5L,
      "retained in model and cleaner figure",
      "retained in model; omitted only from >=5-species figure"
    )
  )
missing_family_audit <- tibble::tibble(
  family_id = "Missing / unknown",
  n_species = dplyr::n_distinct(auc_taxonomy$species_id[!family_known]),
  n_species_experiment_effects = sum(!family_known),
  n_experiments = dplyr::n_distinct(
    auc_taxonomy$experiment_id[!family_known]
  ),
  included_in_model = FALSE,
  included_in_cleaner_figure = FALSE,
  handling = "excluded; no Unknown family was created"
)
family_audit <- dplyr::bind_rows(family_counts, missing_family_audit)
invisible(write_csv_atomic(family_audit, output_files$family_audit))

# Trait fits were already fitted in Step 16 with exactly the requested REML
# M0/M2 formulas, rows, V ordering, scaling, and random structure. Reuse those
# immutable fitted objects after verifying every effect key.
traits <- c("LDMC", "LMA", "vegetative_height")
trait_labels <- c(
  LDMC = "LDMC", LMA = "LMA", vegetative_height = "Vegetative height"
)
trait_model_dir <- file.path(
  "output", "step_by_step_modeling", "species_both_traits"
)

load_trait_pair <- function(trait) {
  prepared <- trait_complete_data(auc_data, trait)
  data <- prepared$data |>
    dplyr::mutate(
      Trait_z = .data$trait_z,
      species_id = factor(as.character(.data$species_id)),
      experiment_id = factor(as.character(.data$experiment_id)),
      TP_z = .data$T_z * .data$P_z
    )
  V <- subset_auc_v(V_auc, data$effect_id)
  model_path <- file.path(trait_model_dir, trait)
  baseline_capture <- readRDS(file.path(model_path, "M0_REML.rds"))
  predictor_capture <- readRDS(file.path(model_path, "M2_REML.rds"))
  for (object in list(baseline_capture, predictor_capture)) {
    if (is.null(object$fit) || object$fit$k != nrow(data) ||
        !identical(as.character(object$fit$data$effect_id), data$effect_id)) {
      stop(trait, ": saved fit does not match the exact trait row set.",
           call. = FALSE)
    }
  }
  if (!identical(rownames(V), data$effect_id)) {
    stop(trait, ": V subset ordering is not exact.", call. = FALSE)
  }
  list(
    data = data, V = V, scaling = prepared$scaling,
    baseline = baseline_capture, predictor = predictor_capture
  )
}

trait_results <- stats::setNames(lapply(traits, load_trait_pair), traits)
detail_rows <- list()
headline_rows <- list()
for (trait in traits) {
  item <- trait_results[[trait]]
  baseline_row <- extract_separate_species(
    item$baseline, trait_labels[[trait]], "matched baseline", item$data
  )
  predictor_row <- extract_separate_species(
    item$predictor, trait_labels[[trait]], "full predictor", item$data
  )
  detail_rows[[paste0(trait, "_baseline")]] <- baseline_row
  detail_rows[[paste0(trait, "_predictor")]] <- predictor_row
  headline_rows[[trait]] <- percent_resolved_row(
    trait_labels[[trait]], baseline_row, predictor_row
  )
}

# Fit graminoid status and woody status as two separate binary predictors.
# The same full angiosperm row set, sampling covariance, climate scaling, and
# random structure are used for each matched comparison.
fit_functional_status <- function(status_column, model_column, status_levels,
                                  formula, label, cache_prefix) {
  data <- auc_taxonomy |>
    dplyr::filter(!is.na(.data[[status_column]])) |>
    dplyr::mutate(
      "{model_column}" := factor(
        .data[[status_column]], levels = status_levels
      )
    )
  V <- subset_auc_v(V_auc, data$effect_id)
  baseline <- load_or_fit_temporary(
    paste0(cache_prefix, "_baseline_v2"), data$effect_id,
    function() fit_species_model(
      data, V, climate_formula, paste(label, "matched M_both baseline")
    )
  )
  predictor <- load_or_fit_temporary(
    paste0(cache_prefix, "_predictor_v2"), data$effect_id,
    function() fit_species_model(
      data, V, formula, paste(label, "full model")
    )
  )
  baseline_row <- extract_separate_species(
    baseline, label, "matched baseline", data
  )
  predictor_row <- extract_separate_species(
    predictor, label, "full predictor", data
  )
  list(
    data = data, V = V, baseline = baseline, predictor = predictor,
    formula = formula, model_column = model_column, levels = status_levels,
    baseline_row = baseline_row, predictor_row = predictor_row
  )
}

graminoid_result <- fit_functional_status(
  "graminoid_status", "GraminoidStatus",
  c("non_graminoid", "graminoid"), graminoid_formula,
  "Graminoid status", "graminoid_status"
)
woody_result <- fit_functional_status(
  "woody_status", "WoodyStatus", c("non_woody", "woody"),
  woody_formula, "Woody status", "woody_status"
)
detail_rows$graminoid_baseline <- graminoid_result$baseline_row
detail_rows$graminoid_predictor <- graminoid_result$predictor_row
detail_rows$woody_baseline <- woody_result$baseline_row
detail_rows$woody_predictor <- woody_result$predictor_row
headline_rows$graminoid <- percent_resolved_row(
  "Graminoid status", graminoid_result$baseline_row,
  graminoid_result$predictor_row
)
headline_rows$woody <- percent_resolved_row(
  "Woody status", woody_result$baseline_row, woody_result$predictor_row
)

# Family: no unknown pseudo-family and no collapsing of small known families.
family_data <- auc_taxonomy |>
  dplyr::filter(.env$family_known)
family_V <- subset_auc_v(V_auc, family_data$effect_id)

existing_both_path <- file.path(
  "output", "step_by_step_modeling", "species_both",
  "auc_change_species_both_model.rds"
)
family_baseline_model <- readRDS(existing_both_path)
if (family_baseline_model$k != nrow(family_data) ||
    !identical(as.character(family_baseline_model$data$effect_id),
               family_data$effect_id)) {
  family_baseline <- fit_species_model(
    family_data, family_V, climate_formula, "Family matched M_both baseline"
  )
} else {
  family_baseline <- list(fit = family_baseline_model, warnings = character(),
                          error = "")
}
family_predictor <- load_or_fit_temporary(
  "family_predictor_bobyqa", family_data$effect_id,
  function() fit_family_model(family_data, family_V)
)
if (!identical(colnames(family_predictor$fit$G),
               c("intrcpt", "W_z", "T_z", "P_z", "TP_z")) ||
    !identical(colnames(family_predictor$fit$H),
               c("intrcpt", "W_z", "T_z", "P_z", "TP_z"))) {
  stop("Family model random coefficient ordering is unexpected.", call. = FALSE)
}
if (!identical(as.numeric(family_predictor$fit$rho[1:7]), rep(0, 7)) ||
    !identical(as.numeric(family_predictor$fit$phi[1:7]), rep(0, 7))) {
  stop("Family model did not retain the requested structural zero correlations.",
       call. = FALSE)
}

family_baseline_row <- extract_separate_species(
  family_baseline, "Family", "matched baseline", family_data
)
family_species_covariance <- as.matrix(family_predictor$fit$H)
family_species_row <- detail_row(
  "Family", "full predictor", "residual species", family_predictor,
  family_species_covariance["intrcpt", "intrcpt"],
  family_species_covariance["W_z", "W_z"],
  family_species_covariance[c("T_z", "P_z", "TP_z"),
                            c("T_z", "P_z", "TP_z")],
  nrow(family_data), dplyr::n_distinct(family_data$species_id)
)
family_covariance <- as.matrix(family_predictor$fit$G)
family_level_row <- detail_row(
  "Family", "full predictor", "family", family_predictor,
  family_covariance["intrcpt", "intrcpt"],
  family_covariance["W_z", "W_z"],
  family_covariance[c("T_z", "P_z", "TP_z"),
                    c("T_z", "P_z", "TP_z")],
  nrow(family_data), dplyr::n_distinct(family_data$family_id)
)
detail_rows$family_baseline <- family_baseline_row
detail_rows$family_species <- family_species_row
detail_rows$family_level <- family_level_row
headline_rows$family <- percent_resolved_row(
  "Family", family_baseline_row, family_species_row
)

detailed_table <- dplyr::bind_rows(detail_rows)
headline_table <- dplyr::bind_rows(headline_rows) |>
  dplyr::mutate(
    `Species information` = factor(
      .data$`Species information`,
      levels = c("LDMC", "LMA", "Vegetative height",
                 "Functional group", "Family")
    )
  ) |>
  dplyr::arrange(.data$`Species information`) |>
  dplyr::mutate(
    `Species information` = as.character(.data$`Species information`),
    dplyr::across(dplyr::where(is.numeric), function(x) round(x, 1))
  )

invisible(write_csv_atomic(headline_table, output_files$headline))
invisible(write_csv_atomic(detailed_table, output_files$detailed))

# Shared supported low/high climate grids. All 42 experiments have equal
# status in the support definition; grid-cell density is only graphical.
support <- build_auc_climate_grid(auc_data, points = 65L)
if (!identical(unique(support$grid$support_method),
               "three-dimensional convex hull")) {
  stop("Response surfaces require the observed three-dimensional climate hull.",
       call. = FALSE)
}
surface_grid <- support$grid |>
  dplyr::filter(
    .data$supported,
    as.character(.data$warming_level) %in%
      c("low (25th percentile)", "high (75th percentile)")
  ) |>
  dplyr::mutate(TP_z = .data$T_z * .data$P_z)

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

warming_values <- experiment_climates |>
  dplyr::mutate(W = .data$W_z * W_scale + W_center) |>
  dplyr::summarise(
    low = stats::quantile(.data$W, 0.25),
    high = stats::quantile(.data$W, 0.75)
  )
warming_labels <- c(
  `low (25th percentile)` = paste0(
    "Low warming (", sprintf("%.2f", warming_values$low), " C)"
  ),
  `high (75th percentile)` = paste0(
    "High warming (", sprintf("%.2f", warming_values$high), " C)"
  )
)
surface_grid <- surface_grid |>
  dplyr::mutate(
    destination_T_JJA = .data$T_z * T_scale + T_center,
    destination_P_JJA = exp(.data$P_z * P_scale + P_center),
    warming_panel = factor(
      unname(.env$warming_labels[as.character(.data$warming_level)]),
      levels = unname(.env$warming_labels)
    )
  )

predict_fixed <- function(fit, formula, newdata) {
  design <- auc_design_matrix(formula, newdata, names(stats::coef(fit)))
  as.numeric(design %*% stats::coef(fit))
}

surface_theme <- function() {
  ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )
}

surface_plot <- function(grid, facet_formula, title, fill_limit) {
  ggplot2::ggplot(
    grid,
    ggplot2::aes(
      x = .data$T_z, y = .data$P_z,
      fill = .data$predicted_delta_auc
    )
  ) +
    ggplot2::geom_raster() +
    ggplot2::facet_grid(facet_formula) +
    ggplot2::scale_fill_gradient2(
      low = "#E66101", mid = "white", high = "#1B9E77",
      midpoint = 0, limits = c(-fill_limit, fill_limit),
      name = "Predicted\nDelta AUC"
    ) +
    ggplot2::scale_x_continuous(
      labels = function(x) sprintf("%.1f", x * T_scale + T_center)
    ) +
    ggplot2::scale_y_continuous(
      labels = function(x) sprintf("%.0f", exp(x * P_scale + P_center))
    ) +
    ggplot2::coord_cartesian(expand = FALSE) +
    surface_theme() +
    ggplot2::labs(
      title = title,
      x = "Destination JJA temperature (C)",
      y = "Destination JJA precipitation"
    )
}

# Three requested 2x2 trait figures: rows low/high W, columns -1/+1 SD trait.
for (trait in traits) {
  item <- trait_results[[trait]]
  trait_grid <- tidyr::crossing(
    surface_grid,
    tibble::tibble(
      Trait_z = c(-1, 1),
      trait_panel = factor(
        c("Low trait (-1 SD)", "High trait (+1 SD)"),
        levels = c("Low trait (-1 SD)", "High trait (+1 SD)")
      )
    )
  )
  trait_grid$predicted_delta_auc <- predict_fixed(
    item$predictor$fit, trait_formula, trait_grid
  )
  fill_limit <- max(abs(trait_grid$predicted_delta_auc), na.rm = TRUE)
  plot <- surface_plot(
    trait_grid, warming_panel ~ trait_panel,
    paste0(trait_labels[[trait]], ": population-level predicted Delta AUC"),
    fill_limit
  )
  figure_file <- switch(
    trait,
    LDMC = output_files$trait_LDMC,
    LMA = output_files$trait_LMA,
    vegetative_height = output_files$trait_height
  )
  ggplot2::ggsave(
    figure_file, plot, width = 9.5, height = 7.2,
    device = grDevices::cairo_pdf
  )
}

# Separate binary functional-status fixed-effect response figures.
save_functional_surface <- function(result, output_path, title) {
  model_column <- result$model_column
  status_grid <- tidyr::crossing(surface_grid, status = result$levels)
  status_grid[[model_column]] <- factor(
    status_grid$status, levels = result$levels
  )
  status_grid$status_panel <- factor(
    gsub("_", " ", status_grid$status, fixed = TRUE),
    levels = gsub("_", " ", result$levels, fixed = TRUE)
  )
  status_grid$predicted_delta_auc <- predict_fixed(
    result$predictor$fit, result$formula, status_grid
  )
  fill_limit <- max(abs(status_grid$predicted_delta_auc), na.rm = TRUE)
  plot <- surface_plot(
    status_grid, warming_panel ~ status_panel, title, fill_limit
  )
  ggplot2::ggsave(
    output_path, plot, width = 9.5, height = 7.2,
    device = grDevices::cairo_pdf
  )
}
save_functional_surface(
  graminoid_result, output_files$graminoid_figure,
  "Graminoid status: population-level predicted Delta AUC"
)
save_functional_surface(
  woody_result, output_files$woody_figure,
  "Woody status: population-level predicted Delta AUC"
)

# Family BLUPs from the family GEN block; species random effects remain zero.
family_fit <- family_predictor$fit
family_X <- as.matrix(family_fit$X)
family_residual <- as.numeric(family_fit$yi - family_X %*% family_fit$beta)
family_weighted_residual <- Matrix::solve(family_fit$M, family_residual)
family_levels <- levels(family_fit$mf.g$outer)
Z_family <- do.call(cbind, lapply(seq_len(ncol(family_fit$G)), function(j) {
  Matrix::Diagonal(x = as.numeric(family_fit$Z.G1[, j])) %*%
    family_fit$Z.G2
}))
D_family <- Matrix::kronecker(
  as.matrix(family_fit$G), Matrix::Diagonal(length(family_levels))
)
family_random_vector <- as.numeric(
  D_family %*% Matrix::crossprod(Z_family, family_weighted_residual)
)
family_random_matrix <- matrix(
  family_random_vector, nrow = length(family_levels),
  ncol = ncol(family_fit$G),
  dimnames = list(family_levels, colnames(family_fit$G))
)

family_grid <- tidyr::crossing(
  surface_grid,
  family_id = family_levels
)
family_grid$fixed_prediction <- predict_fixed(
  family_fit, climate_formula, family_grid
)
family_indices <- match(family_grid$family_id, family_levels)
family_random_design <- cbind(
  intrcpt = 1, W_z = family_grid$W_z, T_z = family_grid$T_z,
  P_z = family_grid$P_z, TP_z = family_grid$TP_z
)
family_grid$family_random_prediction <- rowSums(
  family_random_design * family_random_matrix[family_indices, , drop = FALSE]
)
family_grid$predicted_delta_auc <- family_grid$fixed_prediction +
  family_grid$family_random_prediction
family_grid <- family_grid |>
  dplyr::left_join(
    dplyr::select(
      family_counts, "family_id", "n_species", "n_species_experiment_effects",
      "included_in_cleaner_figure"
    ),
    by = "family_id"
  ) |>
  dplyr::mutate(
    family_panel = paste0(
      .data$family_id, "\n", .data$n_species, " species; ",
      .data$n_species_experiment_effects, " effects"
    )
  )
family_limit <- max(abs(family_grid$predicted_delta_auc), na.rm = TRUE)

save_family_multipage <- function(data, path, families_per_page = 6L) {
  selected_families <- family_counts |>
    dplyr::filter(.data$family_id %in% unique(data$family_id)) |>
    dplyr::arrange(dplyr::desc(.data$n_species), .data$family_id) |>
    dplyr::pull(.data$family_id)
  pages <- split(
    selected_families,
    ceiling(seq_along(selected_families) / families_per_page)
  )
  grDevices::cairo_pdf(path, width = 15, height = 7.5, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (page_families in pages) {
    page <- data |>
      dplyr::filter(.data$family_id %in% .env$page_families) |>
      dplyr::mutate(
        family_panel = factor(
          .data$family_panel,
          levels = unique(.data$family_panel[match(
            .env$page_families, .data$family_id
          )])
        )
      )
    plot <- surface_plot(
      page, warming_panel ~ family_panel,
      "Family-level predicted Delta AUC (fixed effects + family BLUP; species effects zero)",
      family_limit
    ) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
    print(plot)
  }
  invisible(path)
}

save_family_multipage(
  family_grid, output_files$family_complete, families_per_page = 6L
)
save_family_multipage(
  dplyr::filter(family_grid, .data$included_in_cleaner_figure),
  output_files$family_clean, families_per_page = 6L
)

# One requested synthesis figure. Total is shown as a diamond, not a fourth bar.
headline_long <- headline_table |>
  tidyr::pivot_longer(
    cols = c(
      "Overall response heterogeneity resolved (%)",
      "Warming-related response heterogeneity resolved (%)",
      "Background-climate response heterogeneity resolved (%)"
    ),
    names_to = "dimension", values_to = "percent_resolved"
  ) |>
  dplyr::mutate(
    dimension = factor(
      .data$dimension,
      levels = c(
        "Overall response heterogeneity resolved (%)",
        "Warming-related response heterogeneity resolved (%)",
        "Background-climate response heterogeneity resolved (%)"
      ),
      labels = c("Overall response", "Warming-related", "Background climate")
    ),
    `Species information` = factor(
      .data$`Species information`,
      levels = rev(c("LDMC", "LMA", "Vegetative height",
                     "Functional group", "Family"))
    )
  )
headline_total <- headline_table |>
  dplyr::transmute(
    `Species information` = factor(
      .data$`Species information`, levels = levels(headline_long$`Species information`)
    ),
    total = .data$`Total response heterogeneity resolved (%)`
  )
synthesis_plot <- ggplot2::ggplot(
  headline_long,
  ggplot2::aes(
    x = .data$percent_resolved, y = .data$`Species information`,
    fill = .data$dimension
  )
) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey25", linewidth = 0.45) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.72),
                    width = 0.66) +
  ggplot2::geom_point(
    data = headline_total,
    ggplot2::aes(x = .data$total, y = .data$`Species information`),
    inherit.aes = FALSE, shape = 23, size = 3.0, fill = "white",
    colour = "black", stroke = 0.8
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Overall response" = "#4C78A8",
      "Warming-related" = "#F58518",
      "Background climate" = "#54A24B"
    )
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "top", legend.title = ggplot2::element_blank()
  ) +
  ggplot2::labs(
    title = "Species information and residual among-species response heterogeneity",
    subtitle = "Bars show response-scale dimensions; diamonds show total heterogeneity resolved",
    x = "Residual species Delta-AUC heterogeneity resolved (%)",
    y = NULL
  )
ggplot2::ggsave(
  output_files$synthesis_figure, synthesis_plot,
  width = 10.5, height = 5.8, device = grDevices::cairo_pdf
)

# Concise text report, dynamically generated from the fitted results.
format_metric <- function(x) formatC(x, digits = 6, format = "fg")
format_percent <- function(x) sprintf("%.1f", x)
matched_N <- detailed_table |>
  dplyr::filter(
    .data$model_role == "matched baseline",
    .data$heterogeneity_level == "species"
  ) |>
  dplyr::select("species_information", "n_effects", "n_species", "n_experiments")

diagnostic_lines <- detailed_table |>
  dplyr::filter(.data$model_role == "full predictor") |>
  dplyr::mutate(
    line = paste0(
      .data$species_information, " [", .data$heterogeneity_level, "]: optimizer=",
      .data$optimizer, " (code ", .data$optimizer_code, "); finite=",
      .data$finite_fit,
      "; context positive definite=", .data$context_positive_definite,
      "; context min eigenvalue=", format_metric(.data$context_min_eigenvalue),
      "; max |context correlation|=",
      format_metric(.data$max_absolute_context_correlation),
      "; boundary=", .data$any_random_variance_at_boundary,
      "; identifiable=", .data$random_structure_identifiable,
      dplyr::if_else(.data$warnings == "", "", paste0("; warnings=", .data$warnings))
    )
  ) |>
  dplyr::pull(.data$line)

report_lines <- c(
  "SPECIES INFORMATION AND AMONG-SPECIES DELTA-AUC HETEROGENEITY",
  "============================================================",
  "",
  "Question: How much does knowing species information resolve among-species",
  "heterogeneity in predicted Delta AUC response?",
  "",
  "MATCHED BASELINES",
  "-----------------",
  vapply(seq_len(nrow(matched_N)), function(i) {
    paste0(
      matched_N$species_information[[i]], ": ", matched_N$n_effects[[i]],
      " effects; ", matched_N$n_species[[i]], " species; ",
      matched_N$n_experiments[[i]], " experiments"
    )
  }, character(1)),
  "",
  "MODELS",
  "------",
  "",
  "Matched M_both baseline:",
  "M0 <- metafor::rma.mv(",
  "  yi = auc_effect,",
  "  V = V_subset,",
  "  mods = ~ T_z * P_z * W_z,",
  "  random = list(",
  "    ~ 1 | species_id,",
  "    ~ 0 + W_z | species_id,",
  "    ~ 0 + T_z + P_z + TP_z | species_id",
  "  ),",
  "  struct = \"GEN\",",
  "  data = predictor_data,",
  "  method = \"REML\",",
  "  sparse = TRUE",
  ")",
  "",
  "Separate LDMC, LMA, and vegetative-height models:",
  "M_trait <- metafor::rma.mv(",
  "  yi = auc_effect,",
  "  V = V_trait,",
  "  mods = ~ Trait_z * T_z * P_z * W_z,",
  "  random = list(",
  "    ~ 1 | species_id,",
  "    ~ 0 + W_z | species_id,",
  "    ~ 0 + T_z + P_z + TP_z | species_id",
  "  ),",
  "  struct = \"GEN\",",
  "  data = trait_data,",
  "  method = \"REML\",",
  "  sparse = TRUE",
  ")",
  "",
  "Separate graminoid-status model:",
  "M_graminoid <- metafor::rma.mv(",
  "  yi = auc_effect,",
  "  V = V_graminoid,",
  "  mods = ~ GraminoidStatus * T_z * P_z * W_z,",
  "  random = list(",
  "    ~ 1 | species_id,",
  "    ~ 0 + W_z | species_id,",
  "    ~ 0 + T_z + P_z + TP_z | species_id",
  "  ),",
  "  struct = \"GEN\",",
  "  data = graminoid_data,",
  "  method = \"REML\",",
  "  sparse = TRUE",
  ")",
  "",
  "Separate woody-status model:",
  "M_woody <- metafor::rma.mv(",
  "  yi = auc_effect,",
  "  V = V_woody,",
  "  mods = ~ WoodyStatus * T_z * P_z * W_z,",
  "  random = list(",
  "    ~ 1 | species_id,",
  "    ~ 0 + W_z | species_id,",
  "    ~ 0 + T_z + P_z + TP_z | species_id",
  "  ),",
  "  struct = \"GEN\",",
  "  data = woody_data,",
  "  method = \"REML\",",
  "  sparse = TRUE",
  ")",
  "",
  "Hierarchical family model:",
  "M_family <- metafor::rma.mv(",
  "  yi = auc_effect,",
  "  V = V_family,",
  "  mods = ~ T_z * P_z * W_z,",
  "  random = list(",
  "    ~ 1 + W_z + T_z + P_z + TP_z | family_id,",
  "    ~ 1 + W_z + T_z + P_z + TP_z | species_id",
  "  ),",
  "  struct = c(\"GEN\", \"GEN\"),",
  "  rho = c(rep(0, 7), NA, NA, NA),",
  "  phi = c(rep(0, 7), NA, NA, NA),",
  "  data = family_data,",
  "  method = \"REML\",",
  "  sparse = TRUE,",
  "  control = list(optimizer = \"bobyqa\", maxfun = 100000)",
  ")",
  "",
  "RESPONSE-SCALE DEFINITIONS",
  "--------------------------",
  "H_overall = species random-intercept variance.",
  paste0(
    "H_W = mean_e(W_e^2 * sigma_W^2), evaluated over the 42 distinct ",
    "observed experiment climates with equal experiment weight."
  ),
  paste0(
    "H_context = mean_e(x_e' G_context x_e), where x_e = [T_z, P_z, ",
    "T_z*P_z], using the full T/P/T*P covariance matrix and equal ",
    "experiment weight."
  ),
  "H_total = H_overall + H_W + H_context because the three blocks are independent.",
  "Percent resolved = 100 * (H_matched_baseline - H_predictor) / H_matched_baseline.",
  "Negative values are retained and mean residual heterogeneity increased.",
  "These percentages are descriptive and are not formal R-squared values.",
  "",
  "HEADLINE RESULTS (%)",
  "--------------------",
  paste(capture.output(print(
    headline_table |>
      dplyr::mutate(dplyr::across(dplyr::where(is.numeric), format_percent)),
    n = Inf, width = Inf
  )), collapse = "\n"),
  "",
  "FUNCTIONAL-STATUS HANDLING",
  "--------------------------",
  paste0(
    "Graminoid status (graminoid versus non-graminoid) and woody status ",
    "(woody versus non-woody) are independent binary definitions and are ",
    "tested in separate models. They are not collapsed into one mutually ",
    "exclusive forb/graminoid/other factor."
  ),
  paste0(
    "The six prespecified non-angiosperms/conifer are retained in upstream ",
    "winner/loser calculations but excluded before these downstream models: ",
    paste(resident_downstream_model_exclusions(), collapse = "; "), "."
  ),
  "",
  "FAMILY HANDLING",
  "---------------",
  paste0(
    "Genuinely missing families were excluded and no Unknown family was created ",
    "(missing effects/species: ", sum(!family_known), "/",
    dplyr::n_distinct(auc_taxonomy$species_id[!family_known]), ")."
  ),
  paste0(
    "All ", nrow(family_counts), " known families, including small families, ",
    "remained in the hierarchical model and were partially pooled. The >=5 ",
    "species threshold was used only for the cleaner family figure."
  ),
  paste0(
    "The family model used separate five-coefficient family and species GEN ",
    "blocks. Correlations linking intercept or W to other coefficients were ",
    "fixed at zero; only T/P/T*P correlations were estimated, exactly preserving ",
    "the requested independent intercept, independent W, and correlated context blocks."
  ),
  paste0(
    "Optimizer history: the initial nlminb attempt ended with convergence code 1. ",
    "The unchanged model was rerun with bobyqa, which ended normally with ierr=0."
  ),
  "",
  paste0(
    "The structural zeros preserve independent intercept and warming blocks while ",
    "estimating only the T/P/T×P correlations."
  ),
  paste0(
    "The family model converged normally with `bobyqa` after `nlminb` returned ",
    "convergence code 1. However, its family-level covariance is effectively ",
    "non-identifiable: family intercept and warming variances approached zero, ",
    "the maximum context correlation was approximately 0.99998, and one ",
    "eigenvalue was approximately \\(3.1\\times10^{-15}\\). Therefore, the ",
    "family percentages require particular caution."
  ),
  "",
  "CONVERGENCE AND IDENTIFIABILITY",
  "--------------------------------",
  diagnostic_lines,
  "",
  "OUTPUTS",
  "-------",
  vapply(output_files, function(path) normalizePath(path, mustWork = FALSE),
         character(1))
)
writeLines(report_lines, output_files$report, useBytes = TRUE)

# Exact output contract: no model RDS files or unrequested prediction tables.
missing_outputs <- unlist(output_files)[!file.exists(unlist(output_files))]
empty_outputs <- unlist(output_files)[
  file.exists(unlist(output_files)) & file.info(unlist(output_files))$size <= 0
]
if (length(missing_outputs) || length(empty_outputs)) {
  stop(
    "Output contract failed. Missing: ", paste(missing_outputs, collapse = ", "),
    "; empty: ", paste(empty_outputs, collapse = ", "), call. = FALSE
  )
}

cat("\nSpecies-information heterogeneity step completed.\n")
cat("Output directory: ", normalizePath(output_dir), "\n", sep = "")
print(headline_table |>
        dplyr::mutate(dplyr::across(
          dplyr::where(is.numeric), function(x) round(x, 1)
        )))
