# Step 18: family-aggregated resident Delta-AUC analysis ---------------------
#
# Family is the biological analysis unit. Resident species cover is first
# summed within experiment x plot x year x treatment x family; family AUC
# effects are then recalculated from those cover trajectories. Species AUC
# effects are never averaged.

source("config/analysis.R")
source("R/pipeline/utils.R")
source("R/pipeline/winner_loser.R")
source("R/pipeline/frequentist_hierarchical_auc.R")

required_packages <- c(
  "DescTools", "dplyr", "emmeans", "geometry", "ggplot2", "ggrepel",
  "here", "Matrix", "metafor", "readr", "tibble", "tidyr"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "),
       call. = FALSE)
}

output_dir <- file.path(
  "output", "step_by_step_modeling", "family_aggregated_auc"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_files <- list(
  family_auc_data = file.path(output_dir, "family_auc_data.csv"),
  V_family_auc = file.path(output_dir, "V_family_auc.rds"),
  model_summary = file.path(output_dir, "family_model_summary.csv"),
  random_diagnostics = file.path(
    output_dir, "family_random_effect_diagnostics.csv"
  ),
  classification = file.path(output_dir, "family_classification.csv"),
  comparison = file.path(
    output_dir, "species_vs_family_classification_summary.csv"
  ),
  classification_figure = file.path(
    output_dir, "family_minmax_classification.pdf"
  ),
  surface_figure = file.path(
    output_dir, "selected_family_response_surfaces_median_high.pdf"
  ),
  report = file.path(output_dir, "family_aggregated_auc_report.txt")
)

community <- readRDS(file.path("output", "intermediate", "community.rds"))
experiments <- readr::read_csv(
  file.path("output", "tables", "experiments.csv"), show_col_types = FALSE
) |>
  dplyr::filter(.data$included) |>
  dplyr::arrange(.data$experiment_id)
if (nrow(experiments) != 42L) {
  stop("Expected exactly 42 included experiments.", call. = FALSE)
}

# Exact maximum-duration resident species x experiment set that supplied the
# species M_both response. Source species names are retained for plot-level
# aggregation; accepted family assignments come from the audited analysis data.
species_results_upstream <- readr::read_csv(
  file.path("output", "tables", "winner_loser_results.csv"),
  show_col_types = FALSE
) |>
  dplyr::mutate(horizon = as.character(.data$horizon)) |>
  dplyr::filter(
    .data$transform == "raw", .data$contrast == "to_origin",
    .data$pool == "resident", .data$horizon == "max"
  ) |>
  dplyr::transmute(
    experiment_id = as.character(.data$experiment_id),
    source_species = as.character(.data$species)
  ) |>
  dplyr::distinct()

taxonomy <- readr::read_csv(
  file.path("output", "tables", "analysis_dataset.csv"),
  show_col_types = FALSE
) |>
  dplyr::transmute(
    source_species = as.character(.data$species),
    species_id = as.character(.data$accepted_name),
    family_id = trimws(as.character(.data$family))
  ) |>
  dplyr::distinct()
if (anyDuplicated(taxonomy$source_species) ||
    any(is.na(taxonomy$family_id) | taxonomy$family_id == "")) {
  stop("Audited family assignment is missing or non-unique.", call. = FALSE)
}
resident_membership <- species_results_upstream |>
  dplyr::left_join(taxonomy, by = "source_species") |>
  filter_downstream_model_species("species_id")
if (nrow(resident_membership) != 1604L ||
    dplyr::n_distinct(resident_membership$experiment_id) != 42L) {
  stop(
    "The filtered resident downstream species-effect set was not recovered ",
    "(expected 1,604 effects across 42 experiments).", call. = FALSE
  )
}
if (any(is.na(resident_membership$family_id) |
        resident_membership$family_id == "")) {
  stop("At least one resident source species lacks an audited family.",
       call. = FALSE)
}

# Exact family analogue of complete_horizon_community(). Family membership is
# fixed once per experiment from the exact M_both resident species set.
complete_family_community <- function(data, membership) {
  family_universe <- sort(unique(membership$family_id))
  surveyed_plot_years <- data |>
    dplyr::filter(!is.na(.data$relative_year)) |>
    dplyr::distinct(.data$ODT, .data$stable_plot_id, .data$relative_year)
  observed <- data |>
    dplyr::filter(
      .data$species %in% membership$source_species,
      !is.na(.data$relative_year)
    ) |>
    dplyr::left_join(membership, by = c("species" = "source_species")) |>
    dplyr::group_by(
      .data$ODT, .data$stable_plot_id, .data$relative_year, .data$family_id
    ) |>
    dplyr::summarise(
      relative_cover = sum(.data$relative_cover), .groups = "drop"
    )
  surveyed_plot_years |>
    tidyr::crossing(family_id = family_universe) |>
    dplyr::left_join(
      observed,
      by = c("ODT", "stable_plot_id", "relative_year", "family_id")
    ) |>
    dplyr::mutate(relative_cover = dplyr::coalesce(.data$relative_cover, 0)) |>
    dplyr::group_by(.data$ODT) |>
    tidyr::complete(
      stable_plot_id, relative_year, family_id = family_universe
    ) |>
    dplyr::ungroup() |>
    dplyr::group_by(.data$ODT, .data$relative_year, .data$family_id) |>
    dplyr::mutate(
      relative_cover = dplyr::if_else(
        is.na(.data$relative_cover),
        mean(.data$relative_cover, na.rm = TRUE),
        .data$relative_cover
      ),
      relative_cover = dplyr::coalesce(.data$relative_cover, 0)
    ) |>
    dplyr::ungroup()
}

family_auc_by_plot <- function(completed, experiment_id) {
  completed |>
    dplyr::group_by(
      .data$ODT, .data$stable_plot_id, .data$family_id
    ) |>
    dplyr::summarise(
      auc = if (dplyr::n_distinct(.data$relative_year) < 2L) {
        NA_real_
      } else {
        DescTools::AUC(
          .data$relative_year, .data$relative_cover, method = "trapezoid"
        )
      },
      .groups = "drop"
    ) |>
    dplyr::mutate(
      experiment_id = experiment_id,
      horizon = "max",
      transform = "raw"
    )
}

family_cover_audit <- list()
family_auc_plot_blocks <- vector("list", nrow(experiments))
for (i in seq_len(nrow(experiments))) {
  experiment <- experiments[i, ]
  experiment_id <- experiment$experiment_id[[1]]
  data <- experiment_community(community, experiment) |>
    apply_legacy_observation_filter()
  membership <- resident_membership |>
    dplyr::filter(.data$experiment_id == .env$experiment_id)
  if (!nrow(membership)) {
    stop("No resident membership for ", experiment_id, ".", call. = FALSE)
  }
  completed <- complete_family_community(data, membership)
  family_auc_plot_blocks[[i]] <- family_auc_by_plot(completed, experiment_id)
  family_cover_audit[[i]] <- tibble::tibble(
    experiment_id = experiment_id,
    n_resident_species = dplyr::n_distinct(membership$source_species),
    n_resident_families = dplyr::n_distinct(membership$family_id),
    n_completed_family_plot_year_rows = nrow(completed),
    n_origin_plots = dplyr::n_distinct(
      completed$stable_plot_id[completed$ODT == "origin_controls"]
    ),
    n_warmed_plots = dplyr::n_distinct(
      completed$stable_plot_id[completed$ODT == "warmed"]
    )
  )
}
family_auc_plot <- dplyr::bind_rows(family_auc_plot_blocks)
family_cover_audit <- dplyr::bind_rows(family_cover_audit)

# Per-experiment family AUC contrasts and full emmeans covariance blocks.
contrast_blocks <- list()
covariance_blocks <- list()
contrast_warnings <- character()
for (experiment_id in experiments$experiment_id) {
  raw <- family_auc_plot |>
    dplyr::filter(
      .data$experiment_id == .env$experiment_id, !is.na(.data$auc)
    ) |>
    dplyr::mutate(
      ODT = factor(
        as.character(.data$ODT),
        levels = c("destination_controls", "origin_controls", "warmed")
      ),
      family_id = factor(.data$family_id)
    )
  warnings <- character()
  fitted <- tryCatch(
    withCallingHandlers({
      model <- stats::lm(auc ~ ODT * family_id, data = raw)
      means <- emmeans::emmeans(model, ~ ODT | family_id)
      contrast <- emmeans::contrast(
        means, method = list(to_origin = c(0, -1, 1)), adjust = "none"
      )
      list(
        summary = tibble::as_tibble(as.data.frame(summary(contrast))),
        covariance = as.matrix(stats::vcov(contrast))
      )
    }, warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }),
    error = function(condition) condition
  )
  if (inherits(fitted, "condition")) {
    stop("Family AUC contrast failed for ", experiment_id, ": ",
         conditionMessage(fitted), call. = FALSE)
  }
  target_families <- resident_membership |>
    dplyr::filter(.data$experiment_id == .env$experiment_id) |>
    dplyr::distinct(.data$family_id) |>
    dplyr::arrange(.data$family_id) |>
    dplyr::pull(.data$family_id)
  index <- match(target_families, as.character(fitted$summary$family_id))
  if (anyNA(index)) {
    stop("A resident family is missing from contrasts for ", experiment_id,
         ".", call. = FALSE)
  }
  contrast_blocks[[experiment_id]] <- fitted$summary[index, , drop = FALSE] |>
    dplyr::transmute(
      experiment_id = experiment_id,
      family_id = as.character(.data$family_id),
      family_auc_effect = as.numeric(.data$estimate),
      family_auc_se = as.numeric(.data$SE),
      contrast_df = as.numeric(.data$df)
    )
  covariance_blocks[[experiment_id]] <-
    fitted$covariance[index, index, drop = FALSE]
  if (length(warnings)) {
    contrast_warnings <- c(
      contrast_warnings,
      paste(experiment_id, unique(warnings), sep = ": ")
    )
  }
}

family_auc_data <- dplyr::bind_rows(contrast_blocks) |>
  dplyr::arrange(.data$experiment_id, .data$family_id) |>
  dplyr::mutate(
    horizon = "max",
    effect_id = paste("max", .data$experiment_id, .data$family_id, sep = "::")
  )
ordered_block_ids <- unlist(lapply(names(contrast_blocks), function(id) {
  paste("max", id, contrast_blocks[[id]]$family_id, sep = "::")
}), use.names = FALSE)
V_assembled <- as.matrix(Matrix::bdiag(covariance_blocks))
V_index <- match(family_auc_data$effect_id, ordered_block_ids)
if (anyNA(V_index)) stop("Family V ordering could not be aligned.", call. = FALSE)
V_family_auc <- V_assembled[V_index, V_index, drop = FALSE]
dimnames(V_family_auc) <- list(
  family_auc_data$effect_id, family_auc_data$effect_id
)

# Reuse species-analysis climate scaling exactly; never rescale the family set.
species_scaled <- scale_resident_auc_data(
  read_resident_auc_data(here::here(), horizons = "max")
)
climate <- species_scaled$data |>
  dplyr::transmute(
    experiment_id = as.character(.data$experiment_id),
    T = .data$T, P = .data$P, P_raw = .data$P_raw, W = .data$W,
    T_z = .data$T_z, P_z = .data$P_z, W_z = .data$W_z
  ) |>
  dplyr::distinct()
if (nrow(climate) != 42L) {
  stop("Species scaling did not yield 42 experiment climates.", call. = FALSE)
}
family_membership_counts <- resident_membership |>
  dplyr::summarise(
    n_resident_species_in_experiment = dplyr::n_distinct(.data$source_species),
    .by = c("experiment_id", "family_id")
  )
family_auc_data <- family_auc_data |>
  dplyr::left_join(climate, by = "experiment_id") |>
  dplyr::left_join(
    family_membership_counts, by = c("experiment_id", "family_id")
  ) |>
  dplyr::mutate(TP_z = .data$T_z * .data$P_z)
if (any(!stats::complete.cases(family_auc_data)) ||
    any(!is.finite(family_auc_data$family_auc_se) |
          family_auc_data$family_auc_se <= 0) ||
    anyDuplicated(family_auc_data$effect_id)) {
  stop("Family AUC analysis data failed completeness/key/SE validation.",
       call. = FALSE)
}

symmetry_error <- max(abs(V_family_auc - t(V_family_auc)))
diagonal_error <- max(abs(diag(V_family_auc) -
                            family_auc_data$family_auc_se^2))
V_eigenvalues <- eigen(
  V_family_auc, symmetric = TRUE, only.values = TRUE
)$values
V_tolerance <- 1e-8
V_valid <- identical(rownames(V_family_auc), family_auc_data$effect_id) &&
  symmetry_error <= V_tolerance && diagonal_error <= V_tolerance &&
  min(V_eigenvalues) >= -V_tolerance
if (!V_valid) {
  stop("Full family AUC sampling covariance failed validation; it was not repaired.",
       call. = FALSE)
}

invisible(write_csv_atomic(family_auc_data, output_files$family_auc_data))
invisible(write_rds_atomic(
  list(
    V_family_auc = Matrix::Matrix(V_family_auc, sparse = TRUE),
    effect_id = family_auc_data$effect_id,
    validation = tibble::tibble(
      n_effects = nrow(family_auc_data),
      n_experiments = dplyr::n_distinct(family_auc_data$experiment_id),
      dimensions_valid = all(dim(V_family_auc) == nrow(family_auc_data)),
      order_valid = identical(
        rownames(V_family_auc), family_auc_data$effect_id
      ),
      maximum_symmetry_error = symmetry_error,
      maximum_diagonal_variance_error = diagonal_error,
      minimum_eigenvalue = min(V_eigenvalues),
      valid = V_valid
    )
  ),
  output_files$V_family_auc
))

# Exact family analogue of M_both. A second optimizer is allowed only for the
# unchanged model if nlminb genuinely fails; no random component is simplified.
model_data <- family_auc_data |>
  dplyr::mutate(
    family_id = factor(.data$family_id),
    experiment_id = factor(.data$experiment_id)
  )
model_formula <- ~ T_z * P_z * W_z
model_random <- list(
  ~ 1 | family_id,
  ~ 0 + W_z | family_id,
  ~ 0 + T_z + P_z + TP_z | family_id
)

fit_family_model <- function(control = NULL) {
  arguments <- list(
    yi = model_data$family_auc_effect,
    V = V_family_auc,
    mods = model_formula,
    random = model_random,
    struct = "GEN",
    data = model_data,
    method = "REML",
    sparse = TRUE
  )
  if (!is.null(control)) arguments$control <- control
  auc_capture(do.call(metafor::rma.mv, arguments))
}

temporary_cache <- Sys.getenv("TRANSPLANT_FAMILY_AGG_TEMP_CACHE", "")
cache_signature <- list(
  effect_id = model_data$effect_id,
  diagonal = diag(V_family_auc),
  formula = paste(deparse(model_formula), collapse = " ")
)
model_attempts <- character()
fit_capture <- NULL
if (nzchar(temporary_cache) && file.exists(temporary_cache)) {
  cached <- readRDS(temporary_cache)
  if (is.list(cached) && identical(cached$signature, cache_signature)) {
    fit_capture <- cached$fit_capture
    model_attempts <- cached$model_attempts
  }
}
if (is.null(fit_capture)) {
  model_attempts <- "nlminb"
  fit_capture <- fit_family_model()
  if (is.null(fit_capture$fit)) {
    model_attempts <- c(model_attempts, "bobyqa")
    fit_capture <- fit_family_model(
      control = list(optimizer = "bobyqa", maxfun = 100000)
    )
  }
  if (nzchar(temporary_cache)) {
    saveRDS(
      list(
        signature = cache_signature, fit_capture = fit_capture,
        model_attempts = model_attempts
      ),
      temporary_cache
    )
  }
}
if (is.null(fit_capture$fit)) {
  stop("M_family_aggregated failed without simplification: ",
       fit_capture$error, call. = FALSE)
}
M_family_aggregated <- fit_capture$fit
if (!identical(colnames(M_family_aggregated$G), "W_z") ||
    !identical(colnames(M_family_aggregated$H),
               c("T_z", "P_z", "TP_z"))) {
  stop("M_family_aggregated random structure is not the requested structure.",
       call. = FALSE)
}

optimizer <- if (!is.null(M_family_aggregated$opt.res$ierr)) {
  "bobyqa"
} else if (!is.null(M_family_aggregated$opt.res$convergence)) {
  "nlminb"
} else {
  class(M_family_aggregated$opt.res)[[1]]
}
optimizer_code <- auc_value_or(
  M_family_aggregated$opt.res$convergence,
  auc_value_or(M_family_aggregated$opt.res$ierr, NA_integer_)
)
family_intercept_variance <- as.numeric(M_family_aggregated$sigma2[[1]])
family_W_variance <- as.numeric(M_family_aggregated$tau2[[1]])
context_covariance <- as.matrix(M_family_aggregated$H)
context_correlation <- stats::cov2cor(context_covariance)
context_eigenvalues <- eigen(
  context_covariance, symmetric = TRUE, only.values = TRUE
)$values
full_random_covariance <- matrix(
  0, nrow = 5L, ncol = 5L,
  dimnames = list(
    c("intercept", "W_z", "T_z", "P_z", "TP_z"),
    c("intercept", "W_z", "T_z", "P_z", "TP_z")
  )
)
full_random_covariance["intercept", "intercept"] <-
  family_intercept_variance
full_random_covariance["W_z", "W_z"] <- family_W_variance
full_random_covariance[c("T_z", "P_z", "TP_z"),
                       c("T_z", "P_z", "TP_z")] <- context_covariance
full_random_eigenvalues <- eigen(
  full_random_covariance, symmetric = TRUE, only.values = TRUE
)$values
boundary_tolerance <- 1e-8 * max(1, stats::var(model_data$family_auc_effect))
random_variances <- diag(full_random_covariance)
max_absolute_context_correlation <- max(
  abs(context_correlation[upper.tri(context_correlation)])
)
random_identifiable <- all(is.finite(full_random_eigenvalues)) &&
  min(full_random_eigenvalues) > boundary_tolerance &&
  all(random_variances > boundary_tolerance) &&
  max_absolute_context_correlation < 1 - 1e-6

# Family replication and rank over its observed T-P support.
family_support <- model_data |>
  dplyr::distinct(.data$family_id, .data$experiment_id,
                  .data$T_z, .data$P_z, .data$W_z) |>
  dplyr::group_by(.data$family_id) |>
  dplyr::group_modify(function(.x, .y) {
    design <- cbind(
      intercept = 1, T_z = .x$T_z, P_z = .x$P_z,
      TP_z = .x$T_z * .x$P_z
    )
    tp_design <- cbind(intercept = 1, T_z = .x$T_z, P_z = .x$P_z)
    tibble::tibble(
      n_experiments = dplyr::n_distinct(.x$experiment_id),
      tp_design_rank = qr(design)$rank,
      tp_geometric_rank = qr(tp_design)$rank,
      support_class = dplyr::case_when(
        dplyr::n_distinct(.x$experiment_id) == 1L ~ "single-point",
        qr(design)$rank == 4L ~ "full-rank",
        TRUE ~ "partial-rank"
      ),
      T_min = min(.x$T_z), T_max = max(.x$T_z),
      P_min = min(.x$P_z), P_max = max(.x$P_z),
      W_min = min(.x$W_z), W_max = max(.x$W_z)
    )
  }) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    eligible_for_individual_interpretation =
      .data$n_experiments >= 4L & .data$tp_design_rank == 4L
  )

# Conditional random family coefficients (BLUPs).
X_model <- as.matrix(M_family_aggregated$X)
fixed_beta <- stats::coef(M_family_aggregated)
model_residual <- as.numeric(
  M_family_aggregated$yi - X_model %*% M_family_aggregated$beta
)
weighted_residual <- Matrix::solve(M_family_aggregated$M, model_residual)
family_levels <- levels(M_family_aggregated$mf.g$outer)
context_family_levels <- levels(M_family_aggregated$mf.h$outer)
if (!identical(family_levels, context_family_levels)) {
  stop("Family order differs across random blocks.", call. = FALSE)
}
n_families <- length(family_levels)
Z_intercept <- M_family_aggregated$Z.S[[1]]
Z_warming <- Matrix::Diagonal(
  x = as.numeric(M_family_aggregated$Z.G1[, 1])
) %*% M_family_aggregated$Z.G2
Z_context <- do.call(cbind, lapply(seq_len(ncol(context_covariance)), function(j) {
  Matrix::Diagonal(
    x = as.numeric(M_family_aggregated$Z.H1[, j])
  ) %*% M_family_aggregated$Z.H2
}))
random_intercept <- as.numeric(
  family_intercept_variance * Matrix::crossprod(Z_intercept, weighted_residual)
)
random_W <- as.numeric(
  family_W_variance * Matrix::crossprod(Z_warming, weighted_residual)
)
D_context <- Matrix::kronecker(
  context_covariance, Matrix::Diagonal(n_families)
)
random_context_vector <- as.numeric(
  D_context %*% Matrix::crossprod(Z_context, weighted_residual)
)
random_context_matrix <- matrix(
  random_context_vector, nrow = n_families, ncol = 3L,
  dimnames = list(family_levels, c("T_z", "P_z", "TP_z"))
)
family_random <- tibble::tibble(
  family_id = family_levels,
  random_intercept = random_intercept,
  random_W = random_W,
  random_T = random_context_matrix[, "T_z"],
  random_P = random_context_matrix[, "P_z"],
  random_TP = random_context_matrix[, "TP_z"]
) |>
  dplyr::left_join(family_support, by = "family_id")
invisible(write_csv_atomic(family_random, output_files$random_diagnostics))

# Median warming is inherited from the 42 species-analysis climates.
experiment_warming <- climate |>
  dplyr::distinct(.data$experiment_id, .data$W_z, .data$W)
W_median <- stats::median(experiment_warming$W)
W_z_median <- stats::median(experiment_warming$W_z)
W_high <- as.numeric(stats::quantile(experiment_warming$W, 0.75))
W_z_high <- as.numeric(stats::quantile(experiment_warming$W_z, 0.75))
if (round(W_median, 2) != 2.29) {
  stop("Species-analysis median raw warming is not 2.29 C.", call. = FALSE)
}

family_hull_grid <- function(family_name, W_z_value, points = 81L) {
  observed <- model_data |>
    dplyr::filter(as.character(.data$family_id) == .env$family_name) |>
    dplyr::distinct(.data$T_z, .data$P_z)
  if (qr(cbind(1, observed$T_z, observed$P_z))$rank < 3L) {
    stop(family_name, " does not have 2D T-P support.", call. = FALSE)
  }
  candidate <- tidyr::expand_grid(
    T_z = seq(min(observed$T_z), max(observed$T_z), length.out = points),
    P_z = seq(min(observed$P_z), max(observed$P_z), length.out = points)
  ) |>
    dplyr::mutate(
      grid_point_type = "regular_grid",
      tile_width = diff(range(observed$T_z)) / (points - 1L),
      tile_height = diff(range(observed$P_z)) / (points - 1L)
    )
  hull <- geometry::convhulln(as.matrix(observed), options = "Qt")
  candidate |>
    dplyr::filter(geometry::inhulln(
      hull, as.matrix(dplyr::select(candidate, "T_z", "P_z"))
    )) |>
    dplyr::bind_rows(
      observed |>
        dplyr::mutate(
          grid_point_type = "observed_climate",
          tile_width = diff(range(observed$T_z)) / (points - 1L),
          tile_height = diff(range(observed$P_z)) / (points - 1L)
        )
    ) |>
    dplyr::distinct(.data$T_z, .data$P_z, .keep_all = TRUE) |>
    dplyr::mutate(
      W_z = W_z_value,
      TP_z = .data$T_z * .data$P_z
    )
}

# Conditional prediction uncertainty, identical to the species M_both method.
M_inverse_X <- Matrix::solve(M_family_aggregated$M, X_model)
fixed_covariance <- solve(crossprod(X_model, M_inverse_X))
D_family <- full_random_covariance

predict_family_surface <- function(family_name, W_z_value, points = 81L) {
  family_index <- match(family_name, family_levels)
  grid <- family_hull_grid(family_name, W_z_value, points)
  fixed_design <- auc_design_matrix(
    model_formula, grid, names(fixed_beta)
  )
  random_design <- cbind(
    intercept = 1, W_z = grid$W_z, T_z = grid$T_z,
    P_z = grid$P_z, TP_z = grid$TP_z
  )
  random <- family_random |>
    dplyr::filter(.data$family_id == .env$family_name)
  random_vector <- unname(as.numeric(random[
    1, c("random_intercept", "random_W", "random_T", "random_P", "random_TP")
  ]))
  Z_family <- cbind(
    Z_intercept[, family_index, drop = FALSE],
    Z_warming[, family_index, drop = FALSE],
    Z_context[, family_index, drop = FALSE],
    Z_context[, n_families + family_index, drop = FALSE],
    Z_context[, 2L * n_families + family_index, drop = FALSE]
  )
  M_inverse_Z <- Matrix::solve(M_family_aggregated$M, Z_family)
  XtM_inverse_Z <- as.matrix(crossprod(X_model, M_inverse_Z))
  ZtM_inverse_Z <- as.matrix(crossprod(Z_family, M_inverse_Z))
  random_prediction_error_covariance <-
    D_family - D_family %*% ZtM_inverse_Z %*% D_family +
    D_family %*% t(XtM_inverse_Z) %*% fixed_covariance %*%
      XtM_inverse_Z %*% D_family
  fixed_random_error_covariance <-
    -fixed_covariance %*% XtM_inverse_Z %*% D_family
  fixed_variance <- rowSums(
    (fixed_design %*% fixed_covariance) * fixed_design
  )
  random_variance <- rowSums(
    (random_design %*% random_prediction_error_covariance) * random_design
  )
  cross_covariance <- 2 * rowSums(
    (fixed_design %*% fixed_random_error_covariance) * random_design
  )
  prediction_variance <- pmax(
    0, fixed_variance + random_variance + cross_covariance
  )
  estimate <- as.numeric(fixed_design %*% fixed_beta) +
    as.numeric(random_design %*% random_vector)
  grid |>
    dplyr::mutate(
      family_id = family_name,
      estimate = estimate,
      prediction_se = sqrt(prediction_variance),
      conf_low = .data$estimate - stats::qnorm(0.975) * .data$prediction_se,
      conf_high = .data$estimate + stats::qnorm(0.975) * .data$prediction_se
    )
}

eligible_families <- family_support |>
  dplyr::filter(.data$eligible_for_individual_interpretation) |>
  dplyr::pull(.data$family_id) |>
  as.character()
median_predictions <- stats::setNames(lapply(
  eligible_families,
  function(family_name) predict_family_surface(
    family_name, W_z_median, points = 81L
  )
), eligible_families)

classification <- dplyr::bind_rows(lapply(eligible_families, function(family_name) {
  prediction <- median_predictions[[family_name]]
  minimum_index <- which.min(prediction$estimate)
  maximum_index <- which.max(prediction$estimate)
  tibble::tibble(
    family_id = family_name,
    min_fitted_family_delta_auc = prediction$estimate[[minimum_index]],
    max_fitted_family_delta_auc = prediction$estimate[[maximum_index]],
    min_prediction_se = prediction$prediction_se[[minimum_index]],
    min_conf_low = prediction$conf_low[[minimum_index]],
    min_conf_high = prediction$conf_high[[minimum_index]],
    min_T_z = prediction$T_z[[minimum_index]],
    min_P_z = prediction$P_z[[minimum_index]],
    max_prediction_se = prediction$prediction_se[[maximum_index]],
    max_conf_low = prediction$conf_low[[maximum_index]],
    max_conf_high = prediction$conf_high[[maximum_index]],
    max_T_z = prediction$T_z[[maximum_index]],
    max_P_z = prediction$P_z[[maximum_index]]
  )
})) |>
  dplyr::left_join(
    dplyr::select(
      family_random, "family_id", "random_W", "n_experiments",
      "tp_design_rank", "support_class"
    ),
    by = "family_id"
  ) |>
  dplyr::mutate(
    fitted_class = dplyr::case_when(
      .data$min_fitted_family_delta_auc > 0 ~ "Consistent winner",
      .data$max_fitted_family_delta_auc < 0 ~ "Consistent loser",
      TRUE ~ "Context-dependent"
    ),
    fitted_class = factor(
      .data$fitted_class,
      levels = c("Consistent loser", "Context-dependent", "Consistent winner")
    ),
    warming_C = W_median,
    warming_z = W_z_median
  )
if (any(classification$min_fitted_family_delta_auc >
          classification$max_fitted_family_delta_auc)) {
  stop("A family fitted minimum exceeds its maximum.", call. = FALSE)
}

# Two strongest representatives per fitted class where available.
selected_families <- dplyr::bind_rows(
  classification |>
    dplyr::filter(.data$fitted_class == "Consistent winner") |>
    dplyr::slice_max(
      .data$min_fitted_family_delta_auc, n = 2L, with_ties = FALSE
    ),
  classification |>
    dplyr::filter(.data$fitted_class == "Consistent loser") |>
    dplyr::slice_min(
      .data$max_fitted_family_delta_auc, n = 2L, with_ties = FALSE
    ),
  classification |>
    dplyr::filter(.data$fitted_class == "Context-dependent") |>
    dplyr::mutate(
      response_span = .data$max_fitted_family_delta_auc -
        .data$min_fitted_family_delta_auc
    ) |>
    dplyr::slice_max(.data$response_span, n = 2L, with_ties = FALSE)
) |>
  dplyr::distinct(.data$family_id) |>
  dplyr::pull(.data$family_id) |>
  as.character()
classification <- classification |>
  dplyr::mutate(
    selected_for_surface = .data$family_id %in% .env$selected_families,
    fitted_class = as.character(.data$fitted_class)
  )
invisible(write_csv_atomic(classification, output_files$classification))

# Species-versus-family descriptive classification summary.
species_classification <- readr::read_csv(
  file.path(
    "output", "step_by_step_modeling", "species_both",
    "both_full_rank_minmax_classification.csv"
  ),
  show_col_types = FALSE
)
classification_summary_row <- function(data, resolution) {
  tibble::tibble(
    Resolution = resolution,
    `Eligible units` = nrow(data),
    `Consistent winners` = sum(data$fitted_class == "Consistent winner"),
    `Consistent losers` = sum(data$fitted_class == "Consistent loser"),
    `Context-dependent` = sum(data$fitted_class == "Context-dependent"),
    `Percent context-dependent` =
      100 * sum(data$fitted_class == "Context-dependent") / nrow(data)
  )
}
comparison <- dplyr::bind_rows(
  classification_summary_row(species_classification, "Species"),
  classification_summary_row(classification, "Family")
) |>
  dplyr::mutate(
    `Percent context-dependent` = round(.data$`Percent context-dependent`, 1)
  )
invisible(write_csv_atomic(comparison, output_files$comparison))

# Model summary and global random-effect diagnostics.
correlation_value <- function(first, second) {
  context_correlation[first, second]
}
model_summary <- tibble::tibble(
  model = "M_family_aggregated",
  formula = "family_auc_effect ~ T_z * P_z * W_z",
  random_structure = paste(
    "family intercept; independent family W_z slope;",
    "correlated family T_z/P_z/TP_z block"
  ),
  n_effects = nrow(model_data),
  n_families = dplyr::n_distinct(model_data$family_id),
  n_experiments = dplyr::n_distinct(model_data$experiment_id),
  optimizer_attempts = paste(model_attempts, collapse = " -> "),
  optimizer_used = optimizer,
  optimizer_code = optimizer_code,
  warnings = paste(fit_capture$warnings, collapse = " | "),
  finite_coefficients = all(is.finite(stats::coef(M_family_aggregated))),
  finite_fixed_covariance = all(is.finite(stats::vcov(M_family_aggregated))),
  family_intercept_variance = family_intercept_variance,
  family_W_slope_variance = family_W_variance,
  family_T_variance = context_covariance["T_z", "T_z"],
  family_P_variance = context_covariance["P_z", "P_z"],
  family_TP_variance = context_covariance["TP_z", "TP_z"],
  T_P_correlation = correlation_value("T_z", "P_z"),
  T_TP_correlation = correlation_value("T_z", "TP_z"),
  P_TP_correlation = correlation_value("P_z", "TP_z"),
  minimum_context_eigenvalue = min(context_eigenvalues),
  minimum_full_random_eigenvalue = min(full_random_eigenvalues),
  maximum_absolute_context_correlation = max_absolute_context_correlation,
  boundary_tolerance = boundary_tolerance,
  any_variance_on_boundary = any(random_variances <= boundary_tolerance),
  covariance_positive_definite = min(full_random_eigenvalues) > 0,
  random_structure_identifiable = random_identifiable,
  V_source = "full family emmeans covariance reconstructed from family AUC models",
  V_maximum_diagonal_variance_error = diagonal_error,
  V_maximum_symmetry_error = symmetry_error,
  V_minimum_eigenvalue = min(V_eigenvalues),
  V_valid = V_valid,
  contrast_warnings = paste(contrast_warnings, collapse = " | ")
)
invisible(write_csv_atomic(model_summary, output_files$model_summary))

# Family min-versus-max quadrant plot, matching the species visual logic.
class_colours <- c(
  "Consistent loser" = "#E66101",
  "Context-dependent" = "#0072B2",
  "Consistent winner" = "#1B9E77"
)
x_range <- range(c(
  classification$min_conf_low, classification$min_conf_high,
  classification$min_fitted_family_delta_auc, 0
), na.rm = TRUE)
y_range <- range(c(
  classification$max_conf_low, classification$max_conf_high,
  classification$max_fitted_family_delta_auc, 0
), na.rm = TRUE)
x_limits <- x_range + c(-0.08, 0.08) * diff(x_range)
y_limits <- y_range + c(-0.08, 0.08) * diff(y_range)
class_counts <- table(factor(
  classification$fitted_class, levels = names(class_colours)
))
labels <- classification |>
  dplyr::filter(.data$selected_for_surface)
warming_fill_limit <- max(abs(classification$random_W), na.rm = TRUE)

classification_plot <- ggplot2::ggplot(classification) +
  ggplot2::annotate(
    "rect", xmin = x_limits[[1]], xmax = 0,
    ymin = y_limits[[1]], ymax = 0,
    fill = class_colours[["Consistent loser"]], alpha = 0.07
  ) +
  ggplot2::annotate(
    "rect", xmin = x_limits[[1]], xmax = 0,
    ymin = 0, ymax = y_limits[[2]],
    fill = class_colours[["Context-dependent"]], alpha = 0.07
  ) +
  ggplot2::annotate(
    "rect", xmin = 0, xmax = x_limits[[2]],
    ymin = 0, ymax = y_limits[[2]],
    fill = class_colours[["Consistent winner"]], alpha = 0.07
  ) +
  ggplot2::annotate(
    "text", x = x_limits[[1]] + 0.04 * diff(x_limits),
    y = y_limits[[1]] + 0.07 * diff(y_limits), hjust = 0,
    label = paste0("CONSISTENT LOSER\nn = ",
                   class_counts[["Consistent loser"]]),
    colour = class_colours[["Consistent loser"]], fontface = "bold", size = 3.4
  ) +
  ggplot2::annotate(
    "text", x = x_limits[[1]] + 0.04 * diff(x_limits),
    y = y_limits[[2]] - 0.06 * diff(y_limits), hjust = 0,
    label = paste0("CONTEXT-DEPENDENT\nn = ",
                   class_counts[["Context-dependent"]]),
    colour = class_colours[["Context-dependent"]], fontface = "bold", size = 3.4
  ) +
  ggplot2::annotate(
    "text", x = x_limits[[2]] - 0.04 * diff(x_limits),
    y = y_limits[[2]] - 0.06 * diff(y_limits), hjust = 1,
    label = paste0("CONSISTENT WINNER\nn = ",
                   class_counts[["Consistent winner"]]),
    colour = class_colours[["Consistent winner"]], fontface = "bold", size = 3.4
  ) +
  ggplot2::geom_abline(
    intercept = 0, slope = 1, linetype = "dashed",
    colour = "grey55", linewidth = 0.45
  ) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey20", linewidth = 0.55) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey20", linewidth = 0.55) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = .data$min_conf_low, xend = .data$min_conf_high,
      y = .data$max_fitted_family_delta_auc,
      yend = .data$max_fitted_family_delta_auc,
      colour = .data$fitted_class
    ), linewidth = 0.5, alpha = 0.25
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = .data$min_fitted_family_delta_auc,
      xend = .data$min_fitted_family_delta_auc,
      y = .data$max_conf_low, yend = .data$max_conf_high,
      colour = .data$fitted_class
    ), linewidth = 0.5, alpha = 0.25
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      x = .data$min_fitted_family_delta_auc,
      y = .data$max_fitted_family_delta_auc,
      fill = .data$random_W,
      colour = .data$fitted_class
    ), shape = 21, size = 3.5, stroke = 1.0
  ) +
  ggrepel::geom_text_repel(
    data = labels,
    ggplot2::aes(
      x = .data$min_fitted_family_delta_auc,
      y = .data$max_fitted_family_delta_auc,
      label = .data$family_id,
      colour = .data$fitted_class
    ),
    show.legend = FALSE, size = 3.2, min.segment.length = 0,
    segment.colour = "grey55", segment.size = 0.3,
    box.padding = 0.35, point.padding = 0.25, seed = 20260908
  ) +
  ggplot2::scale_colour_manual(values = class_colours, guide = "none") +
  ggplot2::scale_fill_gradient2(
    low = "#E66101", mid = "white", high = "#1B9E77", midpoint = 0,
    limits = c(-warming_fill_limit, warming_fill_limit),
    name = "Conditional family\nrandom W slope"
  ) +
  ggplot2::coord_cartesian(
    xlim = x_limits, ylim = y_limits, expand = FALSE, clip = "off"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "top",
    plot.margin = ggplot2::margin(8, 14, 8, 8)
  ) +
  ggplot2::labs(
    title = "Family resolution: minimum versus maximum fitted Delta AUC",
    subtitle = paste0(
      "Families with >=4 experiments and full-rank T-P support; median warming = ",
      sprintf("%.2f", W_median), " C"
    ),
    x = "Minimum fitted family Delta AUC across observed T-P support",
    y = "Maximum fitted family Delta AUC across observed T-P support",
    caption = paste(
      "Whiskers are pointwise model-based 95% intervals at fitted extrema;",
      "not simultaneous intervals for the whole surface."
    )
  )
ggplot2::ggsave(
  output_files$classification_figure, classification_plot,
  width = 11, height = 8.5, device = grDevices::cairo_pdf
)

# Selected family T-P surfaces: two representatives per class where available,
# common axes and common response scale, rows median/high warming.
warming_lookup <- tibble::tibble(
  warming_panel = factor(
    paste0(
      c("Median warming", "High warming"), "\nDelta T = ",
      sprintf("%.2f", c(W_median, W_high)), " C"
    ),
    levels = paste0(
      c("Median warming", "High warming"), "\nDelta T = ",
      sprintf("%.2f", c(W_median, W_high)), " C"
    )
  ),
  W_z = c(W_z_median, W_z_high)
)
surface_predictions <- dplyr::bind_rows(lapply(selected_families, function(family_name) {
  dplyr::bind_rows(lapply(seq_len(nrow(warming_lookup)), function(i) {
    predict_family_surface(
      family_name, warming_lookup$W_z[[i]], points = 81L
    ) |>
      dplyr::mutate(warming_panel = warming_lookup$warming_panel[[i]])
  }))
})) |>
  dplyr::filter(.data$grid_point_type == "regular_grid") |>
  dplyr::left_join(
    dplyr::select(classification, "family_id", "fitted_class", "n_experiments"),
    by = "family_id"
  ) |>
  dplyr::mutate(
    family_panel = factor(
      paste0(.data$family_id, "\n", .data$fitted_class),
      levels = classification |>
        dplyr::filter(.data$family_id %in% .env$selected_families) |>
        dplyr::arrange(
          factor(.data$fitted_class, levels = names(class_colours)),
          .data$family_id
        ) |>
        dplyr::transmute(
          panel = paste0(.data$family_id, "\n", .data$fitted_class)
        ) |>
        dplyr::pull(.data$panel)
    )
  )
T_scaling <- species_scaled$scaling |>
  dplyr::filter(.data$scaled_field == "T_z")
P_scaling <- species_scaled$scaling |>
  dplyr::filter(.data$scaled_field == "P_z")
raw_T_labels <- function(x) sprintf(
  "%.1f", x * T_scaling$scale[[1]] + T_scaling$center[[1]]
)
raw_P_labels <- function(x) sprintf(
  "%.0f", exp(x * P_scaling$scale[[1]] + P_scaling$center[[1]])
)
observed_family_climates <- model_data |>
  dplyr::filter(as.character(.data$family_id) %in% .env$selected_families) |>
  dplyr::transmute(
    family_id = as.character(.data$family_id),
    T_z = .data$T_z, P_z = .data$P_z
  ) |>
  dplyr::distinct() |>
  dplyr::left_join(
    dplyr::select(classification, "family_id", "fitted_class"),
    by = "family_id"
  ) |>
  dplyr::mutate(
    family_panel = factor(
      paste0(.data$family_id, "\n", .data$fitted_class),
      levels = levels(surface_predictions$family_panel)
    )
  ) |>
  tidyr::crossing(warming_panel = warming_lookup$warming_panel)
surface_limit <- max(abs(surface_predictions$estimate), na.rm = TRUE)
surface_plot <- ggplot2::ggplot(
  surface_predictions,
  ggplot2::aes(x = .data$T_z, y = .data$P_z, fill = .data$estimate)
) +
  ggplot2::geom_tile(
    ggplot2::aes(width = .data$tile_width, height = .data$tile_height)
  ) +
  ggplot2::geom_point(
    data = observed_family_climates,
    ggplot2::aes(x = .data$T_z, y = .data$P_z),
    inherit.aes = FALSE, shape = 21, fill = "white", colour = "grey20",
    size = 1.0, stroke = 0.3
  ) +
  ggplot2::facet_grid(warming_panel ~ family_panel) +
  ggplot2::scale_fill_gradient2(
    low = "#E66101", mid = "white", high = "#1B9E77",
    midpoint = 0, limits = c(-surface_limit, surface_limit),
    name = "Fitted family\nDelta AUC"
  ) +
  ggplot2::scale_x_continuous(labels = raw_T_labels) +
  ggplot2::scale_y_continuous(labels = raw_P_labels) +
  ggplot2::coord_cartesian(expand = FALSE) +
  ggplot2::theme_bw(base_size = 8.5) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold", size = 7.5),
    axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
    legend.position = "right"
  ) +
  ggplot2::labs(
    title = "Selected family response surfaces within family-specific observed T-P support",
    x = "Destination JJA temperature (C)",
    y = "Destination JJA precipitation"
  )
ggplot2::ggsave(
  output_files$surface_figure, surface_plot,
  width = max(11, 2.5 * length(selected_families)), height = 7.0,
  limitsize = FALSE, device = grDevices::cairo_pdf
)

# Concise report with the exact model and requested interpretation scope.
fmt <- function(x) formatC(x, digits = 6, format = "fg")
support_counts <- table(factor(
  family_support$support_class,
  levels = c("full-rank", "partial-rank", "single-point")
))
report_lines <- c(
  "FAMILY-AGGREGATED RESIDENT DELTA-AUC ANALYSIS",
  "============================================",
  "",
  "Scientific question: If biological resolution is changed from species to",
  "family, do we still recover consistent winners, consistent losers, and",
  "background-climate-dependent responders?",
  "",
  "RESPONSE CONSTRUCTION",
  "---------------------",
  paste0(
    "Started from the exact 1,604 maximum-duration resident species x experiment ",
    "effects in the updated downstream analysis set across 42 experiments. The ",
    "six prespecified non-angiosperm/conifer taxa, colonizers, and nonresidents ",
    "were excluded; their upstream winner/loser calculations were retained."
  ),
  paste0(
    "Within each experiment, family relative cover was the sum of member resident ",
    "species relative cover for every treatment x plot x year. Family membership ",
    "was held constant through time; absent family cover was completed with zero ",
    "using the species-pipeline completion procedure."
  ),
  paste0(
    "Family effects were recalculated from plot-level trapezoidal AUC, not by ",
    "averaging species effects. Response = warmed family AUC minus origin-control ",
    "family AUC."
  ),
  "",
  "SAMPLING COVARIANCE",
  "-------------------",
  paste0("Family effects: ", nrow(family_auc_data), "; families: ",
         dplyr::n_distinct(family_auc_data$family_id), "."),
  paste0(
    "The full within-experiment emmeans covariance was reconstructed in exact ",
    "effect_id order. max|diag(V)-SE^2| = ", fmt(diagonal_error),
    "; maximum symmetry error = ", fmt(symmetry_error),
    "; minimum eigenvalue = ", fmt(min(V_eigenvalues)), "."
  ),
  "No eigenvalue repair or diagonal fallback was used.",
  "",
  "CLIMATE SCALING",
  "---------------",
  paste0(
    "T_z, P_z and W_z use the unchanged species-analysis scaling calculated over ",
    "the 42 distinct experiments. TP_z = T_z * P_z."
  ),
  "",
  "MODEL",
  "-----",
  "M_family_aggregated <- metafor::rma.mv(",
  "  yi = family_auc_effect,",
  "  V = V_family_auc,",
  "  mods = ~ T_z * P_z * W_z,",
  "  random = list(",
  "    ~ 1 | family_id,",
  "    ~ 0 + W_z | family_id,",
  "    ~ 0 + T_z + P_z + TP_z | family_id",
  "  ),",
  "  struct = \"GEN\",",
  "  data = family_auc_data,",
  "  method = \"REML\",",
  "  sparse = TRUE",
  ")",
  "",
  "MODEL RESULTS",
  "-------------",
  paste0("Optimizer attempts: ", paste(model_attempts, collapse = " -> "),
         "; retained optimizer: ", optimizer, "; code: ", optimizer_code, "."),
  paste0("Warnings: ", ifelse(length(fit_capture$warnings),
                               paste(fit_capture$warnings, collapse = " | "),
                               "none"), "."),
  paste0("Family intercept variance: ", fmt(family_intercept_variance), "."),
  paste0("Family W-slope variance: ", fmt(family_W_variance), "."),
  paste0(
    "Family T/P/T*P variances: ",
    paste(fmt(diag(context_covariance)), collapse = ", "), "."
  ),
  paste0(
    "T-P, T-(T*P), P-(T*P) correlations: ",
    paste(fmt(context_correlation[upper.tri(context_correlation)]),
          collapse = ", "), "."
  ),
  paste0(
    "Full random-covariance eigenvalues: ",
    paste(fmt(full_random_eigenvalues), collapse = ", "), "."
  ),
  paste0("Random structure identifiable by the prespecified numerical checks: ",
         random_identifiable, "."),
  "The model was not simplified automatically.",
  "",
  "FAMILY REPLICATION AND CLASSIFICATION",
  "-------------------------------------",
  paste0(
    "Support counts: full-rank = ", support_counts[["full-rank"]],
    "; partial-rank = ", support_counts[["partial-rank"]],
    "; single-point = ", support_counts[["single-point"]], "."
  ),
  paste0(
    "Individual-family classification was restricted to ", nrow(classification),
    " families with >=4 experiments and full-rank [1,T,P,T*P] support."
  ),
  paste0("Warming was fixed at the species-analysis median Delta T = ",
         sprintf("%.2f", W_median), " C."),
  paste(capture.output(print(comparison, n = Inf, width = Inf)),
        collapse = "\n"),
  paste0("Selected surface families: ",
         paste(selected_families, collapse = ", "), "."),
  "Predictions were restricted to each family's own observed T-P convex hull.",
  "Classification is descriptive fitted-mean classification; whiskers are",
  "pointwise model-based intervals at the fitted extrema, not simultaneous",
  "confidence regions for the entire family surface.",
  "",
  "OUTPUTS",
  "-------",
  vapply(output_files, function(path) {
    file.path(normalizePath(dirname(path)), basename(path))
  }, character(1))
)
writeLines(report_lines, output_files$report, useBytes = TRUE)

# Exact requested output contract: nine files, all nonempty.
actual_files <- list.files(output_dir, full.names = TRUE)
expected_files <- unname(unlist(output_files))
unexpected <- setdiff(actual_files, expected_files)
missing <- setdiff(expected_files, actual_files)
empty <- expected_files[file.exists(expected_files) &
                          file.info(expected_files)$size <= 0]
if (length(unexpected) || length(missing) || length(empty)) {
  stop(
    "Output contract failed. Unexpected: ", paste(unexpected, collapse = ", "),
    "; missing: ", paste(missing, collapse = ", "),
    "; empty: ", paste(empty, collapse = ", "), call. = FALSE
  )
}

cat("\nFamily-aggregated AUC analysis completed.\n")
cat("Output directory: ", normalizePath(output_dir), "\n", sep = "")
print(comparison)
