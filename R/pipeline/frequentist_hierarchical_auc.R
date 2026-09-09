# Resident continuous-AUC hierarchical meta-regression ----------------------

auc_required_packages <- function() {
  c("clubSandwich", "dplyr", "emmeans", "geometry", "ggplot2", "Matrix",
    "knitr", "metafor", "patchwork", "purrr", "readr", "rmarkdown",
    "sandwich", "scales", "tibble", "tidyr")
}

assert_auc_packages <- function() {
  missing <- auc_required_packages()[
    !vapply(auc_required_packages(), requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing)) {
    stop(
      "Install the resident AUC analysis dependencies before running: ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  invisible(TRUE)
}

auc_capture <- function(expression) {
  warnings <- character()
  error <- ""
  fit <- tryCatch(
    withCallingHandlers(
      expression,
      warning = function(condition) {
        warnings <<- c(warnings, conditionMessage(condition))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(condition) {
      error <<- conditionMessage(condition)
      NULL
    }
  )
  list(fit = fit, warnings = unique(warnings), error = error)
}

auc_formula_text <- function(formula) paste(deparse(formula), collapse = " ")

auc_value_or <- function(x, fallback) {
  if (is.null(x) || !length(x)) fallback else x
}

auc_core_formula <- function() ~ T_z * P_z * W_z

auc_context_formula <- function() ~ climate_context * W_z

auc_reduced_climate_formula <- function() ~ T_z * P_z + W_z

auc_random_structure <- function(kind) {
  switch(
    kind,
    species_slope = list(~ 1 | species_id, ~ 0 + W_z | species_id),
    species_slope_experiment = list(
      ~ 1 | species_id, ~ 0 + W_z | species_id, ~ 1 | experiment_id
    ),
    species_slope_gradient = list(
      ~ 1 | species_id, ~ 0 + W_z | species_id, ~ 1 | gradient_id
    ),
    species_intercept = list(~ 1 | species_id),
    stop("Unknown random structure: ", kind, call. = FALSE)
  )
}

auc_column_map <- function() {
  tibble::tribble(
    ~analysis_field, ~source_field, ~definition,
    "auc_effect", "estimate", "raw warmed AUC minus raw origin-control AUC",
    "auc_se", "SE", "standard error of warmed-minus-origin AUC contrast",
    "species_id", "accepted_name", "reviewed accepted botanical name",
    "experiment_id", "experiment_id", "warmed origin-to-destination comparison",
    "gradient_id", "region", "study/gradient identifier",
    "T", "destination_T_JJA", "destination JJA mean temperature",
    "P", "log(destination_P_JJA)", "log destination JJA precipitation",
    "W", "diff_T_JJA", "destination minus origin JJA temperature",
    "duration", "year_range", "inclusive number of analyzed survey years",
    "integration_span", "year_range - 1", "trapezoid integration time span",
    "plot_area", "plot_size", "experimental plot area",
    "LDMC", "LDMC", "species median; analyzed on log scale",
    "LMA", "LMA", "species median; analyzed on log scale",
    "vegetative_height", "vegetative_height", "species median; analyzed on log scale"
  )
}

read_resident_auc_data <- function(root, horizons = c("max", "3")) {
  result_path <- file.path(root, "output", "tables", "winner_loser_results.csv")
  metadata_path <- file.path(root, "output", "tables", "analysis_dataset.csv")
  context_path <- file.path(root, "output", "tables", "climate_context_experiments.csv")
  for (path in c(result_path, metadata_path, context_path)) {
    if (!file.exists(path)) stop("Missing resident AUC input: ", path, call. = FALSE)
  }

  results <- readr::read_csv(result_path, show_col_types = FALSE) |>
    dplyr::mutate(horizon = as.character(.data$horizon)) |>
    dplyr::filter(
      .data$transform == "raw", .data$contrast == "to_origin",
      .data$pool == "resident", .data$horizon %in% as.character(horizons)
    ) |>
    dplyr::transmute(
      horizon = .data$horizon,
      experiment_id = as.character(.data$experiment_id),
      source_species = as.character(.data$species),
      auc_effect = as.numeric(.data$estimate),
      auc_se = as.numeric(.data$SE)
    )

  metadata_required <- c(
    "experiment_id", "species", "accepted_name", "region", "year_range",
    "plot_size", "destination_T_JJA", "destination_P_JJA", "diff_T_JJA",
    "LDMC", "LMA", "vegetative_height", "climate_complete",
    "climate_summary_status"
  )
  metadata <- readr::read_csv(metadata_path, show_col_types = FALSE)
  assert_columns(metadata, metadata_required, "analysis_dataset.csv")
  metadata <- metadata |>
    dplyr::transmute(
      experiment_id = as.character(.data$experiment_id),
      source_species = as.character(.data$species),
      species_id = as.character(.data$accepted_name),
      gradient_id = as.character(.data$region),
      duration = as.numeric(.data$year_range),
      plot_area = as.numeric(.data$plot_size),
      T = as.numeric(.data$destination_T_JJA),
      P_raw = as.numeric(.data$destination_P_JJA),
      W = as.numeric(.data$diff_T_JJA),
      LDMC = as.numeric(.data$LDMC),
      LMA = as.numeric(.data$LMA),
      vegetative_height = as.numeric(.data$vegetative_height),
      climate_complete = as.logical(.data$climate_complete),
      climate_summary_status = as.character(.data$climate_summary_status)
    ) |>
    dplyr::distinct()
  assert_unique(metadata, c("experiment_id", "source_species"), "resident metadata")

  context <- readr::read_csv(context_path, show_col_types = FALSE) |>
    dplyr::transmute(
      experiment_id = as.character(.data$experiment_id),
      climate_context = factor(
        as.character(.data$climate_context),
        levels = c(
          "cold_dry", "intermediate_dry", "warm_dry",
          "cold_intermediate", "cold_wet"
        )
      )
    )
  assert_unique(context, "experiment_id", "climate context")

  data <- results |>
    dplyr::left_join(metadata, by = c("experiment_id", "source_species")) |>
    dplyr::left_join(context, by = "experiment_id") |>
    dplyr::mutate(
      P = dplyr::if_else(.data$P_raw > 0, log(.data$P_raw), NA_real_),
      log_plot_area = dplyr::if_else(.data$plot_area > 0, log(.data$plot_area), NA_real_),
      integration_span = pmax(.data$duration - 1, 1),
      effect_id = paste(.data$horizon, .data$experiment_id, .data$species_id, sep = "::")
    ) |>
    dplyr::arrange(.data$horizon, .data$experiment_id, .data$species_id)

  # Winner/loser contrasts remain intact upstream. The six declared
  # non-angiosperm taxa are removed only at this downstream model boundary.
  data <- filter_downstream_model_species(data, "species_id")

  required_complete <- c(
    "species_id", "gradient_id", "auc_effect", "auc_se", "T", "P", "W",
    "duration", "log_plot_area"
  )
  if (any(!stats::complete.cases(data[, required_complete]))) {
    stop("Resident AUC rows have missing required response or design fields.", call. = FALSE)
  }
  if (any(!is.finite(data$auc_se) | data$auc_se <= 0)) {
    stop("Resident AUC standard errors must be positive and finite.", call. = FALSE)
  }
  assert_unique(data, "effect_id", "resident AUC data")
  data
}

scale_resident_auc_data <- function(data) {
  experiments <- data |>
    dplyr::select(dplyr::all_of(c(
      "experiment_id", "gradient_id", "T", "P", "W", "duration",
      "log_plot_area"
    ))) |>
    dplyr::distinct()
  assert_unique(experiments, "experiment_id", "experiment-level predictors")
  variables <- c(
    T_z = "T", P_z = "P", W_z = "W", duration_z = "duration",
    log_plot_area_z = "log_plot_area"
  )
  scaling <- purrr::imap_dfr(variables, function(source, scaled) {
    values <- experiments[[source]]
    tibble::tibble(
      scaled_field = scaled, source_field = source,
      center = mean(values), scale = stats::sd(values),
      weighting = "distinct experiments", n = length(values)
    )
  })
  if (any(!is.finite(scaling$scale) | scaling$scale <= 0)) {
    stop("All experiment-level variables need a positive finite scale.", call. = FALSE)
  }
  for (i in seq_len(nrow(scaling))) {
    data[[scaling$scaled_field[[i]]]] <-
      (data[[scaling$source_field[[i]]]] - scaling$center[[i]]) /
      scaling$scale[[i]]
  }
  data$species_id <- factor(data$species_id)
  data$experiment_id <- factor(data$experiment_id)
  data$gradient_id <- factor(data$gradient_id)
  list(data = data, scaling = scaling)
}

trait_complete_data <- function(data, trait) {
  if (!trait %in% c("LDMC", "LMA", "vegetative_height")) {
    stop("Unknown focal trait: ", trait, call. = FALSE)
  }
  species_values <- data |>
    dplyr::transmute(
      species_id = as.character(.data$species_id),
      trait_value = as.numeric(.data[[trait]])
    ) |>
    dplyr::distinct()
  assert_unique(species_values, "species_id", paste(trait, "species values"))
  complete_species <- species_values |>
    dplyr::filter(is.finite(.data$trait_value), .data$trait_value > 0) |>
    dplyr::mutate(log_trait = log(.data$trait_value))
  center <- mean(complete_species$log_trait)
  scale <- stats::sd(complete_species$log_trait)
  if (!is.finite(scale) || scale <= 0) stop(trait, " has no usable scale.", call. = FALSE)
  result <- data |>
    dplyr::filter(as.character(.data$species_id) %in% complete_species$species_id) |>
    dplyr::mutate(
      log_trait = log(.data[[trait]]),
      trait_z = (.data$log_trait - .env$center) / .env$scale
    ) |>
    droplevels()
  list(
    data = result,
    scaling = tibble::tibble(
      trait = trait, source_field = trait, transform = "natural log",
      center = center, scale = scale, weighting = "distinct species",
      n = nrow(complete_species)
    )
  )
}

validate_auc_covariance <- function(V, data, tolerance = 1e-8) {
  dimensions_valid <- is.matrix(V) && identical(dim(V), c(nrow(data), nrow(data)))
  if (!dimensions_valid) {
    return(tibble::tibble(
      dimensions_valid = FALSE, order_valid = FALSE,
      maximum_symmetry_error = NA_real_, minimum_eigenvalue = NA_real_,
      maximum_diagonal_se_error = NA_real_, valid = FALSE
    ))
  }
  order_valid <- !is.null(rownames(V)) && !is.null(colnames(V)) &&
    identical(rownames(V), data$effect_id) && identical(colnames(V), data$effect_id)
  symmetry_error <- max(abs(V - t(V)))
  eigen_min <- if (symmetry_error <= tolerance) {
    min(eigen(V, symmetric = TRUE, only.values = TRUE)$values)
  } else {
    NA_real_
  }
  diagonal_error <- if (all(diag(V) >= 0)) {
    max(abs(sqrt(diag(V)) - data$auc_se))
  } else {
    Inf
  }
  tibble::tibble(
    dimensions_valid = dimensions_valid, order_valid = order_valid,
    maximum_symmetry_error = symmetry_error, minimum_eigenvalue = eigen_min,
    maximum_diagonal_se_error = diagonal_error,
    valid = dimensions_valid && order_valid && symmetry_error <= tolerance &&
      is.finite(eigen_min) && eigen_min >= -tolerance && diagonal_error <= tolerance
  )
}

select_auc_covariance <- function(full_V, data, tolerance = 1e-8,
                                  reconstruction_error = "") {
  validation <- validate_auc_covariance(full_V, data, tolerance)
  if (isTRUE(validation$valid[[1]])) {
    computational_V <- Matrix::drop0(
      Matrix::Matrix(full_V, sparse = TRUE), tol = tolerance
    )
    return(list(
      V = computational_V,
      selected_v_type = paste0(
        "full emmeans covariance reconstructed from per-experiment AUC models; ",
        "entries below ", format(tolerance, scientific = TRUE),
        " stored as exact sparse zeros"
      ),
      fallback_reason = "", validation = validation
    ))
  }
  if (any(!is.finite(data$auc_se) | data$auc_se <= 0)) {
    stop("Full V failed and diagonal sampling variances are invalid.", call. = FALSE)
  }
  diagonal <- diag(data$auc_se^2)
  dimnames(diagonal) <- list(data$effect_id, data$effect_id)
  list(
    V = diagonal, selected_v_type = "diagonal sampling-variance approximation",
    fallback_reason = paste(
      c(reconstruction_error, "full covariance validation failed")[
        nzchar(c(reconstruction_error, "full covariance validation failed"))
      ],
      collapse = " | "
    ),
    validation = validation
  )
}

reconstruct_auc_covariance <- function(data, auc_by_plot, tolerance = 1e-8) {
  horizon_value <- unique(as.character(data$horizon))
  if (length(horizon_value) != 1L) stop("Reconstruct V one horizon at a time.", call. = FALSE)
  blocks <- list()
  ordered_ids <- character()
  block_audit <- list()
  errors <- character()

  for (experiment in levels(droplevels(data$experiment_id))) {
    target <- data |>
      dplyr::filter(as.character(.data$experiment_id) == .env$experiment) |>
      dplyr::arrange(.data$source_species)
    raw <- auc_by_plot |>
      dplyr::filter(
        as.character(.data$experiment_id) == .env$experiment,
        as.character(.data$horizon) == .env$horizon_value,
        as.character(.data$transform) == "raw", !is.na(.data$auc)
      ) |>
      dplyr::mutate(
        ODT = factor(
          as.character(.data$ODT),
          levels = c("destination_controls", "origin_controls", "warmed")
        )
      )
    captured <- auc_capture({
      model <- stats::lm(auc ~ ODT * species, data = raw)
      means <- emmeans::emmeans(model, ~ ODT | species)
      contrast <- emmeans::contrast(
        means, method = list(to_origin = c(0, -1, 1)), adjust = "none"
      )
      summary <- tibble::as_tibble(as.data.frame(summary(contrast)))
      covariance <- as.matrix(stats::vcov(contrast))
      index <- match(target$source_species, as.character(summary$species))
      if (anyNA(index)) stop("Resident species missing from reconstructed contrasts.")
      list(summary = summary[index, , drop = FALSE], covariance = covariance[index, index, drop = FALSE])
    })
    if (is.null(captured$fit)) {
      errors <- c(errors, paste(experiment, captured$error, sep = ": "))
      next
    }
    reconstructed <- captured$fit
    estimate_error <- max(abs(reconstructed$summary$estimate - target$auc_effect))
    se_error <- max(abs(reconstructed$summary$SE - target$auc_se))
    symmetry_error <- max(abs(reconstructed$covariance - t(reconstructed$covariance)))
    eigen_min <- min(eigen(reconstructed$covariance, symmetric = TRUE, only.values = TRUE)$values)
    blocks[[experiment]] <- reconstructed$covariance
    ordered_ids <- c(ordered_ids, target$effect_id)
    block_audit[[experiment]] <- tibble::tibble(
      horizon = horizon_value, experiment_id = experiment, n_effects = nrow(target),
      maximum_estimate_error = estimate_error, maximum_se_error = se_error,
      maximum_symmetry_error = symmetry_error, minimum_eigenvalue = eigen_min,
      valid = estimate_error <= tolerance && se_error <= tolerance &&
        symmetry_error <= tolerance && eigen_min >= -tolerance,
      warning = paste(captured$warnings, collapse = " | "), error = ""
    )
  }

  full_valid <- !length(errors) && length(blocks) == dplyr::n_distinct(data$experiment_id)
  V_full <- NULL
  if (full_valid) {
    assembled <- as.matrix(Matrix::bdiag(blocks))
    index <- match(data$effect_id, ordered_ids)
    full_valid <- !anyNA(index)
    if (full_valid) V_full <- assembled[index, index, drop = FALSE]
  }
  audit_blocks <- dplyr::bind_rows(block_audit)
  full_valid <- full_valid && nrow(audit_blocks) && all(audit_blocks$valid)
  if (full_valid) dimnames(V_full) <- list(data$effect_id, data$effect_id)
  selected <- select_auc_covariance(
    if (full_valid) V_full else matrix(numeric(), 0, 0), data, tolerance,
    paste(errors, collapse = " | ")
  )
  V <- selected$V
  v_type <- selected$selected_v_type
  fallback_reason <- selected$fallback_reason
  validation <- selected$validation
  full_valid <- isTRUE(validation$valid[[1]])
  dim_valid <- validation$dimensions_valid[[1]]
  symmetry_error <- validation$maximum_symmetry_error[[1]]
  eigen_min <- validation$minimum_eigenvalue[[1]]
  diagonal_error <- validation$maximum_diagonal_se_error[[1]]
  off_diagonal <- if (!is.null(V_full) && length(V_full)) {
    V_full[row(V_full) != col(V_full)]
  } else {
    V[row(V) != col(V)]
  }
  audit <- tibble::tibble(
    horizon = horizon_value, n_effects = nrow(data), n_experiments = length(blocks),
    dimensions_valid = dim_valid, maximum_symmetry_error = symmetry_error,
    minimum_eigenvalue = eigen_min, maximum_diagonal_se_error = diagonal_error,
    maximum_absolute_off_diagonal = if (length(off_diagonal)) max(abs(off_diagonal)) else 0,
    full_covariance_valid = full_valid, selected_v_type = v_type,
    fallback_reason = fallback_reason
  )
  list(V = V, audit = audit, block_audit = audit_blocks, selected_v_type = v_type)
}

subset_auc_v <- function(V, effect_ids) {
  if (is.null(rownames(V)) || any(!effect_ids %in% rownames(V))) {
    stop("Cannot align V to requested effect_id values.", call. = FALSE)
  }
  V[effect_ids, effect_ids, drop = FALSE]
}

auc_model_specifications <- function() {
  tibble::tribble(
    ~model_id, ~formula, ~random_kind, ~fallback_level, ~purpose,
    "S1_PRIMARY", auc_core_formula(), "species_slope", 1L,
    "primary biological species hierarchy",
    "S2_EXP_RE", auc_core_formula(), "species_slope_experiment", 3L,
    "experiment random-intercept sensitivity",
    "S3_GRADIENT_RE", auc_core_formula(), "species_slope_gradient", 3L,
    "gradient/study random-intercept sensitivity",
    "S4_DESIGN", ~ T_z * P_z * W_z + duration_z + log_plot_area_z,
    "species_slope", 1L, "duration and plot-area sensitivity",
    "S5_SPECIES_INTERCEPT_ONLY", auc_core_formula(), "species_intercept", 2L,
    "species-slope fallback and comparison"
  )
}

fit_auc_rma <- function(data, V, formula, random_kind, method = "REML",
                        timeout_seconds = 45) {
  design <- stats::model.matrix(formula, data = data)
  design_rank <- qr(design)$rank
  condition_number <- tryCatch(kappa(design), error = function(condition) NA_real_)
  captured <- auc_capture({
    setTimeLimit(elapsed = timeout_seconds, transient = TRUE)
    on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE), add = TRUE)
    arguments <- list(
      yi = data$auc_effect, V = V, mods = formula,
      random = auc_random_structure(random_kind), data = data,
      method = method, sparse = TRUE
    )
    if (grepl("slope", random_kind, fixed = TRUE)) arguments$struct <- "GEN"
    do.call(metafor::rma.mv, arguments)
  })
  fit <- captured$fit
  covariance <- if (is.null(fit)) NULL else
    tryCatch(stats::vcov(fit), error = function(condition) NULL)
  coefficients <- if (is.null(fit)) numeric() else stats::coef(fit)
  finite <- length(coefficients) && all(is.finite(coefficients)) &&
    !is.null(covariance) && all(is.finite(covariance))
  slope_variance <- if (!is.null(fit) && grepl("slope", random_kind, fixed = TRUE) &&
                        length(fit$tau2)) as.numeric(fit$tau2[[1]]) else NA_real_
  boundary_tolerance <- 1e-8 * max(1, stats::var(data$auc_effect))
  slope_identifiable <- !grepl("slope", random_kind, fixed = TRUE) ||
    (is.finite(slope_variance) && slope_variance > boundary_tolerance)
  stable <- !is.null(fit) && !length(captured$warnings) && !nzchar(captured$error) &&
    finite && design_rank == ncol(design) && slope_identifiable
  list(
    fit = fit, data = data, V = V, formula = formula,
    random_kind = random_kind, method = method,
    warnings = captured$warnings, error = captured$error,
    design_rank = design_rank, design_columns = ncol(design),
    condition_number = condition_number, slope_variance = slope_variance,
    boundary_tolerance = boundary_tolerance,
    slope_identifiable = slope_identifiable, stable = stable
  )
}

profile_auc_slope <- function(model_result, steps = 6L) {
  if (is.null(model_result$fit) || !grepl("slope", model_result$random_kind, fixed = TRUE)) {
    return(list(profile = NULL, diagnostic = "no species warming slope"))
  }
  captured <- auc_capture(
    metafor::profile.rma.mv(
      model_result$fit, tau2 = 1, steps = steps,
      progbar = FALSE, plot = FALSE
    )
  )
  diagnostic <- if (is.null(captured$fit)) captured$error else
    paste(captured$warnings, collapse = " | ")
  list(profile = captured$fit, diagnostic = diagnostic)
}

auc_profile_is_informative <- function(profile_result, likelihood_drop = 1.92) {
  profile <- profile_result$profile
  !is.null(profile) && length(profile$ll) >= 3L &&
    all(is.finite(profile$ll)) && diff(range(profile$ll)) >= likelihood_drop
}

auc_variance_components <- function(model_result, model_id) {
  fit <- model_result$fit
  if (is.null(fit)) return(tibble::tibble())
  sigma_labels <- switch(
    model_result$random_kind,
    species_slope = "species_intercept",
    species_slope_experiment = c("species_intercept", "experiment_intercept"),
    species_slope_gradient = c("species_intercept", "gradient_intercept"),
    species_intercept = "species_intercept"
  )
  dplyr::bind_rows(
    tibble::tibble(
      model_id = model_id, component = sigma_labels,
      estimate = as.numeric(fit$sigma2)
    ),
    if (grepl("slope", model_result$random_kind, fixed = TRUE)) {
      tibble::tibble(
        model_id = model_id, component = "species_W_slope",
        estimate = as.numeric(fit$tau2[[1]])
      )
    }
  )
}

tidy_auc_model <- function(model_result, model_id) {
  fit <- model_result$fit
  if (is.null(fit)) return(tibble::tibble())
  estimate <- as.numeric(stats::coef(fit))
  se <- sqrt(diag(stats::vcov(fit)))
  tibble::tibble(
    model_id = model_id, inference = "model_based", term = names(stats::coef(fit)),
    estimate = estimate, standard_error = se,
    statistic = estimate / se, denominator_df = Inf,
    p_value = 2 * stats::pnorm(abs(.data$statistic), lower.tail = FALSE),
    conf_low = estimate - stats::qnorm(0.975) * se,
    conf_high = estimate + stats::qnorm(0.975) * se,
    df_flag = "model-based asymptotic"
  )
}

robust_auc_model <- function(model_result, model_id) {
  if (is.null(model_result$fit)) {
    return(list(model = NULL, table = tibble::tibble(), error = "model unavailable"))
  }
  cr2 <- auc_capture(
    metafor::robust(
      model_result$fit, cluster = model_result$data$experiment_id,
      clubSandwich = TRUE
    )
  )
  selected_method <- "experiment_clustered_CR2"
  captured <- cr2
  cr2_reason <- cr2$error
  if (is.null(cr2$fit)) {
    if (grepl("Could not obtain the cluster-robust", cr2_reason, fixed = TRUE) &&
        grepl("species", model_result$random_kind, fixed = TRUE)) {
      cr2_reason <- paste(
        "clubSandwich CR2 requires random effects nested within the experiment",
        "cluster; species effects cross experiments in this model"
      )
    }
    captured <- auc_capture(
      metafor::robust(
        model_result$fit, cluster = model_result$data$experiment_id,
        adjust = TRUE, clubSandwich = FALSE
      )
    )
    selected_method <- "experiment_clustered_CR1_fallback"
  }
  if (is.null(captured$fit)) {
    return(list(
      model = NULL, table = tibble::tibble(),
      error = paste(
        "CR2:", cr2_reason, "| CR1:", captured$error
      ),
      requested_method = "CR2", selected_method = "none"
    ))
  }
  robust <- captured$fit
  df <- as.numeric(robust$ddf)
  table <- tibble::tibble(
    model_id = model_id, inference = selected_method,
    term = names(stats::coef(model_result$fit)),
    estimate = as.numeric(robust$beta), standard_error = as.numeric(robust$se),
    statistic = as.numeric(robust$zval), denominator_df = df,
    p_value = as.numeric(robust$pval), conf_low = as.numeric(robust$ci.lb),
    conf_high = as.numeric(robust$ci.ub),
    df_flag = dplyr::case_when(
      is.na(.data$denominator_df) ~ "unavailable",
      .data$denominator_df < 4 ~ "extremely weak (<4)",
      .data$denominator_df < 10 ~ "limited (<10)",
      TRUE ~ "adequate (>=10)"
    )
  )
  list(
    model = robust, table = table,
    error = paste(
      c(
        if (is.null(cr2$fit)) paste0("CR2 unavailable: ", cr2_reason) else character(),
        captured$warnings, captured$error
      ),
      collapse = " | "
    ),
    requested_method = "CR2", selected_method = selected_method
  )
}

auc_block_test <- function(model_result, robust_result, terms, block_name) {
  if (is.null(model_result$fit) || is.null(robust_result$model)) {
    return(tibble::tibble(
      block = block_name, terms = paste(terms, collapse = " | "),
      statistic = NA_real_, numerator_df = NA_real_, denominator_df = NA_real_,
      p_value = NA_real_, error = "model or CR2 result unavailable"
    ))
  }
  coefficient_names <- names(stats::coef(model_result$fit))
  indices <- match(terms, coefficient_names)
  if (anyNA(indices)) {
    return(tibble::tibble(
      block = block_name, terms = paste(terms, collapse = " | "),
      statistic = NA_real_, numerator_df = NA_real_, denominator_df = NA_real_,
      p_value = NA_real_, error = "one or more requested terms are absent"
    ))
  }
  if (identical(robust_result$selected_method, "experiment_clustered_CR2")) {
    captured <- auc_capture(
      clubSandwich::Wald_test(
        model_result$fit,
        constraints = clubSandwich::constrain_zero(indices),
        vcov = robust_result$model$vb,
        cluster = model_result$data$experiment_id,
        test = "HTZ"
      )
    )
  } else {
    beta <- as.numeric(stats::coef(model_result$fit))[indices]
    covariance <- as.matrix(robust_result$model$vb)[indices, indices, drop = FALSE]
    captured <- auc_capture({
      statistic <- as.numeric(crossprod(beta, solve(covariance, beta))) / length(indices)
      tibble::tibble(
        Fstat = statistic, df_num = length(indices),
        df_denom = unique(as.numeric(robust_result$model$ddf))[[1]],
        p_val = stats::pf(
          statistic, length(indices),
          unique(as.numeric(robust_result$model$ddf))[[1]], lower.tail = FALSE
        )
      )
    })
  }
  if (is.null(captured$fit)) {
    return(tibble::tibble(
      block = block_name, terms = paste(terms, collapse = " | "),
      statistic = NA_real_, numerator_df = length(indices), denominator_df = NA_real_,
      p_value = NA_real_, error = captured$error
    ))
  }
  result <- captured$fit
  tibble::tibble(
    block = block_name, terms = paste(terms, collapse = " | "),
    statistic = as.numeric(result$Fstat[[1]]),
    numerator_df = as.numeric(result$df_num[[1]]),
    denominator_df = as.numeric(result$df_denom[[1]]),
    p_value = as.numeric(result$p_val[[1]]),
    error = paste(captured$warnings, collapse = " | ")
  )
}

auc_fit_statistics <- function(fit) {
  if (is.null(fit)) {
    return(c(logLik = NA_real_, AIC = NA_real_, AICc = NA_real_, BIC = NA_real_))
  }
  values <- tryCatch(metafor::fitstats(fit), error = function(condition) NULL)
  if (is.null(values)) {
    return(c(logLik = NA_real_, AIC = NA_real_, AICc = NA_real_, BIC = NA_real_))
  }
  as.numeric(values[c("logLik", "AIC", "AICc", "BIC"), 1]) |>
    stats::setNames(c("logLik", "AIC", "AICc", "BIC"))
}

auc_model_ledger_row <- function(model_result, model_id, fallback_level,
                                 v_type, selected = FALSE, reason = "") {
  fit <- model_result$fit
  stats <- auc_fit_statistics(fit)
  variances <- auc_variance_components(model_result, model_id)
  variance <- function(component) {
    value <- variances$estimate[variances$component == component]
    if (length(value)) value[[1]] else NA_real_
  }
  tibble::tibble(
    model_id = model_id, formula = auc_formula_text(model_result$formula),
    random_structure = model_result$random_kind, method = model_result$method,
    fallback_level = fallback_level, V_type = v_type,
    n_observations = nrow(model_result$data),
    n_species = dplyr::n_distinct(model_result$data$species_id),
    n_experiments = dplyr::n_distinct(model_result$data$experiment_id),
    design_rank = model_result$design_rank,
    design_columns = model_result$design_columns,
    condition_number = model_result$condition_number,
    stable = model_result$stable,
    species_intercept_variance = variance("species_intercept"),
    species_W_slope_variance = variance("species_W_slope"),
    experiment_variance = variance("experiment_intercept"),
    gradient_variance = variance("gradient_intercept"),
    logLik = stats[["logLik"]], AIC = stats[["AIC"]],
    AICc = stats[["AICc"]], BIC = stats[["BIC"]],
    selected = selected, reason = reason,
    warning = paste(model_result$warnings, collapse = " | "),
    error = model_result$error
  )
}

fit_auc_model_sequence <- function(data, V, v_type) {
  specifications <- auc_model_specifications()
  models <- list()
  robust <- list()
  ledger <- list()
  coefficients <- list()
  variances <- list()
  block_tests <- list()
  for (i in seq_len(nrow(specifications))) {
    specification <- specifications[i, ]
    id <- specification$model_id[[1]]
    message("  Fitting ", id, " ...")
    timeout <- if (id %in% c("S2_EXP_RE", "S3_GRADIENT_RE")) 20 else 45
    model <- fit_auc_rma(
      data, V, specification$formula[[1]], specification$random_kind[[1]], "REML",
      timeout_seconds = timeout
    )
    message("  ", id, ": ", if (model$stable) "stable" else "not retained")
    models[[id]] <- model
    robust[[id]] <- robust_auc_model(model, id)
    coefficients[[id]] <- dplyr::bind_rows(
      tidy_auc_model(model, id), robust[[id]]$table
    )
    variances[[id]] <- auc_variance_components(model, id)
  }

  selected_id <- if (isTRUE(models$S1_PRIMARY$stable)) {
    "S1_PRIMARY"
  } else if (isTRUE(models$S5_SPECIES_INTERCEPT_ONLY$stable)) {
    "S5_SPECIES_INTERCEPT_ONLY"
  } else {
    NA_character_
  }

  fallback_models <- list()
  if (is.na(selected_id)) {
    context <- fit_auc_rma(data, V, auc_context_formula(), "species_intercept", "REML")
    fallback_models$F1_CLIMATE_CONTEXT <- context
    models$F1_CLIMATE_CONTEXT <- context
    robust$F1_CLIMATE_CONTEXT <- robust_auc_model(context, "F1_CLIMATE_CONTEXT")
    coefficients$F1_CLIMATE_CONTEXT <- dplyr::bind_rows(
      tidy_auc_model(context, "F1_CLIMATE_CONTEXT"),
      robust$F1_CLIMATE_CONTEXT$table
    )
    variances$F1_CLIMATE_CONTEXT <- auc_variance_components(context, "F1_CLIMATE_CONTEXT")
    if (isTRUE(context$stable)) selected_id <- "F1_CLIMATE_CONTEXT"
  }
  wls_fallback <- NULL
  if (is.na(selected_id)) {
    wls_fallback <- fit_auc_wls_fallback(data, auc_core_formula())
    if (isTRUE(wls_fallback$stable)) selected_id <- "F2_WLS_TWO_WAY_CLUSTER"
  }
  if (is.na(selected_id)) {
    reduced <- fit_auc_rma(data, V, auc_reduced_climate_formula(), "species_intercept", "REML")
    fallback_models$F3_REDUCED_DESCRIPTIVE <- reduced
    models$F3_REDUCED_DESCRIPTIVE <- reduced
    robust$F3_REDUCED_DESCRIPTIVE <- robust_auc_model(reduced, "F3_REDUCED_DESCRIPTIVE")
    coefficients$F3_REDUCED_DESCRIPTIVE <- dplyr::bind_rows(
      tidy_auc_model(reduced, "F3_REDUCED_DESCRIPTIVE"),
      robust$F3_REDUCED_DESCRIPTIVE$table
    )
    variances$F3_REDUCED_DESCRIPTIVE <- auc_variance_components(reduced, "F3_REDUCED_DESCRIPTIVE")
    if (isTRUE(reduced$stable)) selected_id <- "F3_REDUCED_DESCRIPTIVE"
  }

  for (i in seq_len(nrow(specifications))) {
    specification <- specifications[i, ]
    id <- specification$model_id[[1]]
    reason <- if (id == selected_id) {
      "selected by prespecified hierarchy"
    } else if (!models[[id]]$stable) {
      "rejected by convergence, rank, covariance, or variance-boundary diagnostics"
    } else {
      specification$purpose[[1]]
    }
    ledger[[id]] <- auc_model_ledger_row(
      models[[id]], id, specification$fallback_level[[1]], v_type,
      selected = id == selected_id, reason = reason
    )
  }
  if (length(fallback_models)) {
    fallback_levels <- c(F1_CLIMATE_CONTEXT = 4L, F3_REDUCED_DESCRIPTIVE = 8L)
    for (id in names(fallback_models)) {
      ledger[[id]] <- auc_model_ledger_row(
        models[[id]], id, fallback_levels[[id]], v_type,
        selected = id == selected_id,
        reason = if (id == selected_id) "selected fallback" else "fallback unsuccessful"
      )
    }
  }
  if (!is.null(wls_fallback)) {
    wls_stats <- if (is.null(wls_fallback$fit)) {
      c(logLik = NA_real_, AIC = NA_real_, BIC = NA_real_)
    } else {
      c(
        logLik = as.numeric(stats::logLik(wls_fallback$fit)),
        AIC = stats::AIC(wls_fallback$fit), BIC = stats::BIC(wls_fallback$fit)
      )
    }
    ledger$F2_WLS_TWO_WAY_CLUSTER <- tibble::tibble(
      model_id = "F2_WLS_TWO_WAY_CLUSTER", formula = auc_formula_text(auc_core_formula()),
      random_structure = "none; two-way clustered covariance",
      method = "inverse-variance weighted least squares", fallback_level = 6L,
      V_type = "1 / auc_se^2 weights", n_observations = nrow(data),
      n_species = dplyr::n_distinct(data$species_id),
      n_experiments = dplyr::n_distinct(data$experiment_id),
      design_rank = qr(stats::model.matrix(auc_core_formula(), data))$rank,
      design_columns = ncol(stats::model.matrix(auc_core_formula(), data)),
      condition_number = kappa(stats::model.matrix(auc_core_formula(), data)),
      stable = wls_fallback$stable, species_intercept_variance = NA_real_,
      species_W_slope_variance = NA_real_, experiment_variance = NA_real_,
      gradient_variance = NA_real_, logLik = wls_stats[["logLik"]],
      AIC = wls_stats[["AIC"]], AICc = NA_real_, BIC = wls_stats[["BIC"]],
      selected = selected_id == "F2_WLS_TWO_WAY_CLUSTER",
      reason = if (selected_id == "F2_WLS_TWO_WAY_CLUSTER")
        "selected last-resort fixed-effect fallback" else "fallback unsuccessful",
      warning = paste(wls_fallback$warnings, collapse = " | "),
      error = wls_fallback$error
    )
    coefficients$F2_WLS_TWO_WAY_CLUSTER <- wls_fallback$table |>
      dplyr::mutate(model_id = "F2_WLS_TWO_WAY_CLUSTER", .before = 1)
  }
  selected <- if (!is.na(selected_id) && selected_id %in% names(models)) {
    models[[selected_id]]
  } else {
    NULL
  }
  if (!is.null(selected) &&
      identical(auc_formula_text(selected$formula), auc_formula_text(auc_core_formula()))) {
    block_tests$climate_warming_modification <- auc_block_test(
      selected, robust[[selected_id]],
      c("T_z:W_z", "P_z:W_z", "T_z:P_z:W_z"),
      "background climate modifies warming response"
    )
  }
  ledger_table <- dplyr::bind_rows(ledger)
  coefficient_table <- dplyr::bind_rows(coefficients)
  robust_summary <- if ("denominator_df" %in% names(coefficient_table)) {
    coefficient_table |>
      dplyr::filter(startsWith(.data$inference, "experiment_clustered_")) |>
      dplyr::summarise(
        minimum_robust_df = min(.data$denominator_df, na.rm = TRUE),
        central_TPW_estimate = .data$estimate[.data$term == "T_z:P_z:W_z"][1],
        .by = "model_id"
      )
  } else {
    tibble::tibble(
      model_id = character(), minimum_robust_df = numeric(),
      central_TPW_estimate = numeric()
    )
  }
  selected_central <- robust_summary$central_TPW_estimate[
    robust_summary$model_id == selected_id
  ]
  ledger_table <- ledger_table |>
    dplyr::left_join(robust_summary, by = "model_id") |>
    dplyr::mutate(
      robust_method = vapply(
        .data$model_id,
        function(id) auc_value_or(robust[[id]]$selected_method, "none"), character(1)
      ),
      robust_note = vapply(
        .data$model_id,
        function(id) auc_value_or(robust[[id]]$error, ""), character(1)
      )
    ) |>
    dplyr::mutate(
      TPW_conclusion_changes = dplyr::if_else(
        length(.env$selected_central) == 1L & is.finite(.env$selected_central) &
          is.finite(.data$central_TPW_estimate),
        sign(.data$central_TPW_estimate) != sign(.env$selected_central), NA
      )
    )
  list(
    models = models, robust = robust, selected_id = selected_id,
    selected = selected, ledger = ledger_table,
    coefficients = coefficient_table,
    variance_components = dplyr::bind_rows(variances),
    block_tests = dplyr::bind_rows(block_tests), wls_fallback = wls_fallback,
    experiment_level_check = experiment_level_auc_check(data)
  )
}

auc_trait_formulas <- function() {
  list(
    H0_TRAIT_ONLY = ~ trait_z,
    H1_CLIMATE = auc_core_formula(),
    H2_ADDITIVE = ~ trait_z + T_z * P_z * W_z,
    H3_FULL = ~ trait_z * T_z * P_z * W_z,
    H3_FALLBACK_A_CONTEXT = ~ trait_z * climate_context * W_z,
    H3_FALLBACK_B_SIMPLE = ~ T_z * P_z * W_z + trait_z +
      trait_z:T_z + trait_z:P_z + trait_z:W_z,
    H3_FALLBACK_C_ADDITIVE = ~ trait_z + T_z * P_z * W_z,
    H3_FALLBACK_D_BENCHMARK = ~ trait_z
  )
}

compare_auc_ml_models <- function(first, second, comparison) {
  first_stats <- auc_fit_statistics(first$fit)
  second_stats <- auc_fit_statistics(second$fit)
  captured <- auc_capture(stats::anova(first$fit, second$fit))
  table <- captured$fit
  lrt <- p <- df <- NA_real_
  if (!is.null(table)) {
    candidates <- names(table)
    lrt_column <- intersect(c("LRT", "Chisq", "QM"), candidates)
    p_column <- intersect(c("pval", "Pr(>Chisq)", "QMp"), candidates)
    df_column <- intersect(c("df", "Df"), candidates)
    if (length(lrt_column)) lrt <- as.numeric(tail(table[[lrt_column[[1]]]], 1))
    if (length(p_column)) p <- as.numeric(tail(table[[p_column[[1]]]], 1))
    if (length(df_column)) df <- as.numeric(tail(table[[df_column[[1]]]], 1))
  }
  tibble::tibble(
    comparison = comparison,
    first_logLik = first_stats[["logLik"]], second_logLik = second_stats[["logLik"]],
    first_AIC = first_stats[["AIC"]], second_AIC = second_stats[["AIC"]],
    first_AICc = first_stats[["AICc"]], second_AICc = second_stats[["AICc"]],
    first_BIC = first_stats[["BIC"]], second_BIC = second_stats[["BIC"]],
    likelihood_ratio = lrt, numerator_df = df, p_value = p,
    valid = !is.null(table),
    warning = paste(c(captured$warnings, captured$error), collapse = " | ")
  )
}

fit_auc_trait_sequence <- function(data, V, trait, v_type) {
  message("Fitting trait sequence: ", trait)
  prepared <- trait_complete_data(data, trait)
  trait_data <- prepared$data
  trait_V <- subset_auc_v(V, trait_data$effect_id)
  formulas <- auc_trait_formulas()

  slope_probe <- fit_auc_rma(
    trait_data, trait_V, formulas$H1_CLIMATE, "species_slope", "REML"
  )
  random_kind <- if (isTRUE(slope_probe$stable)) "species_slope" else "species_intercept"
  models <- list(
    H0_TRAIT_ONLY = fit_auc_rma(
      trait_data, trait_V, formulas$H0_TRAIT_ONLY, "species_intercept", "REML"
    ),
    H1_CLIMATE = if (random_kind == "species_slope") slope_probe else
      fit_auc_rma(trait_data, trait_V, formulas$H1_CLIMATE, random_kind, "REML"),
    H2_ADDITIVE = fit_auc_rma(
      trait_data, trait_V, formulas$H2_ADDITIVE, random_kind, "REML"
    ),
    H3_FULL = fit_auc_rma(
      trait_data, trait_V, formulas$H3_FULL, random_kind, "REML"
    ),
    H3_FIXED_CONTEXT = fit_auc_rma(
      trait_data, trait_V, formulas$H3_FULL, "species_intercept", "REML"
    )
  )
  message("  ", trait, " REML H0-H3 fits completed")
  if (random_kind == "species_slope") {
    models$H3_HIERARCHICAL <- models$H3_FULL
  } else {
    models$H3_HIERARCHICAL <- fit_auc_rma(
      trait_data, trait_V, formulas$H3_FULL, "species_slope", "REML"
    )
  }

  selected_id <- if (isTRUE(models$H3_FULL$stable)) "H3_FULL" else NA_character_
  fallback_order <- c(
    "H3_FALLBACK_A_CONTEXT", "H3_FALLBACK_B_SIMPLE",
    "H3_FALLBACK_C_ADDITIVE", "H3_FALLBACK_D_BENCHMARK"
  )
  if (is.na(selected_id)) {
    for (id in fallback_order) {
      candidate_random <- if (id == "H3_FALLBACK_D_BENCHMARK") {
        "species_intercept"
      } else {
        random_kind
      }
      models[[id]] <- fit_auc_rma(
        trait_data, trait_V, formulas[[id]], candidate_random, "REML"
      )
      if (isTRUE(models[[id]]$stable)) {
        selected_id <- id
        break
      }
    }
  }

  ml <- purrr::map(
    c("H1_CLIMATE", "H2_ADDITIVE", "H3_FULL"),
    function(id) fit_auc_rma(
      trait_data, trait_V, formulas[[id]], random_kind, "ML"
    )
  ) |>
    stats::setNames(c("H1_CLIMATE", "H2_ADDITIVE", "H3_FULL"))
  message("  ", trait, " ML comparison fits completed")
  comparisons <- dplyr::bind_rows(
    compare_auc_ml_models(ml$H1_CLIMATE, ml$H2_ADDITIVE, "H1 vs H2"),
    compare_auc_ml_models(ml$H2_ADDITIVE, ml$H3_FULL, "H2 vs H3")
  ) |>
    dplyr::mutate(trait = trait, .before = 1)

  robust <- purrr::imap(models, robust_auc_model)
  coefficients <- purrr::imap_dfr(models, function(model, id) {
    dplyr::bind_rows(tidy_auc_model(model, id), robust[[id]]$table) |>
      dplyr::mutate(trait = trait, .before = 1)
  })
  variances <- purrr::imap_dfr(models, auc_variance_components) |>
    dplyr::mutate(trait = trait, .before = 1)
  ledger <- purrr::imap_dfr(models, function(model, id) {
    fallback_level <- if (grepl("FALLBACK", id)) 5L else 0L
    auc_model_ledger_row(
      model, id, fallback_level, v_type, selected = id == selected_id,
      reason = dplyr::case_when(
        id == selected_id ~ "selected by prespecified trait hierarchy",
        !model$stable ~ "failed convergence, rank, covariance, or variance diagnostics",
        id == "H3_FIXED_CONTEXT" ~ "fixed-context comparison",
        id == "H3_HIERARCHICAL" ~ "hierarchical varying-coefficient comparison",
        TRUE ~ "prespecified benchmark/comparison"
      )
    ) |>
      dplyr::mutate(trait = trait, .before = 1)
  })
  trait_robust_summary <- coefficients |>
    dplyr::filter(startsWith(.data$inference, "experiment_clustered_")) |>
    dplyr::summarise(
      minimum_robust_df = min(.data$denominator_df, na.rm = TRUE),
      .by = "model_id"
    )
  ledger <- ledger |>
    dplyr::left_join(trait_robust_summary, by = "model_id") |>
    dplyr::mutate(
      robust_method = vapply(
        .data$model_id,
        function(id) auc_value_or(robust[[id]]$selected_method, "none"), character(1)
      ),
      robust_note = vapply(
        .data$model_id,
        function(id) auc_value_or(robust[[id]]$error, ""), character(1)
      )
    )

  selected <- if (!is.na(selected_id)) models[[selected_id]] else NULL
  block_tests <- list()
  if (!is.null(models$H2_ADDITIVE)) {
    block_tests$trait_additive <- auc_block_test(
      models$H2_ADDITIVE, robust$H2_ADDITIVE, "trait_z",
      "trait adds information beyond climate"
    )
  }
  if (!is.null(models$H3_FULL)) {
    block_tests$trait_context <- auc_block_test(
      models$H3_FULL, robust$H3_FULL,
      c(
        "trait_z:W_z", "trait_z:T_z:W_z", "trait_z:P_z:W_z",
        "trait_z:T_z:P_z:W_z"
      ),
      "trait modifies warming across joint climate context"
    )
    added_terms <- setdiff(
      names(stats::coef(models$H3_FULL$fit)),
      names(stats::coef(models$H2_ADDITIVE$fit))
    )
    block_tests$all_added_context <- auc_block_test(
      models$H3_FULL, robust$H3_FULL, added_terms,
      "all context-dependent trait terms added in H3"
    )
  }
  block_tests <- dplyr::bind_rows(block_tests) |>
    dplyr::mutate(trait = trait, .before = 1)

  h1_slope <- variances$estimate[
    variances$model_id == "H1_CLIMATE" & variances$component == "species_W_slope"
  ]
  h3_slope <- variances$estimate[
    variances$model_id == "H3_HIERARCHICAL" & variances$component == "species_W_slope"
  ]
  reduction <- tibble::tibble(
    trait = trait,
    H1_species_W_slope_variance = if (length(h1_slope)) h1_slope[[1]] else NA_real_,
    H3_species_W_slope_variance = if (length(h3_slope)) h3_slope[[1]] else NA_real_
  ) |>
    dplyr::mutate(
      percent_reduction_residual_species_W_slope_variance = dplyr::if_else(
        is.finite(.data$H1_species_W_slope_variance) &
          .data$H1_species_W_slope_variance > 0 &
          is.finite(.data$H3_species_W_slope_variance) &
          models$H1_CLIMATE$stable & models$H3_HIERARCHICAL$stable,
        100 * (.data$H1_species_W_slope_variance -
                 .data$H3_species_W_slope_variance) /
          .data$H1_species_W_slope_variance,
        NA_real_
      ),
      interpretation = dplyr::if_else(
        is.finite(.data$percent_reduction_residual_species_W_slope_variance),
        "descriptive percent reduction; not a formal R-squared",
        "not estimable because both random slopes were not identifiable"
      )
    )

  list(
    trait = trait, data = trait_data, V = trait_V, scaling = prepared$scaling,
    random_kind = random_kind, models = models, robust = robust,
    selected_id = selected_id, selected = selected, coefficients = coefficients,
    variance_components = variances, ledger = ledger, comparisons = comparisons,
    block_tests = block_tests, variance_reduction = reduction, ml_models = ml
  )
}

fit_auc_wls_fallback <- function(data, formula) {
  captured <- auc_capture(
    stats::lm(
      stats::update.formula(formula, auc_effect ~ .), data = data,
      weights = 1 / auc_se^2
    )
  )
  if (is.null(captured$fit)) {
    return(list(fit = NULL, covariance = NULL, table = tibble::tibble(),
                warnings = captured$warnings, error = captured$error, stable = FALSE))
  }
  covariance <- tryCatch(
    sandwich::vcovCL(
      captured$fit,
      cluster = data.frame(
        experiment_id = data$experiment_id, species_id = data$species_id
      ),
      type = "HC1", fix = FALSE
    ),
    error = function(condition) NULL
  )
  beta <- stats::coef(captured$fit)
  stable <- !is.null(covariance) && all(is.finite(beta)) && all(is.finite(covariance))
  se <- if (stable) sqrt(diag(covariance)) else rep(NA_real_, length(beta))
  table <- tibble::tibble(
    inference = "inverse_variance_WLS_two_way_clustered",
    term = names(beta), estimate = as.numeric(beta), standard_error = se,
    statistic = as.numeric(beta) / se,
    p_value = 2 * stats::pnorm(abs(.data$statistic), lower.tail = FALSE),
    conf_low = as.numeric(beta) - stats::qnorm(0.975) * se,
    conf_high = as.numeric(beta) + stats::qnorm(0.975) * se
  )
  list(
    fit = captured$fit, covariance = covariance, table = table,
    warnings = captured$warnings, error = captured$error, stable = stable
  )
}

experiment_level_auc_check <- function(data) {
  summaries <- data |>
    dplyr::mutate(weight = 1 / .data$auc_se^2) |>
    dplyr::summarise(
      experiment_mean_auc = stats::weighted.mean(.data$auc_effect, .data$weight),
      experiment_mean_se = sqrt(1 / sum(.data$weight)),
      T_z = dplyr::first(.data$T_z), P_z = dplyr::first(.data$P_z),
      W_z = dplyr::first(.data$W_z), n_species = dplyr::n(),
      .by = "experiment_id"
    )
  fit <- stats::lm(
    experiment_mean_auc ~ T_z * P_z * W_z, data = summaries,
    weights = 1 / experiment_mean_se^2
  )
  covariance <- sandwich::vcovHC(fit, type = "HC3")
  beta <- stats::coef(fit)
  se <- sqrt(diag(covariance))
  table <- tibble::tibble(
    term = names(beta), estimate = as.numeric(beta), standard_error = se,
    statistic = as.numeric(beta) / se,
    denominator_df = stats::df.residual(fit),
    p_value = 2 * stats::pt(abs(.data$statistic), df = .data$denominator_df,
                           lower.tail = FALSE),
    conf_low = as.numeric(beta) - stats::qt(0.975, .data$denominator_df) * se,
    conf_high = as.numeric(beta) + stats::qt(0.975, .data$denominator_df) * se
  )
  list(data = summaries, fit = fit, covariance = covariance, table = table)
}

auc_design_matrix <- function(formula, newdata, coefficient_names) {
  design <- stats::model.matrix(formula, data = newdata)
  if ("intrcpt" %in% coefficient_names && "(Intercept)" %in% colnames(design)) {
    colnames(design)[colnames(design) == "(Intercept)"] <- "intrcpt"
  }
  missing <- setdiff(coefficient_names, colnames(design))
  if (length(missing)) {
    stop("Prediction matrix is missing coefficients: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  design[, coefficient_names, drop = FALSE]
}

predict_auc_fixed <- function(model_result, robust_result, newdata) {
  fit <- model_result$fit
  coefficient_names <- names(stats::coef(fit))
  design <- auc_design_matrix(model_result$formula, newdata, coefficient_names)
  beta <- as.numeric(stats::coef(fit))
  model_covariance <- as.matrix(stats::vcov(fit))
  estimate <- as.numeric(design %*% beta)
  model_se <- sqrt(pmax(0, rowSums((design %*% model_covariance) * design)))
  robust_covariance <- if (is.null(robust_result$model)) NULL else
    as.matrix(robust_result$model$vb)
  robust_se <- if (is.null(robust_covariance)) rep(NA_real_, nrow(design)) else
    sqrt(pmax(0, rowSums((design %*% robust_covariance) * design)))
  dplyr::bind_cols(
    newdata,
    tibble::tibble(
      estimate = estimate, model_se = model_se,
      model_conf_low = estimate - stats::qnorm(0.975) * model_se,
      model_conf_high = estimate + stats::qnorm(0.975) * model_se,
      robust_se = robust_se,
      robust_conf_low = estimate - stats::qnorm(0.975) * robust_se,
      robust_conf_high = estimate + stats::qnorm(0.975) * robust_se
    )
  )
}

predict_auc_derivative <- function(model_result, robust_result, newdata,
                                   variable = "W_z", epsilon = 1e-6) {
  upper <- newdata
  lower <- newdata
  upper[[variable]] <- upper[[variable]] + epsilon
  lower[[variable]] <- lower[[variable]] - epsilon
  coefficient_names <- names(stats::coef(model_result$fit))
  design <- (
    auc_design_matrix(model_result$formula, upper, coefficient_names) -
      auc_design_matrix(model_result$formula, lower, coefficient_names)
  ) / (2 * epsilon)
  beta <- as.numeric(stats::coef(model_result$fit))
  estimate <- as.numeric(design %*% beta)
  model_covariance <- as.matrix(stats::vcov(model_result$fit))
  model_se <- sqrt(pmax(0, rowSums((design %*% model_covariance) * design)))
  robust_covariance <- if (is.null(robust_result$model)) NULL else
    as.matrix(robust_result$model$vb)
  robust_se <- if (is.null(robust_covariance)) rep(NA_real_, nrow(design)) else
    sqrt(pmax(0, rowSums((design %*% robust_covariance) * design)))
  dplyr::bind_cols(
    newdata,
    tibble::tibble(
      derivative = variable, estimate = estimate, model_se = model_se,
      model_conf_low = estimate - stats::qnorm(0.975) * model_se,
      model_conf_high = estimate + stats::qnorm(0.975) * model_se,
      robust_se = robust_se,
      robust_conf_low = estimate - stats::qnorm(0.975) * robust_se,
      robust_conf_high = estimate + stats::qnorm(0.975) * robust_se
    )
  )
}

auc_support_mask <- function(experiments, grid) {
  observed <- as.matrix(experiments[, c("T_z", "P_z", "W_z")])
  candidate <- as.matrix(grid[, c("T_z", "P_z", "W_z")])
  captured <- auc_capture({
    hull <- geometry::convhulln(observed, options = "Qt")
    geometry::inhulln(hull, candidate)
  })
  if (is.null(captured$fit)) {
    return(list(
      grid = experiments |>
        dplyr::distinct(.data$experiment_id, .data$T_z, .data$P_z, .data$W_z) |>
        dplyr::mutate(supported = TRUE, support_method = "observed experiment points"),
      diagnostic = captured$error
    ))
  }
  list(
    grid = grid |>
      dplyr::mutate(
        supported = as.logical(captured$fit),
        support_method = "three-dimensional convex hull"
      ),
    diagnostic = paste(captured$warnings, collapse = " | ")
  )
}

build_auc_climate_grid <- function(data, points = 60L) {
  experiments <- data |>
    dplyr::distinct(.data$experiment_id, .data$T_z, .data$P_z, .data$W_z)
  warming <- stats::quantile(
    experiments$W_z, probs = c(0.25, 0.50, 0.75), names = FALSE
  )
  warming_labels <- c("low (25th percentile)", "median", "high (75th percentile)")
  grid <- tidyr::expand_grid(
    T_z = seq(min(experiments$T_z), max(experiments$T_z), length.out = points),
    P_z = seq(min(experiments$P_z), max(experiments$P_z), length.out = points),
    W_z = warming
  ) |>
    dplyr::mutate(
      warming_level = factor(
        .env$warming_labels[match(.data$W_z, .env$warming)],
        levels = .env$warming_labels
      ),
      duration_z = 0, log_plot_area_z = 0,
      climate_context = factor("cold_dry", levels = levels(data$climate_context))
    )
  auc_support_mask(experiments, grid)
}

auc_prediction_outputs <- function(model_result, robust_result, data, points = 60L) {
  support <- build_auc_climate_grid(data, points)
  grid <- support$grid
  if (!"warming_level" %in% names(grid)) {
    grid <- grid |>
      dplyr::mutate(
        warming_level = "observed",
        duration_z = 0, log_plot_area_z = 0,
        climate_context = factor("cold_dry", levels = levels(data$climate_context))
      )
  }
  supported <- grid |>
    dplyr::filter(.data$supported)
  list(
    predictions = predict_auc_fixed(model_result, robust_result, supported),
    warming_slopes = predict_auc_derivative(
      model_result, robust_result, supported, variable = "W_z"
    ),
    support = grid,
    support_diagnostic = support$diagnostic
  )
}

auc_trait_prediction_outputs <- function(trait_result, points = 45L) {
  model <- trait_result$selected
  if (is.null(model) || is.null(model$fit)) return(NULL)
  support <- build_auc_climate_grid(trait_result$data, points)
  grid <- support$grid
  if (!"warming_level" %in% names(grid)) grid$warming_level <- "observed"
  grid <- grid |>
    dplyr::filter(.data$supported) |>
    tidyr::crossing(
      trait_z = c(-1, 0, 1),
      trait_level = factor(
        c("-1 SD", "mean", "+1 SD"), levels = c("-1 SD", "mean", "+1 SD")
      )
    ) |>
    dplyr::mutate(
      duration_z = dplyr::coalesce(.data$duration_z, 0),
      log_plot_area_z = dplyr::coalesce(.data$log_plot_area_z, 0),
      climate_context = factor(
        as.character(.data$climate_context), levels = levels(trait_result$data$climate_context)
      )
    )
  robust <- trait_result$robust[[trait_result$selected_id]]
  predictions <- predict_auc_fixed(model, robust, grid)
  slopes <- predict_auc_derivative(model, robust, grid, "W_z")
  reference <- grid |>
    dplyr::filter(.data$trait_z == 0) |>
    dplyr::select(-.data$trait_z, -.data$trait_level)
  plus_one <- reference |>
    dplyr::mutate(trait_z = 1)
  coefficient_names <- names(stats::coef(model$fit))
  contrast_design <-
    auc_design_matrix(model$formula, plus_one, coefficient_names) -
    auc_design_matrix(
      model$formula, dplyr::mutate(reference, trait_z = 0), coefficient_names
    )
  beta <- as.numeric(stats::coef(model$fit))
  covariance <- if (is.null(robust$model)) as.matrix(stats::vcov(model$fit)) else
    as.matrix(robust$model$vb)
  trait_effect <- as.numeric(contrast_design %*% beta)
  trait_se <- sqrt(pmax(0, rowSums((contrast_design %*% covariance) * contrast_design)))
  effects <- dplyr::bind_cols(
    reference,
    tibble::tibble(
      estimate = trait_effect, robust_se = trait_se,
      robust_conf_low = trait_effect - stats::qnorm(0.975) * trait_se,
      robust_conf_high = trait_effect + stats::qnorm(0.975) * trait_se
    )
  )
  list(
    predictions = predictions, warming_slopes = slopes, trait_effect = effects,
    support = support$grid, support_diagnostic = support$diagnostic
  )
}

manual_species_random_effects <- function(model_result) {
  fit <- model_result$fit
  if (is.null(fit) || model_result$random_kind != "species_slope") {
    return(tibble::tibble())
  }
  W_inverse <- as.matrix(Matrix::solve(fit$M))
  X <- as.matrix(fit$X)
  fixed_covariance <- solve(crossprod(X, W_inverse %*% X))
  hat <- X %*% fixed_covariance %*% crossprod(X, W_inverse)
  identity <- diag(nrow(X))
  residual <- as.numeric(fit$yi - X %*% fit$beta)

  calculate <- function(Z, variance, species, component) {
    Z <- as.matrix(Z)
    D <- diag(variance, ncol(Z))
    DZtW <- D %*% crossprod(Z, W_inverse)
    prediction <- as.numeric(DZtW %*% residual)
    prediction_covariance <- D - DZtW %*% (identity - hat) %*% Z %*% D
    se <- sqrt(pmax(0, diag(prediction_covariance)))
    tibble::tibble(
      species_id = species, component = component, estimate = prediction,
      standard_error = se,
      conf_low = prediction - stats::qnorm(0.975) * se,
      conf_high = prediction + stats::qnorm(0.975) * se
    )
  }
  intercept <- calculate(
    fit$Z.S[[1]], fit$sigma2[[1]], fit$s.levels[[1]], "species_intercept"
  )
  slope_Z <- Matrix::Diagonal(x = as.numeric(fit$Z.G1[, 1])) %*% fit$Z.G2
  slope <- calculate(
    slope_Z, fit$tau2[[1]], levels(fit$mf.g$outer), "species_W_slope"
  )
  dplyr::bind_rows(intercept, slope)
}

species_auc_outputs <- function(model_result, data, display_min = 3L) {
  effects <- manual_species_random_effects(model_result)
  support <- data |>
    dplyr::summarise(
      n_experiments = dplyr::n_distinct(.data$experiment_id),
      T_min = min(.data$T), T_max = max(.data$T),
      P_min = min(.data$P_raw), P_max = max(.data$P_raw),
      W_min = min(.data$W), W_max = max(.data$W),
      .by = "species_id"
    ) |>
    dplyr::mutate(display_eligible = .data$n_experiments >= display_min)
  if (!nrow(effects)) {
    return(list(random_effects = effects, support = support,
                conditional_predictions = tibble::tibble()))
  }
  wide <- effects |>
    dplyr::select(.data$species_id, .data$component, .data$estimate) |>
    tidyr::pivot_wider(names_from = .data$component, values_from = .data$estimate)
  fixed <- predict_auc_fixed(
    model_result, list(model = NULL), data |>
      dplyr::select(dplyr::everything())
  )
  conditional <- fixed |>
    dplyr::left_join(wide, by = "species_id") |>
    dplyr::mutate(
      conditional_estimate = .data$estimate + .data$species_intercept +
        .data$species_W_slope * .data$W_z
    ) |>
    dplyr::select(
      .data$effect_id, .data$species_id, .data$experiment_id,
      .data$T, .data$P_raw, .data$W, .data$W_z,
      fixed_estimate = .data$estimate, .data$conditional_estimate
    ) |>
    dplyr::left_join(support, by = "species_id") |>
    dplyr::mutate(
      descriptive_response = dplyr::case_when(
        !.data$display_eligible ~ "below display threshold",
        .data$conditional_estimate > 0 ~ "positive at this observed condition",
        .data$conditional_estimate < 0 ~ "negative at this observed condition",
        TRUE ~ "zero at this observed condition"
      )
    )
  classification <- conditional |>
    dplyr::filter(.data$display_eligible) |>
    dplyr::summarise(
      supported_response = dplyr::case_when(
        all(.data$conditional_estimate > 0) ~ "consistent positive / winner-like",
        all(.data$conditional_estimate < 0) ~ "consistent negative / loser-like",
        TRUE ~ "context-dependent / switcher"
      ),
      .by = "species_id"
    )
  support <- support |>
    dplyr::left_join(classification, by = "species_id")
  list(random_effects = effects, support = support,
       conditional_predictions = conditional)
}

auc_residual_diagnostics <- function(model_result) {
  fit <- model_result$fit
  if (is.null(fit)) return(list(rows = tibble::tibble(), experiments = tibble::tibble()))
  standardized <- tryCatch(metafor::rstandard.rma.mv(fit)$z, error = function(condition) {
    rep(NA_real_, nrow(model_result$data))
  })
  leverage <- tryCatch(stats::hatvalues(fit), error = function(condition) {
    rep(NA_real_, nrow(model_result$data))
  })
  cooks <- (as.numeric(standardized)^2 * as.numeric(leverage)) /
    (ncol(fit$X) * pmax(1 - as.numeric(leverage), 1e-8))
  rows <- model_result$data |>
    dplyr::transmute(
      .data$effect_id, .data$species_id, .data$experiment_id,
      residual = as.numeric(stats::residuals(fit)),
      standardized_residual = as.numeric(standardized),
      leverage = as.numeric(leverage),
      cooks_distance_approximation = as.numeric(cooks)
    )
  experiments <- rows |>
    dplyr::summarise(
      n_effects = dplyr::n(), maximum_leverage = max(.data$leverage, na.rm = TRUE),
      maximum_cooks_distance = max(.data$cooks_distance_approximation, na.rm = TRUE),
      sum_cooks_distance = sum(.data$cooks_distance_approximation, na.rm = TRUE),
      .by = "experiment_id"
    ) |>
    dplyr::arrange(dplyr::desc(.data$sum_cooks_distance))
  list(rows = rows, experiments = experiments)
}

leave_one_experiment_out_auc <- function(model_result, prediction_grid = NULL) {
  if (is.null(model_result) || is.null(model_result$fit)) {
    return(list(
      diagnostics = tibble::tibble(), coefficients = tibble::tibble(),
      prediction_summary = tibble::tibble()
    ))
  }
  data <- model_result$data
  V <- model_result$V
  experiments <- levels(droplevels(data$experiment_id))
  diagnostics <- list()
  coefficients <- list()
  predictions <- list()
  if (!is.null(prediction_grid)) prediction_grid$grid_id <- seq_len(nrow(prediction_grid))

  for (experiment_index in seq_along(experiments)) {
    experiment <- experiments[[experiment_index]]
    if (experiment_index == 1L || experiment_index %% 10L == 0L ||
        experiment_index == length(experiments)) {
      message(
        "  Leave-one-experiment-out refit ", experiment_index, "/",
        length(experiments)
      )
    }
    keep <- as.character(data$experiment_id) != experiment
    subset_data <- droplevels(data[keep, , drop = FALSE])
    subset_V <- subset_auc_v(V, subset_data$effect_id)
    refit <- fit_auc_rma(
      subset_data, subset_V, model_result$formula,
      model_result$random_kind, model_result$method
    )
    diagnostics[[experiment]] <- tibble::tibble(
      omitted_experiment = experiment, stable = refit$stable,
      n_observations = nrow(subset_data),
      warning = paste(refit$warnings, collapse = " | "), error = refit$error
    )
    if (!is.null(refit$fit)) {
      coefficients[[experiment]] <- tibble::tibble(
        omitted_experiment = experiment,
        term = names(stats::coef(refit$fit)),
        estimate = as.numeric(stats::coef(refit$fit))
      )
      if (!is.null(prediction_grid)) {
        predicted <- predict_auc_fixed(
          refit, list(model = NULL), prediction_grid
        )
        predictions[[experiment]] <- predicted |>
          dplyr::transmute(
            .data$grid_id, omitted_experiment = experiment,
            loo_estimate = .data$estimate
          )
      }
    }
  }
  prediction_rows <- dplyr::bind_rows(predictions)
  prediction_summary <- if (nrow(prediction_rows)) {
    prediction_rows |>
      dplyr::summarise(
        n_successful_refits = dplyr::n(),
        loo_mean = mean(.data$loo_estimate), loo_sd = stats::sd(.data$loo_estimate),
        loo_min = min(.data$loo_estimate), loo_max = max(.data$loo_estimate),
        .by = "grid_id"
      )
  } else {
    tibble::tibble()
  }
  if (!is.null(prediction_grid) && nrow(prediction_summary)) {
    prediction_summary <- prediction_grid |>
      dplyr::left_join(prediction_summary, by = "grid_id")
  }
  list(
    diagnostics = dplyr::bind_rows(diagnostics),
    coefficients = dplyr::bind_rows(coefficients),
    prediction_summary = prediction_summary
  )
}

validate_auc_pre_model_data <- function(data, horizons = c("max", "3")) {
  required <- c(
    "horizon", "effect_id", "experiment_id", "species_id", "gradient_id",
    "auc_effect", "auc_se", "T", "P_raw", "W", "duration",
    "integration_span", "plot_area", "LDMC", "LMA", "vegetative_height"
  )
  missing_columns <- setdiff(required, names(data))
  if (length(missing_columns)) {
    stop(
      "Pre-model resident AUC audit is missing required columns: ",
      paste(missing_columns, collapse = ", "), call. = FALSE
    )
  }
  requested <- unique(as.character(horizons))
  observed <- unique(as.character(data$horizon))
  missing_horizons <- setdiff(requested, observed)
  if (length(missing_horizons)) {
    stop(
      "Pre-model resident AUC audit is missing horizons: ",
      paste(missing_horizons, collapse = ", "), call. = FALSE
    )
  }
  focal <- data[as.character(data$horizon) %in% requested, , drop = FALSE]
  character_fields <- c("effect_id", "experiment_id", "species_id", "gradient_id")
  for (field in character_fields) {
    values <- as.character(focal[[field]])
    if (anyNA(values) || any(!nzchar(trimws(values)))) {
      stop("Pre-model resident AUC audit found missing ", field, ".", call. = FALSE)
    }
  }
  if (anyDuplicated(focal$effect_id)) {
    stop("Pre-model resident AUC audit found duplicate effect_id values.", call. = FALSE)
  }
  keys <- paste(focal$horizon, focal$experiment_id, focal$species_id, sep = "::")
  if (anyDuplicated(keys)) {
    stop(
      "Pre-model resident AUC audit found duplicate horizon-species-experiment keys.",
      call. = FALSE
    )
  }
  if (any(!is.finite(focal$auc_effect))) {
    stop("Pre-model resident AUC audit found non-finite AUC effects.", call. = FALSE)
  }
  if (any(!is.finite(focal$auc_se) | focal$auc_se <= 0)) {
    stop("Pre-model resident AUC audit found non-positive or non-finite SEs.", call. = FALSE)
  }
  core_numeric <- c("T", "P_raw", "W", "duration", "integration_span", "plot_area")
  for (field in core_numeric) {
    if (any(!is.finite(focal[[field]]))) {
      stop("Pre-model resident AUC audit found non-finite ", field, ".", call. = FALSE)
    }
  }
  if (any(focal$P_raw <= 0) || any(focal$plot_area <= 0) ||
      any(focal$duration <= 0) || any(focal$integration_span < 0)) {
    stop(
      "Pre-model resident AUC audit requires positive precipitation, plot area, and duration, and a non-negative integration span.",
      call. = FALSE
    )
  }
  experiment_consistency <- focal |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(c("gradient_id", "T", "P_raw", "W", "duration", "integration_span", "plot_area")),
        dplyr::n_distinct
      ),
      .by = c("horizon", "experiment_id")
    )
  inconsistent <- experiment_consistency |>
    dplyr::filter(dplyr::if_any(-dplyr::all_of(c("horizon", "experiment_id")), ~ .x != 1L))
  if (nrow(inconsistent)) {
    stop(
      "Pre-model resident AUC audit found inconsistent experiment-level covariates.",
      call. = FALSE
    )
  }
  tibble::tibble(
    check = c(
      "required_columns", "requested_horizons", "complete_identifiers",
      "unique_effect_id", "unique_horizon_species_experiment",
      "finite_auc_effect", "positive_finite_auc_se",
      "finite_core_covariates", "consistent_experiment_covariates"
    ),
    passed = TRUE,
    detail = c(
      paste(length(required), "required columns present"),
      paste(requested, collapse = ", "),
      paste(nrow(focal), "rows checked"),
      paste(dplyr::n_distinct(focal$effect_id), "unique effect IDs"),
      paste(dplyr::n_distinct(keys), "unique keys"),
      paste(nrow(focal), "finite effects"),
      paste(nrow(focal), "positive finite SEs"),
      "T, precipitation, warming, duration, integration span, and plot area valid",
      paste(nrow(experiment_consistency), "horizon-experiment groups checked")
    )
  )
}

auc_plot_replication <- function(auc_by_plot, horizon) {
  horizon_value <- as.character(horizon)
  auc_by_plot |>
    dplyr::filter(
      as.character(.data$horizon) == .env$horizon_value,
      .data$transform == "raw", .data$ODT %in% c("origin_controls", "warmed")
    ) |>
    dplyr::summarise(
      origin_plots = dplyr::n_distinct(
        .data$stable_plot_id[.data$ODT == "origin_controls"]
      ),
      warmed_plots = dplyr::n_distinct(
        .data$stable_plot_id[.data$ODT == "warmed"]
      ),
      .by = "experiment_id"
    ) |>
    dplyr::mutate(experiment_id = as.character(.data$experiment_id))
}

auc_experiment_audit <- function(data, horizon, auc_by_plot = NULL) {
  horizon_value <- as.character(horizon)
  focal <- data |>
    dplyr::filter(as.character(.data$horizon) == .env$horizon_value)
  result <- focal |>
    dplyr::summarise(
      gradient_id = as.character(dplyr::first(.data$gradient_id)),
      T = dplyr::first(.data$T), P_raw = dplyr::first(.data$P_raw),
      W = dplyr::first(.data$W), duration = dplyr::first(.data$duration),
      integration_span = dplyr::first(.data$integration_span),
      plot_area = dplyr::first(.data$plot_area),
      n_species = dplyr::n_distinct(.data$species_id),
      .by = "experiment_id"
    ) |>
    dplyr::mutate(horizon = .env$horizon_value, .before = 1)
  if (!is.null(auc_by_plot)) {
    result <- result |>
      dplyr::left_join(auc_plot_replication(auc_by_plot, horizon_value), by = "experiment_id")
  } else {
    result <- result |>
      dplyr::mutate(origin_plots = NA_integer_, warmed_plots = NA_integer_)
  }
  result |>
    dplyr::arrange(.data$experiment_id)
}

auc_trait_coverage <- function(data, horizons = unique(data$horizon)) {
  purrr::map_dfr(as.character(horizons), function(horizon_value) {
    focal <- data |>
      dplyr::filter(as.character(.data$horizon) == .env$horizon_value)
    purrr::map_dfr(c("LDMC", "LMA", "vegetative_height"), function(trait) {
      available <- is.finite(focal[[trait]]) & focal[[trait]] > 0
      available_species <- dplyr::n_distinct(focal$species_id[available])
      total_species <- dplyr::n_distinct(focal$species_id)
      tibble::tibble(
        horizon = horizon_value, trait = trait,
        available_effects = sum(available), total_effects = nrow(focal),
        effect_coverage = mean(available),
        available_species = available_species,
        total_species = total_species,
        species_coverage = available_species / total_species
      )
    })
  })
}

auc_data_audit <- function(data, horizon, auc_by_plot = NULL) {
  horizon_value <- as.character(horizon)
  focal <- data |>
    dplyr::filter(as.character(.data$horizon) == .env$horizon_value)
  replication <- focal |>
    dplyr::summarise(
      n_experiments = dplyr::n_distinct(.data$experiment_id),
      .by = "species_id"
    )
  numeric_summary <- function(values) {
    c(
      minimum = min(values, na.rm = TRUE), median = stats::median(values, na.rm = TRUE),
      mean = mean(values, na.rm = TRUE), maximum = max(values, na.rm = TRUE),
      missing = sum(is.na(values))
    )
  }
  experiment_audit <- auc_experiment_audit(focal, horizon_value, auc_by_plot)
  keys <- paste(focal$experiment_id, focal$species_id, sep = "::")
  trait_coverage <- auc_trait_coverage(focal, horizon_value)
  metrics <- c(
    total_effect_rows = nrow(focal),
    species_experiment_combinations = dplyr::n_distinct(keys),
    unique_species = dplyr::n_distinct(focal$species_id),
    experiments = dplyr::n_distinct(focal$experiment_id),
    gradients_studies = dplyr::n_distinct(focal$gradient_id),
    negative_auc_effects = sum(focal$auc_effect < 0),
    zero_auc_effects = sum(focal$auc_effect == 0),
    positive_auc_effects = sum(focal$auc_effect > 0),
    positive_finite_SE = sum(is.finite(focal$auc_se) & focal$auc_se > 0),
    species_one_experiment = sum(replication$n_experiments == 1),
    species_two_or_more_experiments = sum(replication$n_experiments >= 2),
    species_three_or_more_experiments = sum(replication$n_experiments >= 3),
    incomplete_climate_rows = sum(!focal$climate_complete),
    LDMC_available = trait_coverage$available_effects[trait_coverage$trait == "LDMC"],
    LMA_available = trait_coverage$available_effects[trait_coverage$trait == "LMA"],
    vegetative_height_available = trait_coverage$available_effects[
      trait_coverage$trait == "vegetative_height"
    ],
    LDMC_species_available = trait_coverage$available_species[trait_coverage$trait == "LDMC"],
    LMA_species_available = trait_coverage$available_species[trait_coverage$trait == "LMA"],
    vegetative_height_species_available = trait_coverage$available_species[
      trait_coverage$trait == "vegetative_height"
    ]
  )
  ranges <- c()
  for (variable in c("auc_effect", "auc_se")) {
    summary <- numeric_summary(focal[[variable]])
    ranges <- c(ranges, stats::setNames(summary, paste(variable, names(summary), sep = "_")))
  }
  for (variable in c("T", "P_raw", "W", "duration", "integration_span", "plot_area")) {
    summary <- numeric_summary(experiment_audit[[variable]])
    ranges <- c(ranges, stats::setNames(summary, paste(variable, names(summary), sep = "_")))
  }
  plot_metrics <- numeric()
  if (!is.null(auc_by_plot)) {
    plot_metrics <- c(
      origin_plots_minimum = min(experiment_audit$origin_plots),
      origin_plots_median = stats::median(experiment_audit$origin_plots),
      origin_plots_maximum = max(experiment_audit$origin_plots),
      warmed_plots_minimum = min(experiment_audit$warmed_plots),
      warmed_plots_median = stats::median(experiment_audit$warmed_plots),
      warmed_plots_maximum = max(experiment_audit$warmed_plots)
    )
  }
  tibble::tibble(
    horizon = horizon_value, metric = names(c(metrics, ranges, plot_metrics)),
    value = as.numeric(c(metrics, ranges, plot_metrics))
  )
}

auc_missingness_audit <- function(data, horizons = unique(data$horizon)) {
  purrr::map_dfr(as.character(horizons), function(horizon_value) {
    focal <- data |>
      dplyr::filter(as.character(.data$horizon) == .env$horizon_value)
    focal |>
      dplyr::summarise(dplyr::across(dplyr::everything(), ~ sum(is.na(.x)))) |>
      tidyr::pivot_longer(
        dplyr::everything(), names_to = "column", values_to = "n_missing"
      ) |>
      dplyr::mutate(
        horizon = .env$horizon_value, .before = 1,
        n_rows = nrow(focal), proportion_missing = .data$n_missing / .data$n_rows
      )
  })
}

auc_species_replication <- function(data) {
  data |>
    dplyr::summarise(
      n_experiments = dplyr::n_distinct(.data$experiment_id),
      T_min = min(.data$T), T_max = max(.data$T),
      P_min = min(.data$P_raw), P_max = max(.data$P_raw),
      W_min = min(.data$W), W_max = max(.data$W),
      .by = c("horizon", "species_id")
    )
}

auc_context_support <- function(data) {
  data |>
    dplyr::distinct(
      .data$horizon, .data$experiment_id, .data$climate_context
    ) |>
    dplyr::count(.data$horizon, .data$climate_context, name = "n_experiments") |>
    tidyr::complete(
      horizon, climate_context,
      fill = list(n_experiments = 0L)
    )
}

auc_horizon_comparison <- function(data) {
  wide <- data |>
    dplyr::select(
      .data$experiment_id, .data$species_id, .data$horizon,
      .data$auc_effect, .data$auc_se, .data$duration
    ) |>
    tidyr::pivot_wider(
      names_from = .data$horizon,
      values_from = c(.data$auc_effect, .data$auc_se),
      names_sep = "_"
    )
  if (!all(c("auc_effect_max", "auc_effect_3") %in% names(wide))) {
    return(list(paired = wide, summary = tibble::tibble()))
  }
  wide <- wide |>
    dplyr::mutate(
      difference_max_minus_3 = .data$auc_effect_max - .data$auc_effect_3,
      sign_changed = sign(.data$auc_effect_max) != sign(.data$auc_effect_3),
      max_auc_per_inclusive_year = .data$auc_effect_max / .data$duration,
      max_auc_per_integration_span = .data$auc_effect_max / pmax(.data$duration - 1, 1)
    )
  summary <- tibble::tibble(
    n_paired = nrow(wide),
    correlation = stats::cor(wide$auc_effect_max, wide$auc_effect_3),
    mean_absolute_difference = mean(abs(wide$difference_max_minus_3)),
    n_sign_changed = sum(wide$sign_changed),
    correlation_max_with_duration = stats::cor(wide$auc_effect_max, wide$duration),
    correlation_abs_max_with_duration = stats::cor(abs(wide$auc_effect_max), wide$duration)
  )
  list(paired = wide, summary = summary)
}

save_auc_plot <- function(plot, path, width = 10, height = 7) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot = plot, width = width, height = height,
                  device = grDevices::cairo_pdf)
  path
}

plot_auc_climate_surfaces <- function(predictions, path) {
  plot <- ggplot2::ggplot(
    predictions, ggplot2::aes(x = .data$T_z, y = .data$P_z, fill = .data$estimate)
  ) +
    ggplot2::geom_raster() +
    ggplot2::geom_contour(
      ggplot2::aes(z = .data$estimate), breaks = 0,
      colour = "black", linewidth = 0.35
    ) +
    ggplot2::facet_wrap(~ warming_level) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0
    ) +
    ggplot2::coord_equal() + ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Destination JJA temperature (SD)",
      y = "Log destination JJA precipitation (SD)",
      fill = "Predicted\nAUC effect",
      title = "Resident warmed-minus-origin AUC within observed climate support"
    )
  save_auc_plot(plot, path, 11, 4.5)
}

plot_auc_warming_slopes <- function(slopes, path) {
  plot <- ggplot2::ggplot(
    slopes, ggplot2::aes(x = .data$T_z, y = .data$P_z, fill = .data$estimate)
  ) +
    ggplot2::geom_raster() +
    ggplot2::geom_contour(
      ggplot2::aes(z = .data$estimate), breaks = 0,
      colour = "black", linewidth = 0.35
    ) +
    ggplot2::facet_wrap(~ warming_level) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0
    ) +
    ggplot2::coord_equal() + ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Destination JJA temperature (SD)",
      y = "Log destination JJA precipitation (SD)",
      fill = "AUC change per\n1 SD warming",
      title = "Climate-dependent experimental-warming slope"
    )
  save_auc_plot(plot, path, 11, 4.5)
}

plot_auc_uncertainty <- function(predictions, path) {
  plot <- ggplot2::ggplot(
    predictions, ggplot2::aes(x = .data$T_z, y = .data$P_z, fill = .data$robust_se)
  ) +
    ggplot2::geom_raster() + ggplot2::facet_wrap(~ warming_level) +
    ggplot2::scale_fill_viridis_c() + ggplot2::coord_equal() + ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Destination JJA temperature (SD)",
      y = "Log destination JJA precipitation (SD)",
      fill = "Cluster-robust SE",
      title = "Uncertainty of the resident climate surface"
    )
  save_auc_plot(plot, path, 11, 4.5)
}

plot_auc_trait_surfaces <- function(predictions, trait, path) {
  plot <- ggplot2::ggplot(
    predictions, ggplot2::aes(x = .data$T_z, y = .data$P_z, fill = .data$estimate)
  ) +
    ggplot2::geom_raster() +
    ggplot2::geom_contour(
      ggplot2::aes(z = .data$estimate), breaks = 0,
      colour = "black", linewidth = 0.25
    ) +
    ggplot2::facet_grid(trait_level ~ warming_level) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0
    ) +
    ggplot2::coord_equal() + ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Destination JJA temperature (SD)",
      y = "Log destination JJA precipitation (SD)",
      fill = "Predicted\nAUC effect", title = paste(trait, "context-dependent response")
    )
  save_auc_plot(plot, path, 12, 8)
}

plot_auc_trait_effect <- function(effects, trait, path) {
  plot <- ggplot2::ggplot(
    effects, ggplot2::aes(x = .data$T_z, y = .data$P_z, fill = .data$estimate)
  ) +
    ggplot2::geom_raster() +
    ggplot2::geom_contour(
      ggplot2::aes(z = .data$estimate), breaks = 0,
      colour = "black", linewidth = 0.35
    ) +
    ggplot2::facet_wrap(~ warming_level) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0
    ) +
    ggplot2::coord_equal() + ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Destination JJA temperature (SD)",
      y = "Log destination JJA precipitation (SD)",
      fill = "Effect of +1 SD",
      title = paste("Effect of a +1 SD", trait, "change on resident AUC")
    )
  save_auc_plot(plot, path, 11, 4.5)
}

plot_auc_species_effects <- function(random_effects, support, path) {
  display <- random_effects |>
    dplyr::left_join(support, by = "species_id") |>
    dplyr::filter(.data$display_eligible) |>
    dplyr::group_by(.data$component) |>
    dplyr::slice_max(abs(.data$estimate), n = 30, with_ties = FALSE) |>
    dplyr::ungroup()
  plot <- ggplot2::ggplot(
    display,
    ggplot2::aes(
      x = stats::reorder(.data$species_id, .data$estimate),
      y = .data$estimate, ymin = .data$conf_low, ymax = .data$conf_high
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey45") +
    ggplot2::geom_pointrange() + ggplot2::coord_flip() +
    ggplot2::facet_wrap(~ component, scales = "free_y") + ggplot2::theme_bw() +
    ggplot2::labs(
      x = NULL, y = "Species deviation (95% interval)",
      title = "Resident species random intercepts and warming slopes"
    )
  save_auc_plot(plot, path, 11, 10)
}

plot_auc_horizon_comparison <- function(paired, path) {
  plot <- ggplot2::ggplot(
    paired, ggplot2::aes(x = .data$auc_effect_3, y = .data$auc_effect_max)
  ) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey70") +
    ggplot2::geom_vline(xintercept = 0, colour = "grey70") +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
    ggplot2::geom_point(alpha = 0.35, size = 1) + ggplot2::coord_equal() +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "First-3-year AUC effect", y = "Maximum-duration AUC effect",
      title = "Resident AUC horizon comparison"
    )
  save_auc_plot(plot, path, 7, 7)
}

plot_auc_duration <- function(data, path) {
  plot <- ggplot2::ggplot(
    data, ggplot2::aes(x = .data$duration, y = .data$auc_effect)
  ) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey70") +
    ggplot2::geom_jitter(width = 0.08, alpha = 0.25, size = 1) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, colour = "black") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Inclusive experimental duration (years)",
      y = "Maximum-duration warmed-minus-origin AUC",
      title = "Mechanical duration sensitivity of raw maximum AUC"
    )
  save_auc_plot(plot, path, 7, 6)
}

write_auc_markdown_report <- function(path, primary, comparison, v_audit,
                                      horizon_comparison, traits, scaling) {
  selected <- primary$sequence$selected_id
  selected_model <- primary$sequence$selected
  selected_formula <- if (is.null(selected_model)) auc_formula_text(auc_core_formula()) else
    auc_formula_text(selected_model$formula)
  selected_random <- if (is.null(selected_model))
    "No stable rma.mv model; inspect the fallback ledger." else selected_model$random_kind
  selected_robust <- primary$sequence$robust[[selected]]
  robust_method <- if (is.null(selected_robust)) "none" else selected_robust$selected_method
  robust_note <- if (identical(robust_method, "experiment_clustered_CR2")) {
    "Experiment-clustered CR2 with Satterthwaite degrees of freedom."
  } else {
    paste0(
      "Experiment-clustered CR1 was used because CR2 was not supported for the ",
      "crossed species/experiment random structure. CR2 attempt: ",
      if (is.null(selected_robust)) "not available" else selected_robust$error
    )
  }
  slope_supported <- !is.null(selected_model) &&
    identical(selected_model$random_kind, "species_slope") && selected_model$stable
  selected_coefficients <- primary$sequence$coefficients |>
    dplyr::filter(
      .data$model_id == .env$selected,
      startsWith(.data$inference, "experiment_clustered_")
    )
  central <- selected_coefficients |>
    dplyr::filter(.data$term == "T_z:P_z:W_z")
  climate_block <- primary$sequence$block_tests |>
    dplyr::filter(.data$block == "background climate modifies warming response")
  slope_variance <- primary$sequence$variance_components |>
    dplyr::filter(
      .data$model_id == .env$selected, .data$component == "species_W_slope"
    ) |>
    dplyr::pull(.data$estimate)
  nuisance <- primary$sequence$ledger |>
    dplyr::filter(.data$model_id %in% c("S2_EXP_RE", "S3_GRADIENT_RE"))
  climate_claim <- if (nrow(climate_block) && is.finite(climate_block$p_value[[1]]) &&
                       climate_block$p_value[[1]] < 0.05) {
    "The experiment-clustered joint test supports climate-dependent modification of the warming response."
  } else {
    "The experiment-clustered joint test did not provide clear evidence that background climate modifies the warming response."
  }
  trait_lines <- vapply(traits, function(result) {
    robust_block <- result$block_tests |>
      dplyr::filter(.data$block == "trait modifies warming across joint climate context")
    comparison <- result$comparisons |>
      dplyr::filter(.data$comparison == "H2 vs H3")
    reduction <- result$variance_reduction$percent_reduction_residual_species_W_slope_variance[[1]]
    paste0(
      "- ", result$trait, ": selected `", result$selected_id,
      "`; random structure `", result$random_kind, "`; H2-vs-H3 ML LRT p = ",
      formatC(comparison$p_value[[1]], digits = 3, format = "g"),
      "; cluster-robust trait–warming-context block p = ",
      formatC(robust_block$p_value[[1]], digits = 3, format = "g"),
      "; residual species-slope variance change = ",
      if (is.finite(reduction)) paste0(formatC(reduction, digits = 3, format = "f"), "%") else "not estimable",
      "."
    )
  }, character(1))
  lines <- c(
    "# Resident continuous-AUC hierarchical analysis",
    "",
    "## Analysis definition",
    "",
    "- Response: continuous raw relative-cover AUC in warmed plots minus origin-control AUC.",
    "- Primary horizon: maximum available experimental duration (`max`).",
    "- Comparison horizon: common first three years (`3`).",
    "- Population: resident species only. Winner/loser significance categories did not define or filter the response.",
    "- Colonizer/nonresident inference is deferred by design.",
    paste0("- Sampling uncertainty: `", v_audit$selected_v_type[[1]], "`."),
    "- Precipitation and all focal traits were natural-log transformed before scaling.",
    "- Climate scaling used distinct experiments; trait scaling used distinct species.",
    "",
    "## Selected climate model",
    "",
    paste0("- Model ID: `", selected, "`."),
    paste0("- Fixed formula: `", selected_formula, "`."),
    paste0("- Random structure: `", selected_random, "`."),
    paste0("- Dependence correction: ", robust_note),
    paste0("- Effective climate replication: ",
           dplyr::n_distinct(primary$data$experiment_id), " experiments, not ",
           nrow(primary$data), " species–experiment rows."),
    "- `T_z:P_z:W_z` and the joint warming-modification block express whether warming response depends on the full temperature–precipitation background.",
    paste0(
      "- Central T×P×W estimate = ", formatC(central$estimate[[1]], digits = 4, format = "g"),
      " (cluster-robust p = ", formatC(central$p_value[[1]], digits = 3, format = "g"),
      "; df = ", formatC(central$denominator_df[[1]], digits = 3, format = "g"), ")."
    ),
    paste0(
      "- Joint {T×W, P×W, T×P×W} block: F(", climate_block$numerator_df[[1]],
      ", ", formatC(climate_block$denominator_df[[1]], digits = 3, format = "g"),
      ") = ", formatC(climate_block$statistic[[1]], digits = 3, format = "f"),
      ", p = ", formatC(climate_block$p_value[[1]], digits = 3, format = "g"), "."
    ),
    paste0("- Interpretation: ", climate_claim),
    paste0(
      "- Nuisance sensitivities: ",
      paste0(nuisance$model_id, " = ", ifelse(nuisance$stable, "stable", "not retained"),
             collapse = "; "), "."
    ),
    "",
    "## Species hierarchy",
    "",
    if (slope_supported) {
      paste0(
        "The species warming-slope variance was estimable (", formatC(slope_variance[[1]], digits = 4, format = "g"),
        "). Species slopes are residual deviations around the climate-predicted warming slope; species displays are restricted to observed species support."
      )
    } else {
      "The species warming slope was not retained. Species-specific warming sensitivity and Pollock/Jamil/Martinez-style residual slope claims are not supported."
    },
    "",
    "## Trait models",
    "",
    trait_lines,
    "",
    "Trait variance reductions are descriptive and are not labelled R-squared. Full context-dependent trait claims require the H3 hierarchy and adequate cluster-robust degrees of freedom.",
    "The H3 models were retained because they were stable in the prespecified hierarchy, not because individual interaction p-values were selected. Cluster-robust trait–warming-context block tests are the primary inferential checks.",
    "",
    "## Horizon and design sensitivity",
    "",
    paste0("- Paired max/year-3 AUC correlation: ",
           formatC(horizon_comparison$correlation[[1]], digits = 3, format = "f"), "."),
    paste0("- Sign changes between horizons: ", horizon_comparison$n_sign_changed[[1]], "."),
    paste0("- Correlation of absolute maximum AUC with duration: ",
           formatC(horizon_comparison$correlation_abs_max_with_duration[[1]],
                   digits = 3, format = "f"), "."),
    "- S4, year-3, and duration-normalized results are sensitivities; they do not replace the primary maximum-duration response.",
    paste0("- Experiments with partial climate coverage: ",
           dplyr::n_distinct(primary$data$experiment_id[!primary$data$climate_complete]), "."),
    "",
    "## Limitations and interpretation",
    "",
    "Use the model ledger, cluster-robust tables, block tests, variance profiles, and leave-one-experiment-out results together. High-order terms with denominator df below 10 are flagged as limited and below 4 as extremely weak. Prediction surfaces are confidence surfaces for the mean and are masked to observed three-dimensional climate support.",
    "",
    "The existing per-experiment AUC model implies essentially diagonal cross-species sampling covariance. Experiment clustering addresses shared experiment dependence not represented by those near-zero off-diagonal sampling covariances."
  )
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path, useBytes = TRUE)
  path
}

write_auc_table <- function(x, output_dir, subdirectory, filename) {
  write_csv_atomic(x, file.path(output_dir, subdirectory, filename))
}

run_frequentist_hierarchical_auc <- function(
    primary_horizon = "max", comparison_horizon = "3",
    species_display_min = 3L,
    output_dir = file.path("output", "frequentist_hierarchical_auc"),
    root = here::here(), run_leave_one_out = TRUE,
    render_visual_report = TRUE) {
  assert_auc_packages()
  if (!exists("write_auc_pre_model_audit", mode = "function", inherits = TRUE)) {
    reporting_path <- file.path(
      root, "R", "pipeline", "frequentist_hierarchical_auc_reporting.R"
    )
    if (!file.exists(reporting_path)) {
      stop("Missing resident AUC reporting module: ", reporting_path, call. = FALSE)
    }
    source(reporting_path)
  }
  primary_horizon <- as.character(primary_horizon)
  comparison_horizon <- as.character(comparison_horizon)
  if (primary_horizon == comparison_horizon) {
    stop("Primary and comparison horizons must differ.", call. = FALSE)
  }
  if (length(species_display_min) != 1L || !is.finite(species_display_min) ||
      species_display_min < 1) {
    stop("species_display_min must be one positive integer.", call. = FALSE)
  }
  output_dir <- if (grepl("^/", output_dir)) output_dir else file.path(root, output_dir)
  purrr::walk(
    c("tables", "models", "predictions", "diagnostics", "figures"),
    ~ dir.create(file.path(output_dir, .x), recursive = TRUE, showWarnings = FALSE)
  )
  obsolete_tables <- file.path(
    output_dir, "tables",
    c("climate_CR2_block_tests.csv", "trait_CR2_block_tests.csv")
  )
  unlink(obsolete_tables[file.exists(obsolete_tables)])

  horizons <- c(primary_horizon, comparison_horizon)
  raw_data <- read_resident_auc_data(root, horizons)
  scaled <- scale_resident_auc_data(raw_data)
  data <- scaled$data
  auc_path <- file.path(root, "output", "intermediate", "auc_by_plot.rds")
  if (!file.exists(auc_path)) stop("Missing plot-level AUC input: ", auc_path, call. = FALSE)
  auc_by_plot <- readRDS(auc_path)

  pre_model_audit <- write_auc_pre_model_audit(
    data = data, horizons = horizons, auc_by_plot = auc_by_plot,
    scaled = scaled, output_dir = output_dir, root = root, render_html = TRUE
  )
  message("Pre-model audit completed: ", pre_model_audit$report)

  covariance <- list()
  horizon_results <- list()
  predictions <- list()
  profiles <- list()
  diagnostics <- list()
  for (horizon in horizons) {
    message("Preparing horizon ", horizon, " ...")
    horizon_data <- droplevels(data[as.character(data$horizon) == horizon, , drop = FALSE])
    reconstructed <- reconstruct_auc_covariance(horizon_data, auc_by_plot)
    covariance[[horizon]] <- reconstructed
    sequence <- fit_auc_model_sequence(
      horizon_data, reconstructed$V, reconstructed$selected_v_type
    )
    message("Profiling selected horizon ", horizon, " model ...")
    selected_profile <- if (!is.null(sequence$selected) &&
                            sequence$selected$random_kind == "species_slope") {
      profile_auc_slope(sequence$selected)
    } else {
      list(profile = NULL, diagnostic = "selected model has no species warming slope")
    }
    if (!is.null(sequence$selected) && sequence$selected$random_kind == "species_slope" &&
        !auc_profile_is_informative(selected_profile) &&
        isTRUE(sequence$models$S5_SPECIES_INTERCEPT_ONLY$stable)) {
      sequence$ledger <- sequence$ledger |>
        dplyr::mutate(
          selected = .data$model_id == "S5_SPECIES_INTERCEPT_ONLY",
          reason = dplyr::case_when(
            .data$model_id == "S1_PRIMARY" ~
              "species warming-slope profile was flat or lacked a finite profile range",
            .data$model_id == "S5_SPECIES_INTERCEPT_ONLY" ~
              "selected after species warming-slope profile failure",
            TRUE ~ .data$reason
          )
        )
      sequence$selected_id <- "S5_SPECIES_INTERCEPT_ONLY"
      sequence$selected <- sequence$models$S5_SPECIES_INTERCEPT_ONLY
    }
    profiles[[horizon]] <- selected_profile
    message("Generating horizon ", horizon, " predictions and residual diagnostics ...")
    prediction <- if (!is.null(sequence$selected)) {
      auc_prediction_outputs(
        sequence$selected, sequence$robust[[sequence$selected_id]], horizon_data
      )
    } else {
      NULL
    }
    predictions[[horizon]] <- prediction
    residual <- if (!is.null(sequence$selected))
      auc_residual_diagnostics(sequence$selected) else
      list(rows = tibble::tibble(), experiments = tibble::tibble())
    diagnostics[[horizon]] <- residual
    horizon_results[[horizon]] <- list(
      data = horizon_data, covariance = reconstructed,
      sequence = sequence, prediction = prediction,
      profile = selected_profile, residual = residual
    )
  }

  primary <- horizon_results[[primary_horizon]]
  comparison <- horizon_results[[comparison_horizon]]
  if (is.null(primary$sequence$selected)) {
    stop(
      "No stable hierarchical/metafor resident model was available. ",
      "The WLS fallback is saved in the model ledger, but downstream hierarchical ",
      "predictions cannot be generated.", call. = FALSE
    )
  }

  trait_results <- purrr::map(
    c("LDMC", "LMA", "vegetative_height"),
    ~ fit_auc_trait_sequence(
      primary$data, primary$covariance$V, .x,
      primary$covariance$selected_v_type
    )
  ) |>
    stats::setNames(c("LDMC", "LMA", "vegetative_height"))
  message("Profiling trait warming-slope variance components ...")
  trait_profiles <- purrr::imap(trait_results, function(result, trait) {
    list(
      H1_CLIMATE = profile_auc_slope(result$models$H1_CLIMATE, steps = 6L),
      H3_HIERARCHICAL = profile_auc_slope(
        result$models$H3_HIERARCHICAL, steps = 6L
      )
    )
  })
  trait_predictions <- purrr::map(trait_results, auc_trait_prediction_outputs)

  species <- species_auc_outputs(
    primary$sequence$selected, primary$data, as.integer(species_display_min)
  )
  message("Species random-effect outputs completed")
  horizon_comparison <- auc_horizon_comparison(data)

  normalized <- list()
  for (denominator in c("inclusive_years", "integration_span")) {
    message("Fitting duration-normalized sensitivity: ", denominator)
    divisor <- if (denominator == "inclusive_years") {
      primary$data$duration
    } else {
      pmax(primary$data$duration - 1, 1)
    }
    normalized_data <- primary$data
    normalized_data$auc_effect <- normalized_data$auc_effect / divisor
    normalized_data$auc_se <- normalized_data$auc_se / divisor
    inverse_divisor <- Matrix::Diagonal(x = 1 / divisor)
    normalized_V <- inverse_divisor %*% primary$covariance$V %*% inverse_divisor
    dimnames(normalized_V) <- dimnames(primary$covariance$V)
    model <- fit_auc_rma(
      normalized_data, normalized_V, auc_core_formula(),
      primary$sequence$selected$random_kind, "REML"
    )
    normalized[[denominator]] <- list(
      model = model, robust = robust_auc_model(model, paste0("DURATION_", denominator))
    )
  }

  loo <- list(climate = NULL, traits = list())
  if (isTRUE(run_leave_one_out)) {
    climate_grid <- primary$prediction$predictions |>
      dplyr::select(
        dplyr::any_of(c(
          "T_z", "P_z", "W_z", "warming_level", "duration_z",
          "log_plot_area_z", "climate_context"
        ))
      )
    loo$climate <- leave_one_experiment_out_auc(
      primary$sequence$selected, climate_grid
    )
    loo$traits <- purrr::imap(trait_results, function(result, trait) {
      grid <- trait_predictions[[trait]]$predictions |>
        dplyr::filter(.data$trait_z == 0) |>
        dplyr::select(
          dplyr::any_of(c(
            "T_z", "P_z", "W_z", "warming_level", "duration_z",
            "log_plot_area_z", "climate_context", "trait_z"
          ))
        )
      leave_one_experiment_out_auc(result$selected, grid)
    })
  }

  message("Writing resident AUC outputs ...")
  write_auc_table(
    dplyr::bind_rows(purrr::map(covariance, "audit")),
    output_dir, "tables", "auc_uncertainty_audit.csv"
  )
  write_auc_table(
    dplyr::bind_rows(purrr::map(covariance, "block_audit")),
    output_dir, "tables", "auc_uncertainty_blocks.csv"
  )
  write_auc_table(horizon_comparison$paired, output_dir, "tables", "horizon_paired_effects.csv")
  write_auc_table(horizon_comparison$summary, output_dir, "tables", "horizon_comparison.csv")

  model_ledger <- purrr::imap_dfr(horizon_results, function(result, horizon) {
    result$sequence$ledger |>
      dplyr::mutate(horizon = horizon, analysis = "climate", .before = 1)
  })
  trait_ledger <- purrr::imap_dfr(trait_results, function(result, trait) {
    result$ledger |>
      dplyr::mutate(horizon = primary_horizon, analysis = paste0("trait_", trait), .before = 1)
  })
  write_auc_table(
    dplyr::bind_rows(model_ledger, trait_ledger),
    output_dir, "tables", "model_robustness_ledger.csv"
  )
  write_auc_table(
    purrr::imap_dfr(horizon_results, function(result, horizon) {
      result$sequence$coefficients |>
        dplyr::mutate(horizon = horizon, .before = 1)
    }), output_dir, "tables", "climate_coefficients.csv"
  )
  write_auc_table(
    purrr::imap_dfr(horizon_results, function(result, horizon) {
      result$sequence$variance_components |>
        dplyr::mutate(horizon = horizon, .before = 1)
    }), output_dir, "tables", "climate_variance_components.csv"
  )
  write_auc_table(
    purrr::imap_dfr(horizon_results, function(result, horizon) {
      result$sequence$block_tests |>
        dplyr::mutate(horizon = horizon, .before = 1)
    }), output_dir, "tables", "climate_cluster_robust_block_tests.csv"
  )
  write_auc_table(
    purrr::map_dfr(trait_results, "coefficients"),
    output_dir, "tables", "trait_coefficients.csv"
  )
  write_auc_table(
    purrr::map_dfr(trait_results, "variance_components"),
    output_dir, "tables", "trait_variance_components.csv"
  )
  write_auc_table(
    purrr::map_dfr(trait_results, "comparisons"),
    output_dir, "tables", "trait_ML_comparisons.csv"
  )
  write_auc_table(
    purrr::map_dfr(trait_results, "block_tests"),
    output_dir, "tables", "trait_cluster_robust_block_tests.csv"
  )
  write_auc_table(
    purrr::map_dfr(trait_results, "variance_reduction"),
    output_dir, "tables", "trait_species_slope_variance_reduction.csv"
  )
  profile_diagnostics <- dplyr::bind_rows(
    purrr::imap_dfr(profiles, function(result, horizon) {
      tibble::tibble(
        analysis = paste0("climate_horizon_", horizon), model_id = "selected",
        informative = auc_profile_is_informative(result), diagnostic = result$diagnostic
      )
    }),
    purrr::imap_dfr(trait_profiles, function(result, trait) {
      purrr::imap_dfr(result, function(profile, model_id) {
        tibble::tibble(
          analysis = paste0("trait_", trait), model_id = model_id,
          informative = auc_profile_is_informative(profile),
          diagnostic = profile$diagnostic
        )
      })
    })
  )
  write_auc_table(
    profile_diagnostics, output_dir, "diagnostics", "variance_profile_diagnostics.csv"
  )

  write_auc_table(species$random_effects, output_dir, "tables", "species_random_effects.csv")
  write_auc_table(species$support, output_dir, "tables", "species_support.csv")
  write_auc_table(
    species$conditional_predictions, output_dir, "predictions",
    "species_conditional_predictions_at_observed_conditions.csv"
  )
  purrr::iwalk(predictions, function(result, horizon) {
    if (is.null(result)) return(invisible(NULL))
    write_auc_table(result$predictions, output_dir, "predictions",
                    paste0("climate_surface_horizon_", horizon, ".csv"))
    write_auc_table(result$warming_slopes, output_dir, "predictions",
                    paste0("warming_slope_surface_horizon_", horizon, ".csv"))
    write_auc_table(result$support, output_dir, "predictions",
                    paste0("climate_support_mask_horizon_", horizon, ".csv"))
  })
  purrr::iwalk(trait_predictions, function(result, trait) {
    if (is.null(result)) return(invisible(NULL))
    write_auc_table(result$predictions, output_dir, "predictions",
                    paste0("trait_surface_", trait, ".csv"))
    write_auc_table(result$warming_slopes, output_dir, "predictions",
                    paste0("trait_warming_slopes_", trait, ".csv"))
    write_auc_table(result$trait_effect, output_dir, "predictions",
                    paste0("trait_plus_1SD_effect_", trait, ".csv"))
  })

  write_auc_table(primary$residual$rows, output_dir, "diagnostics", "residuals.csv")
  write_auc_table(primary$residual$experiments, output_dir, "diagnostics",
                  "influential_experiments.csv")
  write_auc_table(
    primary$sequence$experiment_level_check$data,
    output_dir, "diagnostics", "experiment_level_weighted_effects.csv"
  )
  write_auc_table(
    primary$sequence$experiment_level_check$table,
    output_dir, "diagnostics", "experiment_level_climate_check.csv"
  )
  if (isTRUE(run_leave_one_out)) {
    write_auc_table(loo$climate$diagnostics, output_dir, "diagnostics", "climate_LOEO_diagnostics.csv")
    write_auc_table(loo$climate$coefficients, output_dir, "diagnostics", "climate_LOEO_coefficients.csv")
    write_auc_table(loo$climate$prediction_summary, output_dir, "diagnostics",
                    "climate_LOEO_prediction_envelope.csv")
    purrr::iwalk(loo$traits, function(result, trait) {
      write_auc_table(result$diagnostics, output_dir, "diagnostics",
                      paste0("trait_LOEO_diagnostics_", trait, ".csv"))
      write_auc_table(result$coefficients, output_dir, "diagnostics",
                      paste0("trait_LOEO_coefficients_", trait, ".csv"))
      write_auc_table(result$prediction_summary, output_dir, "diagnostics",
                      paste0("trait_LOEO_prediction_envelope_", trait, ".csv"))
    })
  }

  normalized_coefficients <- purrr::imap_dfr(normalized, function(result, denominator) {
    dplyr::bind_rows(
      tidy_auc_model(result$model, paste0("DURATION_", denominator)),
      result$robust$table
    ) |>
      dplyr::mutate(denominator = denominator, .before = 1)
  })
  write_auc_table(normalized_coefficients, output_dir, "tables",
                  "duration_normalized_sensitivity.csv")

  write_rds_atomic(
    purrr::map(covariance, ~ .x$V),
    file.path(output_dir, "models", "sampling_covariance_matrices.rds")
  )
  write_rds_atomic(
    purrr::map(horizon_results, ~ purrr::map(.x$sequence$models, "fit")),
    file.path(output_dir, "models", "climate_models.rds")
  )
  write_rds_atomic(
    purrr::map(trait_results, ~ purrr::map(.x$models, "fit")),
    file.path(output_dir, "models", "trait_models.rds")
  )
  write_rds_atomic(
    list(climate = profiles, traits = trait_profiles),
    file.path(output_dir, "models", "variance_profiles.rds")
  )
  write_rds_atomic(
    list(
      climate = purrr::map(horizon_results, ~ purrr::map(.x$sequence$robust, "model")),
      traits = purrr::map(trait_results, ~ purrr::map(.x$robust, "model"))
    ),
    file.path(output_dir, "models", "cluster_robust_models.rds")
  )

  plot_auc_climate_surfaces(
    primary$prediction$predictions,
    file.path(output_dir, "figures", "climate_response_surface.pdf")
  )
  plot_auc_warming_slopes(
    primary$prediction$warming_slopes,
    file.path(output_dir, "figures", "climate_warming_slope_surface.pdf")
  )
  plot_auc_uncertainty(
    primary$prediction$predictions,
    file.path(output_dir, "figures", "climate_response_uncertainty.pdf")
  )
  plot_auc_horizon_comparison(
    horizon_comparison$paired,
    file.path(output_dir, "figures", "max_vs_year3_AUC.pdf")
  )
  plot_auc_duration(
    primary$data, file.path(output_dir, "figures", "maximum_AUC_duration.pdf")
  )
  if (nrow(species$random_effects)) {
    plot_auc_species_effects(
      species$random_effects, species$support,
      file.path(output_dir, "figures", "species_random_effects.pdf")
    )
  }
  purrr::iwalk(trait_predictions, function(result, trait) {
    if (is.null(result)) return(invisible(NULL))
    plot_auc_trait_surfaces(
      result$predictions, trait,
      file.path(output_dir, "figures", paste0("trait_context_", trait, ".pdf"))
    )
    plot_auc_trait_effect(
      result$trait_effect, trait,
      file.path(output_dir, "figures", paste0("trait_plus_1SD_effect_", trait, ".pdf"))
    )
  })

  report_path <- write_auc_markdown_report(
    file.path(output_dir, "report.md"), primary, comparison,
    primary$covariance$audit, horizon_comparison$summary,
    trait_results, scaled$scaling
  )
  visual_report_path <- if (isTRUE(render_visual_report)) {
    message("Rendering the resident AUC visual report ...")
    tryCatch(
      render_frequentist_hierarchical_auc_visual_report(output_dir, root = root),
      error = function(condition) {
        stop(
          "Model outputs were completed, but visual-report rendering failed: ",
          conditionMessage(condition),
          ". Rerun with Rscript R/14_render_frequentist_hierarchical_auc_report.R",
          call. = FALSE
        )
      }
    )
  } else {
    NA_character_
  }
  invisible(list(
    data = data, horizons = horizon_results, traits = trait_results,
    trait_predictions = trait_predictions, species = species,
    normalized = normalized, leave_one_out = loo, audit = pre_model_audit,
    report = report_path, visual_report = visual_report_path, output_dir = output_dir
  ))
}
