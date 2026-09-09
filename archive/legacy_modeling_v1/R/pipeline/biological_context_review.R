biological_context_climate_terms <- function() {
  c(
    "z_destination_T", "z_log_destination_P", "z_warming",
    "z_destination_T:z_log_destination_P",
    "z_destination_T:z_warming",
    "z_log_destination_P:z_warming"
  )
}

biological_context_axes <- function() {
  c(
    temperature = "z_destination_T",
    precipitation = "z_log_destination_P",
    warming = "z_warming"
  )
}

biological_context_base_rhs <- function() {
  paste(biological_context_climate_terms(), collapse = " + ")
}

biological_context_formula <- function(response, added_terms = character(),
                                       mixed_binomial = FALSE) {
  rhs <- paste(c(biological_context_base_rhs(), added_terms), collapse = " + ")
  if (mixed_binomial) {
    rhs <- paste(
      rhs,
      "+ (1 | accepted_species) + (1 | experiment_id)"
    )
  }
  stats::as.formula(paste(response, "~", rhs), env = baseenv())
}

biological_context_specs <- function(kind, trait_column = NULL) {
  axes <- biological_context_axes()
  if (identical(kind, "family")) {
    additions <- c(
      climate = "",
      main = "family_model",
      stats::setNames(
        paste0("family_model + family_model:", axes), names(axes)
      )
    )
  } else if (identical(kind, "functional_group")) {
    additions <- c(
      climate = "",
      main = "functional_group_model",
      interaction = paste(
        c(
          "functional_group_model",
          paste0("functional_group_model:", axes)
        ),
        collapse = " + "
      )
    )
  } else if (identical(kind, "trait")) {
    if (is.null(trait_column) || length(trait_column) != 1L) {
      stop("A single standardized trait column is required.", call. = FALSE)
    }
    additions <- c(
      climate = "",
      main = trait_column,
      interaction = paste(
        c(trait_column, paste0(trait_column, ":", axes)),
        collapse = " + "
      )
    )
  } else {
    stop("Unknown biological-context model kind: ", kind, call. = FALSE)
  }
  lapply(additions, function(addition) {
    terms <- if (nzchar(addition)) addition else character()
    biological_context_formula("resident_outcome", terms)
  })
}

prepare_biological_context_review <- function(data, family_minimum = 50L) {
  required <- c(
    "transform", "horizon", "accepted_name", "family",
    "functional_group", "experiment_id",
    "pool", "analysis_category", "destination_T_JJA",
    "destination_P_JJA", "diff_T_JJA", "LDMC", "LMA",
    "vegetative_height"
  )
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(
      "Biological-context review data are missing: ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  if (length(family_minimum) != 1L || is.na(family_minimum) ||
      family_minimum < 1L) {
    stop("family_minimum must be one positive integer.", call. = FALSE)
  }
  review_input <- data |>
    dplyr::filter(
      as.character(.data$transform) == "raw",
      as.character(.data$horizon) == "max"
    )
  review_input <- filter_downstream_model_species(
    review_input, "accepted_name"
  )
  if (!nrow(review_input)) {
    stop("No raw, maximum-duration records are available for review.",
         call. = FALSE)
  }
  pool_data <- prepare_pool_climate_response_data(review_input)
  species_lookup <- review_input |>
    dplyr::transmute(
      accepted_species = as.character(.data$accepted_name),
      family = trimws(as.character(.data$family)),
      functional_group = as.character(.data$functional_group),
      LDMC = .data$LDMC,
      LMA = .data$LMA,
      vegetative_height = .data$vegetative_height
    ) |>
    dplyr::distinct()
  lookup_duplicates <- species_lookup |>
    dplyr::count(.data$accepted_species) |>
    dplyr::filter(.data$n != 1L)
  if (nrow(lookup_duplicates)) {
    stop(
      "Every accepted species must have one family/group/trait record.",
      call. = FALSE
    )
  }
  species_lookup <- species_lookup |>
    dplyr::mutate(
      family = dplyr::if_else(
        is.na(.data$family) | .data$family == "", "Unresolved", .data$family
      )
    )
  joined <- pool_data$data |>
    dplyr::mutate(accepted_species = as.character(.data$accepted_species)) |>
    dplyr::left_join(species_lookup, by = "accepted_species")
  if (anyNA(joined$family) || anyNA(joined$functional_group)) {
    stop("Family or functional group is missing after the species join.", call. = FALSE)
  }
  family_counts <- joined |>
    dplyr::count(.data$family, name = "n_species_experiment") |>
    dplyr::mutate(
      common_family = .data$n_species_experiment >= .env$family_minimum,
      family_model = dplyr::if_else(
        .data$common_family, .data$family, "Other"
      )
    ) |>
    dplyr::arrange(dplyr::desc(.data$n_species_experiment), .data$family)
  common_families <- family_counts |>
    dplyr::filter(.data$common_family) |>
    dplyr::pull(.data$family)
  family_levels <- c("Other", sort(common_families))
  joined <- joined |>
    dplyr::mutate(
      family_model = factor(
        dplyr::if_else(.data$family %in% .env$common_families,
                       .data$family, "Other"),
        levels = .env$family_levels
      ),
      functional_group_model = factor(
        .data$functional_group, levels = c("forb", "graminoid")
      ),
      accepted_species = factor(.data$accepted_species),
      experiment_id = factor(.data$experiment_id)
    )
  traits <- c("LDMC", "LMA", "vegetative_height")
  trait_scaling <- purrr::map_dfr(traits, function(trait) {
    values <- species_lookup |>
      dplyr::filter(is.finite(.data[[trait]]), .data[[trait]] > 0) |>
      dplyr::distinct(.data$accepted_species, value = .data[[trait]]) |>
      dplyr::mutate(log_value = log(.data$value))
    tibble::tibble(
      trait = trait,
      source_column = trait,
      transformed_column = paste0("log_", trait),
      standardized_column = paste0("z_trait_", trait),
      center = mean(values$log_value),
      scale = stats::sd(values$log_value),
      n_species = nrow(values),
      weighting = "distinct species"
    )
  })
  if (any(!is.finite(trait_scaling$scale)) || any(trait_scaling$scale <= 0)) {
    stop("Every focal trait must have a positive finite species-level scale.",
         call. = FALSE)
  }
  for (index in seq_len(nrow(trait_scaling))) {
    trait <- trait_scaling$trait[[index]]
    joined[[paste0("log_", trait)]] <- dplyr::if_else(
      is.finite(joined[[trait]]) & joined[[trait]] > 0,
      log(joined[[trait]]), NA_real_
    )
    joined[[paste0("z_trait_", trait)]] <-
      (joined[[paste0("log_", trait)]] - trait_scaling$center[[index]]) /
      trait_scaling$scale[[index]]
  }
  resident <- joined |>
    dplyr::filter(.data$pool == "resident") |>
    dplyr::mutate(
      resident_outcome = factor(
        dplyr::case_when(
          .data$analysis_category == "loser" ~ "loser",
          .data$analysis_category == "winner" ~ "winner",
          TRUE ~ "non_significant"
        ),
        levels = c("non_significant", "loser", "winner")
      ),
      loser_event = .data$resident_outcome == "loser",
      winner_event = .data$resident_outcome == "winner"
    )
  nonresident <- joined |>
    dplyr::filter(.data$pool %in% c("colonizer", "outsider")) |>
    dplyr::mutate(
      colonizer_event = .data$analysis_category == "colonizer"
    )
  key_check <- joined |>
    dplyr::count(.data$accepted_species, .data$experiment_id) |>
    dplyr::filter(.data$n != 1L)
  if (nrow(key_check)) {
    stop("Biological-context data must remain one row per species and experiment.",
         call. = FALSE)
  }
  family_event_support <- dplyr::bind_rows(
    resident |>
      dplyr::count(.data$family_model, outcome = .data$resident_outcome,
                   name = "n") |>
      tidyr::complete(
        family_model, outcome = levels(resident$resident_outcome),
        fill = list(n = 0L)
      ) |>
      dplyr::mutate(pool = "resident"),
    nonresident |>
      dplyr::count(
        .data$family_model,
        outcome = dplyr::if_else(
          .data$colonizer_event, "colonizer", "non_significant"
        ),
        name = "n"
      ) |>
      tidyr::complete(
        family_model,
        outcome = c("non_significant", "colonizer"),
        fill = list(n = 0L)
      ) |>
      dplyr::mutate(pool = "nonresident")
  ) |>
    dplyr::mutate(family_model = as.character(.data$family_model))
  functional_group_counts <- joined |>
    dplyr::mutate(
      outcome = dplyr::case_when(
        .data$analysis_category %in% c("loser", "winner", "colonizer") ~
          as.character(.data$analysis_category),
        TRUE ~ "non_significant"
      )
    ) |>
    dplyr::count(.data$functional_group, .data$pool, .data$outcome, name = "n")
  snapshot <- tibble::tibble(
    n_species_experiment = nrow(joined),
    n_resident = nrow(resident),
    n_nonresident = nrow(nonresident),
    n_common_families = length(common_families),
    family_minimum = as.integer(family_minimum),
    n_functional_group_excluded = sum(
      joined$functional_group == "other_unresolved"
    ),
    n_species = dplyr::n_distinct(joined$accepted_species),
    n_experiments = dplyr::n_distinct(joined$experiment_id)
  )
  list(
    data = joined, resident = resident, nonresident = nonresident,
    experiments = pool_data$experiments, climate_scaling = pool_data$scaling,
    trait_scaling = trait_scaling, family_selection = family_counts,
    family_event_support = family_event_support,
    functional_group_counts = functional_group_counts,
    snapshot = snapshot
  )
}

review_fit_subprocess <- function(engine, formula, data, timeout = 45) {
  started <- Sys.time()
  result <- tryCatch(
    callr::r(
      function(engine, formula, data) {
        warnings <- character()
        fit <- tryCatch(
          withCallingHandlers(
            {
              if (engine == "mixed_multinomial") {
                mclogit::mblogit(
                  formula,
                  random = list(
                    ~ 1 | accepted_species,
                    ~ 1 | experiment_id
                  ),
                  data = data, estimator = "ML", trace = FALSE,
                  control = mclogit::mmclogit.control(
                    maxit = 30, maxit.inner = 300, trace = FALSE,
                    trace.inner = FALSE
                  )
                )
              } else if (engine == "binomial_glmm") {
                lme4::glmer(
                  formula, data = data, family = stats::binomial(),
                  control = lme4::glmerControl(
                    optimizer = "bobyqa", optCtrl = list(maxfun = 200000)
                  )
                )
              } else {
                stop("Unknown review fit engine: ", engine)
              }
            },
            warning = function(condition) {
              warnings <<- c(warnings, conditionMessage(condition))
              invokeRestart("muffleWarning")
            }
          ),
          error = function(condition) condition
        )
        if (inherits(fit, "condition")) {
          list(fit = NULL, warnings = unique(warnings),
               error = conditionMessage(fit))
        } else {
          list(fit = fit, warnings = unique(warnings), error = "")
        }
      },
      args = list(engine = engine, formula = formula, data = data),
      timeout = timeout, spinner = FALSE, show = FALSE
    ),
    error = function(condition) {
      list(
        fit = NULL, warnings = character(),
        error = conditionMessage(condition)
      )
    }
  )
  result$elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  result$timed_out <- grepl("timed out|timeout", result$error, ignore.case = TRUE)
  result
}

review_zero_event_levels <- function(data, engine, outcome,
                                     categorical_predictor = NULL) {
  if (is.null(categorical_predictor)) return(character())
  predictor <- as.character(data[[categorical_predictor]])
  if (engine == "mixed_multinomial") {
    response <- as.character(data$resident_outcome)
    cells <- tidyr::expand_grid(
      predictor = unique(predictor),
      response = levels(data$resident_outcome)
    ) |>
      dplyr::left_join(
        tibble::tibble(predictor = predictor, response = response) |>
          dplyr::count(.data$predictor, .data$response, name = "n"),
        by = c("predictor", "response")
      ) |>
      dplyr::mutate(n = dplyr::coalesce(.data$n, 0L)) |>
      dplyr::filter(.data$n == 0L)
    if (!nrow(cells)) return(character())
    return(paste0(cells$predictor, "::", cells$response))
  }
  event <- data$binary_event
  support <- tibble::tibble(predictor = predictor, event = event) |>
    dplyr::summarise(
      n = dplyr::n(), events = sum(.data$event), .by = "predictor"
    ) |>
    dplyr::filter(.data$events == 0L | .data$events == .data$n)
  if (!nrow(support)) return(character())
  paste0(support$predictor, "::", support$events, "/", support$n)
}

review_model_diagnostic <- function(captured, engine, subsystem, outcome,
                                    stage, data, zero_event_levels = character(),
                                    selected = FALSE) {
  fit <- captured$fit
  base <- tibble::tibble(
    subsystem = subsystem, outcome = outcome, stage = stage, engine = engine,
    selected = selected, n = nrow(data),
    events = if (engine == "mixed_multinomial") {
      sum(data$resident_outcome != "non_significant")
    } else {
      sum(data$binary_event)
    },
    n_species = dplyr::n_distinct(data$accepted_species),
    n_experiments = dplyr::n_distinct(data$experiment_id),
    elapsed_seconds = captured$elapsed_seconds,
    timed_out = captured$timed_out,
    error = captured$error,
    warning = paste(unique(captured$warnings), collapse = " | "),
    zero_event_levels = paste(zero_event_levels, collapse = " | ")
  )
  if (is.null(fit)) {
    return(base |>
      dplyr::mutate(
        converged = FALSE, singular = NA, covariance_finite = FALSE,
        covariance_minimum_eigenvalue = NA_real_,
        hessian_minimum_eigenvalue = NA_real_,
        random_covariance_minimum_eigenvalue = NA_real_,
        species_variance = NA_real_, experiment_variance = NA_real_,
        max_absolute_coefficient = NA_real_, max_standard_error = NA_real_,
        separated = length(zero_event_levels) > 0L,
        clean = FALSE,
        status = dplyr::if_else(.data$timed_out, "timeout", "error")
      ))
  }
  coefficients <- if (engine == "mixed_multinomial") {
    stats::coef(fit)
  } else {
    lme4::fixef(fit)
  }
  covariance <- tryCatch(
    as.matrix(stats::vcov(fit)), error = function(condition) NULL
  )
  covariance_finite <- !is.null(covariance) && all(is.finite(covariance))
  covariance_eigen <- if (covariance_finite) {
    eigen(covariance, symmetric = TRUE, only.values = TRUE)$values
  } else {
    NA_real_
  }
  standard_errors <- if (covariance_finite) sqrt(diag(covariance)) else NA_real_
  if (engine == "binomial_glmm") {
    optimizer_messages <- unlist(
      fit@optinfo$conv$lme4$messages, use.names = FALSE
    )
    optimizer_messages <- optimizer_messages[!is.na(optimizer_messages)]
    hessian <- fit@optinfo$derivs$Hessian
    hessian_eigen <- if (!is.null(hessian) && all(is.finite(hessian))) {
      eigen(hessian, symmetric = TRUE, only.values = TRUE)$values
    } else {
      NA_real_
    }
    random <- as.data.frame(lme4::VarCorr(fit))
    species_values <- random$vcov[random$grp == "accepted_species"]
    experiment_values <- random$vcov[random$grp == "experiment_id"]
    species_variance <- if (length(species_values)) species_values[[1]] else NA_real_
    experiment_variance <- if (length(experiment_values)) {
      experiment_values[[1]]
    } else {
      NA_real_
    }
    random_min <- if (any(is.finite(random$vcov))) {
      min(random$vcov[is.finite(random$vcov)])
    } else {
      NA_real_
    }
    converged <- isTRUE(fit@optinfo$conv$opt == 0) &&
      length(optimizer_messages) == 0L
    singular <- lme4::isSingular(fit, tol = 1e-6)
    extra_warning <- paste(optimizer_messages, collapse = " | ")
  } else {
    hessian_eigen <- if (covariance_finite) 1 / rev(covariance_eigen) else NA_real_
    random_eigen <- tryCatch(
      unlist(lapply(
        fit$VarCov,
        function(matrix) eigen(
          matrix, symmetric = TRUE, only.values = TRUE
        )$values
      )),
      error = function(condition) NA_real_
    )
    species_variance <- NA_real_
    experiment_variance <- NA_real_
    random_min <- if (any(is.finite(random_eigen))) {
      min(random_eigen[is.finite(random_eigen)])
    } else {
      NA_real_
    }
    converged <- isTRUE(fit$converged)
    singular <- !all(is.finite(random_eigen)) || random_min <= 1e-6
    extra_warning <- ""
  }
  warning_text <- paste(
    unique(c(captured$warnings, extra_warning)[nzchar(c(captured$warnings, extra_warning))]),
    collapse = " | "
  )
  max_beta <- if (any(is.finite(coefficients))) {
    max(abs(coefficients[is.finite(coefficients)]))
  } else {
    NA_real_
  }
  max_se <- if (any(is.finite(standard_errors))) {
    max(standard_errors[is.finite(standard_errors)])
  } else {
    NA_real_
  }
  separated <- length(zero_event_levels) > 0L || max_beta > 20 || max_se > 100
  hessian_min <- if (any(is.finite(hessian_eigen))) {
    min(hessian_eigen[is.finite(hessian_eigen)])
  } else {
    NA_real_
  }
  covariance_min <- if (any(is.finite(covariance_eigen))) {
    min(covariance_eigen[is.finite(covariance_eigen)])
  } else {
    NA_real_
  }
  clean <- converged && !singular && covariance_finite &&
    is.finite(covariance_min) && covariance_min > 0 &&
    is.finite(hessian_min) && hessian_min > 0 &&
    !separated && !nzchar(captured$error) && !nzchar(warning_text)
  status <- dplyr::case_when(
    captured$timed_out ~ "timeout",
    nzchar(captured$error) ~ "error",
    nzchar(warning_text) ~ "optimizer_warning",
    !converged ~ "not_converged",
    separated ~ "separated",
    singular ~ "singular",
    !covariance_finite | !is.finite(covariance_min) | covariance_min <= 0 ~
      "invalid_covariance",
    !is.finite(hessian_min) | hessian_min <= 0 ~ "nonpositive_hessian",
    TRUE ~ "clean"
  )
  base |>
    dplyr::mutate(
      warning = warning_text, converged = converged, singular = singular,
      covariance_finite = covariance_finite,
      covariance_minimum_eigenvalue = covariance_min,
      hessian_minimum_eigenvalue = hessian_min,
      random_covariance_minimum_eigenvalue = random_min,
      species_variance = species_variance,
      experiment_variance = experiment_variance,
      max_absolute_coefficient = max_beta,
      max_standard_error = max_se,
      separated = separated, clean = clean, status = status
    )
}

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x)) y else x
}

review_skipped_capture <- function(reason) {
  list(
    fit = NULL, warnings = character(), error = reason,
    elapsed_seconds = 0, timed_out = FALSE
  )
}

fit_review_model_set <- function(data, specs, subsystem,
                                 categorical_predictor = NULL,
                                 mixed_timeout = 35, binary_timeout = 90) {
  resident <- data$resident
  nonresident <- data$nonresident
  mixed_models <- list()
  mixed_diagnostics <- list()
  mixed_ok <- TRUE
  for (stage in names(specs)) {
    if (!mixed_ok) break
    formula <- specs[[stage]]
    zero <- review_zero_event_levels(
      resident, "mixed_multinomial", "resident_joint",
      if (stage == "climate") NULL else categorical_predictor
    )
    captured <- review_fit_subprocess(
      "mixed_multinomial", formula, resident, timeout = mixed_timeout
    )
    diagnostic <- review_model_diagnostic(
      captured, "mixed_multinomial", subsystem, "resident_joint", stage,
      resident, zero_event_levels = zero
    )
    mixed_models[[stage]] <- captured$fit
    mixed_diagnostics[[stage]] <- diagnostic
    mixed_ok <- isTRUE(diagnostic$clean[[1]])
  }
  selected_engine <- if (mixed_ok && length(mixed_models) == length(specs)) {
    "mixed_multinomial"
  } else {
    "binomial_glmm"
  }
  models <- list()
  diagnostics <- list()
  if (selected_engine == "mixed_multinomial") {
    models$resident_joint <- mixed_models
    diagnostics <- mixed_diagnostics
  } else {
    diagnostics <- mixed_diagnostics
    resident_events <- c(loser = "loser_event", winner = "winner_event")
    for (outcome in names(resident_events)) {
      frame <- resident |>
        dplyr::mutate(binary_event = .data[[resident_events[[outcome]]]]) |>
        droplevels()
      models[[outcome]] <- list()
      can_continue <- TRUE
      for (stage in names(specs)) {
        fixed_terms <- attr(stats::terms(specs[[stage]]), "term.labels")
        formula <- biological_context_formula(
          "binary_event", fixed_terms, mixed_binomial = TRUE
        )
        zero <- review_zero_event_levels(
          frame, "binomial_glmm", outcome,
          if (stage == "climate") NULL else categorical_predictor
        )
        captured <- if (can_continue) {
          review_fit_subprocess(
            "binomial_glmm", formula, frame, timeout = binary_timeout
          )
        } else {
          review_skipped_capture(
            "Not attempted because an earlier nested model failed diagnostics."
          )
        }
        diagnostic <- review_model_diagnostic(
          captured, "binomial_glmm", subsystem, outcome, stage, frame,
          zero_event_levels = zero, selected = TRUE
        )
        models[[outcome]][[stage]] <- captured$fit
        diagnostics[[paste(outcome, stage, sep = "__")]] <- diagnostic
        can_continue <- can_continue && isTRUE(diagnostic$clean[[1]])
      }
    }
  }
  colonizer_frame <- nonresident |>
    dplyr::mutate(binary_event = .data$colonizer_event) |>
    droplevels()
  models$colonizer <- list()
  colonizer_can_continue <- TRUE
  for (stage in names(specs)) {
    fixed_terms <- attr(stats::terms(specs[[stage]]), "term.labels")
    formula <- biological_context_formula(
      "binary_event", fixed_terms, mixed_binomial = TRUE
    )
    zero <- review_zero_event_levels(
      colonizer_frame, "binomial_glmm", "colonizer",
      if (stage == "climate") NULL else categorical_predictor
    )
    captured <- if (colonizer_can_continue) {
      review_fit_subprocess(
        "binomial_glmm", formula, colonizer_frame, timeout = binary_timeout
      )
    } else {
      review_skipped_capture(
        "Not attempted because an earlier nested model failed diagnostics."
      )
    }
    diagnostic <- review_model_diagnostic(
      captured, "binomial_glmm", subsystem, "colonizer", stage,
      colonizer_frame, zero_event_levels = zero, selected = TRUE
    )
    models$colonizer[[stage]] <- captured$fit
    diagnostics[[paste("colonizer", stage, sep = "__")]] <- diagnostic
    colonizer_can_continue <- colonizer_can_continue &&
      isTRUE(diagnostic$clean[[1]])
  }
  diagnostic_table <- dplyr::bind_rows(diagnostics) |>
    dplyr::mutate(
      selected = dplyr::if_else(
        .data$engine == selected_engine | .data$outcome == "colonizer",
        TRUE, .data$selected
      ),
      selected_resident_engine = selected_engine
    )
  list(
    models = models, attempted_mixed_models = mixed_models,
    diagnostics = diagnostic_table,
    selected_resident_engine = selected_engine,
    resident_data = resident, nonresident_data = nonresident,
    specs = specs, subsystem = subsystem
  )
}

review_subsystem_data <- function(prepared, kind, trait = NULL) {
  if (kind == "family") {
    resident <- prepared$resident |> droplevels()
    nonresident <- prepared$nonresident |> droplevels()
    predictor <- "family_model"
  } else if (kind == "functional_group") {
    resident <- prepared$resident |>
      dplyr::filter(!is.na(.data$functional_group_model)) |>
      droplevels()
    nonresident <- prepared$nonresident |>
      dplyr::filter(!is.na(.data$functional_group_model)) |>
      droplevels()
    predictor <- "functional_group_model"
  } else if (kind == "trait") {
    column <- paste0("z_trait_", trait)
    resident <- prepared$resident |>
      dplyr::filter(is.finite(.data[[column]])) |>
      droplevels()
    nonresident <- prepared$nonresident |>
      dplyr::filter(is.finite(.data[[column]])) |>
      droplevels()
    predictor <- NULL
  } else {
    stop("Unknown subsystem kind.", call. = FALSE)
  }
  list(resident = resident, nonresident = nonresident,
       categorical_predictor = predictor)
}

fit_biological_context_review <- function(prepared, mixed_timeout = 35,
                                          binary_timeout = 90) {
  family_data <- review_subsystem_data(prepared, "family")
  family <- fit_review_model_set(
    family_data, biological_context_specs("family"), "family",
    categorical_predictor = family_data$categorical_predictor,
    mixed_timeout = mixed_timeout, binary_timeout = binary_timeout
  )
  group_data <- review_subsystem_data(prepared, "functional_group")
  functional_group <- fit_review_model_set(
    group_data, biological_context_specs("functional_group"),
    "functional_group",
    categorical_predictor = group_data$categorical_predictor,
    mixed_timeout = mixed_timeout, binary_timeout = binary_timeout
  )
  traits <- purrr::map(
    c("LDMC", "LMA", "vegetative_height"),
    function(trait) {
      trait_data <- review_subsystem_data(prepared, "trait", trait)
      fit_review_model_set(
        trait_data,
        biological_context_specs("trait", paste0("z_trait_", trait)),
        paste0("trait_", trait), categorical_predictor = NULL,
        mixed_timeout = mixed_timeout, binary_timeout = binary_timeout
      )
    }
  )
  names(traits) <- c("LDMC", "LMA", "vegetative_height")
  list(family = family, functional_group = functional_group, traits = traits)
}

review_selected_model <- function(bundle, outcome, stage) {
  if (outcome %in% c("loser", "winner") &&
      identical(bundle$selected_resident_engine, "mixed_multinomial")) {
    return(bundle$models$resident_joint[[stage]])
  }
  bundle$models[[outcome]][[stage]]
}

review_selected_engine <- function(bundle, outcome) {
  if (outcome %in% c("loser", "winner")) {
    bundle$selected_resident_engine
  } else {
    "binomial_glmm"
  }
}

review_selected_diagnostic <- function(bundle, outcome, stage) {
  diagnostic_outcome <- if (
    outcome %in% c("loser", "winner") &&
      identical(bundle$selected_resident_engine, "mixed_multinomial")
  ) {
    "resident_joint"
  } else {
    outcome
  }
  bundle$diagnostics |>
    dplyr::filter(
      .data$outcome == .env$diagnostic_outcome,
      .data$stage == .env$stage,
      .data$selected
    ) |>
    dplyr::slice_head(n = 1L)
}

tidy_review_fit <- function(fit, engine, subsystem, outcome, stage,
                            valid_for_inference) {
  if (is.null(fit)) return(tibble::tibble())
  coefficients <- if (engine == "mixed_multinomial") {
    stats::coef(fit)
  } else {
    lme4::fixef(fit)
  }
  covariance <- as.matrix(stats::vcov(fit))
  if (engine == "mixed_multinomial") {
    prefix <- paste0(outcome, "~")
    keep <- startsWith(names(coefficients), prefix)
    coefficients <- coefficients[keep]
    coefficient_names <- names(coefficients)
    terms <- sub(paste0("^", outcome, "~"), "", coefficient_names)
  } else {
    coefficient_names <- names(coefficients)
    terms <- coefficient_names
  }
  standard_errors <- sqrt(diag(covariance)[coefficient_names])
  statistic <- coefficients / standard_errors
  p_value <- 2 * stats::pnorm(abs(statistic), lower.tail = FALSE)
  tibble::tibble(
    subsystem = subsystem, outcome = outcome, stage = stage, engine = engine,
    term = terms, estimate = unname(coefficients),
    standard_error = unname(standard_errors), statistic = unname(statistic),
    p_value = unname(p_value), odds_ratio = exp(estimate),
    odds_low = exp(estimate - stats::qnorm(0.975) * standard_error),
    odds_high = exp(estimate + stats::qnorm(0.975) * standard_error),
    valid_for_inference = valid_for_inference
  )
}

collect_review_coefficients_one <- function(bundle) {
  purrr::map_dfr(c("loser", "winner", "colonizer"), function(outcome) {
    engine <- review_selected_engine(bundle, outcome)
    purrr::map_dfr(names(bundle$specs), function(stage) {
      diagnostic <- review_selected_diagnostic(bundle, outcome, stage)
      valid <- nrow(diagnostic) == 1L && isTRUE(diagnostic$clean[[1]])
      tidy_review_fit(
        review_selected_model(bundle, outcome, stage), engine,
        bundle$subsystem, outcome, stage, valid
      )
    })
  }) |>
    dplyr::mutate(
      p_holm = NA_real_,
      p_adjustment_note = paste(
        "No multiplicity adjustment is applied to individual coefficients;",
        "trait Holm adjustment is applied to interaction-block LRTs."
      )
    )
}

collect_biological_context_coefficients <- function(fits) {
  dplyr::bind_rows(
    collect_review_coefficients_one(fits$family),
    collect_review_coefficients_one(fits$functional_group),
    purrr::map_dfr(fits$traits, collect_review_coefficients_one)
  )
}

review_likelihood_ratio <- function(bundle, outcome, reduced_stage, full_stage,
                                    comparison) {
  reduced <- review_selected_model(bundle, outcome, reduced_stage)
  full <- review_selected_model(bundle, outcome, full_stage)
  reduced_diagnostic <- review_selected_diagnostic(bundle, outcome, reduced_stage)
  full_diagnostic <- review_selected_diagnostic(bundle, outcome, full_stage)
  valid <- nrow(reduced_diagnostic) == 1L && nrow(full_diagnostic) == 1L &&
    isTRUE(reduced_diagnostic$clean[[1]]) &&
    isTRUE(full_diagnostic$clean[[1]]) &&
    !is.null(reduced) && !is.null(full)
  if (!valid) {
    return(tibble::tibble(
      subsystem = bundle$subsystem, outcome = outcome,
      comparison = comparison, reduced_stage = reduced_stage,
      full_stage = full_stage, statistic = NA_real_, df = NA_real_,
      p_value = NA_real_, p_holm = NA_real_, valid_for_inference = FALSE,
      status = "not tested because one or both nested models failed diagnostics"
    ))
  }
  reduced_loglik <- stats::logLik(reduced)
  full_loglik <- stats::logLik(full)
  degrees <- attr(full_loglik, "df") - attr(reduced_loglik, "df")
  statistic <- 2 * as.numeric(full_loglik - reduced_loglik)
  tibble::tibble(
    subsystem = bundle$subsystem, outcome = outcome,
    comparison = comparison, reduced_stage = reduced_stage,
    full_stage = full_stage, statistic = statistic, df = degrees,
    p_value = stats::pchisq(statistic, df = degrees, lower.tail = FALSE),
    p_holm = NA_real_, valid_for_inference = TRUE, status = "clean"
  )
}

collect_review_likelihood_ratios_one <- function(bundle) {
  outcomes <- if (identical(bundle$selected_resident_engine,
                            "mixed_multinomial")) {
    c("resident_joint", "colonizer")
  } else {
    c("loser", "winner", "colonizer")
  }
  model_outcome <- function(outcome) {
    if (outcome == "resident_joint") "loser" else outcome
  }
  purrr::map_dfr(outcomes, function(outcome_label) {
    outcome <- model_outcome(outcome_label)
    if (bundle$subsystem == "family") {
      comparisons <- tibble::tribble(
        ~reduced, ~full, ~comparison,
        "climate", "main", "family main effect",
        "main", "temperature", "family by temperature",
        "main", "precipitation", "family by precipitation",
        "main", "warming", "family by warming"
      )
    } else {
      comparisons <- tibble::tribble(
        ~reduced, ~full, ~comparison,
        "climate", "main", "biological main effect",
        "main", "interaction", "biological context interactions"
      )
    }
    purrr::pmap_dfr(comparisons, function(reduced, full, comparison) {
      review_likelihood_ratio(bundle, outcome, reduced, full, comparison) |>
        dplyr::mutate(outcome = outcome_label)
    })
  })
}

collect_biological_context_likelihood_ratios <- function(fits) {
  tests <- dplyr::bind_rows(
    collect_review_likelihood_ratios_one(fits$family),
    collect_review_likelihood_ratios_one(fits$functional_group),
    purrr::map_dfr(fits$traits, collect_review_likelihood_ratios_one)
  )
  trait_interactions <- stringr::str_starts(tests$subsystem, "trait_") &
    tests$comparison == "biological context interactions" &
    is.finite(tests$p_value)
  tests$p_holm[trait_interactions] <- ave(
    tests$p_value[trait_interactions],
    tests$outcome[trait_interactions],
    FUN = function(value) stats::p.adjust(value, method = "holm")
  )
  tests |>
    dplyr::mutate(
      holm_family = dplyr::if_else(
        stringr::str_starts(.data$subsystem, "trait_") &
          .data$comparison == "biological context interactions",
        paste0("three trait interaction blocks within ", .data$outcome),
        NA_character_
      ),
      holm_planned_tests = dplyr::if_else(
        !is.na(.data$holm_family), 3L, NA_integer_
      )
    ) |>
    dplyr::group_by(.data$outcome) |>
    dplyr::mutate(
      holm_estimable_tests = dplyr::if_else(
        !is.na(.data$holm_family),
        sum(
          stringr::str_starts(.data$subsystem, "trait_") &
            .data$comparison == "biological context interactions" &
            is.finite(.data$p_value)
        ),
        NA_integer_
      )
    ) |>
    dplyr::ungroup()
}

collect_biological_context_diagnostics <- function(fits) {
  dplyr::bind_rows(
    fits$family$diagnostics,
    fits$functional_group$diagnostics,
    purrr::map_dfr(fits$traits, ~ .x$diagnostics)
  )
}

review_fixed_components <- function(fit, engine, newdata) {
  if (engine == "mixed_multinomial") {
    coefficients <- stats::coef(fit)
    covariance <- as.matrix(stats::vcov(fit))
    rhs <- stats::delete.response(stats::terms(fit))
    design <- stats::model.matrix(
      rhs, data = newdata, contrasts.arg = fit$contrasts,
      xlev = fit$xlevels
    )
  } else {
    coefficients <- lme4::fixef(fit)
    covariance <- as.matrix(stats::vcov(fit))
    fixed_formula <- reformulas::nobars(stats::formula(fit))
    rhs <- stats::delete.response(stats::terms(fixed_formula))
    design <- stats::model.matrix(rhs, data = newdata)
  }
  list(coefficients = coefficients, covariance = covariance, design = design)
}

review_probability_matrix <- function(fit, engine, outcome, newdata,
                                      coefficient_draws) {
  components <- review_fixed_components(fit, engine, newdata)
  design <- components$design
  if (engine == "mixed_multinomial") {
    loser_columns <- match(
      paste0("loser~", colnames(design)), colnames(coefficient_draws)
    )
    winner_columns <- match(
      paste0("winner~", colnames(design)), colnames(coefficient_draws)
    )
    if (anyNA(c(loser_columns, winner_columns))) {
      stop("Mixed multinomial coefficients do not match the prediction design.",
           call. = FALSE)
    }
    loser_eta <- design %*% t(
      coefficient_draws[, loser_columns, drop = FALSE]
    )
    winner_eta <- design %*% t(
      coefficient_draws[, winner_columns, drop = FALSE]
    )
    maximum_eta <- pmax(0, loser_eta, winner_eta)
    reference_exp <- exp(-maximum_eta)
    loser_exp <- exp(loser_eta - maximum_eta)
    winner_exp <- exp(winner_eta - maximum_eta)
    denominator <- reference_exp + loser_exp + winner_exp
    if (outcome == "loser") loser_exp / denominator else winner_exp / denominator
  } else {
    selected <- match(colnames(design), colnames(coefficient_draws))
    if (anyNA(selected)) {
      stop("Binomial coefficients do not match the prediction design.",
           call. = FALSE)
    }
    stats::plogis(
      design %*% t(coefficient_draws[, selected, drop = FALSE])
    )
  }
}

review_probability_draws <- function(fit, engine, outcome, newdata,
                                     n_draws = 2000L, seed = 47163L) {
  components <- review_fixed_components(fit, engine, newdata)
  set.seed(seed)
  draws <- pool_climate_mvn_draws(
    components$coefficients, components$covariance, n_draws
  )
  central <- matrix(components$coefficients, nrow = 1L)
  colnames(central) <- names(components$coefficients)
  list(
    draws = review_probability_matrix(
      fit, engine, outcome, newdata, draws
    ),
    central = as.numeric(review_probability_matrix(
      fit, engine, outcome, newdata, central
    ))
  )
}

review_probability_summary <- function(probabilities) {
  intervals <- matrixStats::rowQuantiles(
    probabilities$draws, probs = c(0.025, 0.975), drop = FALSE
  )
  tibble::tibble(
    estimate = probabilities$central,
    conf_low = intervals[, 1], conf_high = intervals[, 2],
    interval_width = intervals[, 2] - intervals[, 1]
  )
}

review_axis_grid <- function(prepared, axis, length_out = 60L) {
  variable <- biological_context_axes()[[axis]]
  values <- prepared$experiments[[variable]]
  scaling_variable <- c(
    temperature = "destination_T",
    precipitation = "log_destination_P",
    warming = "warming"
  )[[axis]]
  scaling <- prepared$climate_scaling |>
    dplyr::filter(.data$variable == .env$scaling_variable)
  standardized <- seq(min(values), max(values), length.out = length_out)
  transformed <- scaling$center[[1]] + standardized * scaling$scale[[1]]
  raw <- if (axis == "precipitation") exp(transformed) else transformed
  tibble::tibble(
    axis = axis, axis_variable = variable,
    axis_value = standardized, raw_axis_value = raw
  )
}

review_baseline_frame <- function(n) {
  tibble::tibble(
    z_destination_T = rep(0, n),
    z_log_destination_P = rep(0, n),
    z_warming = rep(0, n)
  )
}

review_family_context_predictions <- function(fit, prepared,
                                              n_draws = 2000L,
                                              seed = 47163L) {
  families <- levels(prepared$resident$family_model)
  outcomes <- c("loser", "winner", "colonizer")
  axes <- names(biological_context_axes())
  purrr::map_dfr(seq_along(outcomes), function(outcome_index) {
    outcome <- outcomes[[outcome_index]]
    engine <- review_selected_engine(fit, outcome)
    pool_data <- if (outcome == "colonizer") {
      fit$nonresident_data
    } else {
      fit$resident_data
    }
    purrr::map_dfr(seq_along(axes), function(axis_index) {
      axis <- axes[[axis_index]]
      stage <- axis
      diagnostic <- review_selected_diagnostic(fit, outcome, stage)
      valid <- nrow(diagnostic) == 1L && isTRUE(diagnostic$clean[[1]])
      model <- review_selected_model(fit, outcome, stage)
      if (!valid || is.null(model)) {
        return(tibble::tibble(
          family_model = families, outcome = outcome, axis = axis,
          estimate = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
          supported = FALSE,
          status = if (nrow(diagnostic)) diagnostic$status[[1]] else "missing"
        ))
      }
      axis_variable <- biological_context_axes()[[axis]]
      endpoints <- stats::quantile(
        prepared$experiments[[axis_variable]], c(0.25, 0.75), names = FALSE
      )
      newdata <- dplyr::bind_rows(
        review_baseline_frame(length(families)) |>
          dplyr::mutate(endpoint = "low"),
        review_baseline_frame(length(families)) |>
          dplyr::mutate(endpoint = "high")
      ) |>
        dplyr::mutate(
          family_model = factor(
            rep(families, times = 2L), levels = families
          ),
          !!axis_variable := rep(endpoints, each = length(families))
        )
      probabilities <- review_probability_draws(
        model, engine, outcome, newdata, n_draws,
        seed + outcome_index * 100L + axis_index
      )
      delta_draws <- probabilities$draws[
        seq_along(families) + length(families), , drop = FALSE
      ] - probabilities$draws[seq_along(families), , drop = FALSE]
      delta_intervals <- matrixStats::rowQuantiles(
        delta_draws, probs = c(0.025, 0.975), drop = FALSE
      )
      delta_central <- probabilities$central[
        seq_along(families) + length(families)
      ] - probabilities$central[seq_along(families)]
      support <- pool_data |>
        dplyr::summarise(
          minimum = min(.data[[axis_variable]]),
          maximum = max(.data[[axis_variable]]),
          .by = "family_model"
        ) |>
        dplyr::mutate(
          family_model = as.character(.data$family_model),
          supported = .data$minimum <= endpoints[[1]] &
            .data$maximum >= endpoints[[2]]
        )
      tibble::tibble(
        family_model = families, outcome = outcome, axis = axis,
        estimate = delta_central,
        conf_low = delta_intervals[, 1], conf_high = delta_intervals[, 2],
        status = "clean"
      ) |>
        dplyr::left_join(
          dplyr::select(support, "family_model", "supported"),
          by = "family_model"
        ) |>
        dplyr::mutate(supported = dplyr::coalesce(.data$supported, FALSE))
    })
  })
}

review_curve_predictions <- function(fit, prepared, kind,
                                     n_draws = 2000L, seed = 47163L,
                                     grid_size = 60L, trait = NULL) {
  outcomes <- c("loser", "winner", "colonizer")
  axes <- names(biological_context_axes())
  stage <- "interaction"
  purrr::map_dfr(seq_along(outcomes), function(outcome_index) {
    outcome <- outcomes[[outcome_index]]
    engine <- review_selected_engine(fit, outcome)
    diagnostic <- review_selected_diagnostic(fit, outcome, stage)
    valid <- nrow(diagnostic) == 1L && isTRUE(diagnostic$clean[[1]])
    model <- review_selected_model(fit, outcome, stage)
    purrr::map_dfr(seq_along(axes), function(axis_index) {
      axis <- axes[[axis_index]]
      axis_grid <- review_axis_grid(prepared, axis, grid_size)
      if (kind == "functional_group") {
        focal_levels <- c("forb", "graminoid")
        level_column <- "functional_group_model"
        level_label <- focal_levels
      } else {
        focal_levels <- c(-1, 0, 1)
        level_column <- paste0("z_trait_", trait)
        level_label <- c("Low", "Mean", "High")
      }
      base <- tidyr::crossing(
        axis_grid,
        level_index = seq_along(focal_levels)
      ) |>
        dplyr::bind_cols(review_baseline_frame(nrow(axis_grid) * length(focal_levels)))
      base[[biological_context_axes()[[axis]]]] <- base$axis_value
      if (kind == "functional_group") {
        base[[level_column]] <- factor(
          focal_levels[base$level_index], levels = focal_levels
        )
      } else {
        base[[level_column]] <- focal_levels[base$level_index]
      }
      base$level <- level_label[base$level_index]
      if (!valid || is.null(model)) {
        return(base |>
          dplyr::transmute(
            outcome = .env$outcome, axis = .data$axis,
            axis_variable = .data$axis_variable,
            axis_value = .data$axis_value,
            raw_axis_value = .data$raw_axis_value, level = .data$level,
            estimate = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
            supported = FALSE,
            status = if (nrow(diagnostic)) diagnostic$status[[1]] else "missing"
          ))
      }
      probabilities <- review_probability_draws(
        model, engine, outcome, base, n_draws,
        seed + outcome_index * 100L + axis_index
      )
      summary <- review_probability_summary(probabilities)
      dplyr::bind_cols(base, summary) |>
        dplyr::transmute(
          outcome = .env$outcome, axis = .data$axis,
          axis_variable = .data$axis_variable,
          axis_value = .data$axis_value,
          raw_axis_value = .data$raw_axis_value, level = .data$level,
          estimate = .data$estimate, conf_low = .data$conf_low,
          conf_high = .data$conf_high, supported = TRUE, status = "clean"
        )
    })
  }) |>
    dplyr::mutate(
      model_kind = kind,
      trait = if (kind == "trait") trait else NA_character_
    )
}

review_all_predictions <- function(fits, prepared, n_draws = 2000L,
                                   seed = 47163L, grid_size = 60L) {
  list(
    family = review_family_context_predictions(
      fits$family, prepared, n_draws, seed
    ),
    functional_group = review_curve_predictions(
      fits$functional_group, prepared, "functional_group",
      n_draws, seed + 1000L, grid_size
    ),
    traits = purrr::imap_dfr(fits$traits, function(fit, trait) {
      review_curve_predictions(
        fit, prepared, "trait", n_draws,
        seed + 2000L + match(trait, names(fits$traits)) * 100L,
        grid_size, trait
      )
    })
  )
}

validate_biological_context_review <- function(prepared, fits, predictions) {
  key_duplicates <- prepared$data |>
    dplyr::count(.data$accepted_species, .data$experiment_id) |>
    dplyr::filter(.data$n != 1L)
  common_valid <- all(
    prepared$family_selection$n_species_experiment[
      prepared$family_selection$common_family
    ] >= prepared$snapshot$family_minimum[[1]]
  ) && all(
    prepared$family_selection$n_species_experiment[
      !prepared$family_selection$common_family
    ] < prepared$snapshot$family_minimum[[1]]
  )
  formula_text <- unlist(purrr::map(
    c(
      biological_context_specs("family"),
      biological_context_specs("functional_group"),
      biological_context_specs("trait", "z_trait_LDMC")
    ),
    ~ paste(deparse(.x), collapse = " ")
  ))
  diagnostics <- collect_biological_context_diagnostics(fits)
  fallback_valid <- purrr::every(
    list(fits$family, fits$functional_group, fits$traits[[1]],
         fits$traits[[2]], fits$traits[[3]]),
    function(bundle) {
      if (bundle$selected_resident_engine == "mixed_multinomial") return(TRUE)
      mixed <- bundle$diagnostics |>
        dplyr::filter(.data$engine == "mixed_multinomial")
      nrow(mixed) > 0L && any(!mixed$clean)
    }
  )
  probability_values <- c(
    predictions$functional_group$estimate,
    predictions$functional_group$conf_low,
    predictions$functional_group$conf_high,
    predictions$traits$estimate,
    predictions$traits$conf_low,
    predictions$traits$conf_high
  )
  probability_values <- probability_values[is.finite(probability_values)]
  family_values <- c(
    predictions$family$estimate,
    predictions$family$conf_low,
    predictions$family$conf_high
  )
  family_values <- family_values[is.finite(family_values)]
  checks <- tibble::tribble(
    ~check, ~passed, ~detail,
    "resident denominator", nrow(prepared$resident) == 1604L,
      as.character(nrow(prepared$resident)),
    "nonresident denominator", nrow(prepared$nonresident) == 611L,
      as.character(nrow(prepared$nonresident)),
    "one row per species and experiment", nrow(key_duplicates) == 0L,
      as.character(nrow(key_duplicates)),
    "13 denominator-selected families", prepared$snapshot$n_common_families[[1]] == 13L,
      as.character(prepared$snapshot$n_common_families[[1]]),
    "family threshold is reversible", common_valid,
      paste0("minimum = ", prepared$snapshot$family_minimum[[1]]),
    "forb and graminoid inference filter",
      all(stats::na.omit(unique(as.character(
        fits$functional_group$resident_data$functional_group_model
      ))) %in% c("forb", "graminoid")),
      paste(sort(unique(as.character(
        fits$functional_group$resident_data$functional_group_model
      ))), collapse = ", "),
    "distinct-species trait scaling",
      all(prepared$trait_scaling$weighting == "distinct species") &&
        all(is.finite(prepared$trait_scaling$scale)) &&
        all(prepared$trait_scaling$scale > 0),
      paste(prepared$trait_scaling$n_species, collapse = ", "),
    "no precipitation change or three-way term",
      !any(grepl("diff_P|z_destination_T:z_log_destination_P:z_warming",
                 formula_text, fixed = FALSE)),
      paste(formula_text, collapse = " | "),
    "only prescribed model engines",
      all(diagnostics$engine %in% c("mixed_multinomial", "binomial_glmm")),
      paste(unique(diagnostics$engine), collapse = ", "),
    "multinomial fallback follows failed diagnostics", fallback_valid,
      paste(c(
        fits$family$selected_resident_engine,
        fits$functional_group$selected_resident_engine,
        purrr::map_chr(fits$traits, "selected_resident_engine")
      ), collapse = ", "),
    "probability predictions bounded",
      !length(probability_values) ||
        all(probability_values >= 0 & probability_values <= 1),
      as.character(length(probability_values)),
    "family probability changes bounded",
      !length(family_values) || all(family_values >= -1 & family_values <= 1),
      as.character(length(family_values))
  )
  if (any(!checks$passed)) {
    failed <- checks$check[!checks$passed]
    stop(
      "Biological-context review validation failed: ",
      paste(failed, collapse = ", "), call. = FALSE
    )
  }
  checks
}
