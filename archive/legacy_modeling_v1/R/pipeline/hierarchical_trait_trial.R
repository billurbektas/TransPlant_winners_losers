hierarchical_ldmc_fixed_formula <- function(response = "response") {
  stats::as.formula(paste(
    response, "~ temp_z + precip_z + warming_z +",
    "temp_z:LDMC_z + precip_z:LDMC_z + warming_z:LDMC_z"
  ))
}

hierarchical_ldmc_random_formulas <- function() {
  list(
    stats::as.formula("~ 1 | species"),
    stats::as.formula("~ 0 + temp_z | species"),
    stats::as.formula("~ 0 + precip_z | species"),
    stats::as.formula("~ 0 + warming_z | species"),
    stats::as.formula("~ 1 | experiment"),
    stats::as.formula("~ 1 | gradient")
  )
}

hierarchical_ldmc_binomial_formula <- function(response = "binary_response") {
  fixed <- paste(deparse(hierarchical_ldmc_fixed_formula(response)), collapse = " ")
  stats::as.formula(paste(
    fixed,
    "+ (1 | species) + (0 + temp_z | species) +",
    "(0 + precip_z | species) + (0 + warming_z | species) +",
    "(1 | experiment) + (1 | gradient)"
  ))
}

prepare_hierarchical_ldmc_trial <- function(data) {
  required <- c(
    "transform", "horizon", "accepted_name", "experiment_id", "region",
    "pool", "category", "analysis_category", "destination_T_JJA",
    "destination_P_JJA", "diff_T_JJA", "LDMC"
  )
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(
      "Hierarchical LDMC data are missing columns: ",
      paste(missing, collapse = ", "), call. = FALSE
    )
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
    stop("No raw, maximum-duration records are available.", call. = FALSE)
  }

  species_traits <- review_input |>
    dplyr::transmute(
      species = as.character(.data$accepted_name),
      LDMC = as.numeric(.data$LDMC)
    ) |>
    dplyr::distinct()
  duplicate_traits <- species_traits |>
    dplyr::count(.data$species) |>
    dplyr::filter(.data$n != 1L)
  if (nrow(duplicate_traits)) {
    stop("Every species must have exactly one project LDMC value.", call. = FALSE)
  }
  valid_trait_values <- species_traits |>
    dplyr::filter(is.finite(.data$LDMC), .data$LDMC > 0) |>
    dplyr::mutate(log_LDMC = log(.data$LDMC))
  trait_center <- mean(valid_trait_values$log_LDMC)
  trait_scale <- stats::sd(valid_trait_values$log_LDMC)
  if (!is.finite(trait_scale) || trait_scale <= 0) {
    stop("LDMC must have a positive finite distinct-species scale.", call. = FALSE)
  }
  species_traits <- species_traits |>
    dplyr::mutate(
      log_LDMC = dplyr::if_else(
        is.finite(.data$LDMC) & .data$LDMC > 0,
        log(.data$LDMC), NA_real_
      ),
      LDMC_z = (.data$log_LDMC - .env$trait_center) / .env$trait_scale
    )

  experiment_gradients <- review_input |>
    dplyr::transmute(
      experiment = as.character(.data$experiment_id),
      gradient = trimws(as.character(.data$region))
    ) |>
    dplyr::distinct()
  duplicate_gradients <- experiment_gradients |>
    dplyr::count(.data$experiment) |>
    dplyr::filter(.data$n != 1L)
  if (nrow(duplicate_gradients)) {
    stop("Every experiment must map to exactly one gradient/region.", call. = FALSE)
  }
  if (anyNA(experiment_gradients$gradient) ||
      any(experiment_gradients$gradient == "")) {
    stop("Every experiment needs a non-empty gradient/region.", call. = FALSE)
  }

  pool_data <- prepare_pool_climate_response_data(review_input)
  joined <- pool_data$data |>
    dplyr::transmute(
      experiment = as.character(.data$experiment_id),
      species = as.character(.data$accepted_species),
      pool = as.character(.data$pool),
      category = as.character(.data$category),
      analysis_category = as.character(.data$analysis_category),
      temp_z = as.numeric(.data$z_destination_T),
      precip_z = as.numeric(.data$z_log_destination_P),
      warming_z = as.numeric(.data$z_warming)
    ) |>
    dplyr::left_join(species_traits, by = "species") |>
    dplyr::left_join(experiment_gradients, by = "experiment")
  if (anyNA(joined$gradient)) {
    stop("Gradient join left unmatched experiments.", call. = FALSE)
  }
  duplicate_keys <- joined |>
    dplyr::count(.data$species, .data$experiment) |>
    dplyr::filter(.data$n != 1L)
  if (nrow(duplicate_keys)) {
    stop("Trial data must contain one row per species and experiment.", call. = FALSE)
  }

  joined <- joined |>
    dplyr::mutate(
      species = factor(.data$species),
      experiment = factor(.data$experiment),
      gradient = factor(.data$gradient)
    )
  resident_all <- joined |>
    dplyr::filter(.data$pool == "resident") |>
    dplyr::mutate(
      response = factor(
        dplyr::case_when(
          .data$category == "loser" ~ "loser",
          .data$category == "winner" ~ "winner",
          TRUE ~ "non_significant"
        ),
        levels = c("non_significant", "loser", "winner")
      )
    )
  nonresident_all <- joined |>
    dplyr::filter(.data$pool %in% c("colonizer", "outsider")) |>
    dplyr::mutate(
      response = factor(
        dplyr::if_else(
          .data$analysis_category == "colonizer",
          "colonizer", "non_significant"
        ),
        levels = c("non_significant", "colonizer")
      )
    )

  model_columns <- c(
    "response", "temp_z", "precip_z", "warming_z", "LDMC_z",
    "species", "experiment", "gradient"
  )
  resident <- resident_all |>
    dplyr::filter(stats::complete.cases(dplyr::pick(dplyr::all_of(model_columns)))) |>
    droplevels()
  nonresident <- nonresident_all |>
    dplyr::filter(stats::complete.cases(dplyr::pick(dplyr::all_of(model_columns)))) |>
    droplevels()

  if (!identical(levels(resident$response), c("non_significant", "loser", "winner"))) {
    stop("Resident response reference/order is invalid.", call. = FALSE)
  }
  if (!identical(levels(nonresident$response), c("non_significant", "colonizer"))) {
    stop("Nonresident response reference/order is invalid.", call. = FALSE)
  }

  replication <- dplyr::bind_rows(
    resident |>
      dplyr::summarise(
        n_records = dplyr::n(),
        n_experiments = dplyr::n_distinct(.data$experiment),
        within_temp_sd = stats::sd(.data$temp_z),
        within_precip_sd = stats::sd(.data$precip_z),
        within_warming_sd = stats::sd(.data$warming_z),
        .by = "species"
      ) |>
      dplyr::mutate(pool = "resident"),
    nonresident |>
      dplyr::summarise(
        n_records = dplyr::n(),
        n_experiments = dplyr::n_distinct(.data$experiment),
        within_temp_sd = stats::sd(.data$temp_z),
        within_precip_sd = stats::sd(.data$precip_z),
        within_warming_sd = stats::sd(.data$warming_z),
        .by = "species"
      ) |>
      dplyr::mutate(pool = "nonresident")
  ) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::starts_with("within_"),
        ~ dplyr::coalesce(.x, 0)
      )
    )

  counts <- dplyr::bind_rows(
    resident |>
      dplyr::count(outcome = as.character(.data$response), name = "n") |>
      dplyr::mutate(pool = "resident", denominator = nrow(resident)),
    nonresident |>
      dplyr::count(outcome = as.character(.data$response), name = "n") |>
      dplyr::mutate(pool = "nonresident", denominator = nrow(nonresident))
  ) |>
    dplyr::mutate(proportion = .data$n / .data$denominator)

  snapshot <- tibble::tibble(
    n_input_species_experiment = nrow(joined),
    n_resident_complete = nrow(resident),
    n_nonresident_complete = nrow(nonresident),
    n_species_resident = dplyr::n_distinct(resident$species),
    n_species_nonresident = dplyr::n_distinct(nonresident$species),
    n_experiments_resident = dplyr::n_distinct(resident$experiment),
    n_experiments_nonresident = dplyr::n_distinct(nonresident$experiment),
    n_gradients_resident = dplyr::n_distinct(resident$gradient),
    n_gradients_nonresident = dplyr::n_distinct(nonresident$gradient),
    n_loser = sum(resident$response == "loser"),
    n_winner = sum(resident$response == "winner"),
    n_colonizer = sum(nonresident$response == "colonizer"),
    resident_species_single_record = sum(
      replication$pool == "resident" & replication$n_records == 1L
    ),
    nonresident_species_single_record = sum(
      replication$pool == "nonresident" & replication$n_records == 1L
    ),
    resident_species_without_within_temp = sum(
      replication$pool == "resident" & replication$within_temp_sd == 0
    ),
    resident_species_without_within_precip = sum(
      replication$pool == "resident" & replication$within_precip_sd == 0
    ),
    resident_species_without_within_warming = sum(
      replication$pool == "resident" & replication$within_warming_sd == 0
    ),
    nonresident_species_without_within_temp = sum(
      replication$pool == "nonresident" & replication$within_temp_sd == 0
    ),
    nonresident_species_without_within_precip = sum(
      replication$pool == "nonresident" & replication$within_precip_sd == 0
    ),
    nonresident_species_without_within_warming = sum(
      replication$pool == "nonresident" & replication$within_warming_sd == 0
    )
  )

  list(
    resident = resident,
    nonresident = nonresident,
    climate_scaling = pool_data$scaling,
    trait_scaling = tibble::tibble(
      trait = "LDMC", transformation = "natural log",
      center = trait_center, scale = trait_scale,
      weighting = "distinct species",
      n_species = nrow(valid_trait_values)
    ),
    experiment_gradient_map = experiment_gradients,
    counts = counts,
    replication = replication,
    snapshot = snapshot
  )
}

capture_hierarchical_trial_fit <- function(expression) {
  started <- Sys.time()
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
  list(
    fit = fit, warnings = unique(warnings), error = error,
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    timed_out = FALSE
  )
}

fit_hierarchical_mblogit_subprocess <- function(data, timeout = 180) {
  started <- Sys.time()
  result <- tryCatch(
    callr::r(
      function(data) {
        warnings <- character()
        fit <- tryCatch(
          withCallingHandlers(
            mclogit::mblogit(
              response ~ temp_z + precip_z + warming_z +
                temp_z:LDMC_z + precip_z:LDMC_z + warming_z:LDMC_z,
              random = list(
                ~ 1 | species,
                ~ 0 + temp_z | species,
                ~ 0 + precip_z | species,
                ~ 0 + warming_z | species,
                ~ 1 | experiment,
                ~ 1 | gradient
              ),
              data = data,
              method = "PQL", estimator = "ML", catCov = "diagonal",
              trace = FALSE,
              control = mclogit::mmclogit.control(
                maxit = 100, maxit.inner = 1000,
                trace = FALSE, trace.inner = FALSE
              )
            ),
            warning = function(condition) {
              warnings <<- c(warnings, conditionMessage(condition))
              invokeRestart("muffleWarning")
            }
          ),
          error = function(condition) condition
        )
        if (inherits(fit, "condition")) {
          list(
            fit = NULL, warnings = unique(warnings),
            error = conditionMessage(fit)
          )
        } else {
          list(fit = fit, warnings = unique(warnings), error = "")
        }
      },
      args = list(data = data), timeout = timeout,
      spinner = FALSE, show = FALSE
    ),
    error = function(condition) {
      list(
        fit = NULL, warnings = character(),
        error = conditionMessage(condition)
      )
    }
  )
  result$elapsed_seconds <- as.numeric(
    difftime(Sys.time(), started, units = "secs")
  )
  result$timed_out <- grepl(
    "timed out|timeout", result$error, ignore.case = TRUE
  )
  result
}

flatten_mblogit_random_covariance <- function(fit) {
  if (is.null(fit) || is.null(fit$VarCov)) return(tibble::tibble())
  purrr::imap_dfr(fit$VarCov, function(component, component_name) {
    matrix <- as.matrix(component)
    tibble::as_tibble(as.table(matrix), .name_repair = "minimal") |>
      rlang::set_names(c("term_1", "term_2", "covariance")) |>
      dplyr::mutate(
        component = as.character(component_name),
        variance = .data$term_1 == .data$term_2,
        .before = 1
      )
  })
}

hierarchical_trial_fit_diagnostic <- function(captured, model_name, engine, data) {
  fit <- captured$fit
  base <- tibble::tibble(
    model = model_name,
    engine = engine,
    n = nrow(data),
    n_species = dplyr::n_distinct(data$species),
    n_experiments = dplyr::n_distinct(data$experiment),
    n_gradients = dplyr::n_distinct(data$gradient),
    elapsed_seconds = captured$elapsed_seconds,
    timed_out = captured$timed_out,
    warning = paste(captured$warnings, collapse = " | "),
    error = captured$error
  )
  if (is.null(fit)) {
    return(base |>
      dplyr::mutate(
        fitted = FALSE, converged = FALSE, singular = NA,
        finite_coefficients = FALSE, finite_covariance = FALSE,
        minimum_random_variance = NA_real_, clean_hierarchical_fit = FALSE
      ))
  }
  covariance <- tryCatch(as.matrix(stats::vcov(fit)), error = function(e) NULL)
  coefficient_values <- if (identical(engine, "mblogit_PQL_diagonal")) {
    stats::coef(fit)
  } else {
    lme4::fixef(fit)
  }
  finite_coefficients <- all(is.finite(coefficient_values))
  finite_covariance <- !is.null(covariance) && all(is.finite(covariance))

  if (identical(engine, "mblogit_PQL_diagonal")) {
    random <- flatten_mblogit_random_covariance(fit)
    variances <- random$covariance[random$variance]
    minimum_random_variance <- if (length(variances)) min(variances) else NA_real_
    random_eigenvalues <- tryCatch(
      unlist(lapply(fit$VarCov, function(x) {
        eigen(as.matrix(x), symmetric = TRUE, only.values = TRUE)$values
      })),
      error = function(e) NA_real_
    )
    singular <- !length(random_eigenvalues) ||
      any(!is.finite(random_eigenvalues)) || min(random_eigenvalues) <= 1e-6
    converged <- isTRUE(fit$converged)
  } else {
    random <- as.data.frame(lme4::VarCorr(fit))
    variances <- random$vcov[is.na(random$var2) | random$var1 == random$var2]
    minimum_random_variance <- if (length(variances)) min(variances) else NA_real_
    messages <- unlist(fit@optinfo$conv$lme4$messages, use.names = FALSE)
    messages <- messages[!is.na(messages)]
    if (length(messages)) {
      captured$warnings <- unique(c(captured$warnings, messages))
    }
    converged <- isTRUE(fit@optinfo$conv$opt == 0) && !length(messages)
    singular <- lme4::isSingular(fit, tol = 1e-6)
  }
  warning_text <- paste(captured$warnings, collapse = " | ")
  base$warning <- warning_text
  base |>
    dplyr::mutate(
      fitted = TRUE,
      converged = .env$converged,
      singular = .env$singular,
      finite_coefficients = .env$finite_coefficients,
      finite_covariance = .env$finite_covariance,
      minimum_random_variance = .env$minimum_random_variance,
      clean_hierarchical_fit =
        .env$converged && !.env$singular && .env$finite_coefficients &&
        .env$finite_covariance && !nzchar(.data$warning) && !nzchar(.data$error)
    )
}

fit_hierarchical_ldmc_trial <- function(prepared, multinomial_timeout = 180) {
  resident_multinomial <- fit_hierarchical_mblogit_subprocess(
    prepared$resident, timeout = multinomial_timeout
  )
  multinomial_diagnostic <- hierarchical_trial_fit_diagnostic(
    resident_multinomial, "resident_multinomial", "mblogit_PQL_diagonal",
    prepared$resident
  )

  binary_inputs <- list()
  binary_fits <- list()
  binary_diagnostics <- list()
  if (!isTRUE(multinomial_diagnostic$clean_hierarchical_fit[[1]])) {
    binary_inputs$loser <- prepared$resident |>
      dplyr::filter(.data$response != "winner") |>
      dplyr::mutate(binary_response = .data$response == "loser") |>
      droplevels()
    binary_inputs$winner <- prepared$resident |>
      dplyr::filter(.data$response != "loser") |>
      dplyr::mutate(binary_response = .data$response == "winner") |>
      droplevels()
  }
  binary_inputs$colonizer <- prepared$nonresident |>
    dplyr::mutate(binary_response = .data$response == "colonizer") |>
    droplevels()

  for (model_name in names(binary_inputs)) {
    captured <- capture_hierarchical_trial_fit(
      lme4::glmer(
        hierarchical_ldmc_binomial_formula("binary_response"),
        data = binary_inputs[[model_name]], family = stats::binomial(),
        control = lme4::glmerControl(
          optimizer = "bobyqa", optCtrl = list(maxfun = 200000)
        )
      )
    )
    binary_fits[[model_name]] <- captured
    binary_diagnostics[[model_name]] <- hierarchical_trial_fit_diagnostic(
      captured, paste0(model_name, "_binomial"), "glmer_binomial",
      binary_inputs[[model_name]]
    )
  }

  list(
    resident_multinomial = resident_multinomial,
    binary = binary_fits,
    binary_data = binary_inputs,
    diagnostics = dplyr::bind_rows(
      multinomial_diagnostic,
      dplyr::bind_rows(binary_diagnostics)
    )
  )
}

summarize_hierarchical_trial_event_support <- function(prepared) {
  inputs <- list(
    loser = prepared$resident |>
      dplyr::filter(.data$response != "winner") |>
      dplyr::mutate(event = .data$response == "loser") |>
      droplevels(),
    winner = prepared$resident |>
      dplyr::filter(.data$response != "loser") |>
      dplyr::mutate(event = .data$response == "winner") |>
      droplevels(),
    colonizer = prepared$nonresident |>
      dplyr::mutate(event = .data$response == "colonizer") |>
      droplevels()
  )
  purrr::imap_dfr(inputs, function(data, outcome) {
    purrr::map_dfr(c("species", "experiment", "gradient"), function(grouping) {
      support <- data |>
        dplyr::summarise(
          n_records = dplyr::n(),
          n_events = sum(.data$event),
          has_both = any(.data$event) && any(!.data$event),
          all_event = all(.data$event),
          no_event = !any(.data$event),
          .by = dplyr::all_of(grouping)
        )
      tibble::tibble(
        outcome = outcome,
        grouping = grouping,
        n_groups = nrow(support),
        groups_with_event = sum(support$n_events > 0),
        groups_with_both_outcomes = sum(support$has_both),
        groups_all_event = sum(support$all_event),
        groups_no_event = sum(support$no_event),
        groups_single_record = sum(support$n_records == 1L),
        median_records = stats::median(support$n_records),
        maximum_records = max(support$n_records)
      )
    })
  })
}

tidy_hierarchical_trial_coefficients <- function(fits) {
  tidy_one <- function(fit, model_name, engine) {
    if (is.null(fit)) return(tibble::tibble())
    estimates <- if (identical(engine, "mblogit_PQL_diagonal")) {
      stats::coef(fit)
    } else {
      lme4::fixef(fit)
    }
    covariance <- tryCatch(as.matrix(stats::vcov(fit)), error = function(e) NULL)
    if (is.null(covariance) || any(!is.finite(covariance))) {
      return(tibble::tibble())
    }
    standard_error <- sqrt(diag(covariance))
    tibble::tibble(
      model = model_name,
      engine = engine,
      coefficient = names(estimates),
      estimate = as.numeric(estimates),
      standard_error = as.numeric(standard_error)
    ) |>
      dplyr::mutate(
      odds_ratio = exp(.data$estimate),
      conf_low = exp(.data$estimate - stats::qnorm(0.975) * .data$standard_error),
      conf_high = exp(.data$estimate + stats::qnorm(0.975) * .data$standard_error),
      z = .data$estimate / .data$standard_error,
      p_value = 2 * stats::pnorm(abs(.data$z), lower.tail = FALSE)
      )
  }
  output <- list(
    tidy_one(
      fits$resident_multinomial$fit,
      "resident_multinomial", "mblogit_PQL_diagonal"
    )
  )
  for (model_name in names(fits$binary)) {
    output[[length(output) + 1L]] <- tidy_one(
      fits$binary[[model_name]]$fit,
      paste0(model_name, "_binomial"), "glmer_binomial"
    )
  }
  dplyr::bind_rows(output) |>
    dplyr::left_join(
      fits$diagnostics |>
        dplyr::select(
          "model", usable_for_inference = "clean_hierarchical_fit"
        ),
      by = "model"
    )
}

collect_hierarchical_trial_random_effects <- function(fits) {
  output <- list()
  multinomial <- fits$resident_multinomial$fit
  if (!is.null(multinomial)) {
    output[["resident_multinomial"]] <-
      flatten_mblogit_random_covariance(multinomial) |>
      dplyr::mutate(
        model = "resident_multinomial", engine = "mblogit_PQL_diagonal",
        group = .data$component,
        term = paste(.data$term_1, .data$term_2, sep = " :: "),
        value = .data$covariance
      ) |>
      dplyr::select("model", "engine", "group", "term", "value", "variance")
  }
  for (model_name in names(fits$binary)) {
    fit <- fits$binary[[model_name]]$fit
    if (is.null(fit)) next
    output[[model_name]] <- as.data.frame(lme4::VarCorr(fit)) |>
      dplyr::transmute(
        model = paste0(.env$model_name, "_binomial"),
        engine = "glmer_binomial",
        group = .data$grp,
        term = dplyr::if_else(
          is.na(.data$var2), as.character(.data$var1),
          paste(.data$var1, .data$var2, sep = " :: ")
        ),
        value = .data$vcov,
        variance = is.na(.data$var2) | .data$var1 == .data$var2
      )
  }
  dplyr::bind_rows(output)
}

validate_hierarchical_ldmc_trial <- function(prepared, fits) {
  checks <- c(
    resident_reference = identical(
      levels(prepared$resident$response),
      c("non_significant", "loser", "winner")
    ),
    nonresident_reference = identical(
      levels(prepared$nonresident$response),
      c("non_significant", "colonizer")
    ),
    experiment_maps_once = all(
      dplyr::count(prepared$experiment_gradient_map, .data$experiment)$n == 1L
    ),
    gradient_is_region = dplyr::n_distinct(
      prepared$experiment_gradient_map$gradient
    ) == 20L,
    no_duplicate_resident_keys = all(
      dplyr::count(prepared$resident, .data$species, .data$experiment)$n == 1L
    ),
    no_duplicate_nonresident_keys = all(
      dplyr::count(prepared$nonresident, .data$species, .data$experiment)$n == 1L
    ),
    exact_fixed_formula = identical(
      attr(stats::terms(hierarchical_ldmc_fixed_formula()), "term.labels"),
      c(
        "temp_z", "precip_z", "warming_z", "temp_z:LDMC_z",
        "precip_z:LDMC_z", "warming_z:LDMC_z"
      )
    ),
    exact_random_component_count =
      length(hierarchical_ldmc_random_formulas()) == 6L,
    colonizer_binomial_attempted = "colonizer" %in% names(fits$binary),
    resident_binomial_rule = if (
      isTRUE(fits$diagnostics$clean_hierarchical_fit[
        fits$diagnostics$model == "resident_multinomial"
      ])
    ) {
      !any(c("loser", "winner") %in% names(fits$binary))
    } else {
      all(c("loser", "winner") %in% names(fits$binary))
    }
  )
  tibble::tibble(check = names(checks), passed = unname(checks))
}

hierarchical_ldmc_binomial_inputs <- function(prepared) {
  list(
    loser = prepared$resident |>
      dplyr::filter(.data$response != "winner") |>
      dplyr::mutate(binary_response = .data$response == "loser") |>
      droplevels(),
    winner = prepared$resident |>
      dplyr::filter(.data$response != "loser") |>
      dplyr::mutate(binary_response = .data$response == "winner") |>
      droplevels(),
    colonizer = prepared$nonresident |>
      dplyr::mutate(binary_response = .data$response == "colonizer") |>
      droplevels()
  )
}

hierarchical_ldmc_binomial_reduction_specs <- function() {
  climate_main <- c("temp_z", "precip_z", "warming_z")
  interactions <- c(
    "temp_z:LDMC_z", "precip_z:LDMC_z", "warming_z:LDMC_z"
  )
  species_intercept <- "(1 | species)"
  all_species_slopes <- paste(
    species_intercept,
    "+ (0 + temp_z | species) + (0 + precip_z | species) +",
    "(0 + warming_z | species)"
  )
  make_spec <- function(name, fixed, random, reduction, random_structure) {
    list(
      name = name,
      fixed_terms = fixed,
      random_terms = random,
      reduction = reduction,
      random_structure = random_structure
    )
  }
  list(
    make_spec(
      "full_species_all_slopes", c(climate_main, interactions),
      all_species_slopes, "none", "species intercept + three slopes"
    ),
    make_spec(
      "full_species_intercept", c(climate_main, interactions),
      species_intercept, "remove all species slopes",
      "species intercept"
    ),
    make_spec(
      "drop_temp_interaction", c(climate_main, interactions[-1]),
      species_intercept, "remove temperature x LDMC",
      "species intercept"
    ),
    make_spec(
      "drop_precip_interaction", c(climate_main, interactions[-2]),
      species_intercept, "remove precipitation x LDMC",
      "species intercept"
    ),
    make_spec(
      "drop_warming_interaction", c(climate_main, interactions[-3]),
      species_intercept, "remove warming x LDMC",
      "species intercept"
    ),
    make_spec(
      "climate_main_only", climate_main,
      species_intercept, "remove all LDMC interactions",
      "species intercept"
    ),
    make_spec(
      "temperature_axis", c("temp_z", "temp_z:LDMC_z"),
      paste(species_intercept, "+ (0 + temp_z | species)"),
      "retain temperature axis only",
      "species intercept + temperature slope"
    ),
    make_spec(
      "precipitation_axis", c("precip_z", "precip_z:LDMC_z"),
      paste(species_intercept, "+ (0 + precip_z | species)"),
      "retain precipitation axis only",
      "species intercept + precipitation slope"
    ),
    make_spec(
      "warming_axis", c("warming_z", "warming_z:LDMC_z"),
      paste(species_intercept, "+ (0 + warming_z | species)"),
      "retain warming axis only",
      "species intercept + warming slope"
    )
  )
}

hierarchical_ldmc_reduction_formula <- function(spec) {
  stats::as.formula(paste(
    "binary_response ~", paste(spec$fixed_terms, collapse = " + "),
    "+", spec$random_terms
  ))
}

fit_hierarchical_ldmc_binomial_reductions <- function(prepared) {
  inputs <- hierarchical_ldmc_binomial_inputs(prepared)
  specs <- hierarchical_ldmc_binomial_reduction_specs()
  models <- list()
  diagnostics <- list()

  for (outcome in names(inputs)) {
    models[[outcome]] <- list()
    for (spec in specs) {
      formula <- hierarchical_ldmc_reduction_formula(spec)
      captured <- capture_hierarchical_trial_fit(
        suppressMessages(
          lme4::glmer(
            formula, data = inputs[[outcome]], family = stats::binomial(),
            control = lme4::glmerControl(
              optimizer = "bobyqa", optCtrl = list(maxfun = 200000)
            )
          )
        )
      )
      models[[outcome]][[spec$name]] <- captured
      diagnostic <- hierarchical_trial_fit_diagnostic(
        captured, paste(outcome, spec$name, sep = "__"),
        "glmer_binomial", inputs[[outcome]]
      ) |>
        dplyr::mutate(
          outcome = .env$outcome,
          specification = spec$name,
          reduction = spec$reduction,
          random_structure = spec$random_structure,
          formula = paste(deparse(.env$formula), collapse = " "),
          event_count = sum(inputs[[outcome]]$binary_response),
          non_event_count = sum(!inputs[[outcome]]$binary_response),
          aic = if (is.null(captured$fit)) NA_real_ else stats::AIC(captured$fit),
          bic = if (is.null(captured$fit)) NA_real_ else stats::BIC(captured$fit),
          log_likelihood = if (is.null(captured$fit)) {
            NA_real_
          } else {
            as.numeric(stats::logLik(captured$fit))
          },
          .before = 1
        )
      diagnostics[[paste(outcome, spec$name, sep = "__")]] <- diagnostic
    }
  }
  list(
    inputs = inputs,
    specs = specs,
    models = models,
    diagnostics = dplyr::bind_rows(diagnostics)
  )
}

tidy_hierarchical_ldmc_reduction_coefficients <- function(fits) {
  output <- list()
  for (outcome in names(fits$models)) {
    for (specification in names(fits$models[[outcome]])) {
      fit <- fits$models[[outcome]][[specification]]$fit
      if (is.null(fit)) next
      estimates <- lme4::fixef(fit)
      covariance <- tryCatch(
        as.matrix(stats::vcov(fit)), error = function(condition) NULL
      )
      if (is.null(covariance) || any(!is.finite(covariance))) next
      standard_error <- sqrt(diag(covariance))
      output[[paste(outcome, specification, sep = "__")]] <- tibble::tibble(
        outcome = outcome,
        specification = specification,
        coefficient = names(estimates),
        estimate = as.numeric(estimates),
        standard_error = as.numeric(standard_error)
      ) |>
        dplyr::mutate(
          odds_ratio = exp(.data$estimate),
          conf_low = exp(
            .data$estimate - stats::qnorm(0.975) * .data$standard_error
          ),
          conf_high = exp(
            .data$estimate + stats::qnorm(0.975) * .data$standard_error
          ),
          z = .data$estimate / .data$standard_error,
          p_value = 2 * stats::pnorm(abs(.data$z), lower.tail = FALSE)
        )
    }
  }
  dplyr::bind_rows(output) |>
    dplyr::left_join(
      fits$diagnostics |>
        dplyr::select(
          "outcome", "specification",
          clean_numerical_fit = "clean_hierarchical_fit"
        ),
      by = c("outcome", "specification")
    )
}

collect_hierarchical_ldmc_reduction_random_effects <- function(fits) {
  output <- list()
  for (outcome in names(fits$models)) {
    for (specification in names(fits$models[[outcome]])) {
      fit <- fits$models[[outcome]][[specification]]$fit
      if (is.null(fit)) next
      output[[paste(outcome, specification, sep = "__")]] <-
        as.data.frame(lme4::VarCorr(fit)) |>
        dplyr::transmute(
          outcome = .env$outcome,
          specification = .env$specification,
          group = .data$grp,
          term = dplyr::if_else(
            is.na(.data$var2), as.character(.data$var1),
            paste(.data$var1, .data$var2, sep = " :: ")
          ),
          variance = .data$vcov,
          standard_deviation = .data$sdcor
        )
    }
  }
  dplyr::bind_rows(output)
}

hierarchical_ldmc_reduction_plan <- function() {
  purrr::map_dfr(hierarchical_ldmc_binomial_reduction_specs(), function(spec) {
    formula <- hierarchical_ldmc_reduction_formula(spec)
    tibble::tibble(
      specification = spec$name,
      reduction = spec$reduction,
      random_structure = spec$random_structure,
      formula = paste(deparse(formula), collapse = " "),
      includes_experiment_random = grepl("experiment", deparse(formula)),
      includes_gradient_random = grepl("gradient", deparse(formula))
    )
  })
}

validate_hierarchical_ldmc_binomial_reductions <- function(fits) {
  plan <- hierarchical_ldmc_reduction_plan()
  checks <- c(
    three_outcomes = identical(
      names(fits$models), c("loser", "winner", "colonizer")
    ),
    nine_specs_per_outcome = all(
      lengths(fits$models) == nrow(plan) & lengths(fits$models) == 9L
    ),
    no_experiment_random = !any(plan$includes_experiment_random),
    no_gradient_random = !any(plan$includes_gradient_random),
    all_models_reported = nrow(fits$diagnostics) == 27L,
    all_models_have_species_random = all(
      grepl("species", plan$formula, fixed = TRUE)
    ),
    all_interactions_retain_climate_main = all(
      !grepl(":LDMC_z", plan$formula, fixed = TRUE) |
        grepl("_z", plan$formula, fixed = TRUE)
    )
  )
  tibble::tibble(check = names(checks), passed = unname(checks))
}

plot_hierarchical_ldmc_resident_odds_ratios <- function(
    coefficients, random_effects, path) {
  required_coefficients <- c(
    "outcome", "specification", "coefficient", "odds_ratio",
    "conf_low", "conf_high", "p_value", "clean_numerical_fit"
  )
  missing <- setdiff(required_coefficients, names(coefficients))
  if (length(missing)) {
    stop(
      "Resident odds-ratio data are missing: ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  plot_data <- coefficients |>
    dplyr::filter(
      .data$outcome %in% c("loser", "winner"),
      .data$specification == "full_species_intercept",
      .data$coefficient != "(Intercept)",
      .data$clean_numerical_fit
    )
  expected_terms <- c(
    "temp_z", "precip_z", "warming_z", "temp_z:LDMC_z",
    "precip_z:LDMC_z", "warming_z:LDMC_z"
  )
  if (nrow(plot_data) != 12L ||
      !setequal(unique(plot_data$coefficient), expected_terms)) {
    stop(
      "The plot requires six clean coefficients for each resident outcome.",
      call. = FALSE
    )
  }
  species_sd <- random_effects |>
    dplyr::filter(
      .data$outcome %in% c("loser", "winner"),
      .data$specification == "full_species_intercept",
      .data$group == "species", .data$term == "(Intercept)"
    ) |>
    dplyr::select("outcome", species_sd = "standard_deviation")
  if (nrow(species_sd) != 2L) {
    stop("Both resident species-intercept SDs are required.", call. = FALSE)
  }
  term_labels <- c(
    temp_z = "Destination JJA temperature",
    precip_z = "Destination JJA precipitation",
    warming_z = "Experimental JJA warming",
    `temp_z:LDMC_z` = "Temperature × LDMC",
    `precip_z:LDMC_z` = "Precipitation × LDMC",
    `warming_z:LDMC_z` = "Warming × LDMC"
  )
  term_order <- rev(unname(term_labels[expected_terms]))
  colours <- c(loser = "#E66101", winner = "#1B9E77")
  plot_data <- plot_data |>
    dplyr::left_join(species_sd, by = "outcome") |>
    dplyr::mutate(
      term_label = factor(
        unname(.env$term_labels[.data$coefficient]),
        levels = .env$term_order
      ),
      outcome_label = dplyr::if_else(
        .data$outcome == "loser",
        sprintf(
          "Loser vs non-significant resident\nSpecies intercept SD = %.2f",
          .data$species_sd
        ),
        sprintf(
          "Winner vs non-significant resident\nSpecies intercept SD = %.2f",
          .data$species_sd
        )
      ),
      p_label = dplyr::if_else(
        .data$p_value < 0.001,
        "p < 0.001", sprintf("p = %.3f", .data$p_value)
      ),
      result_label = sprintf(
        "%.2f [%.2f, %.2f]\n%s",
        .data$odds_ratio, .data$conf_low, .data$conf_high, .data$p_label
      )
    )
  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data$odds_ratio, y = .data$term_label,
      colour = .data$outcome
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 1, colour = "grey55", linewidth = 0.7,
      linetype = "dashed"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$conf_low, xmax = .data$conf_high),
      orientation = "y", width = 0, linewidth = 1
    ) +
    ggplot2::geom_point(size = 4.2) +
    ggplot2::geom_text(
      ggplot2::aes(x = 3.05, label = .data$result_label),
      colour = "grey15", hjust = 0, lineheight = 0.95, size = 4.2
    ) +
    ggplot2::facet_wrap(~ outcome_label, ncol = 2) +
    ggplot2::scale_colour_manual(values = colours, guide = "none") +
    ggplot2::scale_x_log10(
      breaks = c(0.5, 0.75, 1, 1.5, 2, 3),
      labels = c("0.5", "0.75", "1", "1.5", "2", "3"),
      limits = NULL
    ) +
    ggplot2::coord_cartesian(xlim = c(0.4, 5.4)) +
    transplant_figure_theme(base_size = 16) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(
        colour = "grey90", linewidth = 0.35
      ),
      strip.text = ggplot2::element_text(size = 16, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 14),
      panel.spacing = grid::unit(0.8, "cm"),
      plot.margin = ggplot2::margin(12, 18, 10, 10)
    ) +
    ggplot2::labs(
      x = "Odds ratio (95% confidence interval)", y = NULL
    )
  save_plot_atomic(plot, path, width = 14, height = 7.5)
}

extract_hierarchical_ldmc_repeated_species_effects <- function(
    models, prepared, signed_species) {
  required_signed <- c(
    "species", "consistency_group", "species_position", "loser", "winner",
    "colonizer", "significant_total", "neutral_total"
  )
  missing_signed <- setdiff(required_signed, names(signed_species))
  if (length(missing_signed)) {
    stop(
      "Repeated-species data are missing: ",
      paste(missing_signed, collapse = ", "), call. = FALSE
    )
  }
  if (anyDuplicated(signed_species$species)) {
    stop("Repeated-species names must be unique.", call. = FALSE)
  }

  binary_inputs <- hierarchical_ldmc_binomial_inputs(prepared)
  effects <- purrr::map_dfr(c("loser", "winner"), function(outcome) {
    captured <- models[[outcome]][["full_species_intercept"]]
    fit <- captured$fit
    if (is.null(fit) || !inherits(fit, "glmerMod")) {
      stop(
        "The clean full-fixed species-intercept model is unavailable for ",
        outcome, ".", call. = FALSE
      )
    }
    random <- lme4::ranef(fit, condVar = TRUE)$species
    conditional_variance <- attr(random, "postVar")
    if (is.null(conditional_variance) ||
        length(dim(conditional_variance)) != 3L ||
        dim(conditional_variance)[[3]] != nrow(random)) {
      stop(
        "Conditional species variances are unavailable for ", outcome, ".",
        call. = FALSE
      )
    }
    conditional_sd <- sqrt(vapply(
      seq_len(nrow(random)),
      function(index) conditional_variance[1, 1, index],
      numeric(1)
    ))
    model_counts <- binary_inputs[[outcome]] |>
      dplyr::mutate(species = as.character(.data$species)) |>
      dplyr::summarise(
        model_records = dplyr::n(),
        event_records = sum(.data$binary_response),
        non_event_records = sum(!.data$binary_response),
        model_experiments = dplyr::n_distinct(.data$experiment),
        .by = "species"
      )
    tibble::tibble(
      species = rownames(random),
      outcome = outcome,
      conditional_mode_log_odds = as.numeric(random[["(Intercept)"]]),
      conditional_sd_log_odds = conditional_sd
    ) |>
      dplyr::mutate(
        conditional_low_log_odds =
          .data$conditional_mode_log_odds - 1.96 * .data$conditional_sd_log_odds,
        conditional_high_log_odds =
          .data$conditional_mode_log_odds + 1.96 * .data$conditional_sd_log_odds,
        odds_multiplier = exp(.data$conditional_mode_log_odds),
        odds_multiplier_low = exp(.data$conditional_low_log_odds),
        odds_multiplier_high = exp(.data$conditional_high_log_odds)
      ) |>
      dplyr::left_join(model_counts, by = "species")
  })

  repeated_grid <- tidyr::crossing(
    signed_species |>
      dplyr::select(dplyr::all_of(required_signed)),
    outcome = c("loser", "winner")
  )
  repeated_grid |>
    dplyr::left_join(effects, by = c("species", "outcome")) |>
    dplyr::mutate(
      included_in_model = is.finite(.data$conditional_mode_log_odds),
      outcome = factor(.data$outcome, levels = c("loser", "winner")),
      outcome_label = factor(
        dplyr::recode(
          as.character(.data$outcome),
          loser = "Loser vs non-significant resident",
          winner = "Winner vs non-significant resident"
        ),
        levels = c(
          "Loser vs non-significant resident",
          "Winner vs non-significant resident"
        )
      )
    ) |>
    dplyr::arrange(.data$outcome, .data$species_position)
}

validate_hierarchical_ldmc_repeated_species_effects <- function(
    effects, signed_species) {
  expected_rows <- 2L * nrow(signed_species)
  included <- effects |>
    dplyr::filter(.data$included_in_model)
  checks <- c(
    exact_species_set = setequal(effects$species, signed_species$species),
    one_row_per_species_outcome =
      nrow(effects) == expected_rows &&
      !anyDuplicated(effects[c("species", "outcome")]),
    outcomes_are_loser_and_winner = setequal(
      as.character(unique(effects$outcome)), c("loser", "winner")
    ),
    finite_included_effects = all(is.finite(unlist(included[c(
      "conditional_mode_log_odds", "conditional_sd_log_odds",
      "conditional_low_log_odds", "conditional_high_log_odds",
      "odds_multiplier", "odds_multiplier_low", "odds_multiplier_high"
    )]))),
    positive_conditional_uncertainty = all(
      included$conditional_sd_log_odds > 0 &
        included$conditional_high_log_odds >= included$conditional_low_log_odds
    ),
    model_counts_are_conserved = all(
      included$model_records ==
        included$event_records + included$non_event_records
    )
  )
  tibble::tibble(check = names(checks), passed = unname(checks))
}

plot_hierarchical_ldmc_repeated_species_effects <- function(effects, path) {
  required <- c(
    "species", "consistency_group", "species_position", "outcome",
    "outcome_label", "conditional_mode_log_odds", "conditional_low_log_odds",
    "conditional_high_log_odds", "included_in_model"
  )
  missing <- setdiff(required, names(effects))
  if (length(missing)) {
    stop(
      "Species-effect plot data are missing: ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  species_order <- effects |>
    dplyr::distinct(.data$species, .data$species_position) |>
    dplyr::arrange(.data$species_position) |>
    dplyr::pull(.data$species)
  plot_data <- effects |>
    dplyr::mutate(
      species_label = factor(.data$species, levels = .env$species_order)
    )
  included <- plot_data |>
    dplyr::filter(.data$included_in_model)
  colours <- c(loser = "#E66101", winner = "#1B9E77")

  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$species_label)
  ) +
    ggplot2::geom_blank(ggplot2::aes(y = 0)) +
    ggplot2::geom_hline(
      yintercept = 0, colour = "grey55", linewidth = 0.55,
      linetype = "dashed"
    ) +
    ggplot2::geom_errorbar(
      data = included,
      mapping = ggplot2::aes(
        ymin = .data$conditional_low_log_odds,
        ymax = .data$conditional_high_log_odds,
        colour = .data$outcome
      ),
      width = 0, linewidth = 0.55, alpha = 0.72
    ) +
    ggplot2::geom_point(
      data = included,
      mapping = ggplot2::aes(
        y = .data$conditional_mode_log_odds, colour = .data$outcome
      ),
      size = 2.1
    ) +
    ggplot2::facet_grid(
      outcome_label ~ consistency_group,
      scales = "free_x", space = "free_x", drop = FALSE
    ) +
    ggplot2::scale_colour_manual(values = colours, guide = "none") +
    ggplot2::scale_y_continuous(
      breaks = seq(-4, 6, by = 2),
      expand = ggplot2::expansion(mult = c(0.03, 0.04))
    ) +
    transplant_figure_theme(base_size = 14) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(
        colour = "grey92", linewidth = 0.3
      ),
      strip.text.x = ggplot2::element_text(size = 14, face = "bold"),
      strip.text.y = ggplot2::element_text(size = 12, face = "bold"),
      axis.text.x = ggplot2::element_text(
        size = 8.5, face = "italic", angle = 60, hjust = 1
      ),
      panel.spacing.x = grid::unit(0.35, "cm"),
      panel.spacing.y = grid::unit(0.7, "cm"),
      plot.margin = ggplot2::margin(10, 12, 10, 10)
    ) +
    ggplot2::labs(
      x = NULL, y = "Species-specific intercept deviation (log-odds scale)"
    )
  save_plot_atomic(plot, path, width = 16, height = 10)
}
