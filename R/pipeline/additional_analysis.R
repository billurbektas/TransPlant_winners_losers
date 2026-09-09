tnrs_family_lookup <- function(accepted_names) {
  response <- TNRS::TNRS(unique(accepted_names[!is.na(accepted_names)]))
  family_col <- intersect(c("Accepted_family", "Name_matched_accepted_family", "Family_matched"), names(response))[[1]]
  tibble::tibble(
    accepted_name = response$Accepted_name,
    family = response[[family_col]],
    family_match_score = response$Overall_score
  ) |>
    dplyr::distinct(.data$accepted_name, .keep_all = TRUE)
}

derive_functional_groups <- function(taxonomy_family, overrides) {
  graminoid_families <- c("Poaceae", "Cyperaceae", "Juncaceae")
  fern_families <- c(
    "Aspleniaceae", "Athyriaceae", "Blechnaceae", "Cystopteridaceae", "Dennstaedtiaceae",
    "Dryopteridaceae", "Equisetaceae", "Lycopodiaceae", "Ophioglossaceae", "Polypodiaceae",
    "Pteridaceae", "Selaginellaceae", "Thelypteridaceae", "Woodsiaceae"
  )
  woody_genera <- c(
    "Arctostaphylos", "Berberis", "Betula", "Calluna", "Cassiope",
    "Daphne", "Dasiphora", "Dryas", "Empetrum", "Helianthemum",
    "Hippophae", "Juniperus", "Kalmia", "Kleinia", "Loiseleuria",
    "Phyllodoce", "Rhododendron", "Rosa", "Salix", "Sorbus", "Spiraea",
    "Symphoricarpos", "Thymus", "Vaccinium"
  )
  required_override_columns <- c(
    "accepted_name", "graminoid_status", "woody_status", "note"
  )
  missing_override_columns <- setdiff(required_override_columns, names(overrides))
  if (length(missing_override_columns)) {
    stop(
      "Functional-status overrides are missing: ",
      paste(missing_override_columns, collapse = ", "), call. = FALSE
    )
  }
  invalid_graminoid <- setdiff(
    unique(trimws(stats::na.omit(overrides$graminoid_status))),
    c("", "graminoid", "non_graminoid")
  )
  invalid_woody <- setdiff(
    unique(trimws(stats::na.omit(overrides$woody_status))),
    c("", "woody", "non_woody")
  )
  if (length(invalid_graminoid) || length(invalid_woody)) {
    stop(
      "Invalid functional-status override values. graminoid_status: ",
      paste(invalid_graminoid, collapse = ", "), "; woody_status: ",
      paste(invalid_woody, collapse = ", "), call. = FALSE
    )
  }
  taxonomy_family |>
    dplyr::mutate(genus = stringr::word(.data$accepted_name, 1)) |>
    dplyr::left_join(overrides, by = "accepted_name", suffix = c("", "_override")) |>
    dplyr::mutate(
      graminoid_override_supplied = !is.na(.data$graminoid_status) &
        nzchar(trimws(.data$graminoid_status)),
      woody_override_supplied = !is.na(.data$woody_status) &
        nzchar(trimws(.data$woody_status)),
      graminoid_status_auto = dplyr::case_when(
        is.na(.data$family) | is.na(.data$accepted_name) ~ NA_character_,
        .data$family %in% graminoid_families ~ "graminoid",
        TRUE ~ "non_graminoid"
      ),
      woody_status_auto = dplyr::case_when(
        is.na(.data$family) | is.na(.data$accepted_name) ~ NA_character_,
        .data$genus %in% woody_genera ~ "woody",
        TRUE ~ "non_woody"
      ),
      graminoid_status = dplyr::coalesce(
        .data$graminoid_status, .data$graminoid_status_auto
      ),
      woody_status = dplyr::coalesce(
        .data$woody_status, .data$woody_status_auto
      ),
      graminoid_review_basis = dplyr::case_when(
        .data$graminoid_override_supplied ~
          paste0("override: ", .data$note),
        .data$family %in% graminoid_families ~
          "reviewed graminoid-family rule",
        !is.na(.data$family) ~ "reviewed non-graminoid family",
        TRUE ~ "unresolved family"
      ),
      woody_review_basis = dplyr::case_when(
        .data$woody_override_supplied ~
          paste0("override: ", .data$note),
        .data$genus %in% woody_genera ~ "reviewed woody-genus rule",
        !is.na(.data$genus) ~ "reviewed non-woody genus",
        TRUE ~ "unresolved genus"
      ),
      # Deprecated compatibility label for old descriptive scripts only. New
      # inference must use graminoid_status and woody_status separately.
      functional_group = dplyr::case_when(
        .data$family %in% fern_families ~ "other_unresolved",
        .data$woody_status == "woody" ~ "other_unresolved",
        .data$graminoid_status == "graminoid" ~ "graminoid",
        is.na(.data$graminoid_status) | is.na(.data$woody_status) ~
          "other_unresolved",
        TRUE ~ "forb"
      ),
      review_basis = dplyr::case_when(
        !is.na(.data$note) ~ paste0("override: ", .data$note),
        .data$family %in% graminoid_families ~ "reviewed graminoid-family rule",
        .data$family %in% fern_families ~ "fern-family exclusion",
        .data$genus %in% woody_genera ~ "woody-genus exclusion",
        is.na(.data$family) ~ "unresolved family",
        TRUE ~ "herbaceous non-graminoid candidate"
      )
    ) |>
    dplyr::select(
      accepted_name, family, genus, graminoid_status, woody_status,
      graminoid_review_basis, woody_review_basis,
      functional_group, review_basis
    )
}

join_experiment_climate <- function(data, experiment_climate) {
  duplicated_metadata <- setdiff(
    intersect(names(data), names(experiment_climate)),
    "experiment_id"
  )
  data |>
    dplyr::left_join(
      experiment_climate |>
        dplyr::select(-dplyr::all_of(duplicated_metadata)),
      by = "experiment_id"
    )
}

reviewed_five_climate_contexts <- function(data) {
  temperature_cutoff_low <- 12.86686795491147
  temperature_cutoff_high <- 16.94267310789054
  precipitation_cutoff_dry <- 550
  precipitation_cutoff_wet <- 768
  context_levels <- c(
    "cold_dry", "intermediate_dry", "warm_dry",
    "cold_intermediate", "cold_wet"
  )
  climate <- data |>
    dplyr::distinct(
      .data$experiment_id, .data$destination_T_JJA, .data$destination_P_JJA,
      .data$climate_summary_status
    ) |>
    dplyr::filter(
      is.finite(.data$destination_T_JJA),
      is.finite(.data$destination_P_JJA)
    )
  if (nrow(climate) == 0L) stop("No finite destination climates are available.", call. = FALSE)
  duplicates <- climate |>
    dplyr::count(.data$experiment_id) |>
    dplyr::filter(.data$n != 1L)
  if (nrow(duplicates)) {
    stop("Every experiment must have one unique destination climate.", call. = FALSE)
  }
  classified <- climate |>
    dplyr::mutate(
      climate_context = dplyr::case_when(
        .data$destination_P_JJA < .env$precipitation_cutoff_dry &
          .data$destination_T_JJA < .env$temperature_cutoff_low ~ "cold_dry",
        .data$destination_P_JJA < .env$precipitation_cutoff_dry &
          .data$destination_T_JJA < .env$temperature_cutoff_high ~ "intermediate_dry",
        .data$destination_P_JJA < .env$precipitation_cutoff_dry ~ "warm_dry",
        .data$destination_P_JJA < .env$precipitation_cutoff_wet &
          .data$destination_T_JJA < .env$temperature_cutoff_high ~ "cold_intermediate",
        .data$destination_P_JJA >= .env$precipitation_cutoff_wet &
          .data$destination_T_JJA < .env$temperature_cutoff_high ~ "cold_wet",
        TRUE ~ NA_character_
      ),
      climate_context = factor(.data$climate_context, levels = .env$context_levels),
      climate_context_label = factor(
        stringr::str_to_title(stringr::str_replace(.data$climate_context, "_", "/")),
        levels = c(
          "Cold/Dry", "Intermediate/Dry", "Warm/Dry",
          "Cold/Intermediate", "Cold/Wet"
        )
      ),
      temperature_cutoff_low = .env$temperature_cutoff_low,
      temperature_cutoff_high = .env$temperature_cutoff_high,
      precipitation_cutoff_dry = .env$precipitation_cutoff_dry,
      precipitation_cutoff_wet = .env$precipitation_cutoff_wet,
      context_review_status = dplyr::if_else(
        is.na(.data$climate_context), "review_required", "reviewed"
      ),
      context_review_reason = dplyr::if_else(
        is.na(.data$climate_context),
        "Warm non-dry destination climate is outside the reviewed five-group scheme.",
        NA_character_
      )
    ) |>
    dplyr::arrange(.data$climate_context, .data$experiment_id)
  list(
    assigned = dplyr::filter(classified, !is.na(.data$climate_context)),
    review = dplyr::filter(classified, is.na(.data$climate_context))
  )
}

validate_five_climate_contexts <- function(context_result) {
  if (nrow(context_result$review)) {
    stop(
      "Destination climate contexts require review; see climate_context_review.csv.",
      call. = FALSE
    )
  }
  context_result$assigned
}

summarise_species_climate_context <- function(data, climate_contexts, minimum_experiments = 5L) {
  category_column <- if ("analysis_category" %in% names(data)) "analysis_category" else "category"
  observations <- data |>
    dplyr::transmute(
      species = dplyr::coalesce(.data$accepted_name, .data$species),
      experiment_id = .data$experiment_id,
      category = as.character(.data[[category_column]])
    ) |>
    dplyr::distinct()
  duplicated_cells <- observations |>
    dplyr::count(.data$species, .data$experiment_id) |>
    dplyr::filter(.data$n != 1L)
  if (nrow(duplicated_cells) > 0L) {
    stop("Accepted names do not uniquely identify species x experiment cells.", call. = FALSE)
  }
  eligible_species <- observations |>
    dplyr::summarise(n_experiments_total = dplyr::n_distinct(.data$experiment_id), .by = "species") |>
    dplyr::filter(.data$n_experiments_total >= .env$minimum_experiments)
  context_summary <- observations |>
    dplyr::semi_join(eligible_species, by = "species") |>
    dplyr::left_join(climate_contexts, by = "experiment_id") |>
    dplyr::mutate(
      response_score = dplyr::case_when(
        .data$category %in% c("winner", "colonizer") ~ 1,
        .data$category == "loser" ~ -1,
        TRUE ~ 0
      )
    ) |>
    dplyr::summarise(
      n_experiments = dplyr::n(),
      n_loser = sum(.data$category == "loser"),
      n_winner = sum(.data$category == "winner"),
      n_colonizer = sum(.data$category == "colonizer"),
      n_non_significant_loser = sum(.data$category == "non_significant_loser"),
      n_non_significant_winner = sum(.data$category == "non_significant_winner"),
      n_non_significant_colonizer = sum(.data$category == "non_significant_colonizer"),
      n_no_change = sum(.data$category == "no_change"),
      n_non_significant = sum(!.data$category %in% c("loser", "winner", "colonizer")),
      response_balance = mean(.data$response_score),
      .by = c(
        "species", "temperature_context", "precipitation_context", "climate_context",
        "temperature_bin_lower", "temperature_bin_upper",
        "precipitation_bin_lower", "precipitation_bin_upper"
      )
    ) |>
    dplyr::mutate(
      n_positive = .data$n_winner + .data$n_colonizer,
      within_context_reversal = .data$n_positive > 0L & .data$n_loser > 0L,
      dominant_response = dplyr::case_when(
        .data$within_context_reversal ~ "mixed significant reversal",
        .data$response_balance > 0 ~ "winner/colonizer dominated",
        .data$response_balance < 0 ~ "loser dominated",
        TRUE ~ "non-significant or balanced"
      )
    )
  species_patterns <- context_summary |>
    dplyr::summarise(
      n_experiments_total = sum(.data$n_experiments),
      n_contexts = dplyr::n(),
      overall_response_balance = stats::weighted.mean(.data$response_balance, .data$n_experiments),
      climate_response_range = max(.data$response_balance) - min(.data$response_balance),
      cross_context_reversal = any(.data$response_balance > 0) & any(.data$response_balance < 0),
      .by = "species"
    )
  context_summary |>
    dplyr::left_join(species_patterns, by = "species") |>
    dplyr::arrange(
      dplyr::desc(.data$cross_context_reversal),
      dplyr::desc(.data$climate_response_range),
      .data$species, .data$precipitation_context, .data$temperature_context
    )
}

analysis_dataset <- function(classified, community_taxonomy, trait_summary,
                             functional_groups, experiments, experiment_climate = NULL) {
  taxonomy_lookup <- community_taxonomy |>
    dplyr::select(species = name_original, accepted_name = accepted_name_reviewed) |>
    dplyr::distinct()
  traits <- trait_summary |>
    dplyr::select(accepted_name, trait, value = median) |>
    tidyr::pivot_wider(names_from = trait, values_from = value)
  result <- classified |>
    dplyr::filter(.data$transform == "raw", .data$horizon == "max") |>
    dplyr::mutate(
      calculation_category = .data$category,
      category = droplevels(.data$analysis_category),
      model_category = dplyr::case_when(
        .data$category == "loser" ~ "loser",
        .data$category == "winner" ~ "winner",
        .data$category == "colonizer" ~ "colonizer",
        stringr::str_starts(as.character(.data$category), "non_significant_") |
          .data$category == "no_change" ~ "no_change",
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::left_join(taxonomy_lookup, by = "species") |>
    dplyr::left_join(traits, by = "accepted_name") |>
    dplyr::left_join(functional_groups, by = "accepted_name") |>
    dplyr::left_join(
      experiments |>
        dplyr::select(experiment_id, plot_size, year_range),
      by = "experiment_id"
    )
  if (!is.null(experiment_climate)) {
    result <- join_experiment_climate(result, experiment_climate)
  }
  result
}

transform_model_predictors <- function(data) {
  positive <- c("LDMC", "LMA", "vegetative_height", "destination_P_JJA", "plot_size")
  missing <- setdiff(c(positive, "diff_T_JJA", "destination_T_JJA", "year_range"), names(data))
  if (length(missing)) stop("Model data are missing predictors: ", paste(missing, collapse = ", "), call. = FALSE)
  data |>
    dplyr::mutate(
      dplyr::across(dplyr::all_of(positive), ~ dplyr::if_else(.x > 0, log(.x), NA_real_), .names = "log_{.col}"),
      z_LDMC = safe_scale(.data$log_LDMC),
      z_LMA = safe_scale(.data$log_LMA),
      z_vegetative_height = safe_scale(.data$log_vegetative_height),
      z_diff_T_JJA = safe_scale(.data$diff_T_JJA),
      z_destination_T_JJA = safe_scale(.data$destination_T_JJA),
      z_destination_P_JJA = safe_scale(.data$log_destination_P_JJA),
      z_plot_size = safe_scale(.data$log_plot_size),
      z_year_range = safe_scale(.data$year_range)
    )
}

model_formula_text <- function(reduced = FALSE) {
  traits <- c("z_LDMC", "z_LMA", "z_vegetative_height")
  climate <- c("z_diff_T_JJA", "z_destination_T_JJA", "z_destination_P_JJA")
  main <- c(traits, climate, "z_plot_size", "z_year_range")
  interactions <- if (reduced) {
    paste(traits, "z_diff_T_JJA", sep = ":")
  } else {
    as.vector(outer(traits, climate, paste, sep = ":"))
  }
  paste(c(main, interactions), collapse = " + ")
}

fit_multinomial_model <- function(data) {
  complete <- data |>
    dplyr::filter(.data$model_category %in% c("loser", "no_change", "winner", "colonizer")) |>
    dplyr::filter(stats::complete.cases(dplyr::pick(
      model_category, species, experiment_id,
      dplyr::all_of(c("z_LDMC", "z_LMA", "z_vegetative_height", "z_diff_T_JJA",
                      "z_destination_T_JJA", "z_destination_P_JJA", "z_plot_size", "z_year_range"))
    ))) |>
    dplyr::mutate(model_category = stats::relevel(factor(.data$model_category), ref = "no_change"))
  fit <- function(reduced) {
    formula <- stats::as.formula(paste("model_category ~", model_formula_text(reduced)))
    mclogit::mblogit(
      formula,
      random = list(~ 1 | species, ~ 1 | experiment_id),
      data = complete,
      estimator = "ML",
      trace = FALSE
    )
  }
  primary <- tryCatch(fit(FALSE), error = identity)
  primary_ok <- !inherits(primary, "error") && isTRUE(primary$converged)
  fallback <- if (primary_ok) NULL else tryCatch(fit(TRUE), error = identity)
  list(
    primary = primary,
    fallback = fallback,
    selected = if (primary_ok) primary else fallback,
    used_fallback = !primary_ok,
    n_complete = nrow(complete),
    data = complete,
    random_effect_note = "crossed species and experiment random intercepts"
  )
}

fit_strength_meta_model <- function(data, reduced = FALSE) {
  variables <- c("estimate", "SE", "species", "experiment_id", "z_LDMC", "z_LMA",
                 "z_vegetative_height", "z_diff_T_JJA", "z_destination_T_JJA",
                 "z_destination_P_JJA", "z_plot_size", "z_year_range")
  complete <- data |>
    dplyr::filter(stats::complete.cases(dplyr::pick(dplyr::all_of(variables))), .data$SE > 0)
  model <- metafor::rma.mv(
    yi = estimate,
    V = SE^2,
    mods = stats::as.formula(paste("~", model_formula_text(reduced))),
    random = list(~ 1 | species, ~ 1 | experiment_id),
    data = complete,
    method = "REML"
  )
  robust <- metafor::robust(model, cluster = complete$experiment_id, clubSandwich = TRUE)
  list(model = model, robust = robust, data = complete, n_complete = nrow(complete))
}

fit_within_between_climate <- function(data) {
  repeated <- data |>
    dplyr::group_by(.data$species) |>
    dplyr::filter(dplyr::n_distinct(.data$experiment_id) > 3L) |>
    dplyr::mutate(
      species_mean_T = mean(.data$destination_T_JJA, na.rm = TRUE),
      species_mean_P = mean(log(.data$destination_P_JJA), na.rm = TRUE),
      within_T = .data$destination_T_JJA - .data$species_mean_T,
      within_P = log(.data$destination_P_JJA) - .data$species_mean_P
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(stats::complete.cases(.data$estimate, .data$within_T, .data$within_P,
                                       .data$species_mean_T, .data$species_mean_P))
  model <- lme4::lmer(
    estimate ~ within_T + within_P + species_mean_T + species_mean_P +
      (1 | species) + (1 | experiment_id),
    data = repeated,
    REML = FALSE
  )
  list(model = model, data = repeated)
}

model_diagnostics <- function(data) {
  predictors <- data |>
    dplyr::select(dplyr::starts_with("z_"))
  correlations <- stats::cor(predictors, use = "pairwise.complete.obs")
  tibble::tibble(
    n_rows = nrow(data),
    n_complete = sum(stats::complete.cases(data)),
    max_absolute_predictor_correlation = max(abs(correlations[upper.tri(correlations)]), na.rm = TRUE),
    n_experiments = dplyr::n_distinct(data$experiment_id),
    n_species = dplyr::n_distinct(data$species)
  )
}

tidy_multinomial_model <- function(model_result) {
  model <- model_result$selected
  if (is.null(model) || inherits(model, "error")) {
    diagnostic <- if (inherits(model, "error")) conditionMessage(model) else "no selected multinomial model"
    return(tibble::tibble(term = NA_character_, estimate = NA_real_, std_error = NA_real_,
                          conf_low = NA_real_, conf_high = NA_real_, odds_ratio = NA_real_,
                          odds_low = NA_real_, odds_high = NA_real_, diagnostic = diagnostic))
  }
  coefficients <- stats::coef(model)
  estimates <- if (is.matrix(coefficients)) as.vector(coefficients) else as.numeric(coefficients)
  coefficient_names <- if (is.matrix(coefficients)) {
    as.vector(outer(rownames(coefficients), colnames(coefficients), paste, sep = ":"))
  } else {
    names(coefficients)
  }
  covariance <- stats::vcov(model)
  standard_errors <- sqrt(diag(covariance))
  if (length(standard_errors) != length(estimates)) standard_errors <- rep(NA_real_, length(estimates))
  conf_low <- estimates - 1.96 * standard_errors
  conf_high <- estimates + 1.96 * standard_errors
  tibble::tibble(
    term = coefficient_names,
    estimate = estimates,
    std_error = standard_errors,
    conf_low = conf_low,
    conf_high = conf_high,
    odds_ratio = exp(estimates),
    odds_low = exp(conf_low),
    odds_high = exp(conf_high),
    diagnostic = NA_character_
  )
}

tidy_strength_model <- function(model_result) {
  robust <- model_result$robust
  estimates <- as.numeric(robust$b)
  standard_errors <- as.numeric(robust$se)
  terms <- rownames(robust$b)
  conf_low <- estimates - 1.96 * standard_errors
  conf_high <- estimates + 1.96 * standard_errors
  tibble::tibble(
    term = terms,
    estimate = estimates,
    std_error = standard_errors,
    conf_low = conf_low,
    conf_high = conf_high,
    z_value = as.numeric(robust$zval),
    p_value = as.numeric(robust$pval)
  )
}

multinomial_prediction_grid <- function(model_result, points = 31L) {
  model <- model_result$selected
  if (is.null(model) || inherits(model, "error")) return(tibble::tibble())
  traits <- c("z_LDMC", "z_LMA", "z_vegetative_height")
  climates <- c("z_diff_T_JJA", "z_destination_T_JJA", "z_destination_P_JJA")
  reference <- model_result$data[1, , drop = FALSE]
  grids <- purrr::map_dfr(traits, function(trait) {
    purrr::map_dfr(climates, function(climate) {
      grid <- reference[rep(1, points * 3L), , drop = FALSE]
      for (column in c(traits, climates, "z_plot_size", "z_year_range")) grid[[column]] <- 0
      grid[[trait]] <- rep(c(-1, 0, 1), each = points)
      grid[[climate]] <- rep(seq(-2, 2, length.out = points), times = 3)
      grid$focal_trait <- trait
      grid$focal_climate <- climate
      grid$trait_level <- factor(grid[[trait]], levels = c(-1, 0, 1), labels = c("-1 SD", "mean", "+1 SD"))
      probabilities <- tryCatch(stats::predict(model, newdata = grid, type = "response"), error = function(e) NULL)
      if (is.null(probabilities)) return(tibble::tibble())
      if (is.vector(probabilities)) probabilities <- matrix(probabilities, ncol = 1)
      probability_table <- tibble::as_tibble(probabilities, .name_repair = "unique")
      names(probability_table) <- colnames(probabilities) %||% paste0("category_", seq_len(ncol(probabilities)))
      dplyr::bind_cols(
        grid |>
          dplyr::transmute(.data$focal_trait, .data$focal_climate, .data$trait_level,
                           climate_value = .data[[climate]]),
        probability_table
      ) |>
        tidyr::pivot_longer(-c(.data$focal_trait, .data$focal_climate, .data$trait_level, .data$climate_value),
                            names_to = "category", values_to = "probability")
    })
  })
  grids
}

leave_one_experiment_out_strength <- function(data) {
  experiments <- unique(data$experiment_id)
  full <- fit_strength_meta_model(data)$robust$b[, 1]
  purrr::map_dfr(experiments, function(held_out) {
    result <- tryCatch(fit_strength_meta_model(dplyr::filter(data, .data$experiment_id != .env$held_out)), error = identity)
    if (inherits(result, "error")) {
      return(tibble::tibble(held_out_experiment = held_out, term = NA_character_, estimate = NA_real_,
                            full_estimate = NA_real_, absolute_change = NA_real_, error = conditionMessage(result)))
    }
    estimate <- result$robust$b[, 1]
    full_estimate <- as.numeric(full[names(estimate)])
    tibble::tibble(
      held_out_experiment = held_out,
      term = names(estimate),
      estimate = as.numeric(estimate),
      full_estimate = full_estimate,
      absolute_change = abs(as.numeric(estimate) - full_estimate),
      error = NA_character_
    )
  })
}

`%||%` <- function(x, y) if (is.null(x)) y else x
