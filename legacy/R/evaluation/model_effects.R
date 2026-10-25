# Module: R/evaluation/model_effects
# Responsibility: Summarize fixed coefficients and random-component uncertainty for the appendix.
# Inputs: Fitted resident and predictor models plus profile settings.
# Outputs: Model-effect tables, profile results, and figure data.
# Side effects: Profile calculations fit evaluation models; no files written here.

# Fixed-effect coefficients and species-level random-effect heterogeneity for
# the seven registered resident and predictor models.

#' Model Effect Predictor Key.
#'
#' @return An R object described by the function name; its structure is validated before later stages use it.
model_effect_predictor_key <- function() {
  c(
    "resident_climate" = "Climate",
    "LDMC" = "LDMC",
    "LMA" = "LMA",
    "Height" = "Height",
    "Graminoid status" = "Graminoid",
    "Woodiness" = "Woodiness",
    "Phylogenetic position (PEV1-PEV3)" = "Phylogeny"
  )
}

#' Model Effect Pretty Term.
#'
#' @param term A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param model_key A fitted model object with the structure checked by this function.
#' @return An R object described by the function name; its structure is validated before later stages use it.
model_effect_pretty_term <- function(term, model_key) {
  trait_label <- c(
    LDMC = "LDMC", LMA = "LMA", Height = "Vegetative height"
  )
  token_label <- c(
    intrcpt = "Intercept",
    T_z = "Destination temperature",
    P_z = "Destination precipitation",
    W_z = "Warming",
    TP_z = "Temperature × precipitation",
    GraminoidStatusgraminoid = "Graminoid",
    Woodinesswoody = "Woody",
    PEV1_z = "PEV1",
    PEV2_z = "PEV2",
    PEV3_z = "PEV3"
  )
  pieces <- strsplit(as.character(term), ":", fixed = TRUE)
  vapply(seq_along(pieces), function(index) {
    tokens <- pieces[[index]]
    labels <- vapply(tokens, function(token) {
      if (identical(token, "Trait_z")) {
        label <- trait_label[[as.character(model_key[[index]])]]
        if (is.null(label)) "Standardized trait" else unname(label)
      } else if (token %in% names(token_label)) {
        unname(token_label[[token]])
      } else {
        token
      }
    }, character(1))
    paste(labels, collapse = " × ")
  }, character(1))
}

#' Assemble Model Fixed Effects.
#'
#' @param climate_cr1 A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param predictor_cr1 A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param fit_r2 A fitted model object with the structure checked by this function.
#' @return An R object described by the function name; its structure is validated before later stages use it.
assemble_model_fixed_effects <- function(
    climate_cr1, predictor_cr1, fit_r2) {
  key <- model_performance_key()
  predictor_key <- model_effect_predictor_key()
  fixed <- dplyr::bind_rows(climate_cr1, predictor_cr1) |>
    dplyr::mutate(
      model_key = unname(.env$predictor_key[as.character(.data$predictor)]),
      model = model_performance_display(.data$model_key),
      model_order = .env$key$model_order[
        match(.data$model_key, .env$key$model_key)
      ]
    )
  if (anyNA(fixed$model_key) || anyNA(fixed$model_order)) {
    stop("One or more fixed-effect models could not be identified.",
         call. = FALSE)
  }
  samples <- fit_r2 |>
    dplyr::filter(
      .data$weighting == "equal_gradient",
      .data$candidate_role %in% c("standalone", "full")
    ) |>
    dplyr::select("model_key", "n_effects", "n_species", "n_gradients")
  if (anyDuplicated(samples$model_key) || nrow(samples) != nrow(key)) {
    stop("The fixed-effect figure needs one sample-size row per model.",
         call. = FALSE)
  }
  fixed |>
    dplyr::left_join(samples, by = "model_key") |>
    dplyr::mutate(
      term_order = dplyr::row_number(),
      .by = "model_key"
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      term_label = model_effect_pretty_term(.data$term, .data$model_key),
      supported = .data$CR1_CI_low > 0 | .data$CR1_CI_high < 0,
      inference = dplyr::if_else(
        .data$supported,
        "CR1 95% CI excludes zero",
        "CR1 95% CI includes zero"
      ),
      inference_method = paste(
        "experiment-clustered CR1 coefficient interval;",
        "two-sided 95% confidence level"
      )
    ) |>
    dplyr::select(
      "model_key", "model", "model_order", "n_effects", "n_species",
      "n_gradients", "term", "term_label", "term_order", "estimate",
      "CR1_SE", "CR1_CI_low", "CR1_CI_high", "CR1_p_value",
      "CR1_denominator_df", "supported", "inference", "inference_method"
    )
}

#' Build Predictor Variance Profile Jobs.
#'
#' @param ... Additional arguments forwarded to the delegated operation.
#' @return An R object described by the function name; its structure is validated before later stages use it.
build_predictor_variance_profile_jobs <- function(...) {
  sequences <- list(...)
  specifications <- variance_profile_specifications()
  jobs <- list()
  for (sequence in sequences) {
    for (component in names(specifications)) {
      job_name <- paste(sequence$key, component, sep = "__")
      jobs[[job_name]] <- list(
        model_key = sequence$key,
        model = model_performance_display(sequence$key),
        captured = sequence$fits$full_REML,
        specification = specifications[[component]]
      )
    }
  }
  jobs
}

#' Profile Predictor Variance Component.
#'
#' @param job A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
profile_predictor_variance_component <- function(job) {
  profile <- profile_model_variance(
    job$captured, job$specification
  )
  profile |>
    dplyr::mutate(
      model_key = job$model_key,
      model = job$model,
      .before = 1L
    )
}

#' Assemble Model Random Effect Profiles.
#'
#' @param climate_profiles A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param predictor_profiles A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
assemble_model_random_effect_profiles <- function(
    climate_profiles, predictor_profiles) {
  key <- model_performance_key()
  climate <- climate_profiles |>
    dplyr::mutate(
      model_key = "Climate",
      model = model_performance_display("Climate"),
      .before = 1L
    )
  profiles <- dplyr::bind_rows(
    climate,
    dplyr::bind_rows(predictor_profiles)
  ) |>
    dplyr::mutate(
      model_order = .env$key$model_order[
        match(.data$model_key, .env$key$model_key)
      ],
      component_order = match(
        .data$component,
        c(
          "species_intercept", "species_W", "species_T",
          "species_P", "species_TP"
        )
      ),
      component_label = dplyr::recode(
        .data$component,
        species_intercept = "Species intercept",
        species_W = "Warming slope",
        species_T = "Temperature slope",
        species_P = "Precipitation slope",
        species_TP = "Temperature × precipitation slope"
      ),
      SD_estimate = sqrt(pmax(.data$estimate, 0)),
      SD_CI_low = sqrt(pmax(.data$CI_low, 0)),
      SD_CI_high = sqrt(pmax(.data$CI_high, 0)),
      supported = dplyr::case_when(
        is.finite(.data$CI_low) & is.finite(.data$CI_high) ~
          .data$CI_low > 0,
        TRUE ~ NA
      ),
      inference = dplyr::case_when(
        is.na(.data$supported) ~
          "lower profile-likelihood bound unavailable; support undetermined",
        .data$supported ~ "profile-likelihood 95% CI excludes zero",
        TRUE ~ "profile-likelihood 95% CI includes zero"
      ),
      inference_method = paste(
        "metafor profile-likelihood interval for the variance component;",
        "displayed after square-root transformation to a standard deviation"
      ),
      warnings = dplyr::na_if(as.character(.data$warnings), ""),
      error = dplyr::na_if(as.character(.data$error), "")
    ) |>
    dplyr::arrange(.data$model_order, .data$component_order)
  profiles
}

#' Model Random Diagnostic Row.
#'
#' @param captured A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param model_key A fitted model object with the structure checked by this function.
#' @return An R object described by the function name; its structure is validated before later stages use it.
model_random_diagnostic_row <- function(captured, model_key) {
  diagnostics <- random_effect_diagnostics(
    captured,
    model_performance_display(model_key),
    if (identical(model_key, "Climate")) "resident_climate" else "full_REML"
  )
  diagnostics |>
    dplyr::mutate(
      model_key = model_key,
      model = model_performance_display(model_key),
      .before = 1L
    )
}

#' Assemble Model Random Correlations.
#'
#' @param climate_model A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param ... Additional arguments forwarded to the delegated operation.
#' @return An R object described by the function name; its structure is validated before later stages use it.
assemble_model_random_correlations <- function(climate_model, ...) {
  sequences <- list(...)
  rows <- list(model_random_diagnostic_row(climate_model, "Climate"))
  for (sequence in sequences) {
    rows[[length(rows) + 1L]] <- model_random_diagnostic_row(
      sequence$fits$full_REML, sequence$key
    )
  }
  key <- model_performance_key()
  dplyr::bind_rows(rows) |>
    dplyr::mutate(
      model_order = .env$key$model_order[
        match(.data$model_key, .env$key$model_key)
      ]
    ) |>
    dplyr::arrange(.data$model_order)
}

#' Validate Model Effects.
#'
#' @param fixed A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param random_profiles A model formula or random-effects specification.
#' @param random_diagnostics A model formula or random-effects specification.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
validate_model_effects <- function(fixed, random_profiles,
                                       random_diagnostics) {
  key <- model_performance_key()
  expected_components <- c(
    "species_intercept", "species_W", "species_T", "species_P", "species_TP"
  )
  variance_columns <- c(
    species_intercept = "species_intercept_variance",
    species_W = "species_W_variance",
    species_T = "species_T_variance",
    species_P = "species_P_variance",
    species_TP = "species_TP_variance"
  )
  diagnostic_long <- random_diagnostics |>
    tidyr::pivot_longer(
      dplyr::all_of(unname(variance_columns)),
      names_to = "variance_column", values_to = "diagnostic_estimate"
    ) |>
    dplyr::mutate(
      component = names(variance_columns)[
        match(.data$variance_column, unname(variance_columns))
      ]
    ) |>
    dplyr::select("model_key", "component", "diagnostic_estimate")
  profile_check <- random_profiles |>
    dplyr::left_join(diagnostic_long, by = c("model_key", "component"))
  maximum_profile_error <- max(abs(
    profile_check$estimate - profile_check$diagnostic_estimate
  ))
  correlation_values <- unlist(
    random_diagnostics[c(
      "T_P_correlation", "T_TP_correlation", "P_TP_correlation"
    )],
    use.names = FALSE
  )
  fixed_interval_support <- fixed$CR1_CI_low > 0 | fixed$CR1_CI_high < 0
  fixed_p_support <- fixed$CR1_p_value < 0.05
  complete_profile <- with(random_profiles,
    is.finite(estimate) & is.finite(CI_low) & is.finite(CI_high) &
      CI_low >= 0 & CI_low <= estimate & estimate <= CI_high
  )
  # The October 2026 rebuild has one unresolved lower profile bound. Keep its
  # estimate and warning visible, but do not interpret the incomplete interval.
  documented_incomplete_profile <- with(random_profiles,
    model_key == "LMA" & component == "species_P" &
      is.finite(estimate) & is.na(CI_low) & is.finite(CI_high) &
      estimate <= CI_high & is.na(error) &
      !is.na(warnings) &
      grepl(
        "Cannot obtain lower bound of profile likelihood CI due to convergence problems.",
        warnings, fixed = TRUE
      ) & is.na(supported) &
      inference ==
        "lower profile-likelihood bound unavailable; support undetermined"
  )
  profile_intervals_acceptable <-
    all(is.na(random_profiles$error)) &&
    all(complete_profile | documented_incomplete_profile) &&
    sum(documented_incomplete_profile) <= 1L
  validation <- tibble::tribble(
    ~check, ~pass, ~detail,
    "Seven current species models",
    identical(sort(unique(fixed$model_key)), sort(key$model_key)) &&
      identical(sort(unique(random_profiles$model_key)), sort(key$model_key)),
    paste(nrow(key), "models represented in fixed and random summaries"),
    "All fixed coefficients retained",
    nrow(fixed) == 111L && !anyDuplicated(fixed[c("model_key", "term")]),
    paste(nrow(fixed), "unique fixed coefficients"),
    "Finite CR1 coefficient inference",
    all(is.finite(fixed$estimate)) && all(is.finite(fixed$CR1_SE)) &&
      all(is.finite(fixed$CR1_CI_low)) &&
      all(is.finite(fixed$CR1_CI_high)) &&
      all(is.finite(fixed$CR1_p_value)),
    "Every coefficient has a finite estimate, CR1 SE, interval, and p-value",
    "Fixed-effect support coding",
    identical(fixed$supported, fixed_interval_support) &&
      identical(fixed_interval_support, fixed_p_support),
    paste(sum(fixed$supported), "CR1 intervals exclude zero"),
    "Five random components per model",
    nrow(random_profiles) == 35L &&
      !anyDuplicated(random_profiles[c("model_key", "component")]) &&
      all(table(random_profiles$model_key) == length(expected_components)) &&
      all(expected_components %in% random_profiles$component),
    paste(nrow(random_profiles), "model-component profiles"),
    "Profile intervals complete except documented LMA lower bound",
    profile_intervals_acceptable,
    paste(
      sum(complete_profile), "complete intervals;",
      sum(documented_incomplete_profile),
      "LMA precipitation-slope lower bound unavailable (support undetermined)"
    ),
    "Profile estimates match fitted variances",
    is.finite(maximum_profile_error) && maximum_profile_error <= 1e-10,
    paste("maximum absolute variance difference",
          format(maximum_profile_error, scientific = TRUE)),
    "Valid random-effect correlations",
    all(is.finite(correlation_values)) && all(abs(correlation_values) <= 1),
    "All T/P/TP random-slope correlations are finite and within [-1, 1]"
  )
  validation
}

#' Assert Model Effects Validation.
#'
#' @param validation A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
assert_model_effects_validation <- function(validation) {
  failed <- is.na(validation$pass) | !validation$pass
  if (any(failed)) {
    stop(
      "Model-effect validation failed: ",
      paste(validation$check[failed], collapse = ", "),
      call. = FALSE
    )
  }
  validation
}

#' Model Effect Panel Data.
#'
#' @param fixed A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
model_effect_panel_data <- function(fixed) {
  key <- model_performance_key()
  model_labels <- fixed |>
    dplyr::distinct(
      .data$model_key, .data$model, .data$model_order,
      .data$n_effects, .data$n_species
    ) |>
    dplyr::mutate(
      model_label = paste0(
        .data$model, "\n(n = ", .data$n_effects,
        " effects; ", .data$n_species, " species)"
      )
    )
  model_label_lookup <- stats::setNames(
    model_labels$model_label, model_labels$model_key
  )
  ordered <- fixed |>
    dplyr::arrange(.data$model_order, .data$term_order) |>
    dplyr::mutate(
      plot_id = paste0(
        .data$model_key, "__", sprintf("%02d", .data$term_order),
        "__", .data$term
      )
    )
  plot_levels <- rev(ordered$plot_id)
  term_labels <- stats::setNames(ordered$term_label, ordered$plot_id)
  ordered |>
    dplyr::mutate(
      plot_id = factor(.data$plot_id, levels = .env$plot_levels),
      model_label = factor(
        unname(.env$model_label_lookup[.data$model_key]),
        levels = unname(.env$model_label_lookup[.env$key$model_key])
      )
    ) |>
    structure(term_labels = term_labels)
}

#' Build Model Effects Figure.
#'
#' @param fixed A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param random_profiles A model formula or random-effects specification.
#' @return A plot object, or the declared figure path when rendering to disk.
build_model_effects_figure <- function(fixed, random_profiles) {
  key <- model_performance_key()
  fixed_plot_data <- model_effect_panel_data(fixed)
  term_labels <- attr(fixed_plot_data, "term_labels")
  support_labels <- c(
    `FALSE` = "95% interval includes zero",
    `TRUE` = "95% interval excludes zero"
  )
  support_fill <- c(`FALSE` = "white", `TRUE` = "#2166AC")

  fixed_plot <- ggplot2::ggplot(
    fixed_plot_data,
    ggplot2::aes(x = .data$estimate, y = .data$plot_id)
  ) +
    ggplot2::geom_vline(
      xintercept = 0, colour = "grey55", linewidth = 0.55
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$CR1_CI_low, xmax = .data$CR1_CI_high),
      width = 0.18, linewidth = 0.65, orientation = "y"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = factor(.data$supported)),
      shape = 21, colour = "#202020", stroke = 0.65, size = 3.2
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(.data$model_label),
      scales = "free", ncol = 2
    ) +
    ggplot2::scale_y_discrete(labels = term_labels) +
    ggplot2::scale_fill_manual(
      values = support_fill, labels = support_labels,
      name = "Fixed-effect CR1 interval"
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.08, 0.08))
    ) +
    ggplot2::labs(
      x = "Fixed-effect coefficient (change in ΔAUC)", y = NULL
    ) +
    ggplot2::theme_bw(base_size = 17) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.spacing = grid::unit(1.25, "lines"),
      strip.text = ggplot2::element_text(size = 18, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 13.5, colour = "black"),
      axis.text.y = ggplot2::element_text(size = 12.2, colour = "black"),
      axis.title.x = ggplot2::element_text(size = 17),
      legend.position = "top",
      legend.title = ggplot2::element_text(size = 15, face = "bold"),
      legend.text = ggplot2::element_text(size = 14),
      plot.margin = ggplot2::margin(10, 14, 10, 10)
    )

  model_labels <- fixed |>
    dplyr::distinct(
      .data$model_key, .data$model, .data$model_order,
      .data$n_effects, .data$n_species
    ) |>
    dplyr::mutate(
      model_label = paste0(
        .data$model, "\n(n = ", .data$n_effects,
        " effects; ", .data$n_species, " species)"
      )
    )
  model_label_lookup <- stats::setNames(
    model_labels$model_label, model_labels$model_key
  )
  component_levels <- rev(c(
    "Species intercept", "Warming slope", "Temperature slope",
    "Precipitation slope", "Temperature × precipitation slope"
  ))
  random_plot_data <- random_profiles |>
    dplyr::mutate(
      model_label = factor(
        unname(.env$model_label_lookup[.data$model_key]),
        levels = unname(.env$model_label_lookup[.env$key$model_key])
      ),
      component_label = factor(
        .data$component_label, levels = .env$component_levels
      )
    )
  random_interval_fill <- c(
    "profile-likelihood 95% CI includes zero" = "white",
    "profile-likelihood 95% CI excludes zero" = "#2166AC",
    "lower profile-likelihood bound unavailable; support undetermined" =
      "#E69F00"
  )
  random_plot <- ggplot2::ggplot(
    random_plot_data,
    ggplot2::aes(x = .data$SD_estimate, y = .data$component_label)
  ) +
    ggplot2::geom_vline(
      xintercept = 0, colour = "grey55", linewidth = 0.55
    ) +
    ggplot2::geom_errorbar(
      data = dplyr::filter(
        random_plot_data,
        is.finite(.data$SD_CI_low), is.finite(.data$SD_CI_high)
      ),
      ggplot2::aes(xmin = .data$SD_CI_low, xmax = .data$SD_CI_high),
      width = 0.18, linewidth = 0.7, orientation = "y"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = .data$inference),
      shape = 21, colour = "#202020", stroke = 0.65, size = 3.5
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(.data$model_label), ncol = 4
    ) +
    ggplot2::scale_fill_manual(
      values = random_interval_fill,
      name = "Random-component profile interval"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, NA),
      expand = ggplot2::expansion(mult = c(0.02, 0.08))
    ) +
    ggplot2::labs(
      x = "Species-level random-effect SD",
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 17) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.spacing = grid::unit(1.1, "lines"),
      strip.text = ggplot2::element_text(size = 16, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 13.5, colour = "black"),
      axis.text.y = ggplot2::element_text(size = 12.2, colour = "black"),
      axis.title.x = ggplot2::element_text(size = 17),
      legend.position = "top",
      legend.title = ggplot2::element_text(size = 15, face = "bold"),
      legend.text = ggplot2::element_text(size = 14),
      plot.margin = ggplot2::margin(10, 14, 10, 10)
    )

  patchwork::wrap_plots(
    fixed_plot, random_plot,
    ncol = 1, heights = c(3.5, 1.35)
  ) +
    patchwork::plot_annotation(
      tag_levels = "A",
      theme = ggplot2::theme(
        plot.tag = ggplot2::element_text(size = 24, face = "bold")
      )
    )
}
