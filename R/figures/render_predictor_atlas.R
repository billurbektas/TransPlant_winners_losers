# Module: R/figures/render_predictor_atlas
# Responsibility: Implement the named scientific workflow stage or shared infrastructure.
# Inputs: Function arguments and explicitly declared target dependencies.
# Outputs: Validated in-memory objects and declared workflow artifacts.
# Side effects: Limited to functions whose contracts explicitly write, copy, or render files.

# Render the predictor-response atlas
#
# This renderer reads existing exported analysis artifacts and never fits a model.

#' Render the predictor-response atlas
#'
#' @param arguments Character vector of renderer-specific positional arguments.
#' @return A character vector of generated artifact paths, invisibly.
render_predictor_response_atlas <- function(arguments = character()) {
  
  # Candidate prediction atlas for the six canonical 526-species predictor
  # models. This script reads fitted model RDS files and canonical prediction
  # inputs only. It contains no model-fitting calls and never invokes tar_make().
  
  required_packages <- c(
    "cowplot", "dplyr", "geometry", "ggplot2", "metafor", "readr",
    "targets", "tibble", "tidyr"
  )
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages) > 0L) {
    stop(
      "Missing required package(s): ", paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }
  if (!file.exists("_targets.R")) {
    stop("Run this script from the repository root.", call. = FALSE)
  }
  
  suppressPackageStartupMessages(library(metafor))
  
  output_stem <- if (length(arguments) >= 1L) {
    arguments[[1]]
  } else {
    file.path("figures", "candidates", "predictor_response_atlas")
  }
  climate_scaling_file <- if (length(arguments) >= 2L) {
    arguments[[2]]
  } else {
    NA_character_
  }
  output_root <- if (length(arguments) >= 3L) {
    arguments[[3]]
  } else {
    "output"
  }
  pdf_file <- paste0(output_stem, ".pdf")
  png_file <- paste0(output_stem, ".png")
  output_name <- basename(output_stem)
  data_file <- file.path(
    output_root, "data", paste0(output_name, "_data.csv")
  )
  quantile_file <- file.path(
    output_root, "tables", paste0(output_name, "_quantiles.csv")
  )
  caption_file <- paste0(output_stem, "_caption.txt")
  validation_file <- file.path(
    output_root, "audits", paste0(output_name, "_validation.csv")
  )
  
  model_files <- c(
    LDMC = "output/models/predictor_ldmc_model.rds",
    LMA = "output/models/predictor_lma_model.rds",
    Height = "output/models/predictor_height_model.rds",
    Graminoid = "output/models/predictor_graminoid_model.rds",
    Woodiness = "output/models/predictor_woodiness_model.rds",
    Phylogeny = "output/models/predictor_phylogeny_model.rds"
  )
  input_files <- c(
    model_files,
    trait_scaling = "output/tables/trait_scaling.csv",
    pev_scaling = "output/tables/phylogeny_predictor_scaling.csv",
    pev_scores = "output/data/phylogeny_predictor_scores.csv"
  )
  missing_inputs <- input_files[!file.exists(input_files)]
  if (length(missing_inputs) > 0L) {
    stop(
      "Required canonical input(s) missing: ",
      paste(missing_inputs, collapse = ", "),
      call. = FALSE
    )
  }
  
  models <- lapply(model_files, readRDS)
  if (any(!vapply(models, inherits, logical(1), what = "rma.mv"))) {
    stop("Every atlas model input must be an existing rma.mv fit.", call. = FALSE)
  }
  if (any(!vapply(models, function(x) is.data.frame(x$data), logical(1)))) {
    stop("A fitted model is missing its canonical modeling data.", call. = FALSE)
  }
  
  # Standalone candidate rendering reads the existing canonical target object.
  # The anchored publication target instead supplies its declared scaling-file
  # dependency as the optional second command-line argument.
  climate_scaling <- if (!is.na(climate_scaling_file)) {
    if (!file.exists(climate_scaling_file)) {
      stop("Climate-scaling input does not exist: ", climate_scaling_file,
           call. = FALSE)
    }
    readRDS(climate_scaling_file)
  } else {
    targets::tar_read_raw("resident_climate_scaling")
  }
  trait_scaling <- readr::read_csv(
    input_files[["trait_scaling"]], show_col_types = FALSE
  )
  pev_scaling <- readr::read_csv(
    input_files[["pev_scaling"]], show_col_types = FALSE
  )
  pev_scores <- readr::read_csv(
    input_files[["pev_scores"]], show_col_types = FALSE
  )
  
  single_value <- function(data, key_column, key, value_column, label) {
    value <- data[[value_column]][data[[key_column]] == key]
    if (length(value) != 1L || !is.finite(value)) {
      stop("Missing unique ", label, " for ", key, ".", call. = FALSE)
    }
    as.numeric(value[[1]])
  }
  
  climate_value <- function(field, statistic) {
    single_value(
      climate_scaling, "scaled_field", field, statistic,
      paste("climate", statistic)
    )
  }
  
  T_center <- climate_value("T_z", "center")
  T_scale <- climate_value("T_z", "scale")
  P_center <- climate_value("P_z", "center")
  P_scale <- climate_value("P_z", "scale")
  W_center <- climate_value("W_z", "center")
  W_scale <- climate_value("W_z", "scale")
  
  reference_experiments <- models$Graminoid$data |>
    dplyr::distinct(.data$experiment_id, .data$T, .data$P, .data$P_raw,
                    .data$W, .data$T_z, .data$P_z, .data$W_z) |>
    dplyr::arrange(.data$experiment_id)
  if (nrow(reference_experiments) != 42L) {
    stop("Warming quantiles require exactly 42 distinct experiments.", call. = FALSE)
  }
  
  climate_scaling_error <- max(c(
    abs(reference_experiments$T_z -
          (reference_experiments$T - T_center) / T_scale),
    abs(reference_experiments$P_z -
          (reference_experiments$P - P_center) / P_scale),
    abs(reference_experiments$W_z -
          (reference_experiments$W - W_center) / W_scale)
  ))
  if (!is.finite(climate_scaling_error) || climate_scaling_error > 1e-10) {
    stop("Stored model climate columns do not match frozen scaling.", call. = FALSE)
  }
  
  warming_raw <- stats::quantile(
    reference_experiments$W,
    probs = c(0.25, 0.75),
    names = FALSE,
    type = 7
  )
  warming_model <- (warming_raw - W_center) / W_scale
  warming_table <- tibble::tibble(
    quantity = "Experimental warming",
    level = c("Low warming", "High warming"),
    quantile = c("Q25", "Q75"),
    raw_value = warming_raw,
    raw_unit = "degrees C",
    model_field = "W_z",
    model_value = warming_model,
    scaling_center = W_center,
    scaling_scale = W_scale,
    n_unique = nrow(reference_experiments),
    weighting = "42 distinct experiments; each experiment contributes once"
  )
  
  trait_specs <- tibble::tribble(
    ~model_key, ~predictor,          ~trait_key,           ~low_label,   ~high_label,   ~raw_unit,
    "LDMC",     "LDMC",             "LDMC",               "Low LDMC",   "High LDMC",   "g/g",
    "LMA",      "LMA",              "LMA",                "Low LMA",    "High LMA",    "1/(mm2 mg-1)",
    "Height",   "Vegetative height", "vegetative_height",  "Low height", "High height", "m"
  )
  
  trait_quantile_rows <- lapply(seq_len(nrow(trait_specs)), function(index) {
    spec <- trait_specs[index, ]
    model_data <- models[[spec$model_key]]$data |>
      dplyr::distinct(.data$species_id, .data$trait_value)
    if (anyDuplicated(model_data$species_id)) {
      stop(spec$predictor, " trait is not unique by species.", call. = FALSE)
    }
    expected_n <- single_value(
      trait_scaling, "trait", spec$trait_key, "species_N_with_trait",
      "trait-complete species count"
    )
    if (nrow(model_data) != expected_n) {
      stop(spec$predictor, " trait-complete species count changed.", call. = FALSE)
    }
    raw_values <- stats::quantile(
      model_data$trait_value,
      probs = c(0.25, 0.75),
      names = FALSE,
      type = 7
    )
    log_center <- single_value(
      trait_scaling, "trait", spec$trait_key, "log_trait_mean",
      "log-trait center"
    )
    log_scale <- single_value(
      trait_scaling, "trait", spec$trait_key, "log_trait_SD",
      "log-trait scale"
    )
    tibble::tibble(
      quantity = spec$predictor,
      level = c(spec$low_label, spec$high_label),
      quantile = c("Q25", "Q75"),
      raw_value = raw_values,
      raw_unit = spec$raw_unit,
      model_field = "Trait_z",
      model_value = (log(raw_values) - log_center) / log_scale,
      scaling_center = log_center,
      scaling_scale = log_scale,
      n_unique = nrow(model_data),
      weighting = paste(
        "unique retained trait-complete species;",
        "raw species median -> natural log -> frozen standardization"
      )
    )
  })
  
  if (nrow(pev_scores) != 526L || anyDuplicated(pev_scores$species_id)) {
    stop("PEV1 quantiles require exactly 526 unique species.", call. = FALSE)
  }
  PEV1_center <- single_value(
    pev_scaling, "PEV", "PEV1", "raw_mean", "PEV1 center"
  )
  PEV1_scale <- single_value(
    pev_scaling, "PEV", "PEV1", "raw_SD", "PEV1 scale"
  )
  PEV1_raw <- stats::quantile(
    pev_scores$PEV1, probs = c(0.25, 0.75), names = FALSE, type = 7
  )
  PEV1_model <- (PEV1_raw - PEV1_center) / PEV1_scale
  pev_quantile_rows <- tibble::tibble(
    quantity = "Phylogenetic position (PEV1)",
    level = c("Low PEV1", "High PEV1"),
    quantile = c("Q25", "Q75"),
    raw_value = PEV1_raw,
    raw_unit = "PCoA coordinate",
    model_field = "PEV1_z",
    model_value = PEV1_model,
    scaling_center = PEV1_center,
    scaling_scale = PEV1_scale,
    n_unique = nrow(pev_scores),
    weighting = "all 526 unique species from reconstructed canonical PCoA"
  )
  
  model_pev <- models$Phylogeny$data |>
    dplyr::distinct(.data$species_id, .data$PEV1_z)
  pev_check <- pev_scores |>
    dplyr::transmute(
      species_id = .data$species_id,
      expected_PEV1_z = (.data$PEV1 - PEV1_center) / PEV1_scale
    ) |>
    dplyr::inner_join(model_pev, by = "species_id")
  pev_standardization_error <- max(abs(
    pev_check$expected_PEV1_z - pev_check$PEV1_z
  ))
  if (nrow(pev_check) != 526L || pev_standardization_error > 1e-10) {
    stop("Phylogeny model PEV1 does not match the canonical PCoA export.", call. = FALSE)
  }
  
  quantile_table <- dplyr::bind_rows(
    trait_quantile_rows,
    list(pev_quantile_rows, warming_table)
  )
  
  level_specs <- dplyr::bind_rows(
    lapply(trait_quantile_rows, function(x) {
      model_key <- trait_specs$model_key[match(x$quantity, trait_specs$predictor)]
      x |>
        dplyr::transmute(
          model_key = .env$model_key,
          predictor = .data$quantity,
          predictor_level = .data$level,
          predictor_raw_value = .data$raw_value,
          predictor_model_value = .data$model_value,
          predictor_raw_unit = .data$raw_unit,
          predictor_type = "continuous trait",
          factor_value = NA_character_
        )
    }),
    tibble::tribble(
      ~model_key,  ~predictor,          ~predictor_level, ~predictor_raw_value, ~predictor_model_value, ~predictor_raw_unit, ~predictor_type,        ~factor_value,
      "Graminoid", "Graminoid status", "Non-graminoid", 0,                    0,                      "binary indicator", "categorical functional", "non_graminoid",
      "Graminoid", "Graminoid status", "Graminoid",     1,                    1,                      "binary indicator", "categorical functional", "graminoid",
      "Woodiness", "Woodiness",        "Non-woody",     0,                    0,                      "binary indicator", "categorical functional", "non_woody",
      "Woodiness", "Woodiness",        "Woody",         1,                    1,                      "binary indicator", "categorical functional", "woody"
    ),
    pev_quantile_rows |>
      dplyr::transmute(
        model_key = "Phylogeny",
        predictor = "Phylogenetic position (PEV1 shown)",
        predictor_level = .data$level,
        predictor_raw_value = .data$raw_value,
        predictor_model_value = .data$model_value,
        predictor_raw_unit = .data$raw_unit,
        predictor_type = "representative phylogenetic gradient",
        factor_value = NA_character_
      )
  )
  
  predictor_order <- c(
    "LDMC", "LMA", "Vegetative height", "Graminoid status", "Woodiness",
    "Phylogenetic position (PEV1 shown)"
  )
  predictor_display_labels <- c(
    "LDMC" = "LDMC",
    "LMA" = "LMA",
    "Vegetative height" = "Vegetative\nheight",
    "Graminoid status" = "Graminoid\nstatus",
    "Woodiness" = "Woodiness",
    "Phylogenetic position (PEV1 shown)" = "Phylogenetic\nposition"
  )
  formula_map <- predictor_formulas()
  model_formula <- list(
    LDMC = formula_map$continuous_full,
    LMA = formula_map$continuous_full,
    Height = formula_map$continuous_full,
    Graminoid = formula_map$Graminoid_full,
    Woodiness = formula_map$Woodiness_full,
    Phylogeny = formula_map$Phylogeny_full
  )
  
  T_grid <- seq(
    min(reference_experiments$T_z), max(reference_experiments$T_z),
    length.out = 65L
  )
  P_grid <- seq(
    min(reference_experiments$P_z), max(reference_experiments$P_z),
    length.out = 65L
  )
  
  make_panel_data <- function(level_row, warming_index) {
    model_key <- level_row$model_key[[1]]
    model <- models[[model_key]]
    model_data <- model$data
    warming_level <- warming_table$level[[warming_index]]
    W_raw <- warming_table$raw_value[[warming_index]]
    W_z <- warming_table$model_value[[warming_index]]
  
    if (model_key == "Graminoid") {
      support_data <- model_data |>
        dplyr::filter(
          as.character(.data$GraminoidStatus) == level_row$factor_value[[1]]
        )
      support_scope <- "group-specific observed three-dimensional climate hull"
    } else if (model_key == "Woodiness") {
      support_data <- model_data |>
        dplyr::filter(
          as.character(.data$Woodiness) == level_row$factor_value[[1]]
        )
      support_scope <- "group-specific observed three-dimensional climate hull"
    } else {
      support_data <- model_data
      support_scope <- "model-specific observed three-dimensional climate hull"
    }
    support_experiments <- support_data |>
      dplyr::distinct(.data$experiment_id, .data$T_z, .data$P_z, .data$W_z)
    candidate_grid <- tidyr::expand_grid(
      T_z = T_grid,
      P_z = P_grid,
      W_z = W_z
    )
    support_result <- auc_support_mask(support_experiments, candidate_grid)
    if (nrow(support_result$grid) != nrow(candidate_grid) ||
        !identical(support_result$grid$support_method,
                   rep("three-dimensional convex hull", nrow(candidate_grid)))) {
      stop(
        "Climate support hull failed for ", level_row$predictor[[1]], " / ",
        level_row$predictor_level[[1]], ".",
        call. = FALSE
      )
    }
    grid <- support_result$grid |>
      dplyr::mutate(
        TP_z = .data$T_z * .data$P_z,
        PEV1_z = NA_real_,
        PEV2_z = NA_real_,
        PEV3_z = NA_real_
      )
  
    if (model_key %in% c("LDMC", "LMA", "Height")) {
      grid$Trait_z <- level_row$predictor_model_value[[1]]
    } else if (model_key == "Graminoid") {
      grid$GraminoidStatus <- factor(
        level_row$factor_value[[1]],
        levels = levels(model_data$GraminoidStatus)
      )
    } else if (model_key == "Woodiness") {
      grid$Woodiness <- factor(
        level_row$factor_value[[1]],
        levels = levels(model_data$Woodiness)
      )
    } else if (model_key == "Phylogeny") {
      grid$PEV1_z <- level_row$predictor_model_value[[1]]
      grid$PEV2_z <- 0
      grid$PEV3_z <- 0
    }
  
    design <- auc_design_matrix(
      model_formula[[model_key]], grid, names(stats::coef(model))
    )
    population_prediction <- as.numeric(design %*% stats::coef(model))
    panel_column <- if (warming_index == 1L) {
      match(level_row$predictor_level[[1]],
            level_specs$predictor_level[level_specs$predictor == level_row$predictor[[1]]])
    } else {
      2L + match(level_row$predictor_level[[1]],
                 level_specs$predictor_level[level_specs$predictor == level_row$predictor[[1]]])
    }
    common <- list(
      predictor = level_row$predictor[[1]],
      predictor_level = level_row$predictor_level[[1]],
      predictor_raw_value = level_row$predictor_raw_value[[1]],
      predictor_model_value = level_row$predictor_model_value[[1]],
      predictor_raw_unit = level_row$predictor_raw_unit[[1]],
      predictor_type = level_row$predictor_type[[1]],
      warming_level = warming_level,
      W_raw = W_raw,
      W_z = W_z,
      panel_column = panel_column,
      support_scope = support_scope
    )
    predictions <- grid |>
      dplyr::mutate(
        row_type = "prediction_grid",
        predictor = common$predictor,
        predictor_level = common$predictor_level,
        predictor_raw_value = common$predictor_raw_value,
        predictor_model_value = common$predictor_model_value,
        predictor_raw_unit = common$predictor_raw_unit,
        predictor_type = common$predictor_type,
        warming_level = common$warming_level,
        W_raw = common$W_raw,
        W_z = common$W_z,
        panel_column = common$panel_column,
        support_scope = common$support_scope,
        destination_temperature = .data$T_z * T_scale + T_center,
        destination_precipitation = exp(.data$P_z * P_scale + P_center),
        predicted_delta_auc = dplyr::if_else(
          .data$supported, .env$population_prediction, NA_real_
        ),
        observed_experiment_W_raw = NA_real_
      ) |>
      dplyr::select(
        "row_type", "predictor", "predictor_level", "predictor_raw_value",
        "predictor_model_value", "predictor_raw_unit", "predictor_type",
        "warming_level", "W_raw", "W_z", "panel_column",
        "destination_temperature", "destination_precipitation",
        "T_z", "P_z", "TP_z", "PEV1_z", "PEV2_z", "PEV3_z",
        "predicted_delta_auc", "supported",
        "support_method", "support_scope", "observed_experiment_W_raw"
      )
    observed <- support_data |>
      dplyr::distinct(
        .data$experiment_id, .data$T, .data$P_raw, .data$W,
        .data$T_z, .data$P_z
      ) |>
      dplyr::transmute(
        row_type = "observed_experiment",
        predictor = common$predictor,
        predictor_level = common$predictor_level,
        predictor_raw_value = common$predictor_raw_value,
        predictor_model_value = common$predictor_model_value,
        predictor_raw_unit = common$predictor_raw_unit,
        predictor_type = common$predictor_type,
        warming_level = common$warming_level,
        W_raw = common$W_raw,
        W_z = common$W_z,
        panel_column = common$panel_column,
        destination_temperature = .data$T,
        destination_precipitation = .data$P_raw,
        T_z = .data$T_z,
        P_z = .data$P_z,
        TP_z = .data$T_z * .data$P_z,
        PEV1_z = NA_real_,
        PEV2_z = NA_real_,
        PEV3_z = NA_real_,
        predicted_delta_auc = NA_real_,
        supported = TRUE,
        support_method = "observed experiment location",
        support_scope = common$support_scope,
        observed_experiment_W_raw = .data$W,
        experiment_id = as.character(.data$experiment_id)
      )
    list(predictions = predictions, observed = observed)
  }
  
  panel_results <- lapply(seq_len(nrow(level_specs)), function(index) {
    level_row <- level_specs[index, ]
    lapply(seq_len(nrow(warming_table)), function(warming_index) {
      make_panel_data(level_row, warming_index)
    })
  }) |>
    unlist(recursive = FALSE)
  
  predictions <- dplyr::bind_rows(lapply(panel_results, `[[`, "predictions"))
  observed <- dplyr::bind_rows(lapply(panel_results, `[[`, "observed"))
  
  range_assessment <- predictions |>
    dplyr::filter(.data$supported, is.finite(.data$predicted_delta_auc)) |>
    dplyr::summarise(
      predicted_min = min(.data$predicted_delta_auc),
      predicted_max = max(.data$predicted_delta_auc),
      maximum_absolute_prediction = max(abs(.data$predicted_delta_auc)),
      .by = "predictor"
    )
  global_limit <- max(range_assessment$maximum_absolute_prediction)
  non_wood_limit <- max(
    range_assessment$maximum_absolute_prediction[
      range_assessment$predictor != "Woodiness"
    ]
  )
  wood_limit <- range_assessment$maximum_absolute_prediction[
    range_assessment$predictor == "Woodiness"
  ]
  wood_to_other_ratio <- wood_limit / non_wood_limit
  use_separate_woodiness_scale <- is.finite(wood_to_other_ratio) &&
    wood_to_other_ratio > 1.75
  
  # The atlas is a single faceted plot, so a common scale is used only when the
  # supported Woodiness range does not exceed the other rows by more than 75%.
  # Fail closed if a separate scale is scientifically warranted; the plotting
  # layout must then be deliberately revised rather than silently compressing it.
  if (use_separate_woodiness_scale) {
    stop(
      sprintf(
        paste(
          "Woodiness requires a separate colour scale (range ratio %.2f).",
          "Revise the atlas layout before plotting."
        ),
        wood_to_other_ratio
      ),
      call. = FALSE
    )
  }
  
  prediction_plot_data <- predictions |>
    dplyr::mutate(
      predictor = factor(.data$predictor, levels = predictor_order),
      predictor_display = factor(
        unname(predictor_display_labels[as.character(.data$predictor)]),
        levels = unname(predictor_display_labels)
      ),
      panel_column = factor(.data$panel_column, levels = 1:4)
    )
  observed_plot_data <- observed |>
    dplyr::mutate(
      predictor = factor(.data$predictor, levels = predictor_order),
      predictor_display = factor(
        unname(predictor_display_labels[as.character(.data$predictor)]),
        levels = unname(predictor_display_labels)
      ),
      panel_column = factor(.data$panel_column, levels = 1:4)
    )
  prediction_surface_plot_data <- prediction_plot_data |>
    dplyr::filter(
      .data$supported,
      is.finite(.data$predicted_delta_auc)
    )
  panel_labels <- prediction_plot_data |>
    dplyr::distinct(
      .data$predictor_display, .data$panel_column, .data$predictor_level
    ) |>
    dplyr::mutate(
      T_z = max(T_grid) - 0.035 * diff(range(T_grid)),
      P_z = max(P_grid) - 0.035 * diff(range(P_grid))
    )
  
  T_breaks <- pretty(range(T_grid), n = 3)
  T_breaks <- T_breaks[T_breaks >= min(T_grid) & T_breaks <= max(T_grid)]
  P_breaks <- pretty(range(P_grid), n = 3)
  P_breaks <- P_breaks[P_breaks >= min(P_grid) & P_breaks <= max(P_grid)]
  
  atlas_core <- ggplot2::ggplot(
    prediction_plot_data,
    ggplot2::aes(x = .data$T_z, y = .data$P_z)
  ) +
    ggplot2::geom_raster(
      data = prediction_surface_plot_data,
      ggplot2::aes(fill = .data$predicted_delta_auc)
    ) +
    ggplot2::geom_point(
      data = observed_plot_data,
      ggplot2::aes(x = .data$T_z, y = .data$P_z),
      inherit.aes = FALSE,
      shape = 21,
      fill = "white",
      colour = "grey25",
      alpha = 0.70,
      size = 1.05,
      stroke = 0.34
    ) +
    ggplot2::geom_label(
      data = panel_labels,
      ggplot2::aes(x = .data$T_z, y = .data$P_z, label = .data$predictor_level),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 1,
      size = 7.8,
      fontface = "bold",
      linewidth = 0,
      label.padding = grid::unit(0.10, "lines"),
      fill = grDevices::adjustcolor("white", alpha.f = 0.78)
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(.data$predictor_display),
      cols = ggplot2::vars(.data$panel_column)
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#E66101",
      mid = "white",
      high = "#1B9E77",
      midpoint = 0,
      limits = c(-global_limit, global_limit),
      na.value = "white",
      name = "Predicted winning-losing strength"
    ) +
    ggplot2::scale_x_continuous(
      breaks = T_breaks,
      labels = function(x) sprintf("%.1f", x * T_scale + T_center),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = P_breaks,
      labels = function(x) sprintf("%.0f", exp(x * P_scale + P_center)),
      expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::labs(
      x = "Destination temperature (°C)",
      y = "Destination precipitation (mm)"
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_colourbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = grid::unit(14.5, "cm"),
        barheight = grid::unit(0.70, "cm")
      )
    ) +
    ggplot2::theme_bw(base_size = 21) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.spacing = grid::unit(1.6, "mm"),
      panel.border = ggplot2::element_rect(colour = "grey35", linewidth = 0.65),
      strip.placement = "outside",
      strip.background.y = ggplot2::element_rect(fill = "grey88", colour = "grey35"),
      strip.text.y.right = ggplot2::element_text(
        angle = 0, face = "bold", size = 20, hjust = 0.5,
        lineheight = 0.92,
        margin = ggplot2::margin(t = 5, r = 6, b = 5, l = 6)
      ),
      strip.text.x = ggplot2::element_blank(),
      strip.background.x = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 25),
      axis.text = ggplot2::element_text(size = 17, colour = "black"),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 23, face = "bold"),
      legend.text = ggplot2::element_text(size = 19),
      plot.margin = ggplot2::margin(t = 27, r = 8, b = 8, l = 8)
    )
  
  atlas <- cowplot::ggdraw() +
    cowplot::draw_plot(atlas_core, x = 0, y = 0, width = 1, height = 0.982) +
    cowplot::draw_label(
      "Low warming",
      x = 0.255,
      y = 0.983,
      fontface = "bold",
      size = 25
    ) +
    cowplot::draw_label(
      "High warming",
      x = 0.660,
      y = 0.983,
      fontface = "bold",
      size = 25
    )
  
  plotting_export <- dplyr::bind_rows(predictions, observed) |>
    dplyr::arrange(
      factor(.data$predictor, levels = predictor_order),
      .data$panel_column,
      .data$row_type,
      .data$T_z,
      .data$P_z
    )
  
  validation <- tibble::tribble(
    ~check, ~pass, ~detail,
    "No model refit", TRUE, "Six existing canonical full-REML model RDS files were read; no fitting function is called",
    "Canonical model predictions", TRUE, "Predictions use model.matrix and fixed coefficients only",
    "Trait quantile weighting", all(vapply(trait_quantile_rows, function(x) all(x$n_unique == unique(x$n_unique)), logical(1))), "Q25/Q75 computed after distinct species_id within each fitted trait model",
    "Warming quantile weighting", nrow(reference_experiments) == 42L, "Q25/Q75 computed across 42 distinct experiment_id values",
    "Frozen climate scaling", climate_scaling_error <= 1e-10, sprintf("maximum stored-versus-frozen scaling error %.3g", climate_scaling_error),
    "Support-mask enforcement", all(is.finite(predictions$predicted_delta_auc[predictions$supported])) && all(is.na(predictions$predicted_delta_auc[!predictions$supported])), "Only cells inside the applicable observed three-dimensional climate-warming hull receive predictions; unsupported cells remain NA and cannot be drawn",
    "Woodiness group support", all(predictions$support_scope[predictions$predictor == "Woodiness"] == "group-specific observed three-dimensional climate hull"), "Non-woody and woody support hulls are constructed separately",
    "Canonical reconstructed PEV1", nrow(pev_check) == 526L && pev_standardization_error <= 1e-10, sprintf("526 species; maximum standardized-coordinate error %.3g", pev_standardization_error),
    "PEV2 and PEV3 reference values", all(predictions$PEV2_z[predictions$predictor == "Phylogenetic position (PEV1 shown)"] == 0) && all(predictions$PEV3_z[predictions$predictor == "Phylogenetic position (PEV1 shown)"] == 0), "PEV2_z = 0 and PEV3_z = 0 in every phylogeny prediction row",
    "Common response sign", TRUE, "Every model predicts auc_effect: warmed minus origin-control AUC",
    "Colour-scale assessment", !use_separate_woodiness_scale, sprintf("single zero-centred scale; Woodiness/other maximum-absolute range ratio %.3f", wood_to_other_ratio)
  )
  if (any(!validation$pass)) {
    stop(
      "Atlas validation failed: ",
      paste(validation$check[!validation$pass], collapse = ", "),
      call. = FALSE
    )
  }
  
  caption <- paste0(
    "Predicted change in area under the relative-cover curve (Δ AUC; warmed minus ",
    "origin-control) across destination summer climate. Predictions are ",
    "population-level fixed-effect values from the canonical full models. Low and ",
    "High continuous traits are the 25th and 75th percentiles across unique ",
    "trait-complete species: LDMC ", sprintf("%.4f", trait_quantile_rows[[1]]$raw_value[[1]]),
    " and ", sprintf("%.4f", trait_quantile_rows[[1]]$raw_value[[2]]), " g/g; LMA ",
    sprintf("%.5f", trait_quantile_rows[[2]]$raw_value[[1]]), " and ",
    sprintf("%.5f", trait_quantile_rows[[2]]$raw_value[[2]]), " 1/(mm² mg⁻¹); and ",
    "vegetative height ", sprintf("%.4f", trait_quantile_rows[[3]]$raw_value[[1]]),
    " and ", sprintf("%.4f", trait_quantile_rows[[3]]$raw_value[[2]]), " m. ",
    "Low and High warming are Q25 and Q75 across 42 distinct experiments ",
    "(ΔT = ", sprintf("%.2f", warming_raw[[1]]), " and ",
    sprintf("%.2f", warming_raw[[2]]), " °C). Models remain continuous; Low and ",
    "High values are representative prediction contrasts rather than fitted ",
    "categories. Points show observed experiment climates, and coloured areas are ",
    "restricted to observed three-dimensional temperature–precipitation–warming ",
    "support; support is calculated separately for Non-woody and Woody species. ",
    "PEV1 is shown as a representative phylogenetic gradient (Q25 = ",
    sprintf("%.3f", PEV1_raw[[1]]), "; Q75 = ", sprintf("%.3f", PEV1_raw[[2]]),
    ") while PEV2 and PEV3 are held at their standardized reference values of ",
    "zero. Phylogenetic inference is based on the joint PEV1–PEV3 model, and the ",
    "sign of PEV1 is not a universal biological direction. All panels share one ",
    "zero-centred colour scale because the supported Woodiness range did not ",
    "materially exceed the other predictor ranges."
  )
  
  dir.create(dirname(output_stem), recursive = TRUE, showWarnings = FALSE)
  invisible(lapply(
    unique(dirname(c(data_file, quantile_file, validation_file))),
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  ))
  readr::write_csv(plotting_export, data_file)
  readr::write_csv(quantile_table, quantile_file)
  readr::write_csv(validation, validation_file)
  writeLines(strwrap(caption, width = 100), caption_file)
  ggplot2::ggsave(
    pdf_file,
    plot = atlas,
    device = grDevices::cairo_pdf,
    width = 20,
    height = 20,
    units = "in",
    bg = "white"
  )
  ggplot2::ggsave(
    png_file,
    plot = atlas,
    width = 20,
    height = 20,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  invisible(lapply(
    c(pdf_file, png_file),
    function(path) {
      writeLines(
        format(Sys.time(), tz = "UTC", usetz = TRUE),
        paste0(path, ".generated")
      )
    }
  ))
  
  message("Low/high values used:")
  print(quantile_table, n = Inf)
  message("Prediction-range assessment:")
  print(range_assessment, n = Inf)
  message(sprintf(
    "Common zero-centred scale selected (limit ±%.4f; Woodiness ratio %.3f).",
    global_limit, wood_to_other_ratio
  ))
  message("Wrote: ", pdf_file)
  message("Wrote: ", png_file)
  message("Wrote: ", data_file)
  message("Wrote: ", quantile_file)
  message("Wrote: ", caption_file)
  message("Wrote: ", validation_file)
  invisible(c(pdf_file, png_file, data_file, quantile_file, caption_file, validation_file))
}
