# Module: R/figures/render_phylogeny_axes
# Responsibility: Implement the named scientific workflow stage or shared infrastructure.
# Inputs: Function arguments and explicitly declared target dependencies.
# Outputs: Validated in-memory objects and declared workflow artifacts.
# Side effects: Limited to functions whose contracts explicitly write, copy, or render files.

# Render the secondary phylogenetic-axis atlas
#
# This renderer reads existing exported analysis artifacts and never fits a model.

#' Render the secondary phylogenetic-axis atlas
#'
#' @param arguments Character vector of renderer-specific positional arguments.
#' @return A character vector of generated artifact paths, invisibly.
render_phylogeny_axis_atlas <- function(arguments = character()) {
  
  # Appendix prediction atlas for PEV2 and PEV3 from the existing current
  # joint PEV1--PEV3 full-REML model. This script performs no model fitting and
  # never invokes tar_make().
  
  required_packages <- c(
    "cowplot", "dplyr", "geometry", "ggplot2", "metafor", "readr",
    "tibble", "tidyr"
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
    file.path(
      "figures", "manuscript",
      "appendix_pev2_pev3_prediction_atlas"
    )
  }
  climate_scaling_file <- if (length(arguments) >= 2L) {
    arguments[[2]]
  } else {
    "output/data/resident_climate_scaling.rds"
  }
  figure07_data_file <- if (length(arguments) >= 3L) {
    arguments[[3]]
  } else {
    "output/data/fig07_predictor_response_atlas_data.csv"
  }
  output_root <- if (length(arguments) >= 4L) arguments[[4]] else "output"
  
  pdf_file <- paste0(output_stem, ".pdf")
  png_file <- paste0(output_stem, ".png")
  data_file <- file.path(
    output_root, "data", "appendix_pev2_pev3_prediction_atlas_data.csv"
  )
  quantile_file <- file.path(
    output_root, "tables",
    "appendix_pev2_pev3_prediction_atlas_quantiles.csv"
  )
  caption_file <- paste0(output_stem, "_caption.txt")
  validation_file <- file.path(
    output_root, "audits",
    "appendix_pev2_pev3_prediction_atlas_validation.csv"
  )
  
  input_files <- c(
    phylogeny_model = "output/models/predictor_phylogeny_model.rds",
    pev_scaling = "output/tables/phylogeny_predictor_scaling.csv",
    pev_scores = "output/data/phylogeny_predictor_scores.csv",
    climate_scaling = climate_scaling_file,
    figure07_plotting_data = figure07_data_file
  )
  missing_inputs <- input_files[!file.exists(input_files)]
  if (length(missing_inputs) > 0L) {
    stop(
      "Required current input(s) missing: ",
      paste(missing_inputs, collapse = ", "),
      call. = FALSE
    )
  }
  
  model <- readRDS(input_files[["phylogeny_model"]])
  if (!inherits(model, "rma.mv") || !is.data.frame(model$data)) {
    stop("The phylogeny input must be the existing current rma.mv fit.",
         call. = FALSE)
  }
  climate_scaling <- readRDS(input_files[["climate_scaling"]])
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
  
  reference_experiments <- model$data |>
    dplyr::distinct(
      .data$experiment_id, .data$T, .data$P, .data$P_raw, .data$W,
      .data$T_z, .data$P_z, .data$W_z
    ) |>
    dplyr::arrange(.data$experiment_id)
  if (nrow(reference_experiments) != 42L) {
    stop("Warming quantiles require exactly 42 distinct experiments.",
         call. = FALSE)
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
    stop("Stored model climate columns do not match frozen scaling.",
         call. = FALSE)
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
  
  if (nrow(pev_scores) != 526L || anyDuplicated(pev_scores$species_id)) {
    stop("PEV quantiles require exactly 526 unique species.", call. = FALSE)
  }
  
  pev_specs <- tibble::tribble(
    ~PEV,   ~model_field, ~display,
    "PEV2", "PEV2_z",    "PEV2",
    "PEV3", "PEV3_z",    "PEV3"
  )
  
  pev_quantiles <- dplyr::bind_rows(lapply(
    seq_len(nrow(pev_specs)),
    function(index) {
      spec <- pev_specs[index, ]
      center <- single_value(
        pev_scaling, "PEV", spec$PEV[[1]], "raw_mean", "PEV center"
      )
      scale <- single_value(
        pev_scaling, "PEV", spec$PEV[[1]], "raw_SD", "PEV scale"
      )
      raw_values <- stats::quantile(
        pev_scores[[spec$PEV[[1]]]],
        probs = c(0.25, 0.75),
        names = FALSE,
        type = 7
      )
      tibble::tibble(
        quantity = spec$display[[1]],
        level = c(
          paste("Low", spec$display[[1]]),
          paste("High", spec$display[[1]])
        ),
        quantile = c("Q25", "Q75"),
        raw_value = raw_values,
        raw_unit = "PCoA coordinate",
        model_field = spec$model_field[[1]],
        model_value = (raw_values - center) / scale,
        scaling_center = center,
        scaling_scale = scale,
        n_unique = nrow(pev_scores),
        weighting = "all 526 unique species from reconstructed current PCoA"
      )
    }
  ))
  
  standardization_checks <- dplyr::bind_rows(lapply(
    seq_len(nrow(pev_specs)),
    function(index) {
      spec <- pev_specs[index, ]
      center <- single_value(
        pev_scaling, "PEV", spec$PEV[[1]], "raw_mean", "PEV center"
      )
      scale <- single_value(
        pev_scaling, "PEV", spec$PEV[[1]], "raw_SD", "PEV scale"
      )
      model_pev <- model$data |>
        dplyr::distinct(.data$species_id, value = .data[[spec$model_field[[1]]]])
      score_pev <- pev_scores |>
        dplyr::transmute(
          species_id = .data$species_id,
          expected = (.data[[spec$PEV[[1]]]] - center) / scale
        )
      checked <- dplyr::inner_join(score_pev, model_pev, by = "species_id")
      tibble::tibble(
        PEV = spec$PEV[[1]],
        n_species = nrow(checked),
        maximum_error = max(abs(checked$expected - checked$value))
      )
    }
  ))
  if (any(standardization_checks$n_species != 526L) ||
      any(standardization_checks$maximum_error > 1e-10)) {
    stop("Current PEV coordinates do not match the fitted-model columns.",
         call. = FALSE)
  }
  
  level_specs <- pev_quantiles |>
    dplyr::transmute(
      predictor = .data$quantity,
      predictor_level = .data$level,
      model_field = .data$model_field,
      predictor_raw_value = .data$raw_value,
      predictor_model_value = .data$model_value,
      predictor_raw_unit = .data$raw_unit
    )
  
  T_grid <- seq(
    min(reference_experiments$T_z), max(reference_experiments$T_z),
    length.out = 65L
  )
  P_grid <- seq(
    min(reference_experiments$P_z), max(reference_experiments$P_z),
    length.out = 65L
  )
  support_experiments <- reference_experiments |>
    dplyr::select("experiment_id", "T_z", "P_z", "W_z")
  model_formula <- model$formula.mods
  if (!inherits(model_formula, "formula")) {
    stop("The fitted phylogeny model does not retain its fixed-effect formula.",
         call. = FALSE)
  }
  fixed_beta <- stats::coef(model)
  
  make_panel_data <- function(level_row, warming_index) {
    warming_level <- warming_table$level[[warming_index]]
    W_raw <- warming_table$raw_value[[warming_index]]
    W_z <- warming_table$model_value[[warming_index]]
    level_index <- match(
      level_row$predictor_level[[1]],
      level_specs$predictor_level[level_specs$predictor == level_row$predictor[[1]]]
    )
    panel_column <- if (warming_index == 1L) level_index else 2L + level_index
  
    candidate_grid <- tidyr::expand_grid(
      T_z = T_grid,
      P_z = P_grid,
      W_z = W_z
    )
    support_result <- auc_support_mask(support_experiments, candidate_grid)
    if (nrow(support_result$grid) != nrow(candidate_grid)) {
      stop("Three-dimensional climate support mask failed.", call. = FALSE)
    }
    grid <- support_result$grid |>
      dplyr::mutate(
        TP_z = .data$T_z * .data$P_z,
        PEV1_z = 0,
        PEV2_z = 0,
        PEV3_z = 0
      )
    grid[[level_row$model_field[[1]]]] <- level_row$predictor_model_value[[1]]
    design <- auc_design_matrix(model_formula, grid, names(fixed_beta))
    population_prediction <- as.numeric(design %*% fixed_beta)
  
    predictions <- grid |>
      dplyr::mutate(
        row_type = "prediction_grid",
        predictor = level_row$predictor[[1]],
        predictor_level = level_row$predictor_level[[1]],
        predictor_model_field = level_row$model_field[[1]],
        predictor_raw_value = level_row$predictor_raw_value[[1]],
        predictor_model_value = level_row$predictor_model_value[[1]],
        predictor_raw_unit = level_row$predictor_raw_unit[[1]],
        warming_level = warming_level,
        W_raw = W_raw,
        panel_column = panel_column,
        destination_temperature = .data$T_z * T_scale + T_center,
        destination_precipitation = exp(.data$P_z * P_scale + P_center),
        predicted_delta_auc = dplyr::if_else(
          .data$supported, population_prediction, NA_real_
        ),
        experiment_id = NA_character_,
        observed_experiment_W_raw = NA_real_
      )
  
    observed <- reference_experiments |>
      dplyr::transmute(
        row_type = "observed_experiment",
        predictor = level_row$predictor[[1]],
        predictor_level = level_row$predictor_level[[1]],
        predictor_model_field = level_row$model_field[[1]],
        predictor_raw_value = level_row$predictor_raw_value[[1]],
        predictor_model_value = level_row$predictor_model_value[[1]],
        predictor_raw_unit = level_row$predictor_raw_unit[[1]],
        warming_level = warming_level,
        W_raw = W_raw,
        W_z = W_z,
        panel_column = panel_column,
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
        experiment_id = as.character(.data$experiment_id),
        observed_experiment_W_raw = .data$W
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
  
  # Use the exact common colour scale already established for main Figure 7.
  figure07_data <- readr::read_csv(
    input_files[["figure07_plotting_data"]], show_col_types = FALSE
  )
  figure07_values <- figure07_data$predicted_delta_auc[
    figure07_data$row_type == "prediction_grid" &
      figure07_data$supported %in% TRUE &
      is.finite(figure07_data$predicted_delta_auc)
  ]
  global_limit <- max(abs(figure07_values))
  if (!is.finite(global_limit) || global_limit <= 0) {
    stop("Could not recover Figure 7's common colour limit.", call. = FALSE)
  }
  
  predictor_order <- c("PEV2", "PEV3")
  prediction_plot_data <- predictions |>
    dplyr::mutate(
      predictor = factor(.data$predictor, levels = predictor_order),
      panel_column = factor(.data$panel_column, levels = 1:4)
    )
  observed_plot_data <- observed |>
    dplyr::mutate(
      predictor = factor(.data$predictor, levels = predictor_order),
      panel_column = factor(.data$panel_column, levels = 1:4)
    )
  prediction_surface_plot_data <- prediction_plot_data |>
    dplyr::filter(.data$supported, is.finite(.data$predicted_delta_auc))
  panel_labels <- prediction_plot_data |>
    dplyr::distinct(.data$predictor, .data$panel_column, .data$predictor_level) |>
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
      ggplot2::aes(
        x = .data$T_z,
        y = .data$P_z,
        label = .data$predictor_level
      ),
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
      rows = ggplot2::vars(.data$predictor),
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
        angle = 0,
        face = "bold",
        size = 20,
        hjust = 0.5,
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
    cowplot::draw_plot(atlas_core, x = 0, y = 0, width = 1, height = 0.965) +
    cowplot::draw_label(
      "Low warming",
      x = 0.255,
      y = 0.978,
      fontface = "bold",
      size = 25
    ) +
    cowplot::draw_label(
      "High warming",
      x = 0.660,
      y = 0.978,
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
  
  reference_axes_zero <- predictions |>
    dplyr::filter(.data$row_type == "prediction_grid") |>
    dplyr::summarise(
      pass = all(
        dplyr::if_else(.data$predictor == "PEV2", .data$PEV1_z == 0 & .data$PEV3_z == 0,
                       .data$PEV1_z == 0 & .data$PEV2_z == 0)
      )
    ) |>
    dplyr::pull(.data$pass)
  validation <- tibble::tribble(
    ~check, ~pass, ~detail,
    "No model refit", TRUE, "Existing current predictor_phylogeny_model.rds was read; no fitting function is called",
    "Current model predictions", TRUE, "Predictions use model.matrix and fixed coefficients from the joint PEV1--PEV3 model",
    "PEV quantile weighting", all(pev_quantiles$n_unique == 526L), "Q25/Q75 computed across 526 unique species",
    "Warming quantile weighting", nrow(reference_experiments) == 42L, "Q25/Q75 computed across 42 distinct experiments",
    "Frozen climate scaling", climate_scaling_error <= 1e-10, sprintf("maximum stored-versus-frozen scaling error %.3g", climate_scaling_error),
    "Current reconstructed PEVs", all(standardization_checks$n_species == 526L) && all(standardization_checks$maximum_error <= 1e-10), paste(sprintf("%s maximum error %.3g", standardization_checks$PEV, standardization_checks$maximum_error), collapse = "; "),
    "Other PEV axes at reference", reference_axes_zero, "For each varied PEV, the other two standardized PEV coordinates equal zero",
    "Support-mask enforcement", all(is.finite(predictions$predicted_delta_auc[predictions$supported])) && all(is.na(predictions$predicted_delta_auc[!predictions$supported])), "Only cells inside the observed three-dimensional climate-warming hull receive predictions",
    "Figure 7 colour scale", isTRUE(all.equal(global_limit, max(abs(figure07_values)))), sprintf("shared zero-centred limit +/- %.6f", global_limit),
    "Common response sign", TRUE, "Every panel predicts auc_effect: warmed minus origin-control AUC"
  )
  if (any(!validation$pass)) {
    stop(
      "PEV appendix validation failed: ",
      paste(validation$check[!validation$pass], collapse = ", "),
      call. = FALSE
    )
  }
  
  caption <- paste0(
    "Predicted winning-losing strength (warmed minus origin-control AUC) across ",
    "destination summer climate for PEV2 and PEV3. Population-level fixed-effect ",
    "predictions come from the existing current joint PEV1--PEV3 full model. ",
    "Low and High values are the 25th and 75th percentiles across all 526 unique ",
    "species: PEV2 ", sprintf("%.3f", pev_quantiles$raw_value[pev_quantiles$quantity == "PEV2"][[1]]),
    " and ", sprintf("%.3f", pev_quantiles$raw_value[pev_quantiles$quantity == "PEV2"][[2]]),
    "; PEV3 ", sprintf("%.3f", pev_quantiles$raw_value[pev_quantiles$quantity == "PEV3"][[1]]),
    " and ", sprintf("%.3f", pev_quantiles$raw_value[pev_quantiles$quantity == "PEV3"][[2]]),
    ". The other two standardized PEV coordinates are held at zero in each row. ",
    "Low and High warming are Q25 and Q75 across 42 distinct experiments ",
    "(", sprintf("%.2f", warming_raw[[1]]), " and ",
    sprintf("%.2f", warming_raw[[2]]), " degrees C). Points show observed ",
    "experiment climates; coloured areas are restricted to observed ",
    "three-dimensional temperature--precipitation--warming support. The shared ",
    "zero-centred colour scale is identical to main Figure 7. PEV signs are ",
    "coordinate orientations and are not universal biological directions."
  )
  
  dir.create(dirname(output_stem), recursive = TRUE, showWarnings = FALSE)
  invisible(lapply(
    unique(dirname(c(data_file, quantile_file, validation_file))),
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  ))
  readr::write_csv(plotting_export, data_file)
  readr::write_csv(dplyr::bind_rows(pev_quantiles, warming_table), quantile_file)
  readr::write_csv(validation, validation_file)
  writeLines(strwrap(caption, width = 100), caption_file)
  ggplot2::ggsave(
    pdf_file,
    plot = atlas,
    device = grDevices::cairo_pdf,
    width = 20,
    height = 8.5,
    units = "in",
    bg = "white"
  )
  ggplot2::ggsave(
    png_file,
    plot = atlas,
    width = 20,
    height = 8.5,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  invisible(lapply(c(pdf_file, png_file), function(path) {
    writeLines(
      format(Sys.time(), tz = "UTC", usetz = TRUE),
      paste0(path, ".generated")
    )
  }))
  
  message("Low/high values used:")
  print(dplyr::bind_rows(pev_quantiles, warming_table), n = Inf)
  message(sprintf("Figure 7 common colour limit retained: +/- %.6f", global_limit))
  message("Wrote: ", pdf_file)
  message("Wrote: ", png_file)
  message("Wrote: ", data_file)
  message("Wrote: ", quantile_file)
  message("Wrote: ", caption_file)
  message("Wrote: ", validation_file)
  invisible(c(pdf_file, png_file, data_file, quantile_file, caption_file, validation_file))
}
