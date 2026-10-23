# Module: R/figures/render_family_context
# Responsibility: Render the family-context appendix figure from exports.
# Inputs: Family model summaries and prediction surfaces.
# Outputs: Registered family appendix figure files.
# Side effects: Writes figure files.

# Render the family-context appendix
#
# This renderer reads existing exported analysis artifacts and never fits a model.

#' Render the family-context appendix
#'
#' @param arguments Character vector of renderer-specific positional arguments.
#' @return A character vector of generated artifact paths, invisibly.
render_family_context_appendix <- function(arguments = character()) {
  
  # Figure-only family-resolution analogue of manuscript Figure 4. Existing
  # current model and table exports are read; no model is fitted.
  
  required_packages <- c(
    "cowplot", "dplyr", "geometry", "ggplot2", "metafor", "readr",
    "targets", "tibble", "tidyr"
  )
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages)) {
    stop("Missing required package(s): ", paste(missing_packages, collapse = ", "),
         call. = FALSE)
  }
  if (!file.exists("_targets.R")) {
    stop("Run this script from the repository root.", call. = FALSE)
  }
  
  suppressPackageStartupMessages(library(metafor))
  
  output_stem <- if (length(arguments)) {
    arguments[[1]]
  } else {
    file.path(
      "figures", "manuscript", "appendix_family_minmax_classification"
    )
  }
  climate_scaling_file <- if (length(arguments) >= 2L) {
    arguments[[2]]
  } else {
    "output/data/resident_climate_scaling.rds"
  }
  output_root <- if (length(arguments) >= 3L) arguments[[3]] else "output"
  pdf_file <- paste0(output_stem, ".pdf")
  png_file <- paste0(output_stem, ".png")
  selection_file <- file.path(
    output_root, "tables",
    "appendix_family_minmax_classification_selected_families.csv"
  )
  validation_file <- file.path(
    output_root, "audits",
    "appendix_family_minmax_classification_validation.csv"
  )
  
  input_files <- c(
    model = "output/models/family_context_model.rds",
    data = "output/data/family_auc_data.csv",
    classification = "output/tables/family_context_classification.csv",
    random = "output/tables/family_random_coefficients.csv"
  )
  if (any(!file.exists(input_files))) {
    stop("Missing current family input(s): ",
         paste(input_files[!file.exists(input_files)], collapse = ", "),
         call. = FALSE)
  }
  fit <- readRDS(input_files[["model"]])
  data <- readr::read_csv(input_files[["data"]], show_col_types = FALSE)
  classification <- readr::read_csv(
    input_files[["classification"]], show_col_types = FALSE
  )
  random_coefficients <- readr::read_csv(
    input_files[["random"]], show_col_types = FALSE
  )
  climate_scaling <- if (file.exists(climate_scaling_file)) {
    readRDS(climate_scaling_file)
  } else {
    targets::tar_read_raw("resident_climate_scaling")
  }
  
  selected_families <- c("Fabaceae", "Poaceae", "Rubiaceae", "Rosaceae")
  bundle <- prepare_family_context_callouts(
    fit, data, classification, random_coefficients,
    selected_families = selected_families
  )
  figure <- build_family_context_callout_figure(
    bundle, climate_scaling,
    classification_base_size = 14,
    classification_label_size = 4.8,
    inset_base_size = 14
  )
  
  validation <- tibble::tribble(
    ~check, ~pass, ~detail,
    "No model refit", inherits(fit, "rma.mv"), "Existing current family REML model read from disk; no fitting function is called",
    "Current family scope", nrow(classification) == 33L, "33 eligible full-rank families in current classification export",
    "Current classes", identical(sort(unique(classification$fitted_class)), sort(c("Consistent loser", "Context-dependent", "Consistent winner"))), "All three current family response classes are present",
    "Current callouts", identical(as.character(bundle$selection$family_id), selected_families), "Four existing selected-for-surface families used",
    "Family support masks", all(bundle$surfaces$supported), "Every coloured callout cell lies inside that family's observed three-dimensional climate-warming hull",
    "Finite predictions", all(is.finite(bundle$surfaces$estimate)), "All displayed family predictions are finite"
  )
  if (any(!validation$pass)) {
    stop("Family appendix validation failed: ",
         paste(validation$check[!validation$pass], collapse = ", "),
         call. = FALSE)
  }
  
  dir.create(dirname(output_stem), recursive = TRUE, showWarnings = FALSE)
  invisible(lapply(
    unique(dirname(c(selection_file, validation_file))),
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  ))
  readr::write_csv(bundle$selection, selection_file, na = "")
  readr::write_csv(validation, validation_file, na = "")
  ggplot2::ggsave(
    pdf_file, figure, device = grDevices::cairo_pdf,
    width = 16, height = 12.5, units = "in", bg = "white"
  )
  ggplot2::ggsave(
    png_file, figure, width = 16, height = 12.5, units = "in",
    dpi = 300, bg = "white"
  )
  invisible(lapply(c(pdf_file, png_file), function(path) {
    writeLines(format(Sys.time(), tz = "UTC", usetz = TRUE), paste0(path, ".generated"))
  }))
  
  message("Wrote: ", pdf_file)
  message("Wrote: ", png_file)
  message("Wrote: ", selection_file)
  message("Wrote: ", validation_file)
  invisible(c(pdf_file, png_file, selection_file, validation_file))
}
