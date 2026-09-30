# Module: R/figures/render_predictor_synthesis
# Responsibility: Implement the named scientific workflow stage or shared infrastructure.
# Inputs: Function arguments and explicitly declared target dependencies.
# Outputs: Validated in-memory objects and declared workflow artifacts.
# Side effects: Limited to functions whose contracts explicitly write, copy, or render files.

# Render the predictor-synthesis figure
#
# This renderer reads existing exported analysis artifacts and never fits a model.

#' Render the predictor-synthesis figure
#'
#' @param arguments Character vector of renderer-specific positional arguments.
#' @return A character vector of generated artifact paths, invisibly.
render_predictor_synthesis_figure <- function(arguments = character()) {
  
  # Standalone Figure 6 renderer. This reads exported summary tables
  # only; it does not invoke targets, load model objects, or refit models.
  
  required_packages <- c("cowplot", "dplyr", "ggplot2", "readr", "tidyr")
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
  
  
  output_stem <- if (length(arguments) >= 1L) {
    arguments[[1]]
  } else {
    file.path(
      "figures", "candidates",
      "predictor_synthesis_support_and_heterogeneity"
    )
  }
  data_file <- if (length(arguments) >= 2L) {
    arguments[[2]]
  } else {
    paste0(output_stem, "_data.csv")
  }
  pdf_file <- paste0(output_stem, ".pdf")
  png_file <- paste0(output_stem, ".png")
  fit_file <- if (length(arguments) >= 3L) arguments[[3]] else
    "output/tables/predictor_model_fit_statistics.csv"
  comparison_file <- if (length(arguments) >= 4L) arguments[[4]] else
    "output/tables/predictor_model_comparisons.csv"
  heterogeneity_file <- if (length(arguments) >= 5L) arguments[[5]] else
    "output/tables/heterogeneity_synthesis.csv"
  
  plot_data <- prepare_predictor_synthesis(
    fit_file, comparison_file, heterogeneity_file
  )
  figure <- build_predictor_synthesis_figure(plot_data)
  
  dir.create(dirname(output_stem), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(data_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(plot_data, data_file)
  ggplot2::ggsave(
    filename = pdf_file,
    plot = figure,
    device = grDevices::cairo_pdf,
    width = 16,
    height = 8.5,
    units = "in"
  )
  ggplot2::ggsave(
    filename = png_file,
    plot = figure,
    width = 16,
    height = 8.5,
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
  
  message("Wrote: ", pdf_file)
  message("Wrote: ", png_file)
  message("Wrote: ", data_file)
  invisible(c(pdf_file, png_file, data_file))
}
