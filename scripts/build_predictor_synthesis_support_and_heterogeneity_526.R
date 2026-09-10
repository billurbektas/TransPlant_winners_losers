#!/usr/bin/env Rscript

# Standalone Figure 5 renderer. This reads canonical exported summary tables
# only; it does not invoke targets, load model objects, or refit models.

required_packages <- c("dplyr", "ggplot2", "patchwork", "readr", "tidyr")
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

source("R/functions/manuscript_predictor_synthesis.R")

arguments <- commandArgs(trailingOnly = TRUE)
output_stem <- if (length(arguments) >= 1L) {
  arguments[[1]]
} else {
  file.path(
    "figures", "candidates",
    "predictor_synthesis_support_and_heterogeneity_526"
  )
}
data_file <- if (length(arguments) >= 2L) {
  arguments[[2]]
} else {
  paste0(output_stem, "_data.csv")
}
pdf_file <- paste0(output_stem, ".pdf")
png_file <- paste0(output_stem, ".png")

plot_data <- prepare_predictor_synthesis_526(
  "output/tables/predictor_ML_fit_statistics_526.csv",
  "output/tables/predictor_ML_nested_comparisons_526.csv",
  "output/tables/heterogeneity_resolution_synthesis_526.csv"
)
figure <- build_predictor_synthesis_figure_526(plot_data)

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
