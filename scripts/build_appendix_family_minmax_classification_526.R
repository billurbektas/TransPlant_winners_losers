#!/usr/bin/env Rscript

# Figure-only family-resolution analogue of manuscript Figure 4. Existing
# canonical model and table exports are read; no model is fitted.

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
source("R/functions/auc_model_helpers.R")
source("R/functions/manuscript_family_context.R")

arguments <- commandArgs(trailingOnly = TRUE)
output_stem <- if (length(arguments)) {
  arguments[[1]]
} else {
  file.path(
    "figures", "manuscript", "appendix_family_minmax_classification_526"
  )
}
climate_scaling_file <- if (length(arguments) >= 2L) {
  arguments[[2]]
} else {
  "output/data/downstream_climate_scaling_526.rds"
}
output_root <- if (length(arguments) >= 3L) arguments[[3]] else "output"
pdf_file <- paste0(output_stem, ".pdf")
png_file <- paste0(output_stem, ".png")
selection_file <- file.path(
  output_root, "tables",
  "appendix_family_minmax_classification_526_selected_families.csv"
)
validation_file <- file.path(
  output_root, "audits",
  "appendix_family_minmax_classification_526_validation.csv"
)

input_files <- c(
  model = "output/models/M_family_aggregated_526.rds",
  data = "output/data/family_auc_data_526.csv",
  classification = "output/tables/family_minmax_classification_526.csv",
  random = "output/tables/family_random_coefficients_526.csv"
)
if (any(!file.exists(input_files))) {
  stop("Missing canonical family input(s): ",
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
  targets::tar_read_raw("downstream_climate_scaling_526")
}

selected_families <- c("Fabaceae", "Poaceae", "Rubiaceae", "Rosaceae")
bundle <- prepare_family_context_callouts_526(
  fit, data, classification, random_coefficients,
  selected_families = selected_families
)
figure <- build_family_context_callout_figure_526(
  bundle, climate_scaling,
  classification_base_size = 14,
  classification_label_size = 4.8,
  inset_base_size = 14
)

validation <- tibble::tribble(
  ~check, ~pass, ~detail,
  "No model refit", inherits(fit, "rma.mv"), "Existing canonical family REML model read from disk; no fitting function is called",
  "Canonical family scope", nrow(classification) == 33L, "33 eligible full-rank families in canonical classification export",
  "Canonical classes", identical(sort(unique(classification$fitted_class)), sort(c("Consistent loser", "Context-dependent", "Consistent winner"))), "All three canonical family response classes are present",
  "Canonical callouts", identical(as.character(bundle$selection$family_id), selected_families), "Four existing selected-for-surface families used",
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
