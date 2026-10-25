#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
if (!file.exists("_targets.R")) {
  stop("Run this script from the repository root.", call. = FALSE)
}

targets::tar_source(c(
  "R/infrastructure", "R/responses", "R/models", "R/figures"
))

if (!length(arguments) || identical(arguments[[1]], "--list")) {
  print(figure_renderer_registry(), row.names = FALSE)
  quit(save = "no", status = 0L)
}

render_figure_from_exports(arguments[[1]], arguments[-1L])
