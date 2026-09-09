#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
names <- if (length(args)) args else NULL

if (is.null(names)) {
  targets::tar_make()
} else {
  targets::tar_make(names = tidyselect::all_of(names))
}
