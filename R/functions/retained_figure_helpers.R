# Shared helpers for canonical and retained publication figures.

category_palette <- c(
  loser = "#E66101",
  non_significant_loser = "#FDBB84",
  no_change = "#D9D9D9",
  non_significant_winner = "#99D8C9",
  winner = "#1B9E77",
  colonizer = "#5E3C99",
  non_significant_colonizer = "#B7A6D8"
)

category_display_labels <- c(
  loser = "Loser",
  winner = "Winner",
  colonizer = "Colonizer",
  non_significant_loser = "Loser\n(non-significant)",
  non_significant_winner = "Winner\n(non-significant)",
  non_significant_colonizer = "Colonizer\n(non-significant)",
  no_change = "No significant\nchange"
)

clean_figure_label <- function(x) {
  stringr::str_replace_all(as.character(x), "_", " ")
}

transplant_figure_theme <- function(base_size = 15) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 2),
      strip.text = ggplot2::element_text(size = base_size - 1, face = "bold"),
      legend.text = ggplot2::element_text(size = base_size - 2),
      legend.title = ggplot2::element_text(size = base_size - 1)
    )
}

save_plot_atomic <- function(plot, path, width, height, dpi = 300, bg = NULL) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  extension <- tools::file_ext(path)
  temporary <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path), fileext = paste0(".", extension))
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  ggplot2::ggsave(
    temporary, plot = plot, width = width, height = height, dpi = dpi, bg = bg
  )
  if (!file.rename(temporary, path)) stop("Could not save figure atomically: ", path, call. = FALSE)
  receipt <- paste0(path, ".generated")
  receipt_temporary <- tempfile(pattern = paste0(basename(receipt), "."), tmpdir = dirname(path))
  writeLines(format(Sys.time(), tz = "UTC", usetz = TRUE), receipt_temporary)
  if (!file.rename(receipt_temporary, receipt)) {
    stop("Could not save figure-generation receipt: ", receipt, call. = FALSE)
  }
  path
}

copy_figure_atomic_verified <- function(input, output, expected_sha256) {
  input_hash <- sha256_file(input)
  if (!identical(input_hash, expected_sha256)) {
    stop(
      "Figure source checksum mismatch: ", input,
      "\nExpected: ", expected_sha256,
      "\nObserved: ", input_hash,
      call. = FALSE
    )
  }
  copied <- copy_file_atomic_verified(input, output)
  output_hash <- sha256_file(copied)
  if (!identical(output_hash, expected_sha256)) {
    stop("Published figure checksum differs from its verified source: ", output,
         call. = FALSE)
  }
  receipt <- paste0(output, ".generated")
  receipt_temporary <- tempfile(
    pattern = paste0(basename(receipt), "."), tmpdir = dirname(output)
  )
  writeLines(format(Sys.time(), tz = "UTC", usetz = TRUE), receipt_temporary)
  if (!file.rename(receipt_temporary, receipt)) {
    stop("Could not save figure-generation receipt: ", receipt, call. = FALSE)
  }
  copied
}
