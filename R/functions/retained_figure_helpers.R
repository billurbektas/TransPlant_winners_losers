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

convert_pdf_to_png_atomic <- function(input, output, dpi = 300L, page = 1L) {
  executable <- Sys.which("pdftoppm")
  if (!nzchar(executable)) {
    stop("pdftoppm is required to create the PNG appendix export.", call. = FALSE)
  }
  if (!file.exists(input) || dir.exists(input)) {
    stop("PDF input does not exist: ", input, call. = FALSE)
  }
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  temporary_stem <- tempfile(
    pattern = paste0(tools::file_path_sans_ext(basename(output)), "."),
    tmpdir = dirname(output)
  )
  temporary_output <- paste0(temporary_stem, ".png")
  on.exit(if (file.exists(temporary_output)) unlink(temporary_output), add = TRUE)
  status <- system2(
    executable,
    c(
      "-png", "-f", as.character(page), "-l", as.character(page),
      "-singlefile", "-r", as.character(as.integer(dpi)), input,
      temporary_stem
    ),
    stdout = FALSE,
    stderr = FALSE
  )
  if (!identical(status, 0L) || !file.exists(temporary_output)) {
    stop("Could not convert PDF figure to PNG: ", input, call. = FALSE)
  }
  if (file.exists(output)) unlink(output)
  if (!file.rename(temporary_output, output)) {
    stop("Could not publish converted PNG atomically: ", output, call. = FALSE)
  }
  receipt <- paste0(output, ".generated")
  receipt_temporary <- tempfile(
    pattern = paste0(basename(receipt), "."), tmpdir = dirname(output)
  )
  writeLines(format(Sys.time(), tz = "UTC", usetz = TRUE), receipt_temporary)
  if (!file.rename(receipt_temporary, receipt)) {
    stop("Could not save figure-generation receipt: ", receipt, call. = FALSE)
  }
  output
}
