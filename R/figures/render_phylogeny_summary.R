# Module: R/figures/render_phylogeny_summary
# Responsibility: Implement the named scientific workflow stage or shared infrastructure.
# Inputs: Function arguments and explicitly declared target dependencies.
# Outputs: Validated in-memory objects and declared workflow artifacts.
# Side effects: Limited to functions whose contracts explicitly write, copy, or render files.

# Render the phylogenetic-coordinate summary appendix
#
# This renderer reads existing exported analysis artifacts and never fits a model.

#' Render the phylogenetic-coordinate summary appendix
#'
#' @param arguments Character vector of renderer-specific positional arguments.
#' @return A character vector of generated artifact paths, invisibly.
render_phylogeny_summary_appendix <- function(arguments = character()) {
  
  # Figure-only appendix summary of the canonical resident phylogenetic PCoA.
  # This script reads exported PEV scores and eigenvalues only; it does not
  # reconstruct the phylogeny, rerun PCoA, fit a model, or invoke tar_make().
  
  required_packages <- c("cowplot", "dplyr", "ggplot2", "ggrepel", "readr", "tibble")
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
      "figures", "manuscript", "appendix_phylogeny_pcoa_summary"
    )
  }
  output_root <- if (length(arguments) >= 2L) arguments[[2]] else "output"
  
  pdf_file <- paste0(output_stem, ".pdf")
  png_file <- paste0(output_stem, ".png")
  data_file <- file.path(
    output_root, "data", "appendix_phylogeny_pcoa_summary_data.csv"
  )
  axis_file <- file.path(
    output_root, "tables", "appendix_phylogeny_pcoa_first10_axes.csv"
  )
  validation_file <- file.path(
    output_root, "audits", "appendix_phylogeny_pcoa_summary_validation.csv"
  )
  
  input_files <- c(
    scores = "output/data/phylogeny_predictor_scores.csv",
    eigenvalues = "output/tables/pcoa_eigenvalue_summary.csv"
  )
  missing_inputs <- input_files[!file.exists(input_files)]
  if (length(missing_inputs) > 0L) {
    stop(
      "Required canonical PCoA export(s) missing: ",
      paste(missing_inputs, collapse = ", "),
      call. = FALSE
    )
  }
  
  scores <- readr::read_csv(input_files[["scores"]], show_col_types = FALSE)
  eigenvalues <- readr::read_csv(
    input_files[["eigenvalues"]], show_col_types = FALSE
  )
  required_scores <- c("species_id", "PEV1", "PEV2", "PEV3")
  required_eigenvalues <- c(
    "solution", "axis_number", "PEV", "positive_eigenvalue",
    "percent_positive_variation", "cumulative_percent"
  )
  if (!all(required_scores %in% names(scores)) ||
      !all(required_eigenvalues %in% names(eigenvalues))) {
    stop("Canonical PCoA exports do not contain the required columns.",
         call. = FALSE)
  }
  if (nrow(scores) != 526L || anyDuplicated(scores$species_id) ||
      any(!stats::complete.cases(scores[required_scores]))) {
    stop("Canonical PEV1--PEV3 scores must contain 526 complete unique species.",
         call. = FALSE)
  }
  
  selected_solution <- unique(eigenvalues$solution[
    eigenvalues$positive_eigenvalue %in% TRUE
  ])
  if (length(selected_solution) != 1L) {
    stop("A unique positive-eigenvalue PCoA solution was not found.",
         call. = FALSE)
  }
  first_ten <- eigenvalues |>
    dplyr::filter(
      .data$solution == .env$selected_solution,
      .data$positive_eigenvalue,
      .data$axis_number <= 10L
    ) |>
    dplyr::arrange(.data$axis_number) |>
    dplyr::mutate(
      retained_in_model = .data$axis_number <= 3L,
      axis_label = factor(.data$PEV, levels = .data$PEV)
    )
  if (!identical(as.integer(first_ten$axis_number), 1:10)) {
    stop("The selected PCoA solution does not contain axes 1--10.",
         call. = FALSE)
  }
  
  leverage_labels <- function(data, y_axis, n = 6L) {
    data |>
      dplyr::mutate(
        leverage = abs(.data$PEV1) / max(abs(.data$PEV1)) +
          abs(.data[[y_axis]]) / max(abs(.data[[y_axis]])),
        coordinate_group = paste(
          round(.data$PEV1), round(.data[[y_axis]]), sep = ":"
        )
      ) |>
      dplyr::arrange(dplyr::desc(.data$leverage), .data$species_id) |>
      dplyr::group_by(.data$coordinate_group) |>
      dplyr::slice_head(n = 1L) |>
      dplyr::ungroup() |>
      dplyr::slice_max(.data$leverage, n = n, with_ties = FALSE) |>
      dplyr::select("species_id") |>
      dplyr::pull(.data$species_id)
  }
  
  labels_pev2 <- leverage_labels(scores, "PEV2")
  labels_pev3 <- leverage_labels(scores, "PEV3")
  plotting_data <- scores |>
    dplyr::mutate(
      label_in_pev1_pev2 = .data$species_id %in% .env$labels_pev2,
      label_in_pev1_pev3 = .data$species_id %in% .env$labels_pev3
    )
  
  point_colour <- "#405A78"
  axis_colour <- "#6C5B7B"
  other_axis_colour <- "grey78"
  base_size <- 18
  
  make_scatter <- function(y_axis, label_column) {
    label_data <- plotting_data |>
      dplyr::filter(.data[[label_column]])
    ggplot2::ggplot(
      plotting_data,
      ggplot2::aes(x = .data$PEV1, y = .data[[y_axis]])
    ) +
      ggplot2::geom_hline(
        yintercept = 0, colour = "grey82", linewidth = 0.45
      ) +
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = 0.10)
      ) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = 0.10)
      ) +
      ggplot2::geom_vline(
        xintercept = 0, colour = "grey82", linewidth = 0.45
      ) +
      ggplot2::geom_point(
        colour = point_colour, alpha = 0.72, size = 2.1
      ) +
      ggrepel::geom_text_repel(
        data = label_data,
        ggplot2::aes(label = .data$species_id),
        colour = point_colour,
        fontface = "italic",
        size = 3.8,
        max.overlaps = Inf,
        box.padding = 0.35,
        point.padding = 0.15,
        min.segment.length = 0,
        segment.colour = "grey52",
        segment.size = 0.35,
        seed = 20260824L
      ) +
      ggplot2::labs(x = "PEV1", y = y_axis) +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_line(colour = "grey92"),
        axis.title = ggplot2::element_text(size = base_size + 2),
        axis.text = ggplot2::element_text(size = base_size - 2, colour = "black"),
        plot.margin = ggplot2::margin(8, 10, 8, 8)
      )
  }
  
  pev12_plot <- make_scatter(
    "PEV2", "label_in_pev1_pev2"
  )
  pev13_plot <- make_scatter(
    "PEV3", "label_in_pev1_pev3"
  )
  
  cumulative_first_three <- first_ten$cumulative_percent[
    first_ten$axis_number == 3L
  ]
  variation_plot <- ggplot2::ggplot(
    first_ten,
    ggplot2::aes(
      x = .data$axis_label,
      y = .data$percent_positive_variation,
      fill = .data$retained_in_model
    )
  ) +
    ggplot2::geom_col(width = 0.74, colour = "grey35", linewidth = 0.35) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = sprintf("%.1f%%", .data$percent_positive_variation)
      ),
      vjust = -0.35,
      size = 4.4,
      colour = "grey20"
    ) +
    ggplot2::geom_vline(
      xintercept = 3.5,
      linetype = "dashed",
      linewidth = 0.55,
      colour = axis_colour
    ) +
    ggplot2::annotate(
      "label",
      x = 6.8,
      y = max(first_ten$percent_positive_variation) * 0.93,
      label = sprintf(
        "PEV1–PEV3 cumulative variation: %.1f%%",
        cumulative_first_three
      ),
      size = 4.2,
      linewidth = 0,
      fill = grDevices::adjustcolor("white", alpha.f = 0.85),
      colour = axis_colour,
      fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(
      values = c(`TRUE` = axis_colour, `FALSE` = other_axis_colour),
      breaks = c(TRUE, FALSE),
      labels = c("PEV1–PEV3 used", "Other displayed axes"),
      name = NULL
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, max(first_ten$percent_positive_variation) * 1.12),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(
      x = "Phylogenetic eigenvector axis",
      y = "Phylogenetic variation explained (%)"
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "grey92"),
      axis.title = ggplot2::element_text(size = base_size + 2),
      axis.text = ggplot2::element_text(size = base_size - 2, colour = "black"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = base_size - 2),
      legend.key.size = grid::unit(0.55, "cm"),
      plot.margin = ggplot2::margin(8, 8, 8, 10)
    )
  
  figure <- cowplot::plot_grid(
    pev12_plot,
    pev13_plot,
    variation_plot,
    labels = c("a", "b", "c"),
    label_size = 24,
    label_fontface = "bold",
    label_x = 0.01,
    label_y = 0.995,
    hjust = 0,
    vjust = 1,
    nrow = 1,
    rel_widths = c(1, 1, 1.03),
    align = "hv",
    axis = "tblr"
  )
  
  validation <- tibble::tribble(
    ~check, ~pass, ~detail,
    "No phylogeny or PCoA rerun", TRUE, "Figure reads canonical exported scores and eigenvalues only",
    "Canonical species scope", nrow(scores) == 526L && !anyDuplicated(scores$species_id), "526 unique retained species",
    "Complete first three axes", all(stats::complete.cases(scores[c("PEV1", "PEV2", "PEV3")])), "PEV1–PEV3 are finite and complete for every species",
    "First ten positive axes", identical(as.integer(first_ten$axis_number), 1:10), "Axes 1–10 from the selected uncorrected positive-eigenvalue solution",
    "First three axes marked", identical(first_ten$retained_in_model, c(rep(TRUE, 3), rep(FALSE, 7))), sprintf("PEV1–PEV3 cumulatively explain %.6f%%", cumulative_first_three),
    "Species labels deterministic", length(labels_pev2) == 6L && length(labels_pev3) == 6L, "Six highest-leverage coordinate-distinct representative species labelled per scatter; fixed ggrepel seed"
  )
  if (any(!validation$pass)) {
    stop(
      "PCoA appendix validation failed: ",
      paste(validation$check[!validation$pass], collapse = ", "),
      call. = FALSE
    )
  }
  
  invisible(lapply(
    unique(c(dirname(output_stem), dirname(data_file), dirname(axis_file),
             dirname(validation_file))),
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  ))
  readr::write_csv(plotting_data, data_file)
  readr::write_csv(
    dplyr::select(first_ten, -"axis_label"), axis_file
  )
  readr::write_csv(validation, validation_file)
  ggplot2::ggsave(
    pdf_file,
    plot = figure,
    device = grDevices::cairo_pdf,
    width = 20,
    height = 7.5,
    units = "in",
    bg = "white"
  )
  ggplot2::ggsave(
    png_file,
    plot = figure,
    width = 20,
    height = 7.5,
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
  
  message("Wrote: ", pdf_file)
  message("Wrote: ", png_file)
  message("Wrote: ", data_file)
  message("Wrote: ", axis_file)
  message("Wrote: ", validation_file)
  invisible(c(pdf_file, png_file, data_file, axis_file, validation_file))
}
