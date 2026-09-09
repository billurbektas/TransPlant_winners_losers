# Step 13: neutral linked callouts with low/high warming in rows -------------
#
# This is a separate rendering variant of Step 12. It preserves the existing
# figures, loads the already fitted M_both, uses one neutral connector/frame
# colour, and stacks low and high warming vertically within every callout.

source("R/24_step_12_both_four_species_callout_surfaces.R")

low_high_lookup <- tibble::tibble(
  warming_name = c("Low warming", "High warming"),
  probability = c(0.25, 0.75),
  W_z = as.numeric(stats::quantile(
    experiment_warming$W_z, c(0.25, 0.75)
  )),
  W_C = as.numeric(stats::quantile(
    experiment_warming$W, c(0.25, 0.75)
  ))
)
warming_panel_levels <- paste0(
  low_high_lookup$warming_name, "\nDelta T = ",
  sprintf("%.2f", low_high_lookup$W_C), " °C"
)

surface_low_high <- stats::setNames(lapply(requested_species, function(species_name) {
  dplyr::bind_rows(lapply(seq_len(nrow(low_high_lookup)), function(index) {
    result <- surface_one(
      species_name, W_z_value = low_high_lookup$W_z[[index]]
    )
    result$grid |>
      dplyr::mutate(
        warming_panel = factor(
          warming_panel_levels[[index]], levels = warming_panel_levels
        )
      )
  }))
}), requested_species)

observations_low_high <- stats::setNames(lapply(requested_species, function(species_name) {
  observed <- auc_data |>
    dplyr::filter(.data$species_id == species_name) |>
    dplyr::distinct(.data$experiment_id, .data$T_z, .data$P_z)
  tidyr::crossing(
    observed,
    warming_panel = factor(
      warming_panel_levels, levels = warming_panel_levels
    )
  )
}), requested_species)

surface_limit_low_high <- max(abs(unlist(lapply(
  surface_low_high, function(x) x$estimate
))), na.rm = TRUE)
neutral_colour <- "#3F3F3F"
class_colours[["Context-dependent"]] <- "#0072B2"

# Raise Elymus so its connector exits above, rather than through, the 0.50
# y-axis tick. The previous rendering remains preserved in the Step 11 files.
inset_boxes[["Elymus repens"]][["y"]] <- 0.575

class_counts <- classification |>
  dplyr::count(.data$fitted_class, name = "n_species")
class_n <- function(class_name) {
  class_counts$n_species[as.character(class_counts$fitted_class) == class_name]
}

# Retain the original classification semantics for all ordinary points, but
# give the four highlighted callout points one uniform outline.
classification_neutral <- classification_plot
classification_neutral$layers[[length(classification_neutral$layers)]]$
  aes_params$colour <- neutral_colour
# The second rectangle is the context-dependent background inherited from
# Step 12. Replace its former purple (reserved for colonizers) with blue.
classification_neutral$layers[[2L]]$aes_params$fill <-
  class_colours[["Context-dependent"]]
classification_neutral <- classification_neutral +
  ggplot2::scale_colour_manual(
    values = class_colours, drop = FALSE, guide = "none"
  ) +
  ggplot2::guides(colour = "none") +
  ggplot2::annotate(
    "text",
    x = x_limits[[1]] + 0.16 * diff(x_limits),
    y = y_limits[[1]] + 0.09 * diff(y_limits),
    label = paste0("CONSISTENT LOSER\nn = ", class_n("Consistent loser")),
    colour = class_colours[["Consistent loser"]],
    fontface = "bold", size = 3.1, hjust = 0
  ) +
  ggplot2::annotate(
    "text",
    x = x_limits[[1]] + 0.05 * diff(x_limits),
    y = y_limits[[2]] - 0.07 * diff(y_limits),
    label = paste0(
      "CONTEXT-DEPENDENT\nn = ", class_n("Context-dependent")
    ),
    colour = class_colours[["Context-dependent"]],
    fontface = "bold", size = 3.1, hjust = 0
  ) +
  ggplot2::annotate(
    "text",
    x = x_limits[[2]] - 0.04 * diff(x_limits),
    y = y_limits[[2]] - 0.07 * diff(y_limits),
    label = paste0("CONSISTENT WINNER\nn = ", class_n("Consistent winner")),
    colour = class_colours[["Consistent winner"]],
    fontface = "bold", size = 3.1, hjust = 1
  )

# Connector origins are drawn in the classification plot's own data
# coordinates. This makes their origin mathematically identical to each
# highlighted point center, independent of axes, margins, or output size.
# Exact central-panel bounds on the final cowplot canvas. The lower boundary
# differs from Step 12 because this version removes the class legend/title
# material; using the old bound created visible breaks at the plot frame.
panel_box <- c(left = 0.294, right = 0.734, bottom = 0.231, top = 0.843)
callout_centres_y <- vapply(
  inset_boxes,
  function(box) box[["y"]] + box[["height"]] / 2,
  numeric(1)
)
callout_centres_y <- pmin(
  pmax(callout_centres_y, panel_box[["bottom"]]),
  panel_box[["top"]]
)
canvas_y_to_data <- function(canvas_y) {
  y_limits[[1]] +
    (canvas_y - panel_box[["bottom"]]) /
    (panel_box[["top"]] - panel_box[["bottom"]]) * diff(y_limits)
}
connector_data <- selection |>
  dplyr::transmute(
    species_id = .data$species_id,
    x = .data$min_fitted_delta_auc,
    y = .data$max_fitted_delta_auc,
    xend = ifelse(
      .data$species_id %in% c("Elymus repens", "Silene acaulis"),
      x_limits[[1]], x_limits[[2]]
    ),
    yend = vapply(
      .data$species_id,
      function(species_name) canvas_y_to_data(
        callout_centres_y[[species_name]]
      ),
      numeric(1)
    )
  )
classification_neutral <- classification_neutral +
  ggplot2::geom_segment(
    data = connector_data,
    ggplot2::aes(
      x = .data$x, y = .data$y,
      xend = .data$xend, yend = .data$yend
    ), inherit.aes = FALSE, colour = neutral_colour,
    linewidth = 0.85, lineend = "round"
  ) +
  ggplot2::geom_point(
    data = selection,
    ggplot2::aes(
      x = .data$min_fitted_delta_auc,
      y = .data$max_fitted_delta_auc,
      fill = .data$random_W
    ), inherit.aes = FALSE, shape = 21, size = 4.2,
    stroke = 1.25, colour = neutral_colour, show.legend = FALSE
  )

surface_plot_low_high <- function(species_name, show_legend = FALSE) {
  ggplot2::ggplot(
    surface_low_high[[species_name]],
    ggplot2::aes(x = .data$T_z, y = .data$P_z, fill = .data$estimate)
  ) +
    ggplot2::geom_tile(
      ggplot2::aes(width = .data$tile_width, height = .data$tile_height)
    ) +
    ggplot2::geom_point(
      data = observations_low_high[[species_name]],
      ggplot2::aes(x = .data$T_z, y = .data$P_z), inherit.aes = FALSE,
      shape = 21, fill = "white", colour = "grey15", size = 1.05,
      stroke = 0.3
    ) +
    ggplot2::facet_grid(rows = ggplot2::vars(.data$warming_panel)) +
    ggplot2::scale_fill_gradient2(
      low = "#E66101", mid = "white", high = "#1B9E77", midpoint = 0,
      limits = c(-surface_limit_low_high, surface_limit_low_high),
      guide = if (show_legend) {
        ggplot2::guide_colourbar(
          title.position = "top", direction = "horizontal",
          barwidth = grid::unit(35, "mm"), barheight = grid::unit(3, "mm")
        )
      } else {
        "none"
      }
    ) +
    ggplot2::scale_x_continuous(labels = raw_T_labels) +
    ggplot2::scale_y_continuous(labels = raw_P_labels) +
    ggplot2::theme_bw(base_size = 7.1) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text.y = ggplot2::element_text(face = "bold", size = 6.5),
      plot.title = ggplot2::element_text(
        face = "bold.italic", size = 9, colour = neutral_colour,
        hjust = 0.5
      ),
      plot.background = ggplot2::element_rect(
        colour = neutral_colour, linewidth = 1, fill = "white"
      ),
      plot.margin = ggplot2::margin(4, 4, 3, 4),
      axis.title = ggplot2::element_text(size = 6.7),
      axis.text = ggplot2::element_text(size = 6.1),
      legend.position = if (show_legend) "bottom" else "none"
    ) +
    ggplot2::labs(
      title = paste0(
        species_name, " (",
        dplyr::filter(selection, .data$species_id == species_name)$n_experiments,
        " experiments)"
      ),
      x = "Destination temperature (°C)",
      y = "Destination precipitation",
      fill = "Fitted AUC change"
    )
}

mini_plots <- stats::setNames(
  lapply(
    requested_species,
    function(species_name) surface_plot_low_high(species_name, TRUE) +
      ggplot2::theme(
        legend.title = ggplot2::element_text(size = 6.7),
        legend.text = ggplot2::element_text(size = 6.1),
        legend.margin = ggplot2::margin(-2, 0, -2, 0)
      )
  ),
  requested_species
)

# Use the same boxes and exact registered point centers as Step 12. Lines and
# frames are deliberately neutral: species identity comes from titles/links.
neutral_canvas <- cowplot::ggdraw() +
  cowplot::draw_plot(
    classification_neutral,
    x = class_box[["x"]], y = class_box[["y"]],
    width = class_box[["width"]], height = class_box[["height"]]
  )

for (species_name in names(inset_boxes)) {
  is_left <- species_name %in% c("Elymus repens", "Silene acaulis")
  start <- c(
    x = if (is_left) panel_box[["left"]] else panel_box[["right"]],
    y = callout_centres_y[[species_name]]
  )
  finish <- connector_target(species_name)
  neutral_canvas <- neutral_canvas + cowplot::draw_line(
    x = c(start[["x"]], finish[["x"]]),
    y = c(start[["y"]], finish[["y"]]),
    colour = neutral_colour, linewidth = 0.85, alpha = 0.9
  )
}
for (species_name in names(inset_boxes)) {
  box <- inset_boxes[[species_name]]
  neutral_canvas <- neutral_canvas + cowplot::draw_plot(
    mini_plots[[species_name]],
    x = box[["x"]], y = box[["y"]],
    width = box[["width"]], height = box[["height"]]
  )
}
neutral_canvas <- neutral_canvas

output_stem <- file.path(
  figure_dir,
  "step_12_both_four_species_neutral_linked_callouts_low_high_rows"
)
ggplot2::ggsave(
  paste0(output_stem, ".pdf"), neutral_canvas,
  width = 16, height = 12.5, device = grDevices::cairo_pdf
)
ggplot2::ggsave(
  paste0(output_stem, ".png"), neutral_canvas,
  width = 16, height = 12.5, dpi = 190, bg = "white"
)

cat("M_both was loaded without refitting.\n")
cat("Low/high warming rows: ",
    paste0(low_high_lookup$warming_name, " = ",
           sprintf("%.2f", low_high_lookup$W_C), " deg C",
           collapse = "; "), "\n", sep = "")
cat("Neutral linked-callout figure: ",
    normalizePath(paste0(output_stem, ".pdf")), "\n", sep = "")
