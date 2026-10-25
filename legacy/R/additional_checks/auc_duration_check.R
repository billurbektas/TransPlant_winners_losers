#!/usr/bin/env Rscript

# AUC versus maximum observation duration; an additional descriptive check.
# Run from the repository root: Rscript R/additional_checks/auc_duration_check.R
# AUC accumulates over time: these cross-experiment trends are not causal
# effects of duration and are not within-experiment temporal trajectories.

library(targets)
library(dplyr)
library(ggplot2)

# 1. Resident effects and the original, individual-plot AUC values ----------
effects <- tar_read(resident_auc_data) |>
  transmute(effect_id, experiment_id = as.character(experiment_id),
            source_species = as.character(source_species),
            species_id = as.character(species_id), duration,
            delta_auc = auc_effect)

plot_auc <- tar_read(auc_inputs)$auc |>
  filter(horizon == "max", transform == "raw",
         ODT %in% c("origin_controls", "warmed")) |>
  mutate(experiment_id = as.character(experiment_id),
         source_species = as.character(species), ODT = as.character(ODT))

# Raw AUC names can differ from harmonized species names. Join on the original
# name; use species_id from the resident dataset when grouping species lines.
resident_plots <- plot_auc |>
  inner_join(effects, by = c("experiment_id", "source_species"),
             relationship = "many-to-one")
message("Resident plot AUC rows: ", nrow(resident_plots),
        "; non-finite AUC rows omitted: ", sum(!is.finite(resident_plots$auc)))
resident_plots <- resident_plots |> filter(is.finite(auc))

resident_plots$ODT <- factor(resident_plots$ODT,
  levels = c("origin_controls", "warmed"), labels = c("Origin controls", "Warmed"))
# 2. Signed and absolute Delta AUC together, then individual plot AUC -------
# Each point has equal weight within its page. The individual-plot trends
# therefore give effects with more replicate plots more weight.
# Curves use LOESS with default smoothing settings, without confidence bands.
# Species with two or three distinct durations retain their straight trend.
out_dir <- "output/additional_checks"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
pdf(file.path(out_dir, "auc_duration_check.pdf"), width = 11, height = 7)
for (trend in c("Straight", "Curved")) {
  delta_plots <- list()
  for (view in c("Delta AUC", "Absolute Delta AUC", "Individual plot AUC")) {
    if (view != "Individual plot AUC") {
      # Absolute Delta AUC measures response magnitude, regardless of sign.
      plot_data <- effects |>
        mutate(response = if (view == "Absolute Delta AUC") abs(delta_auc) else delta_auc,
               trend_label = if (trend == "Straight") "Linear" else "LOESS")
      line_data <- plot_data |> group_by(species_id) |>
        filter(n_distinct(duration) >= 2L) |>
        mutate(n_durations = n_distinct(duration)) |> ungroup()
    } else {
      plot_data <- resident_plots |> mutate(response = auc)
      line_data <- plot_data |> group_by(ODT, species_id) |>
        filter(n_distinct(duration) >= 2L) |>
        mutate(n_durations = n_distinct(duration)) |> ungroup()
    }
    p <- ggplot(plot_data, aes(duration, response)) +
      geom_point(alpha = 0.20, size = 1.1, colour = "grey35")
    if (trend == "Straight") {
      p <- p +
        geom_smooth(data = line_data, aes(group = species_id, colour = "Species trends"),
          method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.45) +
        geom_smooth(aes(colour = "Overall trend"),
          method = "lm", formula = y ~ x, se = FALSE, linewidth = 1.1)
      # Ordinary linear regression for the red line, separately by treatment.
      # The slope p-value assumes independent rows; it is not adjusted for
      # shared experiments, plots or species, so this remains descriptive.
      panels <- if (view == "Individual plot AUC") levels(plot_data$ODT) else "Overall"
      annotations <- list()
      for (panel in panels) {
        model_data <- if (view == "Individual plot AUC")
          plot_data[plot_data$ODT == panel, ] else plot_data
        duration_model <- lm(response ~ duration, data = model_data)
        model_summary <- summary(duration_model)
        slope_p <- model_summary$coefficients["duration", "Pr(>|t|)"]
        p_label <- if (slope_p < 0.001) "< 0.001" else
          paste("=", formatC(slope_p, format = "f", digits = 3))
        annotation <- data.frame(label = paste0(
          sprintf("y = %.3f %+.3f duration", coef(duration_model)[1], coef(duration_model)[2]),
          sprintf("\nR\u00b2 = %.3f", model_summary$r.squared),
          "\nDuration p ", p_label))
        if (view == "Individual plot AUC") annotation$ODT <- panel
        annotations[[panel]] <- annotation
      }
      p <- p + geom_text(data = bind_rows(annotations),
        aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
        hjust = -0.03, vjust = 1.15, colour = "#B54A43", size = 5)
    } else {
      p <- p +
        geom_smooth(data = filter(line_data, n_durations < 4L),
          aes(group = species_id, colour = "Species trends"),
          method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.45) +
        geom_smooth(data = filter(line_data, n_durations >= 4L),
          aes(group = species_id, colour = "Species trends"),
          method = "loess", se = FALSE) +
        geom_smooth(aes(colour = "Overall trend"),
          method = "loess", se = FALSE)
    }
    p <- p +
      scale_colour_manual(name = NULL,
        values = c("Overall trend" = "#B54A43", "Species trends" = "grey60"),
        breaks = c("Overall trend", "Species trends")) +
      scale_x_continuous(breaks = sort(unique(effects$duration))) +
      labs(x = "Maximum observation duration in years", y = view) +
      theme_classic(base_size = 16) +
      theme(axis.title = element_text(size = 18),
            axis.text = element_text(size = 15),
            strip.text = element_text(size = 17),
            legend.text = element_text(size = 16),
            legend.position = "bottom")
    if (view == "Individual plot AUC") p <- p + facet_wrap(~ ODT, nrow = 1)
    if (view != "Individual plot AUC") {
      delta_plots[[view]] <- p + facet_wrap(~ trend_label) +
        labs(x = "Maximum observation duration\nin years") +
        guides(colour = guide_legend(ncol = 1))
      if (view == "Absolute Delta AUC") {
        grid::grid.newpage()
        print(delta_plots[["Delta AUC"]], newpage = FALSE,
              vp = grid::viewport(x = 0.25, width = 0.5))
        print(delta_plots[["Absolute Delta AUC"]], newpage = FALSE,
              vp = grid::viewport(x = 0.75, width = 0.5))
      }
    } else {
      print(p)
    }
  }
}
dev.off()
message("Saved ", file.path(out_dir, "auc_duration_check.pdf"), " (four pages).")
