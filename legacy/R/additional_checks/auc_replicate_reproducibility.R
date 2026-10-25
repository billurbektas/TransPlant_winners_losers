#!/usr/bin/env Rscript

# Split replicate plots into two groups and calculate resident Delta AUC twice.
# Run: Rscript R/additional_checks/auc_replicate_reproducibility.R
library(targets)
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Existing observations and resident cohort ----------------------------
# Include the full resident cohort, including US_Arizona, as in the duration
# check. Keep the existing species selection, relative cover and max horizon.
community <- tar_read(analysis_community)
experiments <- tar_read(experiments)
effects <- tar_read(resident_auc_data) |>
  transmute(experiment_id = as.character(experiment_id), effect_id,
    source_species = as.character(source_species),
    species_id = as.character(species_id), delta_auc = auc_effect)
out_dir <- "output/additional_checks/auc_replicate_reproducibility"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20261005)
results <- list()
designs <- list()
assignments <- list()
reference_checks <- list()
experiment_ids <- sort(unique(effects$experiment_id))

for (i in seq_along(experiment_ids)) {
  id <- experiment_ids[[i]]
  message(i, "/", length(experiment_ids), ": ", id)
  experiment <- experiments[experiments$experiment_id == id, ]
  resident <- effects |> filter(experiment_id == id)
  # Reproduce the pipeline's treatment selection and minimum-observation rule.
  raw <- community |>
    filter(region == experiment$region, origin_site == experiment$origin_site) |>
    mutate(ODT = case_when(
      destination_site == experiment$origin_site & treatment == "LocalControl" ~ "origin_controls",
      destination_site == experiment$destination_site & treatment == "Warm" ~ "warmed")) |>
    filter(!is.na(ODT)) |>
    group_by(ODT, species) |> filter(n() >= 2L) |> ungroup() |>
    filter(!is.na(relative_year))

  # 2. Keep recorded blocks together -------------------------------------
  # A shared origin block links controls and transplants where recorded.
  # Destination blocks group plots only within their actual site: matching
  # numbers at different sites do NOT establish a pair. Missing block labels
  # leave plots separate; undocumented pairing cannot be recovered here.
  plots <- raw |>
    distinct(stable_plot_id, ODT, origin_site, destination_site,
             origin_block, destination_block) |>
    arrange(ODT, stable_plot_id)
  stopifnot(!anyDuplicated(plots$stable_plot_id))
  plots$unit <- seq_len(nrow(plots))
  block_keys <- bind_rows(
    plots |> transmute(stable_plot_id, site = origin_site, block = origin_block),
    plots |> transmute(stable_plot_id, site = destination_site, block = destination_block)) |>
    filter(!is.na(block), nzchar(trimws(block))) |>
    mutate(key = paste(site, block, sep = "::"))
  for (key in unique(block_keys$key)) {
    positions <- which(plots$stable_plot_id %in% block_keys$stable_plot_id[block_keys$key == key])
    linked_units <- unique(plots$unit[positions])
    plots$unit[plots$unit %in% linked_units] <- min(linked_units)
  }
  plots$unit <- match(plots$unit, unique(plots$unit))
  n_units <- max(plots$unit)
  unit_counts <- table(factor(plots$unit, levels = seq_len(n_units)),
                       factor(plots$ODT, levels = c("origin_controls", "warmed")))
  n_plots <- colSums(unit_counts)

  # 3. Annual cover matrices ---------------------------------------------
  # Observed absence is zero; an unsurveyed plot-year is NA, not absence.
  # The pipeline imputes missing years using other plots' treatment-year means.
  # Here that mean is calculated separately within each half, avoiding shared
  # imputed information between halves. Mean plot AUC equals trapezoidal AUC of
  # mean cover because integration is linear. No new species pools are selected.
  annual <- list()
  time_weights <- list()
  missing_plot_years <- 0L
  for (arm in c("origin_controls", "warmed")) {
    arm_plots <- plots$stable_plot_id[plots$ODT == arm]
    years <- sort(unique(raw$relative_year[raw$ODT == arm]))
    stopifnot(length(years) >= 2L)
    gaps <- diff(years)
    time_weights[[arm]] <- (c(0, gaps) + c(gaps, 0)) / 2
    annual[[arm]] <- list()
    for (year in years) {
      surveyed <- raw |> filter(ODT == .env$arm, relative_year == .env$year)
      values <- matrix(NA_real_, nrow = length(arm_plots), ncol = nrow(resident),
                       dimnames = list(arm_plots, resident$source_species))
      values[arm_plots %in% surveyed$stable_plot_id, ] <- 0
      observed <- surveyed |> filter(species %in% resident$source_species) |>
        group_by(stable_plot_id, species) |>
        summarise(cover = sum(relative_cover), .groups = "drop")
      values[cbind(match(observed$stable_plot_id, arm_plots),
                    match(observed$species, resident$source_species))] <- observed$cover
      missing_plot_years <- missing_plot_years + sum(is.na(values[, 1]))
      annual[[arm]][[as.character(year)]] <- values
    }
  }
  full_auc <- list()
  for (arm in names(annual)) {
    means <- vapply(annual[[arm]], colMeans, numeric(nrow(resident)), na.rm = TRUE)
    full_auc[[arm]] <- as.vector(means %*% time_weights[[arm]])
  }
  reconstructed <- full_auc$warmed - full_auc$origin_controls
  reference_checks[[id]] <- resident |>
    mutate(reconstructed_delta_auc = reconstructed,
           difference = reconstructed_delta_auc - delta_auc)
  # Essential: verify the same response is being split, not a new definition.
  stopifnot(all(is.finite(reconstructed)), max(abs(reconstructed - resident$delta_auc)) < 1e-8)

  # 4. Balanced partitions of whole blocks/plots --------------------------
  # Fix the first unit in A so exchanging A and B does not duplicate a split.
  # Enumerate small designs; sample larger designs. Retain at most 500 unique
  # partitions, balanced to within one largest block in each treatment.
  n_partitions <- 0L
  if (n_units >= 2L && all(n_plots >= 2L)) {
    if (n_units <= 14L) {
      candidates <- cbind(TRUE, as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), n_units - 1L))))
    } else {
      candidates <- cbind(TRUE, matrix(sample(c(FALSE, TRUE), 50000L * (n_units - 1L),
                                             replace = TRUE), ncol = n_units - 1L))
    }
    candidates <- unique(candidates)
    in_A <- candidates %*% unit_counts
    balanced <- apply(in_A > 0 & sweep(in_A, 2, n_plots, "<"), 1, all) &
      apply(abs(sweep(2 * in_A, 2, n_plots, "-")) <=
              matrix(apply(unit_counts, 2, max), nrow(in_A), 2, byrow = TRUE), 1, all)
    candidates <- candidates[balanced, , drop = FALSE]
    if (nrow(candidates)) {
      candidates <- candidates[sample.int(nrow(candidates), min(500L, nrow(candidates))), , drop = FALSE]
      n_partitions <- nrow(candidates)
    }
  }
  valid <- 0L
  missing_year_splits <- 0L
  for (split in seq_len(n_partitions)) {
    in_A <- as.logical(candidates[split, plots$unit])
    half_auc <- list()
    complete_years <- TRUE
    for (arm in names(annual)) {
      selected <- in_A[plots$ODT == arm]
      for (half in c("A", "B")) {
        keep <- if (half == "A") selected else !selected
        means <- matrix(NA_real_, nrow(resident), length(annual[[arm]]))
        for (year in seq_along(annual[[arm]])) {
          values <- annual[[arm]][[year]][keep, , drop = FALSE]
          if (all(is.na(values[, 1]))) complete_years <- FALSE
          means[, year] <- colMeans(values, na.rm = TRUE)
        }
        half_auc[[paste(arm, half)]] <- as.vector(means %*% time_weights[[arm]])
      }
    }
    if (!complete_years) {
      missing_year_splits <- missing_year_splits + 1L
      next
    }
    valid <- valid + 1L
    results[[paste(id, split)]] <- resident |>
      mutate(split = split,
        delta_A = half_auc[["warmed A"]] - half_auc[["origin_controls A"]],
        delta_B = half_auc[["warmed B"]] - half_auc[["origin_controls B"]])
    assignments[[paste(id, split)]] <- plots |>
      mutate(experiment_id = id, split = split, half = if_else(in_A, "A", "B"))
  }
  designs[[id]] <- data.frame(experiment_id = id,
    n_controls = n_plots[[1]], n_warmed = n_plots[[2]], n_split_units = n_units,
    plots_without_block = sum(!plots$stable_plot_id %in% block_keys$stable_plot_id),
    missing_plot_years = missing_plot_years, candidate_splits = n_partitions,
    valid_splits = valid, rejected_missing_year = missing_year_splits,
    status = if (valid) "Evaluated" else "No valid split: replication, blocks or missing years")
}

# 5. Descriptive agreement, not tests or confidence intervals -------------
splits <- bind_rows(results) |>
  mutate(difference = delta_A - delta_B,
         same_sign = sign(delta_A) == sign(delta_B))
split_scores <- splits |>
  group_by(experiment_id, split) |>
  summarise(n_effects = n(), RMSE = sqrt(mean(difference^2)),
    mean_absolute_difference = mean(abs(difference)),
    sign_agreement = mean(same_sign),
    pearson_r = if (n() > 1L && sd(delta_A) > 0 && sd(delta_B) > 0)
      cor(delta_A, delta_B) else NA_real_, .groups = "drop")
experiment_summary <- split_scores |>
  group_by(experiment_id) |>
  summarise(median_RMSE = median(RMSE),
    RMSE_10 = quantile(RMSE, .1), RMSE_90 = quantile(RMSE, .9),
    mean_sign_agreement = mean(sign_agreement),
    median_pearson_r = if (any(is.finite(pearson_r))) median(pearson_r, na.rm = TRUE) else NA_real_,
    .groups = "drop") |>
  left_join(bind_rows(designs), y = _, by = "experiment_id")
effect_summary <- splits |>
  group_by(experiment_id, effect_id, species_id, delta_auc) |>
  summarise(n_splits = n(), RMSE = sqrt(mean(difference^2)),
    mean_absolute_difference = mean(abs(difference)),
    sign_agreement = mean(same_sign), .groups = "drop")
write.csv(splits, file.path(out_dir, "split_effects.csv"), row.names = FALSE)
write.csv(split_scores, file.path(out_dir, "split_scores.csv"), row.names = FALSE)
write.csv(experiment_summary, file.path(out_dir, "experiment_summary.csv"), row.names = FALSE)
write.csv(effect_summary, file.path(out_dir, "effect_summary.csv"), row.names = FALSE)
write.csv(bind_rows(assignments), file.path(out_dir, "plot_assignments.csv"), row.names = FALSE)
write.csv(bind_rows(reference_checks), file.path(out_dir, "full_data_auc_check.csv"), row.names = FALSE)

# One example split per experiment: never plot averages of A/B across splits,
# since those would tend to agree just because plots repeatedly change halves.
example <- splits |> group_by(experiment_id) |> filter(split == min(split)) |> ungroup()
notes <- paste(
  "A and B use disjoint plots; recorded shared blocks remain together. Line = equal Delta AUC.",
  "Each half uses fewer plots than the full estimate. Missing-year means use only that half.",
  "Splits, species and experiments can be dependent; summaries are descriptive, not confidence intervals.",
  sep = "\n")
style <- theme_classic(base_size = 17) +
  theme(strip.text = element_text(size = 12), plot.caption = element_text(size = 11, hjust = 0),
        plot.tag = element_text(size = 23, face = "bold"), plot.tag.position = "topleft")
experiment_labels <- setNames(gsub("::", "\n", experiment_ids, fixed = TRUE), experiment_ids)
pdf(file.path(out_dir, "auc_replicate_reproducibility.pdf"), width = 16, height = 14)
print(ggplot(example, aes(delta_A, delta_B)) +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_abline(slope = 1, intercept = 0, colour = "#B54A43") +
  geom_point(alpha = .4, size = 1.5) + coord_equal() +
  labs(x = "Delta AUC from half A", y = "Delta AUC from half B", tag = "A",
    caption = paste("First valid split per experiment, chosen without examining agreement.", notes, sep = "\n")) + style)
for (page in split(sort(unique(example$experiment_id)),
                   ceiling(seq_along(unique(example$experiment_id)) / 12))) {
  print(ggplot(filter(example, experiment_id %in% page), aes(delta_A, delta_B)) +
    geom_abline(slope = 1, intercept = 0, colour = "#B54A43") +
    geom_point(alpha = .5, size = 1.2) + facet_wrap(~ experiment_id, ncol = 3,
      labeller = labeller(experiment_id = as_labeller(experiment_labels))) +
    coord_equal() + labs(x = "Delta AUC from half A", y = "Delta AUC from half B", caption = notes) + style)
}
print(ggplot(experiment_summary, aes(median_RMSE, reorder(experiment_id, median_RMSE))) +
  geom_segment(aes(x = RMSE_10, xend = RMSE_90, yend = reorder(experiment_id, median_RMSE)),
               colour = "grey65") + geom_point(size = 2, colour = "#B54A43") +
  labs(x = "Disagreement between halves: RMSE in Delta AUC units", y = NULL,
    caption = "Points = median across splits; lines = 10th to 90th percentiles across splits, not confidence intervals.") + style)
dev.off()
