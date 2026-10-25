#!/usr/bin/env Rscript

# Two checks: abundance alone, and abundance measured in separate control plots.
# Include every resident species, retaining the existing AUCs and trait subsets.
# Run: Rscript R/additional_checks/random_forest_abundance_checks_all_species.R
# Species-inclusive models predict only species represented in training.
# Check 1 retains both weighting and species-ID options.
# Check 2 uses unweighted forests without species ID, as agreed.
# US_Arizona is excluded from all training data and all performance scores.
# Known species: only held-out observations of species present in that fold's
# training data; unseen species are excluded. 

library(targets)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ranger)

# 1. Data and model formulas ----------------------------------------------
# T_z = standardized destination temperature; P_z = standardized log destination precipitation;
# W_z = standardized warming. Trait_z is the existing standardized log trait.
# Keep the same transformations and complete-case rows as the mixed models.
traits <- tar_read(predictor_specifications)
datasets <- list(Climate = tar_read(resident_auc_data), LDMC = traits$LDMC$data,
                 LMA = traits$LMA$data, Height = traits$Height$data)
# Remove only Arizona; do not require a minimum number of occupied gradients.
for (model in names(datasets)) {
  datasets[[model]] <- datasets[[model]] |>
    filter(as.character(gradient_id) != "US_Arizona")
}
# Keep the original inputs for checking the existing climate/trait forests.
baseline_datasets <- datasets
cohort_counts <- list()
for (model in names(datasets)) {
  data <- datasets[[model]]
  cohort_counts[[model]] <- data.frame(model, n_effects = nrow(data),
    n_species = n_distinct(data$species_id), n_experiments = n_distinct(data$experiment_id),
    n_gradients = n_distinct(data$gradient_id))
}
# 2. Focal-species abundance in the first recorded origin-control year -------
# Pick the year for the whole experiment, NOT each species' first occurrence.
# Average over all surveyed origin-control plots, including species absences
# as zero. Missing survey years are not filled. These are first observations,
# not necessarily pre-treatment abundances.
# Use the original species name to match the source of its AUC contrast.
# Existing relative_cover already uses the full-community denominator: do not
# renormalise it for the focal species.
community <- tar_read(analysis_community)
experiments <- tar_read(experiments) |>
  filter(included, experiment_id %in% datasets$Climate$experiment_id)
first_control_plots <- list()
for (i in seq_len(nrow(experiments))) {
  experiment <- experiments[i, ]
  controls <- community |>
    filter(region == experiment$region,
           origin_site == experiment$origin_site,
           destination_site == experiment$origin_site,
           treatment == "LocalControl")
  controls <- controls |> filter(year == min(year))
  plots <- controls |> distinct(plot_year_id, stable_plot_id, year, relative_year)
  focal_species <- datasets$Climate |>
    filter(as.character(experiment_id) == as.character(experiment$experiment_id)) |>
    select(effect_id, experiment_id, species_id, source_species)
  observed_cover <- controls |>
    group_by(plot_year_id, source_species = species) |>
    summarise(relative_abundance = sum(relative_cover), .groups = "drop")
  first_control_plots[[i]] <- focal_species |>
    cross_join(plots) |>
    left_join(observed_cover, by = c("plot_year_id", "source_species"),
              relationship = "many-to-one") |>
    mutate(relative_abundance = coalesce(relative_abundance, 0))
}
first_control_plots <- bind_rows(first_control_plots)
first_control_abundance <- first_control_plots |>
  group_by(effect_id, experiment_id, species_id, source_species) |>
  summarise(first_control_year = first(year),
    first_control_relative_year = first(relative_year),
    n_control_plots = n(),
    n_plots_present = sum(relative_abundance > 0),
    first_origin_abundance = mean(relative_abundance), .groups = "drop")
stopifnot(!anyDuplicated(first_control_abundance$effect_id),
          nrow(first_control_abundance) == nrow(datasets$Climate),
          all(is.finite(first_control_abundance$first_origin_abundance)))
for (model in names(datasets)) {
  datasets[[model]] <- datasets[[model]] |>
    left_join(first_control_abundance |>
                select(effect_id, first_origin_abundance, first_control_year,
                       first_control_relative_year, n_control_plots, n_plots_present),
              by = "effect_id", relationship = "one-to-one")
}

out_dir <- "output/additional_checks/random_forest_abundance_checks_all_species"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
base_seed <- 20260914L
write.csv(first_control_abundance, file.path(out_dir, "first_year_origin_abundance.csv"),
          row.names = FALSE)

write.csv(bind_rows(cohort_counts), file.path(out_dir, "cohort_counts.csv"), row.names = FALSE)

# 3. Abundance alone, climate/trait, and climate/trait plus abundance -------
# Species ID is either included in every model of a comparison, or in none.
# No species-ID forest predicts a species absent from its training data.
formulas <- list(
  Climate = list(
    "Abundance only" = list(
      without_species = auc_effect ~ first_origin_abundance,
      with_species = auc_effect ~ first_origin_abundance + species_id),
    "Climate / trait" = list(
      without_species = auc_effect ~ T_z + P_z + W_z,
      with_species = auc_effect ~ T_z + P_z + W_z + species_id),
    "Climate / trait + abundance" = list(
      without_species = auc_effect ~ T_z + P_z + W_z + first_origin_abundance,
      with_species = auc_effect ~ T_z + P_z + W_z + first_origin_abundance + species_id)),
  LDMC = list(
    "Abundance only" = list(
      without_species = auc_effect ~ first_origin_abundance,
      with_species = auc_effect ~ first_origin_abundance + species_id),
    "Climate / trait" = list(
      without_species = auc_effect ~ T_z + P_z + W_z + Trait_z,
      with_species = auc_effect ~ T_z + P_z + W_z + Trait_z + species_id),
    "Climate / trait + abundance" = list(
      without_species = auc_effect ~ T_z + P_z + W_z + Trait_z + first_origin_abundance,
      with_species = auc_effect ~ T_z + P_z + W_z + Trait_z + first_origin_abundance + species_id)),
  LMA = list(
    "Abundance only" = list(
      without_species = auc_effect ~ first_origin_abundance,
      with_species = auc_effect ~ first_origin_abundance + species_id),
    "Climate / trait" = list(
      without_species = auc_effect ~ T_z + P_z + W_z + Trait_z,
      with_species = auc_effect ~ T_z + P_z + W_z + Trait_z + species_id),
    "Climate / trait + abundance" = list(
      without_species = auc_effect ~ T_z + P_z + W_z + Trait_z + first_origin_abundance,
      with_species = auc_effect ~ T_z + P_z + W_z + Trait_z + first_origin_abundance + species_id)),
  Height = list(
    "Abundance only" = list(
      without_species = auc_effect ~ first_origin_abundance,
      with_species = auc_effect ~ first_origin_abundance + species_id),
    "Climate / trait" = list(
      without_species = auc_effect ~ T_z + P_z + W_z + Trait_z,
      with_species = auc_effect ~ T_z + P_z + W_z + Trait_z + species_id),
    "Climate / trait + abundance" = list(
      without_species = auc_effect ~ T_z + P_z + W_z + Trait_z + first_origin_abundance,
      with_species = auc_effect ~ T_z + P_z + W_z + Trait_z + first_origin_abundance + species_id)))
baseline_formulas <- list()
for (model in names(formulas)) {
  baseline_formulas[[model]] <- formulas[[model]][["Climate / trait"]]
}
baseline_signature <- digest::digest(list(
  data = baseline_datasets, formulas = baseline_formulas, seed = base_seed,
  trees = 500L, mtry = "floor(sqrt(p))", minimum_node_size = 5L,
  replace = TRUE, sample_fraction = 1, weights = NULL, splitrule = "variance",
  unordered_factors = "order", threads = 3L, prediction_rule_version = 1L,
  ranger_version = as.character(packageVersion("ranger"))))
baseline_dir <- "output/additional_checks/random_forest_without_US_Arizona/checkpoints/gradient"
first_signature <- digest::digest(list(datasets, formulas, base_seed,
  trees = 500L, node_size = 5L, mtry = "floor(sqrt(p))",
  ranger_version = as.character(packageVersion("ranger"))))
groups <- sort(unique(as.character(datasets$Climate$gradient_id)))
comparison_predictions <- list()
first_training <- list()
for (weighting in c("Unweighted", "Inverse variance")) {
  checkpoint_dir <- file.path(out_dir, "checkpoints", "abundance_comparison",
                              if (weighting == "Unweighted") "unweighted" else "weighted")
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  baseline_cache <- if (weighting == "Unweighted") baseline_dir else file.path(baseline_dir, "weighted")
  expected_baseline_signature <- if (weighting == "Unweighted") baseline_signature else
    digest::digest(list(baseline_signature, models = names(baseline_datasets),
                        sampling_weights = "1 / auc_se^2"))
  for (g in seq_along(c(groups, "Full_training"))) {
    group <- c(groups, "Full_training")[[g]]
    full_training <- group == "Full_training"
    checkpoint <- file.path(checkpoint_dir, paste0(group, ".rds"))
    if (file.exists(checkpoint)) {
      saved <- readRDS(checkpoint)
      stopifnot(identical(saved$signature, first_signature))
    } else {
      message("Abundance comparison: ", weighting, " / ", group)
      if (!full_training) {
        baseline_saved <- readRDS(file.path(baseline_cache, paste0(group, ".rds")))
        stopifnot(identical(baseline_saved$input_signature, expected_baseline_signature))
      }
      fits <- list()
      predictions <- list()
      training <- list()
      for (m in seq_along(datasets)) {
        model <- names(datasets)[[m]]
        data <- as.data.frame(datasets[[model]])
        data$species_id <- as.character(data$species_id)
        train <- if (full_training) data else data[data$gradient_id != group, ]
        test <- if (full_training) data else data[data$gradient_id == group, ]
        train$species_id <- factor(train$species_id)
        known <- test$species_id %in% levels(train$species_id)
        for (specification in names(formulas[[model]])) {
          for (s in seq_along(formulas[[model]][[specification]])) {
            method <- names(formulas[[model]][[specification]])[[s]]
            formula <- formulas[[model]][[specification]][[s]]
            predictors <- all.vars(formula)[-1L]
            seed <- if (full_training) base_seed + 10000L + 10L * m + s else
              base_seed + 100L * g + 10L * m + s
            key <- paste(model, method, sep = "__")
            if (specification == "Climate / trait") {
              if (full_training) {
                existing <- readRDS(file.path(baseline_cache, "full_data", paste0(key, ".rds")))
                expected <- digest::digest(list(expected_baseline_signature, key, seed,
                                                evaluation = "full_data_training_v1"))
                stopifnot(identical(existing$input_signature, expected))
              } else {
                existing <- baseline_saved$fits[[key]]
              }
              stopifnot(identical(as.character(train$effect_id), existing$training_effect_ids),
                        identical(all.vars(formula), all.vars(existing$formula)))
              fit <- existing$fit
            } else {
              fit <- ranger::ranger(
                formula, data = train[, all.vars(formula), drop = FALSE],
                num.trees = 500L, mtry = floor(sqrt(length(predictors))),
                min.node.size = 5L, splitrule = "variance", replace = TRUE,
                sample.fraction = 1,
                case.weights = if (weighting == "Inverse variance") 1 / train$auc_se^2 else NULL,
                respect.unordered.factors = "order", importance = "none",
                write.forest = TRUE, oob.error = FALSE, na.action = "na.fail",
                seed = seed, num.threads = 3L)
              fits[[paste(specification, key, sep = "__")]] <- list(fit = fit, formula = formula,
                training_effect_ids = as.character(train$effect_id),
                training_species = levels(train$species_id))
            }
            selected <- if (method == "with_species") known else rep(TRUE, nrow(test))
            test_rows <- test[selected, ]
            if (method == "with_species") {
              test_rows$species_id <- factor(test_rows$species_id, levels = levels(train$species_id))
            }
            predicted <- numeric()
            if (nrow(test_rows)) {
              predicted <- predict(fit, data = test_rows[, predictors, drop = FALSE],
                                   num.threads = 3L)$predictions
            }
            rows <- data.frame(specification = rep(specification, nrow(test_rows)),
              model = rep(model, nrow(test_rows)), method = rep(paste0("RF_", method), nrow(test_rows)),
              weighting = rep(weighting, nrow(test_rows)), heldout_group = rep(group, nrow(test_rows)),
              gradient_id = as.character(test_rows$gradient_id), effect_id = test_rows$effect_id,
              species_id = as.character(test_rows$species_id), known_species = known[selected],
              observed = test_rows$auc_effect, predicted = predicted)
            key <- paste(specification, key, sep = "__")
            if (full_training) {
              sizes <- table(rows$gradient_id)
              w <- 1 / (length(sizes) * as.numeric(sizes[rows$gradient_id]))
              training[[key]] <- data.frame(specification, model,
                method = paste0("RF_", method), weighting,
                R2_training_equal_gradient = 1 - sum(w * (rows$predicted - rows$observed)^2) /
                  sum(w * (rows$observed - sum(w * rows$observed))^2))
            } else {
              stopifnot(!any(rows$effect_id %in% train$effect_id), all(is.finite(predicted)))
              predictions[[key]] <- rows
            }
          }
        }
      }
      saved <- list(signature = first_signature, fits = fits,
                    predictions = bind_rows(predictions), training = bind_rows(training))
      saveRDS(saved, paste0(checkpoint, ".tmp"))
      stopifnot(file.rename(paste0(checkpoint, ".tmp"), checkpoint))
    }
    comparison_predictions[[paste(weighting, group)]] <- saved$predictions
    first_training[[paste(weighting, group)]] <- saved$training
  }
}
comparison_predictions <- bind_rows(comparison_predictions)
first_training <- bind_rows(first_training)
# Compare models with/without species ID on the same known-species test rows.
# Separately retain no-ID results for every held-out observation, including
# species absent from training. No fallback or switching between forests.
evaluation_predictions <- bind_rows(
  comparison_predictions |> filter(method == "RF_without_species") |>
    mutate(scope = "Every held-out species"),
  comparison_predictions |> filter(known_species) |>
    mutate(scope = "Species seen in training"))
first_fold_scores <- evaluation_predictions |>
  group_by(specification, weighting, model, method, scope, gradient_id) |>
  summarise(n_effects = n(), n_species = n_distinct(species_id),
    SSE = sum((predicted - observed)^2), SST = sum((observed - mean(observed))^2),
    MSE = mean((predicted - observed)^2), mean_y = mean(observed),
    mean_y2 = mean(observed^2), .groups = "drop") |>
  mutate(R2pred = if_else(SST > 0, 1 - SSE / SST, NA_real_))
# Preserve undefined/empty gradient rows, rather than substituting a zero score.
first_fold_scores <- first_fold_scores |>
  complete(nesting(specification, weighting, model, method, scope), gradient_id = groups,
           fill = list(n_effects = 0L, n_species = 0L))
first_overall <- first_fold_scores |>
  group_by(specification, weighting, model, method, scope) |>
  summarise(n_gradients = sum(n_effects > 0L),
    MSE = mean(MSE[n_effects > 0L]), mean_y = mean(mean_y[n_effects > 0L]),
    mean_y2 = mean(mean_y2[n_effects > 0L]), n_effects = sum(n_effects),
    .groups = "drop") |>
  mutate(R2pred = if_else(mean_y2 - mean_y^2 > 0, 1 - MSE / (mean_y2 - mean_y^2), NA_real_))
first_fold_counts <- comparison_predictions |>
  filter(specification == "Abundance only", method == "RF_without_species",
         weighting == "Unweighted") |>
  group_by(model, gradient_id) |>
  summarise(n_effects = n(), n_species = n_distinct(species_id),
    n_known_effects = sum(known_species),
    n_known_species = n_distinct(species_id[known_species]),
    n_excluded_from_species_id = sum(!known_species), .groups = "drop")
write.csv(first_fold_counts, file.path(out_dir, "fold_counts.csv"), row.names = FALSE)
write.csv(first_overall, file.path(out_dir, "abundance_only_overall.csv"), row.names = FALSE)
write.csv(first_fold_scores, file.path(out_dir, "abundance_only_by_gradient.csv"), row.names = FALSE)
write.csv(first_training, file.path(out_dir, "abundance_only_training.csv"), row.names = FALSE)
saveRDS(comparison_predictions, file.path(out_dir, "abundance_only_predictions.rds"))

# Species ranks within each held-out gradient -----------------------------
# Average observed and predicted Delta AUC across experiments for each species
# within a gradient, so that each species has one rank. Rank signed Delta AUC,
# not its absolute value. Ties receive average ranks.
# Spearman = correlation of the observed ranks with the predicted ranks.
# The summary is the mean correlation across gradients with defined ranks,
# giving those gradients equal weight. It is neither R2 nor a pooled correlation.
first_species_predictions <- evaluation_predictions |>
  group_by(specification, weighting, model, method, scope, gradient_id, species_id) |>
  summarise(n_effects = n(), observed = mean(observed), predicted = mean(predicted),
            .groups = "drop")
first_rank_by_gradient <- first_species_predictions |>
  group_by(specification, weighting, model, method, scope, gradient_id) |>
  summarise(n_species = n(),
    Spearman = if (n() >= 2L && n_distinct(observed) > 1L && n_distinct(predicted) > 1L)
      cor(observed, predicted, method = "spearman") else NA_real_, .groups = "drop")
first_rank_by_gradient <- first_rank_by_gradient |>
  complete(nesting(specification, weighting, model, method, scope), gradient_id = groups,
           fill = list(n_species = 0L))
first_rank_overall <- first_rank_by_gradient |>
  group_by(specification, weighting, model, method, scope) |>
  summarise(n_gradients = n(), n_defined = sum(is.finite(Spearman)),
    Spearman = if (any(is.finite(Spearman))) mean(Spearman, na.rm = TRUE) else NA_real_,
    .groups = "drop")
write.csv(first_species_predictions, file.path(out_dir, "species_rank_predictions.csv"),
          row.names = FALSE)
write.csv(first_rank_by_gradient, file.path(out_dir, "spearman_by_gradient.csv"),
          row.names = FALSE)
write.csv(first_rank_overall, file.path(out_dir, "spearman_overall.csv"), row.names = FALSE)

# 4. Check 2: independent origin-control plots -----------------------------
# The same 20 reproducible splits are used for all four predictor sets.
# A: first-recorded abundance only. B: control AUC in the response.
# Warmed AUC is unchanged. No A observations enter control-B AUC or its
# missing-year means. This avoids shared control observations, not every
# biological association or the mathematical role of control AUC in Delta AUC.
# No old SE weights: they use all controls. This agreed check is unweighted,
# without species ID, comparing climate/trait versus climate/trait + abundance.
split_formulas <- list(
  Climate = list(
    baseline = auc_effect ~ T_z + P_z + W_z,
    abundance = auc_effect ~ T_z + P_z + W_z + first_origin_abundance),
  LDMC = list(
    baseline = auc_effect ~ T_z + P_z + W_z + Trait_z,
    abundance = auc_effect ~ T_z + P_z + W_z + Trait_z + first_origin_abundance),
  LMA = list(
    baseline = auc_effect ~ T_z + P_z + W_z + Trait_z,
    abundance = auc_effect ~ T_z + P_z + W_z + Trait_z + first_origin_abundance),
  Height = list(
    baseline = auc_effect ~ T_z + P_z + W_z + Trait_z,
    abundance = auc_effect ~ T_z + P_z + W_z + Trait_z + first_origin_abundance))
n_splits <- 20L
split_seed <- 20261020L
cached_auc <- tar_read(auc_inputs)$auc |>
  filter(horizon == "max", transform == "raw",
         ODT %in% c("origin_controls", "warmed"))
cached_auc_means <- cached_auc |>
  group_by(experiment_id, source_species = species, ODT) |>
  summarise(mean_auc = mean(auc), .groups = "drop") |>
  pivot_wider(names_from = ODT, values_from = mean_auc)

# Keep shared control sites and all years of a physical plot together.
control_membership <- list()
for (i in seq_len(nrow(experiments))) {
  experiment <- experiments[i, ]
  control_membership[[i]] <- community |>
    filter(region == experiment$region, origin_site == experiment$origin_site,
           destination_site == experiment$origin_site, treatment == "LocalControl") |>
    distinct(stable_plot_id, origin_block, destination_block) |>
    mutate(control_site = paste(experiment$region, experiment$origin_site, sep = "::"))
}
control_plots <- bind_rows(control_membership) |> distinct() |>
  arrange(control_site, stable_plot_id)
stopifnot(!anyDuplicated(control_plots$stable_plot_id))
control_plots$unit <- seq_len(nrow(control_plots))
for (site in unique(control_plots$control_site)) {
  for (column in c("origin_block", "destination_block")) {
    for (block in unique(na.omit(control_plots[[column]][control_plots$control_site == site]))) {
      if (!nzchar(trimws(block))) next
      positions <- which(control_plots$control_site == site & control_plots[[column]] == block)
      linked <- unique(control_plots$unit[positions])
      control_plots$unit[control_plots$unit %in% linked] <- min(linked)
    }
  }
}

# Annual cover arrays reproduce the pipeline's response calculation.
# Surveyed absences are zero. Unsurveyed plot-years remain NA. Integrating
# annual means equals averaging plot AUCs after treatment-year mean imputation.
annual_controls <- list()
first_year_cover <- list()
response_metadata <- list()
full_auc_checks <- list()
experiment_ids <- sort(as.character(experiments$experiment_id))
for (id in experiment_ids) {
  experiment <- experiments[experiments$experiment_id == id, ]
  focal <- datasets$Climate |>
    filter(as.character(experiment_id) == id) |>
    select(effect_id, source_species, species_id, auc_effect) |>
    left_join(filter(cached_auc_means, experiment_id == id), by = "source_species",
              relationship = "one-to-one")
  raw <- community |>
    filter(region == experiment$region, origin_site == experiment$origin_site,
           destination_site == experiment$origin_site, treatment == "LocalControl")
  # Match the original >=2 observed rows per treatment-species rule for AUC.
  response_rows <- raw |> group_by(species) |> filter(n() >= 2L) |> ungroup() |>
    filter(!is.na(relative_year))
  plots <- sort(unique(response_rows$stable_plot_id))
  years <- sort(unique(response_rows$relative_year))
  gaps <- diff(years)
  integration_weights <- (c(0, gaps) + c(gaps, 0)) / 2
  annual <- list()
  for (year in years) {
    surveyed <- filter(response_rows, relative_year == .env$year)
    values <- matrix(NA_real_, length(plots), nrow(focal),
                     dimnames = list(plots, focal$effect_id))
    values[plots %in% surveyed$stable_plot_id, ] <- 0
    observed <- surveyed |> filter(species %in% focal$source_species) |>
      group_by(stable_plot_id, species) |>
      summarise(cover = sum(relative_cover), .groups = "drop")
    values[cbind(match(observed$stable_plot_id, plots),
                 match(observed$species, focal$source_species))] <- observed$cover
    annual[[as.character(year)]] <- values
  }
  full_control_auc <- rep(0, nrow(focal))
  for (j in seq_along(annual)) {
    full_control_auc <- full_control_auc + colMeans(annual[[j]], na.rm = TRUE) * integration_weights[j]
  }
  stopifnot(max(abs(full_control_auc - focal$origin_controls)) < 1e-8,
            max(abs(focal$warmed - full_control_auc - focal$auc_effect)) < 1e-8)
  full_auc_checks[[id]] <- focal |>
    mutate(reconstructed_control_auc = full_control_auc,
           reconstructed_delta_auc = warmed - full_control_auc)

  # Predictor uses the same first-recorded-year definition as Check 1.
  first_rows <- raw |> filter(year == min(year))
  first_plots <- sort(unique(first_rows$stable_plot_id))
  first_values <- matrix(0, length(first_plots), nrow(focal),
                        dimnames = list(first_plots, focal$effect_id))
  observed <- first_rows |> filter(species %in% focal$source_species) |>
    group_by(stable_plot_id, species) |>
    summarise(cover = sum(relative_cover), .groups = "drop")
  first_values[cbind(match(observed$stable_plot_id, first_plots),
                    match(observed$species, focal$source_species))] <- observed$cover
  reference <- first_control_abundance$first_origin_abundance[
    match(focal$effect_id, first_control_abundance$effect_id)]
  stopifnot(max(abs(colMeans(first_values) - reference)) < 1e-12)
  annual_controls[[id]] <- list(annual = annual, weights = integration_weights)
  first_year_cover[[id]] <- first_values
  response_metadata[[id]] <- focal |>
    mutate(control_site = paste(experiment$region, experiment$origin_site, sep = "::"))
}
write.csv(bind_rows(full_auc_checks), file.path(out_dir, "full_auc_reconstruction.csv"),
          row.names = FALSE)

# 5. Draw 20 plot partitions; require B observations at every original year -
# Sample whole recorded blocks. A receives
# floor(n_units/2); B receives the remainder. Reject a split only for missing
# temporal coverage, never for the abundance or Delta AUC values.
split_assignments <- list()
split_responses <- list()
set.seed(split_seed)
for (split in seq_len(n_splits)) {
  assignments <- control_plots
  assignments$half <- "B"
  for (site in unique(assignments$control_site)) {
    site_units <- unique(assignments$unit[assignments$control_site == site])
    stopifnot(length(site_units) >= 2L)
    site_experiments <- as.character(experiments$experiment_id[
      paste(experiments$region, experiments$origin_site, sep = "::") == site])
    valid <- FALSE
    for (attempt in seq_len(1000L)) {
      # sample.int avoids sample(x)'s special case for a single numeric value.
      A_units <- site_units[sample.int(length(site_units), floor(length(site_units) / 2))]
      A_plots <- assignments$stable_plot_id[assignments$unit %in% A_units]
      B_plots <- assignments$stable_plot_id[
        assignments$control_site == site & !assignments$unit %in% A_units]
      valid <- TRUE
      for (id in site_experiments) {
        if (!any(rownames(first_year_cover[[id]]) %in% A_plots)) valid <- FALSE
        for (values in annual_controls[[id]]$annual) {
          subset <- values[rownames(values) %in% B_plots, , drop = FALSE]
          if (!nrow(subset) || all(is.na(subset[, 1]))) valid <- FALSE
        }
      }
      if (valid) break
    }
    stopifnot(valid)
    assignments$half[assignments$stable_plot_id %in% A_plots] <- "A"
  }
  split_assignments[[split]] <- mutate(assignments, split = split)
  for (id in experiment_ids) {
    focal <- response_metadata[[id]]
    A <- assignments$stable_plot_id[assignments$half == "A"]
    B <- assignments$stable_plot_id[assignments$half == "B"]
    first_values <- first_year_cover[[id]]
    abundance_A <- colMeans(first_values[rownames(first_values) %in% A, , drop = FALSE])
    control_B_auc <- rep(0, nrow(focal))
    for (j in seq_along(annual_controls[[id]]$annual)) {
      values <- annual_controls[[id]]$annual[[j]]
      control_B_auc <- control_B_auc +
        colMeans(values[rownames(values) %in% B, , drop = FALSE], na.rm = TRUE) *
        annual_controls[[id]]$weights[j]
    }
    split_responses[[paste(split, id)]] <- focal |>
      transmute(split = split, effect_id, experiment_id,
        first_origin_abundance = abundance_A, control_B_auc = control_B_auc,
        warmed_auc = warmed, auc_effect = warmed - control_B_auc)
  }
}
split_assignments <- bind_rows(split_assignments)
split_responses <- bind_rows(split_responses)
stopifnot(all(is.finite(split_responses$auc_effect)),
          all(is.finite(split_responses$first_origin_abundance)),
          nrow(split_responses) == n_splits * nrow(datasets$Climate))
write.csv(split_assignments, file.path(out_dir, "split_plot_assignments.csv"), row.names = FALSE)
write.csv(split_responses, file.path(out_dir, "split_responses.csv"), row.names = FALSE)
saveRDS(list(datasets = datasets, split_responses = split_responses,
            annual_controls = annual_controls, first_year_cover = first_year_cover),
        file.path(out_dir, "analysis_inputs.rds"))

# 6. Leave each gradient out again, separately for each plot partition ------
# Identical effect IDs, response estimates and seeds for each paired fit.
# Never pool repeated splits as though they were additional independent data.
split_signature <- digest::digest(list(datasets, split_responses, split_formulas,
  split_seed, base_seed, trees = 500L, node_size = 5L, mtry = "floor(sqrt(p))",
  ranger_version = as.character(packageVersion("ranger"))))
split_predictions <- list()
for (split in seq_len(n_splits)) {
  checkpoint_dir <- file.path(out_dir, "checkpoints", sprintf("split_%02d", split))
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  message("Split-control check: partition ", split, "/", n_splits)
  for (g in seq_along(groups)) {
    group <- groups[[g]]
    checkpoint <- file.path(checkpoint_dir, paste0(group, ".rds"))
    if (file.exists(checkpoint)) {
      saved <- readRDS(checkpoint)
      stopifnot(identical(saved$signature, split_signature))
    } else {
      fits <- list()
      predictions <- list()
      for (m in seq_along(datasets)) {
        model <- names(datasets)[[m]]
        data <- datasets[[model]] |>
          select(-auc_effect, -first_origin_abundance, -auc_se) |>
          left_join(split_responses |> filter(.data$split == .env$split) |>
                      select(effect_id, auc_effect, first_origin_abundance),
                    by = "effect_id", relationship = "one-to-one")
        train <- filter(data, as.character(gradient_id) != group)
        test <- filter(data, as.character(gradient_id) == group)
        stopifnot(!any(test$effect_id %in% train$effect_id))
        for (variant in names(split_formulas[[model]])) {
          formula <- split_formulas[[model]][[variant]]
          predictors <- all.vars(formula)[-1L]
          fit <- ranger::ranger(
            formula, data = as.data.frame(train[, all.vars(formula)]),
            num.trees = 500L, mtry = floor(sqrt(length(predictors))),
            min.node.size = 5L, splitrule = "variance", replace = TRUE,
            sample.fraction = 1, case.weights = NULL,
            respect.unordered.factors = "order", importance = "none",
            write.forest = TRUE, oob.error = FALSE, na.action = "na.fail",
            seed = base_seed + 100000L + split * 1000L + 10L * g + m,
            num.threads = 3L)
          key <- paste(model, variant, sep = "__")
          fits[[key]] <- list(fit = fit, formula = formula,
                              training_effect_ids = as.character(train$effect_id))
          predictions[[key]] <- data.frame(split, model, variant, gradient_id = group,
            effect_id = test$effect_id, observed = test$auc_effect,
            predicted = predict(fit, data = as.data.frame(test[, predictors]),
                                num.threads = 3L)$predictions)
        }
      }
      saved <- list(signature = split_signature, fits = fits,
                    predictions = bind_rows(predictions))
      saveRDS(saved, paste0(checkpoint, ".tmp"))
      stopifnot(file.rename(paste0(checkpoint, ".tmp"), checkpoint))
    }
    split_predictions[[paste(split, group)]] <- saved$predictions
  }
}
split_predictions <- bind_rows(split_predictions)
split_fold_scores <- split_predictions |>
  group_by(split, model, variant, gradient_id) |>
  summarise(n_effects = n(), SSE = sum((predicted - observed)^2),
    SST = sum((observed - mean(observed))^2), MSE = mean((predicted - observed)^2),
    mean_y = mean(observed), mean_y2 = mean(observed^2), .groups = "drop") |>
  mutate(R2pred = if_else(SST > 0, 1 - SSE / SST, NA_real_))
split_overall <- split_fold_scores |>
  group_by(split, model, variant) |>
  summarise(n_effects = sum(n_effects), n_gradients = n(),
    MSE = mean(MSE), mean_y = mean(mean_y), mean_y2 = mean(mean_y2),
    .groups = "drop") |>
  mutate(R2pred = 1 - MSE / (mean_y2 - mean_y^2))
split_changes <- split_overall |>
  select(split, model, variant, R2pred) |>
  pivot_wider(names_from = variant, values_from = R2pred) |>
  mutate(R2_change = abundance - baseline)
split_summary <- split_changes |>
  group_by(model) |>
  summarise(n_splits = n(), baseline_median = median(baseline),
    abundance_median = median(abundance), median_change = median(R2_change),
    change_min = min(R2_change), change_max = max(R2_change),
    n_improved = sum(R2_change > 0), .groups = "drop")
write.csv(split_overall, file.path(out_dir, "split_overall.csv"), row.names = FALSE)
write.csv(split_fold_scores, file.path(out_dir, "split_by_gradient.csv"), row.names = FALSE)
write.csv(split_changes, file.path(out_dir, "split_paired_changes.csv"), row.names = FALSE)
write.csv(split_summary, file.path(out_dir, "split_summary.csv"), row.names = FALSE)
saveRDS(split_predictions, file.path(out_dir, "split_predictions.rds"))

# Apply the same species-ranking calculation separately to every plot split.
split_species_predictions <- split_predictions |>
  left_join(datasets$Climate |> select(effect_id, species_id),
            by = "effect_id", relationship = "many-to-one") |>
  group_by(split, model, variant, gradient_id, species_id) |>
  summarise(n_effects = n(), observed = mean(observed), predicted = mean(predicted),
            .groups = "drop")
split_rank_by_gradient <- split_species_predictions |>
  group_by(split, model, variant, gradient_id) |>
  summarise(n_species = n(),
    Spearman = if (n() >= 2L && n_distinct(observed) > 1L && n_distinct(predicted) > 1L)
      cor(observed, predicted, method = "spearman") else NA_real_, .groups = "drop")
# Use the same gradients for the paired baseline/abundance summaries.
# Undefined correlations are not replaced with zero.
split_rank_by_gradient <- split_rank_by_gradient |>
  group_by(split, model, gradient_id) |>
  mutate(both_defined = all(is.finite(Spearman))) |> ungroup()
split_rank_overall <- split_rank_by_gradient |>
  group_by(split, model, variant) |>
  summarise(n_gradients = n(), n_defined = sum(both_defined),
    Spearman = if (any(both_defined)) mean(Spearman[both_defined]) else NA_real_,
    .groups = "drop")
split_rank_changes <- split_rank_overall |>
  select(split, model, variant, Spearman) |>
  pivot_wider(names_from = variant, values_from = Spearman) |>
  mutate(change = abundance - baseline)
split_rank_summary <- split_rank_changes |>
  group_by(model) |>
  summarise(n_splits = n(),
    baseline_median = median(baseline, na.rm = TRUE),
    abundance_median = median(abundance, na.rm = TRUE),
    median_change = median(change, na.rm = TRUE),
    n_improved = sum(change > 0, na.rm = TRUE), .groups = "drop") |>
  left_join(split_rank_overall |> group_by(model) |>
    summarise(n_defined_min = min(n_defined), n_defined_max = max(n_defined),
              .groups = "drop"), by = "model")
write.csv(split_rank_by_gradient, file.path(out_dir, "split_spearman_by_gradient.csv"),
          row.names = FALSE)
write.csv(split_rank_overall, file.path(out_dir, "split_spearman_overall.csv"),
          row.names = FALSE)
write.csv(split_rank_changes, file.path(out_dir, "split_spearman_paired_changes.csv"),
          row.names = FALSE)
write.csv(split_rank_summary, file.path(out_dir, "split_spearman_summary.csv"),
          row.names = FALSE)

# Winner/loser direction on the same species means used for Spearman --------
# Positive mean Delta AUC = winner; negative = loser. This is a sign check,
# NOT the pipeline's significance-based winner/loser classification.
# Use the already held-out predictions: no new forests or thresholds are fitted.
# Keep the original-response and split-control results separate.
category_predictions <- bind_rows(
  first_species_predictions |>
    mutate(check = "Original controls", partition = 0L),
  split_species_predictions |>
    mutate(check = "Split controls", partition = split,
      specification = if_else(variant == "baseline", "Climate / trait",
                              "Climate / trait + abundance"),
      weighting = "Unweighted", method = "RF_without_species",
      scope = "Every held-out species")) |>
  mutate(observed_category = case_when(
           observed > 0 ~ "Winner", observed < 0 ~ "Loser", TRUE ~ "Zero"),
         predicted_category = case_when(
           predicted > 0 ~ "Winner", predicted < 0 ~ "Loser", TRUE ~ "Zero"))

# Count correct and incorrect categories within each held-out gradient.
# An observed zero has neither category and is excluded from binary scoring.
# A predicted zero for an observed winner/loser counts as an incorrect answer.
category_by_gradient <- category_predictions |>
  group_by(check, partition, specification, weighting, model, method, scope, gradient_id) |>
  summarise(n_species = n(),
    observed_winners = sum(observed_category == "Winner"),
    observed_losers = sum(observed_category == "Loser"),
    observed_zero = sum(observed_category == "Zero"),
    predicted_winners = sum(predicted_category == "Winner"),
    predicted_losers = sum(predicted_category == "Loser"),
    predicted_zero = sum(predicted_category == "Zero"),
    correct_winners = sum(observed_category == "Winner" & predicted_category == "Winner"),
    correct_losers = sum(observed_category == "Loser" & predicted_category == "Loser"),
    winners_predicted_loser = sum(observed_category == "Winner" & predicted_category == "Loser"),
    losers_predicted_winner = sum(observed_category == "Loser" & predicted_category == "Winner"),
    winners_predicted_zero = sum(observed_category == "Winner" & predicted_category == "Zero"),
    losers_predicted_zero = sum(observed_category == "Loser" & predicted_category == "Zero"),
    .groups = "drop") |>
  mutate(
    accuracy = if_else(observed_winners + observed_losers > 0,
      (correct_winners + correct_losers) / (observed_winners + observed_losers), NA_real_),
    winner_recall = if_else(observed_winners > 0, correct_winners / observed_winners, NA_real_),
    loser_recall = if_else(observed_losers > 0, correct_losers / observed_losers, NA_real_),
    balanced_accuracy = (winner_recall + loser_recall) / 2,
    both_categories = observed_winners > 0 & observed_losers > 0)

# Accuracy is simply the fraction correctly classified.
# Balanced accuracy gives winner detection and loser detection equal importance:
#   (fraction of observed winners found + fraction of observed losers found) / 2.
# Always predicting one category scores 50% when both observed categories exist.
# Average gradient scores, not pooled counts. A gradient with only one observed
# category still has ordinary accuracy, but its balanced accuracy is undefined.
category_overall <- category_by_gradient |>
  group_by(check, partition, specification, weighting, model, method, scope) |>
  summarise(n_gradients = n(), n_defined = sum(both_categories),
    n_species_gradient_records = sum(n_species),
    observed_winners = sum(observed_winners), observed_losers = sum(observed_losers),
    observed_zero = sum(observed_zero),
    predicted_winners = sum(predicted_winners), predicted_losers = sum(predicted_losers),
    predicted_zero = sum(predicted_zero),
    correct_winners = sum(correct_winners), correct_losers = sum(correct_losers),
    accuracy = if (any(is.finite(accuracy))) mean(accuracy, na.rm = TRUE) else NA_real_,
    winner_recall = if (any(both_categories)) mean(winner_recall[both_categories]) else NA_real_,
    loser_recall = if (any(both_categories)) mean(loser_recall[both_categories]) else NA_real_,
    balanced_accuracy = if (any(both_categories)) mean(balanced_accuracy[both_categories]) else NA_real_,
    .groups = "drop")
first_category_overall <- category_overall |> filter(check == "Original controls")
split_category_summary <- category_overall |>
  filter(check == "Split controls") |>
  group_by(model, specification) |>
  summarise(n_splits = n(),
    accuracy = median(accuracy, na.rm = TRUE),
    winner_recall = median(winner_recall, na.rm = TRUE),
    loser_recall = median(loser_recall, na.rm = TRUE),
    balanced_accuracy = median(balanced_accuracy, na.rm = TRUE),
    n_defined_min = min(n_defined), n_defined_max = max(n_defined), .groups = "drop")
write.csv(category_by_gradient, file.path(out_dir, "winner_loser_by_gradient.csv"),
          row.names = FALSE)
write.csv(category_overall, file.path(out_dir, "winner_loser_overall.csv"),
          row.names = FALSE)
write.csv(split_category_summary, file.path(out_dir, "split_winner_loser_summary.csv"),
          row.names = FALSE)

# 7. One PDF --------------------------------------------------------------
# Training R2 is in-sample. Predictive R2 gives each gradient equal weight.
# Repeated-partition ranges describe sensitivity to the plot split, not CIs.
first_caption <- paste(
  "All resident species; no minimum gradient occurrence; US_Arizona excluded. Trait rows use their own complete observations.",
  "Abundance = first-recorded focal-species mean relative cover in origin controls, with surveyed absences as zero.",
  "Climate / trait = T + P + W, plus the indicated trait. Species-ID models add + Species; their test species must occur in training.",
  "T/P = destination temperature/log precipitation; W = warming. Overall scores give each gradient equal weight.",
  "Climate/trait forests reuse verified identical-row fits. Abundance-only and abundance-added forests are newly fitted.",
  "First-recorded abundance is not necessarily pre-treatment; these full-data responses reuse control observations.",
  sep = "\n")
split_caption <- paste(
  "20 plot partitions; A supplies first-year abundance, B supplies control AUC. All warmed plots are retained.",
  "Forests are unweighted, without species ID. Each pair uses identical responses and held-out rows.",
  "Shared control sites use the same partition; blocks stay together; B-only annual means fill missing control plot-years.",
  "R2 is calculated separately per partition with equal-gradient weighting; partitions are not pooled or treated as independent samples.",
  "Paired improvement is the median of the 20 within-partition differences, not the difference between the two medians.",
  "This removes reused control observations, not all biological dependence or the control term in Delta AUC.",
  sep = "\n")
style <- list(
  geom_tile(colour = "white", linewidth = 1),
  geom_text(aes(label = ifelse(is.na(value), "NA", sprintf("%.3f", value)),
                colour = I(ifelse(!is.na(value) & abs(value) >= 0.65, "white", "black"))),
            size = 6),
  scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC",
    midpoint = 0, limits = c(-1, 1), oob = scales::squish, na.value = "grey80",
    guide = guide_colourbar(barwidth = grid::unit(9, "cm"))),
  scale_x_discrete(position = "top"),
  scale_y_discrete(limits = rev),
  theme_minimal(base_size = 19),
  theme(panel.grid = element_blank(), axis.title.x = element_blank(),
    strip.text = element_text(size = 19, face = "bold"),
    plot.caption = element_text(size = 14, hjust = 0),
    plot.tag = element_text(size = 25, face = "bold"),
    plot.tag.position = "topleft", legend.position = "bottom"))
method_labels <- c(RF_without_species = "Without species ID",
                   RF_with_species = "With species ID")
spec_order <- c("Abundance only", "Climate / trait", "Climate / trait + abundance")
spec_labels <- c("Abundance\nonly", "Climate\n(+ trait)", "Climate (+ trait)\n+ abundance")
pdf(file.path(out_dir, "random_forest_abundance_checks.pdf"), width = 20, height = 15)
rank_caption <- paste(
  "Species are ranked within each held-out gradient by their mean signed Delta AUC across experiments.",
  "Spearman compares observed and predicted ranks; ties receive average ranks. +1 = same order; -1 = reversed order.",
  "Overall = mean of the defined gradient correlations, with equal gradient weight; it is not squared and is not R2.",
  "Fewer than two species or constant values give NA, not zero. n = gradients with a defined correlation.",
  "Cells may summarize different gradients; use the per-gradient tables for comparisons on matched gradients.",
  sep = "\n")
rank_style <- style
rank_style[[2]] <- geom_text(
  aes(label = ifelse(is.na(value), "NA",
                    paste0(sprintf("%.3f", value), "\n(n = ", n_defined, ")")),
      colour = I(ifelse(!is.na(value) & abs(value) >= 0.65, "white", "black"))),
  size = 5.5)
category_caption <- paste(
  "One observed/predicted mean Delta AUC per species per gradient, as for Spearman. Positive = winner; negative = loser.",
  "Winner-loser prediction skill = winner detection rate + loser detection rate - 1 (= 2 x balanced accuracy - 1; the TSS scale).",
  "0 = no-skill reference; 1 = perfect classification; -1 = every category wrong. W/L give the two detection percentages.",
  "Overall = mean gradient prediction skill; n = gradients containing both observed categories.",
  "Single-category gradients retain counts and ordinary accuracy in the CSV, but have undefined prediction skill.",
  "Observed zeros are excluded from binary scoring; predicted zeros count as incorrect for observed winners/losers.",
  "These are signs of estimated responses, NOT significant winner/loser classifications; uncertainty around zero is not tested.",
  sep = "\n")
category_style <- style
category_style[[2]] <- geom_text(
  aes(label = ifelse(is.na(value), "NA", paste0(
    sprintf("%.3f", value),
    "\nW ", sprintf("%.0f%%", 100 * winner_recall),
    " | L ", sprintf("%.0f%%", 100 * loser_recall),
    "\n(n = ", n_defined, ")")),
    colour = I(ifelse(!is.na(value) & abs(value) >= 0.7, "white", "black"))),
  size = 5.2)
category_style[[3]] <- scale_fill_gradient2(
  low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0,
  limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1),
  labels = scales::label_number(accuracy = 0.1), na.value = "grey80",
  guide = guide_colourbar(barwidth = grid::unit(9, "cm")))
for (metric in c("Training", "Predictive", "Spearman", "Category")) {
  shown <- if (metric == "Training")
    first_training |> mutate(value = R2_training_equal_gradient) else
      if (metric == "Predictive") first_overall |> mutate(value = R2pred) else
        if (metric == "Spearman") first_rank_overall |> mutate(value = Spearman) else
          first_category_overall |> mutate(value = 2 * balanced_accuracy - 1)
  shown <- shown |>
    mutate(specification = factor(specification, levels = spec_order, labels = spec_labels),
           method = factor(method, levels = c("RF_without_species", "RF_with_species")),
           weighting = factor(weighting, levels = c("Unweighted", "Inverse variance")))
  if (metric == "Training") {
    shown$view <- factor(as.character(shown$method),
      levels = c("RF_without_species", "RF_with_species"), labels = unname(method_labels))
  } else {
    shown <- shown |> mutate(view = case_when(
      scope == "Every held-out species" ~ "Without species ID\nEvery held-out species",
      as.character(method) == "RF_without_species" ~ "Without species ID\nSpecies seen in training",
      TRUE ~ "With species ID\nSpecies seen in training"),
      view = factor(view, levels = c("Without species ID\nEvery held-out species",
        "Without species ID\nSpecies seen in training", "With species ID\nSpecies seen in training")))
  }
  print(ggplot(shown, aes(specification, model, fill = value)) +
    (if (metric == "Category") category_style else if (metric == "Spearman") rank_style else style) +
    facet_grid(weighting ~ view, labeller = labeller(
      weighting = c(Unweighted = "Unweighted", "Inverse variance" = "SE-weighted"))) +
    labs(y = if (metric == "Category") "Winner-loser prediction skill" else
           if (metric == "Spearman") "Spearman rank correlation" else paste(metric, "R\u00b2"),
         fill = if (metric == "Category") "Winner-loser prediction skill" else
           if (metric == "Spearman") "Spearman" else "R\u00b2",
         caption = if (metric == "Category") paste(
           "All resident species; US_Arizona excluded. Same model specifications and evaluation subsets as the preceding pages.",
           category_caption, sep = "\n") else
           if (metric == "Spearman") paste(first_caption, rank_caption, sep = "\n") else first_caption,
         tag = c(Training = "A", Predictive = "B", Spearman = "C", Category = "D")[[metric]]))
}
shown <- split_summary |>
  select(model, baseline_median, abundance_median, median_change) |>
  pivot_longer(-model, names_to = "measure", values_to = "value") |>
  mutate(measure = factor(measure,
    levels = c("baseline_median", "abundance_median", "median_change"),
    labels = c("Climate / trait\nMedian predictive R\u00b2",
               "+ Independent abundance\nMedian predictive R\u00b2",
               "Paired improvement\nMedian change in R\u00b2")))
print(ggplot(shown, aes(measure, model, fill = value)) + style +
  labs(y = NULL, fill = "R\u00b2 or change", caption = split_caption, tag = "E"))
shown <- split_rank_summary |>
  select(model, baseline_median, abundance_median, median_change,
         n_defined_min, n_defined_max) |>
  mutate(n_defined = ifelse(n_defined_min == n_defined_max, as.character(n_defined_min),
                           paste0(n_defined_min, "-", n_defined_max))) |>
  pivot_longer(c(baseline_median, abundance_median, median_change),
               names_to = "measure", values_to = "value") |>
  mutate(measure = factor(measure,
    levels = c("baseline_median", "abundance_median", "median_change"),
    labels = c("Climate / trait\nMedian Spearman",
               "+ Independent abundance\nMedian Spearman",
               "Paired improvement\nMedian change in Spearman")))
print(ggplot(shown, aes(measure, model, fill = value)) + rank_style +
  labs(y = NULL, fill = "Spearman or change",
    caption = paste(rank_caption,
      "20 control-plot partitions: A supplies abundance, B supplies control AUC; all warmed plots retained.",
      "Unweighted forests, no species ID. Medians summarize partitions, not confidence intervals.",
      "Within each partition, baseline and abundance summaries use the same gradients with defined correlations.",
      sep = "\n"), tag = "F"))
shown <- split_category_summary |>
  mutate(value = 2 * balanced_accuracy - 1,
    n_defined = ifelse(n_defined_min == n_defined_max, as.character(n_defined_min),
                       paste0(n_defined_min, "-", n_defined_max)),
    specification = factor(specification,
      levels = c("Climate / trait", "Climate / trait + abundance"),
      labels = c("Climate / trait", "+ Independent abundance")))
print(ggplot(shown, aes(specification, model, fill = value)) + category_style +
  labs(y = "Winner-loser prediction skill", fill = "Winner-loser prediction skill",
    caption = paste(category_caption,
      "20 control-plot partitions: A supplies abundance, B supplies control AUC; all warmed plots retained.",
      "Unweighted forests, no species ID. Each displayed score or detection percentage is its median across the 20 partitions.",
      sep = "\n"), tag = "G"))
shown <- split_changes |>
  mutate(split = factor(split, levels = rev(seq_len(n_splits))), value = R2_change)
print(ggplot(shown, aes(model, split, fill = value)) +
  style +
  labs(y = "Control-plot partition", fill = "Change in predictive R\u00b2",
    caption = paste(split_caption,
      "Positive values: independent abundance improved prediction in that partition.", sep = "\n"),
    tag = "H"))
dev.off()
message("Saved ", file.path(out_dir, "random_forest_abundance_checks.pdf"))

# 8. R2, species ranks and winner/loser accuracy together -------------------
# Each point is one held-out gradient. Use the same species means for all
# three measures: a species observed in several experiments gets one value.
# This R2 therefore differs from the effect-level R2 in the preceding PDF.
# Only the SE-weighted forests are shown; no models are refitted here.
gradient_metrics <- first_species_predictions |>
  filter(weighting == "Inverse variance") |>
  group_by(specification, weighting, model, method, scope, gradient_id) |>
  summarise(n_species = n(),
    SSE = sum((observed - predicted)^2),
    SST = sum((observed - mean(observed))^2), .groups = "drop") |>
  mutate(R2pred = if_else(SST > 0, 1 - SSE / SST, NA_real_)) |>
  left_join(
    first_rank_by_gradient |>
      select(specification, weighting, model, method, scope, gradient_id, Spearman),
    by = c("specification", "weighting", "model", "method", "scope", "gradient_id"),
    relationship = "one-to-one") |>
  left_join(
    category_by_gradient |>
      filter(check == "Original controls") |>
      select(specification, weighting, model, method, scope, gradient_id,
        observed_winners, observed_losers, observed_zero,
        correct_winners, correct_losers, accuracy, balanced_accuracy),
    by = c("specification", "weighting", "model", "method", "scope", "gradient_id"),
    relationship = "one-to-one") |>
  mutate(winner_loser_skill = 2 * balanced_accuracy - 1,
    view = case_when(
    scope == "Every held-out species" ~ "Without species ID: every held-out species",
    method == "RF_without_species" ~ "Without species ID: species seen in training",
    TRUE ~ "With species ID: species seen in training"),
    status = case_when(
      !is.finite(R2pred) ~ "R2 undefined",
      !is.finite(Spearman) ~ "Spearman undefined",
      !is.finite(balanced_accuracy) ~ "Only one observed category",
      TRUE ~ "OK"))

# Short numeric labels keep 19 gradient names from overlapping in each panel.
# The same gradient has the same number throughout the PDF.
gradient_key <- gradient_metrics |>
  distinct(gradient_id) |>
  arrange(gradient_id) |>
  mutate(gradient_number = row_number())
gradient_metrics <- gradient_metrics |>
  left_join(gradient_key, by = "gradient_id", relationship = "many-to-one")
write.csv(gradient_metrics, file.path(out_dir, "predictive_metrics_by_gradient.csv"),
          row.names = FALSE)
gradient_labels <- paste0(gradient_key$gradient_number, " = ", gradient_key$gradient_id)
gradient_key_lines <- character()
for (first_label in seq(1L, length(gradient_labels), by = 5L)) {
  gradient_key_lines <- c(gradient_key_lines, paste(
    gradient_labels[first_label:min(first_label + 4L, length(gradient_labels))],
    collapse = "     "))
}

metric_views <- c("Without species ID: every held-out species",
                 "Without species ID: species seen in training",
                 "With species ID: species seen in training")
pdf(file.path(out_dir, "predictive_metrics_by_gradient.pdf"), width = 21, height = 18)
for (view_index in seq_along(metric_views)) {
  page_rows <- gradient_metrics |>
    filter(view == metric_views[[view_index]]) |>
    mutate(model = factor(model, levels = c("Climate", "Height", "LDMC", "LMA")),
      specification = factor(specification, levels = spec_order, labels = spec_labels))
  # Undefined ranks cannot be plotted on this axis; do not turn them into zero.
  omitted_points <- page_rows |>
    filter(!is.finite(R2pred) | !is.finite(Spearman)) |>
    count(specification, name = "n_omitted")
  omission_note <- if (nrow(omitted_points) == 0L) "All gradient scores are plotted." else
    paste0("Undefined R2 or Spearman: ", paste(
      paste0(gsub("\n", " ", as.character(omitted_points$specification)),
             " = ", omitted_points$n_omitted), collapse = "; "),
      " gradient/model points omitted (listed in the CSV).")
  plot_rows <- page_rows |> filter(is.finite(R2pred), is.finite(Spearman))
  print(ggplot(plot_rows, aes(R2pred, Spearman)) +
    geom_hline(yintercept = 0, colour = "grey65", linetype = "dashed") +
    geom_vline(xintercept = 0, colour = "grey65", linetype = "dashed") +
    geom_point(aes(fill = winner_loser_skill), shape = 21, colour = "grey20",
               size = 4.3, stroke = 0.6) +
    ggrepel::geom_text_repel(aes(label = gradient_number),
      size = 4.6, seed = 20261020, max.overlaps = Inf,
      min.segment.length = 0, box.padding = 0.3, point.padding = 0.25,
      segment.colour = "grey55", max.time = 3, max.iter = 20000) +
    facet_wrap(vars(model, specification), ncol = 3, scales = "free_x", drop = FALSE,
      labeller = label_wrap_gen(width = 45, multi_line = FALSE)) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 4),
      expand = expansion(mult = c(0.14, 0.2))) +
    scale_y_continuous(limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1),
      expand = expansion(mult = 0.06)) +
    scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC",
      midpoint = 0, limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1),
      labels = scales::label_number(accuracy = 0.1), na.value = "grey65",
      guide = guide_colourbar(barwidth = grid::unit(11, "cm"), title.position = "top")) +
    labs(x = "Predictive R\u00b2", y = "Spearman rank correlation",
      fill = "Winner-loser prediction skill",
      title = paste("SE-weighted forests", metric_views[[view_index]], sep = " - "),
      tag = LETTERS[[view_index]],
      caption = paste(
        "Each point is a held-out gradient. All three scores use one mean observed/predicted Delta AUC per species within that gradient.",
        "R2 = 1 - SSE/SST; Spearman compares ranks. X-axis ranges differ between panels; no negative R2 values are truncated.",
        "Colour = winner detection rate + loser detection rate - 1 (TSS scale): 0 = no-skill reference; 1 = perfect; -1 = every category wrong.",
        "Positive Delta AUC = winner; negative = loser, regardless of significance. Grey = only one observed category; observed zeros are excluded from scoring.",
        "Climate = destination T + log P + warming; trait rows use their own complete observations. US_Arizona is excluded.",
        omission_note, paste(gradient_key_lines, collapse = "\n"), sep = "\n")) +
    theme_bw(base_size = 19) +
    theme(panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey94"),
      strip.text = element_text(size = 17, face = "bold"),
      strip.background = element_rect(fill = "grey96", colour = "grey75"),
      panel.spacing = grid::unit(1, "lines"),
      axis.title = element_text(size = 21),
      axis.text = element_text(size = 15),
      plot.title = element_text(size = 23, face = "bold"),
      plot.tag = element_text(size = 26, face = "bold"),
      plot.tag.position = "topleft",
      plot.caption = element_text(size = 13, hjust = 0),
      legend.position = "bottom",
      plot.margin = margin(15, 20, 15, 20)))
}
dev.off()
message("Saved ", file.path(out_dir, "predictive_metrics_by_gradient.pdf"))

# 9. Distributions of the three scores across gradients -------------------
# One page per species-evaluation category. Overlay the three predictor
# specifications using the same bins, with transparent bars and coloured edges.
# Keep the climate and trait datasets separate: each histogram counts gradients,
# not individual species or repeated scores from several models.
histogram_metrics <- gradient_metrics |>
  select(view, specification, model, gradient_id, R2pred, Spearman, winner_loser_skill) |>
  pivot_longer(c(R2pred, Spearman, winner_loser_skill),
               names_to = "metric", values_to = "value") |>
  mutate(specification = factor(specification,
      levels = c("Climate / trait", "Abundance only", "Climate / trait + abundance")),
    model = factor(model, levels = c("Climate", "Height", "LDMC", "LMA")),
    metric = factor(metric, levels = c("R2pred", "Spearman", "winner_loser_skill"),
      labels = c("Predictive R\u00b2", "Spearman rank correlation",
                 "Winner-loser\nprediction skill")))

# Identical bin widths throughout: R2 = 0.25, Spearman = 0.1, prediction skill = 0.1.
# Include the entire negative R2 range. Within each page, the four dataset rows
# share the same x-axis range for a given metric.
histogram_references <- tibble(
  metric = factor(levels(histogram_metrics$metric), levels = levels(histogram_metrics$metric)),
  reference = c(0, 0, 0))
histogram_page <- 0L
histogram_colours <- c("Climate / trait" = "#0072B2", "Abundance only" = "#D55E00",
                       "Climate / trait + abundance" = "#009E73")
histogram_labels <- c("Climate (+ trait)", "Abundance only", "Climate (+ trait) + abundance")
pdf(file.path(out_dir, "predictive_metric_histograms.pdf"), width = 20, height = 15)
for (current_view in metric_views) {
    histogram_page <- histogram_page + 1L
    page_metrics <- histogram_metrics |>
      filter(view == current_view)
    histogram_counts <- page_metrics |>
      group_by(model, metric, specification) |>
      summarise(n_defined = sum(is.finite(value)), .groups = "drop") |>
      arrange(model, metric, specification) |>
      group_by(model, metric) |>
      summarise(label = paste0("n = ", paste(n_defined, collapse = " | ")),
                .groups = "drop")
    # A missing Spearman score does not discard that gradient's valid R2 or
    # category accuracy. Each metric uses its own available scores.
    page_metrics <- page_metrics |> filter(is.finite(value))
    r2_breaks <- seq(
      min(0, floor(min(page_metrics$value[page_metrics$metric == "Predictive R\u00b2"]) / 0.25) * 0.25),
      1, by = 0.25)
    histogram_plot <- ggplot(page_metrics,
      aes(value, fill = specification, colour = specification)) +
      geom_histogram(data = filter(page_metrics, metric == "Predictive R\u00b2"),
        breaks = r2_breaks, position = "identity", alpha = 0.3, linewidth = 0.6) +
      geom_histogram(data = filter(page_metrics, metric == "Spearman rank correlation"),
        breaks = seq(-1, 1, by = 0.1), position = "identity", alpha = 0.3, linewidth = 0.6) +
      geom_histogram(data = filter(page_metrics, metric == "Winner-loser\nprediction skill"),
        breaks = seq(-1, 1, by = 0.1), position = "identity", alpha = 0.3, linewidth = 0.6) +
      geom_vline(data = histogram_references, aes(xintercept = reference),
        linetype = "dashed", colour = "grey20", linewidth = 0.7) +
      geom_text(data = histogram_counts, aes(label = label),
        x = Inf, y = Inf, hjust = 1.1, vjust = 1.3, size = 5.5, inherit.aes = FALSE) +
      facet_grid(model ~ metric, scales = "free_x", drop = FALSE) +
      scale_fill_manual(values = histogram_colours, labels = histogram_labels, name = NULL) +
      scale_colour_manual(values = histogram_colours, labels = histogram_labels, name = NULL) +
      guides(fill = guide_legend(override.aes = list(alpha = 0.5))) +
      scale_x_continuous(breaks = scales::breaks_pretty(n = 5),
        expand = expansion(mult = c(0.02, 0.02))) +
      scale_y_continuous(breaks = scales::breaks_width(2),
        expand = expansion(mult = c(0, 0.2))) +
      labs(x = NULL, y = "Number of gradients",
        title = "SE-weighted forests",
        subtitle = current_view, tag = LETTERS[[histogram_page]],
        caption = paste(
          "Each histogram shows the scores across held-out gradients, using one observed/predicted mean Delta AUC per species per gradient.",
          "Colours compare the three model specifications within each panel; bars overlap rather than stack, using identical bins.",
          "Dashed lines: zero for all three metrics. n lists valid gradient counts in legend order (19 gradients total).",
          "Undefined scores are omitted separately for each metric, not counted as zero. All negative R2 values are retained; x-axis ranges differ between metrics/pages.",
          "Prediction skill = winner detection rate + loser detection rate - 1 (TSS scale); 0 = no-skill reference, 1 = perfect, -1 = every category wrong.",
          "Winner/loser refers to the sign of Delta AUC, not significance; both observed categories receive equal weight.",
          "Trait rows use their own complete observations, including in the abundance-only models. US_Arizona is excluded.",
          sep = "\n")) +
      theme_bw(base_size = 20) +
      theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92"),
        panel.spacing = grid::unit(1, "lines"),
        strip.text = element_text(size = 19, face = "bold"),
        strip.background = element_rect(fill = "grey96", colour = "grey75"),
        axis.title.y = element_text(size = 22),
        axis.text = element_text(size = 17),
        plot.title = element_text(size = 24, face = "bold"),
        plot.subtitle = element_text(size = 20),
        plot.tag = element_text(size = 26, face = "bold"),
        plot.tag.position = "topleft",
        plot.caption = element_text(size = 13, hjust = 0),
        legend.position = "top", legend.text = element_text(size = 19),
        legend.key.width = grid::unit(1.1, "cm"),
        plot.margin = margin(15, 20, 15, 20))
    print(histogram_plot)
}
dev.off()
message("Saved ", file.path(out_dir, "predictive_metric_histograms.pdf"))

# 10. Density alternative to the histograms -------------------------------
# Smooth the same gradient scores; do not change or pool the model results.
# Each curve is a density, not a count. Keep automatic bandwidth selection.
# Respect the possible score ranges when smoothing near their boundaries.
density_page <- 0L
pdf(file.path(out_dir, "predictive_metric_densities.pdf"), width = 20, height = 15)
for (current_view in metric_views) {
  density_page <- density_page + 1L
  page_metrics <- histogram_metrics |> filter(view == current_view)
  density_counts <- page_metrics |>
    group_by(model, metric, specification) |>
    summarise(n_defined = sum(is.finite(value)), .groups = "drop") |>
    arrange(model, metric, specification) |>
    group_by(model, metric) |>
    summarise(label = paste0("n = ", paste(n_defined, collapse = " | ")),
              .groups = "drop")
  page_metrics <- page_metrics |> filter(is.finite(value))
  density_means <- page_metrics |>
    group_by(model, metric, specification) |>
    summarise(mean_score = mean(value), .groups = "drop") |>
    mutate(label = sprintf("Mean = %.3f", mean_score),
      label_line = 1.2 + (as.integer(specification) - 1L) * 1.4)
  # Keep x ranges consistent down each metric column. Density heights need
  # separate y scales because the three metrics have different numerical units.
  density_ranges <- page_metrics |>
    group_by(metric) |>
    summarise(low = case_when(
        first(metric) == "Predictive R\u00b2" ~ min(value) - 0.1 * diff(range(value)),
        TRUE ~ -1),
      high = 1,
      .groups = "drop") |>
    pivot_longer(c(low, high), values_to = "value") |>
    cross_join(tibble(model = factor(c("Climate", "Height", "LDMC", "LMA"),
      levels = levels(histogram_metrics$model))))
  density_plot <- ggplot(page_metrics,
    aes(value, fill = specification, colour = specification)) +
    geom_blank(data = density_ranges, aes(value, 0), inherit.aes = FALSE) +
    geom_density(data = filter(page_metrics, metric == "Predictive R\u00b2"),
      bounds = c(-Inf, 1), n = 4096, alpha = 0.25, linewidth = 0.9) +
    geom_density(data = filter(page_metrics, metric == "Spearman rank correlation"),
      bounds = c(-1, 1), n = 4096, alpha = 0.25, linewidth = 0.9) +
    geom_density(data = filter(page_metrics, metric == "Winner-loser\nprediction skill"),
      bounds = c(-1, 1), n = 4096, alpha = 0.25, linewidth = 0.9) +
    geom_rug(sides = "b", alpha = 0.5, linewidth = 0.4, show.legend = FALSE) +
    geom_vline(data = histogram_references, aes(xintercept = reference),
      linetype = "dashed", colour = "grey20", linewidth = 0.6) +
    geom_text(data = density_counts, aes(label = label),
      x = Inf, y = Inf, hjust = 1.1, vjust = 1.3, size = 5, inherit.aes = FALSE) +
    geom_text(data = density_means,
      aes(label = label, colour = specification, vjust = label_line),
      x = -Inf, y = Inf, hjust = -0.08, size = 4.8,
      inherit.aes = FALSE, show.legend = FALSE) +
    facet_wrap(vars(model, metric), ncol = 3, scales = "free",
      labeller = label_wrap_gen(width = 50, multi_line = FALSE)) +
    scale_fill_manual(values = histogram_colours, labels = histogram_labels, name = NULL) +
    scale_colour_manual(values = histogram_colours, labels = histogram_labels, name = NULL) +
    guides(fill = guide_legend(override.aes = list(alpha = 0.5))) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 5),
      expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 4),
      expand = expansion(mult = c(0, 0.2))) +
    labs(x = NULL, y = "Density", title = "SE-weighted forests",
      subtitle = current_view, tag = LETTERS[[density_page]],
      caption = paste(
        "Same gradient scores as the histogram PDF. Transparent curves use automatic smoothing; bottom ticks show the actual scores.",
        "Coloured means are arithmetic averages across the available gradients; mean gradient R2 is not the overall predictive R2.",
        "Curves show relative concentration, not gradient counts. Y-axis ranges vary between panels; compare the three colours within a panel.",
        "Dashed lines: zero for all three metrics. n lists valid gradient counts in legend order (19 gradients total).",
        "Prediction skill = winner detection rate + loser detection rate - 1 (TSS scale); 0 = no-skill reference, 1 = perfect, -1 = every category wrong.",
        "Undefined scores are omitted separately for each metric. All negative R2 scores remain visible; smoothing respects the possible score ranges.",
        "Winner/loser refers to the sign of Delta AUC, not significance. Trait rows use matched complete observations; US_Arizona is excluded.",
        sep = "\n")) +
    theme_bw(base_size = 20) +
    theme(panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "grey92"),
      panel.spacing = grid::unit(1, "lines"),
      strip.text = element_text(size = 17, face = "bold"),
      strip.background = element_rect(fill = "grey96", colour = "grey75"),
      axis.title.y = element_text(size = 22), axis.text = element_text(size = 16),
      plot.title = element_text(size = 24, face = "bold"),
      plot.subtitle = element_text(size = 20),
      plot.tag = element_text(size = 26, face = "bold"), plot.tag.position = "topleft",
      plot.caption = element_text(size = 13, hjust = 0),
      legend.position = "top", legend.text = element_text(size = 19),
      legend.key.width = grid::unit(1.1, "cm"), plot.margin = margin(15, 20, 15, 20))
  print(density_plot)
}
dev.off()
message("Saved ", file.path(out_dir, "predictive_metric_densities.pdf"))

