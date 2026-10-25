#!/usr/bin/env Rscript

# Two checks: abundance alone, and abundance measured in separate control plots.
# Keep species recorded in at least six gradients; retain their existing AUCs.
# Run: Rscript R/additional_checks/random_forest_abundance_checks.R
# Species-inclusive models predict only species represented in training.
# Check 1 retains both weighting and species-ID options.
# Check 2 uses unweighted forests without species ID, as agreed.
# US_Arizona is excluded from all training data and all performance scores.
# All species: climate/trait predictions for every held-out observation.
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
# Select the cohort once using occurrence, never the direction or size of Delta AUC.
# Count distinct geographic gradients after excluding Arizona; several experiments
# in one gradient still count as ONE gradient. Apply the same species list to all
# four predictor sets, then retain their existing trait-complete observations.
# This evaluates widespread species, not transfer to species absent from training.
species_selection <- datasets$Climate |>
  filter(as.character(gradient_id) != "US_Arizona") |>
  group_by(species_id) |>
  summarise(n_gradients = n_distinct(gradient_id), n_effects = n(), .groups = "drop") |>
  arrange(desc(n_gradients), desc(n_effects), species_id) |>
  mutate(retained = n_gradients >= 6L)
retained_species <- species_selection |> filter(retained) |> pull(species_id)
cohort_counts <- list()
for (model in names(datasets)) {
  original <- datasets[[model]] |> filter(as.character(gradient_id) != "US_Arizona")
  datasets[[model]] <- original |> filter(species_id %in% retained_species)
  cohort_counts[[model]] <- data.frame(model,
    original_effects = nrow(original), retained_effects = nrow(datasets[[model]]),
    original_species = n_distinct(original$species_id),
    retained_species = n_distinct(datasets[[model]]$species_id),
    retained_gradients = n_distinct(datasets[[model]]$gradient_id))
}
# 2. Focal-species abundance in the first recorded origin-control year -------
# Pick the year for the whole experiment, NOT each species' first occurrence.
# Average over all surveyed origin-control plots, including species absences
# as zero. Missing survey years are not filled. These are first observations,
# not necessarily pre-treatment abundances.
# Use the original species name to match the source of its AUC contrast.
# Existing relative_cover already uses the full-community denominator: do not
# renormalise it after selecting the 32 focal species.
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

out_dir <- "output/additional_checks/random_forest_abundance_checks"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
base_seed <- 20260914L
write.csv(first_control_abundance, file.path(out_dir, "first_year_origin_abundance.csv"),
          row.names = FALSE)

# 3. Check 1: can abundance alone predict almost as well? -------------------
# Repeat abundance-only fits on each trait's own complete observations, so
# comparisons never confuse adding predictors with changing the species/rows.
abundance_formulas <- list(
  without_species = auc_effect ~ first_origin_abundance,
  with_species = auc_effect ~ first_origin_abundance + species_id)
baseline_dir <- "output/additional_checks/random_forest_six_gradients_without_US_Arizona"
augmented_dir <- "output/additional_checks/random_forest_first_origin_abundance_without_US_Arizona"
baseline_predictions <- read.csv(file.path(baseline_dir, "random_forest_predictions.csv"))
augmented_predictions <- read.csv(file.path(augmented_dir, "random_forest_predictions.csv"))
baseline_training <- read.csv(file.path(baseline_dir, "random_forest_training_R2.csv"))
augmented_training <- read.csv(file.path(augmented_dir, "random_forest_training_R2.csv"))
cached_inputs <- read.csv(file.path(augmented_dir, "model_inputs.csv"))
for (model in names(datasets)) {
  old <- cached_inputs |> filter(.data$model == .env$model)
  new <- as.data.frame(datasets[[model]])
  stopifnot(setequal(old$effect_id, new$effect_id))
  old <- old[match(new$effect_id, old$effect_id), ]
  columns <- c("auc_effect", "auc_se", "T_z", "P_z", "W_z", "first_origin_abundance",
               if (model != "Climate") "Trait_z")
  stopifnot(isTRUE(all.equal(old[, columns], new[, columns],
                            check.attributes = FALSE, tolerance = 1e-12)))
}

groups <- sort(unique(as.character(datasets$Climate$gradient_id)))
first_signature <- digest::digest(list(datasets, abundance_formulas, base_seed,
  trees = 500L, node_size = 5L, mtry = "floor(sqrt(p))",
  ranger_version = as.character(packageVersion("ranger"))))
abundance_predictions <- list()
abundance_training <- list()
for (weighting in c("Unweighted", "Inverse variance")) {
  checkpoint_dir <- file.path(out_dir, "checkpoints", "abundance_only",
                              if (weighting == "Unweighted") "unweighted" else "weighted")
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  for (g in seq_along(c(groups, "Full_training"))) {
    group <- c(groups, "Full_training")[[g]]
    full_training <- group == "Full_training"
    checkpoint <- file.path(checkpoint_dir, paste0(group, ".rds"))
    if (file.exists(checkpoint)) {
      saved <- readRDS(checkpoint)
      stopifnot(identical(saved$signature, first_signature))
    } else {
      message("Abundance only: ", weighting, " / ", group)
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
        for (s in seq_along(abundance_formulas)) {
          method <- names(abundance_formulas)[[s]]
          formula <- abundance_formulas[[s]]
          predictors <- all.vars(formula)[-1L]
          seed <- if (full_training) base_seed + 10000L + 10L * m + s else
            base_seed + 100L * g + 10L * m + s
          fit <- ranger::ranger(
            formula, data = train[, all.vars(formula), drop = FALSE],
            num.trees = 500L, mtry = floor(sqrt(length(predictors))),
            min.node.size = 5L, splitrule = "variance", replace = TRUE,
            sample.fraction = 1,
            case.weights = if (weighting == "Inverse variance") 1 / train$auc_se^2 else NULL,
            respect.unordered.factors = "order", importance = "none",
            write.forest = TRUE, oob.error = FALSE, na.action = "na.fail",
            seed = seed, num.threads = 3L)
          selected <- if (method == "with_species") known else rep(TRUE, nrow(test))
          test_rows <- test[selected, ]
          test_rows$species_id <- factor(test_rows$species_id,
                                         levels = levels(train$species_id))
          predicted <- predict(fit, data = test_rows[, predictors, drop = FALSE],
                               num.threads = 3L)$predictions
          key <- paste(model, method, sep = "__")
          fits[[key]] <- list(fit = fit, formula = formula,
            training_effect_ids = as.character(train$effect_id),
            training_species = levels(train$species_id))
          rows <- data.frame(model, method = paste0("RF_", method), weighting,
            heldout_group = group, gradient_id = as.character(test_rows$gradient_id),
            effect_id = test_rows$effect_id, species_id = as.character(test_rows$species_id),
            known_species = known[selected], observed = test_rows$auc_effect,
            predicted = predicted)
          if (full_training) {
            sizes <- table(rows$gradient_id)
            w <- 1 / (length(sizes) * as.numeric(sizes[rows$gradient_id]))
            training[[key]] <- data.frame(model, method = paste0("RF_", method), weighting,
              R2_training_equal_gradient =
                1 - sum(w * (rows$predicted - rows$observed)^2) /
                  sum(w * (rows$observed - sum(w * rows$observed))^2))
          } else {
            stopifnot(!any(rows$effect_id %in% train$effect_id))
            predictions[[key]] <- rows
          }
        }
      }
      saved <- list(signature = first_signature, fits = fits,
        predictions = bind_rows(predictions), training = bind_rows(training))
      saveRDS(saved, paste0(checkpoint, ".tmp"))
      stopifnot(file.rename(paste0(checkpoint, ".tmp"), checkpoint))
    }
    abundance_predictions[[paste(weighting, group)]] <- saved$predictions
    abundance_training[[paste(weighting, group)]] <- saved$training
  }
}
abundance_predictions <- bind_rows(abundance_predictions)
abundance_training <- bind_rows(abundance_training)
comparison_predictions <- bind_rows(
  mutate(abundance_predictions, specification = "Abundance only"),
  mutate(baseline_predictions, specification = "Climate / trait"),
  mutate(augmented_predictions, specification = "Climate / trait + abundance"))
# Compare every specification on the exact same held-out observations.
prediction_keys <- c("model", "method", "weighting", "heldout_group", "effect_id")
for (reference in list(baseline_predictions, augmented_predictions)) {
  paired <- abundance_predictions |>
    inner_join(reference |> select(all_of(prediction_keys), reference_y = observed),
               by = prediction_keys, relationship = "one-to-one")
  stopifnot(nrow(paired) == nrow(abundance_predictions),
            nrow(paired) == nrow(reference),
            max(abs(paired$observed - paired$reference_y)) < 1e-12)
}
# All held-out species in this widespread-species cohort are known to training.
stopifnot(all(comparison_predictions$known_species))
first_fold_scores <- comparison_predictions |>
  group_by(specification, weighting, model, method, gradient_id) |>
  summarise(n_effects = n(), SSE = sum((predicted - observed)^2),
    SST = sum((observed - mean(observed))^2), MSE = mean((predicted - observed)^2),
    mean_y = mean(observed), mean_y2 = mean(observed^2), .groups = "drop") |>
  mutate(R2pred = if_else(SST > 0, 1 - SSE / SST, NA_real_))
first_overall <- first_fold_scores |>
  group_by(specification, weighting, model, method) |>
  summarise(n_effects = sum(n_effects), n_gradients = n(),
    MSE = mean(MSE), mean_y = mean(mean_y), mean_y2 = mean(mean_y2),
    .groups = "drop") |>
  mutate(R2pred = 1 - MSE / (mean_y2 - mean_y^2))
first_training <- bind_rows(
  mutate(abundance_training, specification = "Abundance only"),
  mutate(baseline_training, specification = "Climate / trait"),
  mutate(augmented_training, specification = "Climate / trait + abundance"))
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
first_species_predictions <- comparison_predictions |>
  group_by(specification, weighting, model, method, gradient_id, species_id) |>
  summarise(n_effects = n(), observed = mean(observed), predicted = mean(predicted),
            .groups = "drop")
first_rank_by_gradient <- first_species_predictions |>
  group_by(specification, weighting, model, method, gradient_id) |>
  summarise(n_species = n(),
    Spearman = if (n() >= 2L && n_distinct(observed) > 1L && n_distinct(predicted) > 1L)
      cor(observed, predicted, method = "spearman") else NA_real_, .groups = "drop")
first_rank_overall <- first_rank_by_gradient |>
  group_by(specification, weighting, model, method) |>
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
# Sample whole blocks (one plot per recorded block in this cohort). A receives
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

# 7. One PDF --------------------------------------------------------------
# Training R2 is in-sample. Predictive R2 gives each gradient equal weight.
# Repeated-partition ranges describe sensitivity to the plot split, not CIs.
first_caption <- paste(
  "Same >=6-gradient species cohort; US_Arizona excluded. Trait rows use their own complete observations.",
  "Abundance = first-recorded focal-species mean relative cover in origin controls, with surveyed absences as zero.",
  "Climate / trait = T + P + W, plus the indicated trait. Species-ID panels add + Species to every model.",
  "T/P = destination temperature/log precipitation; W = warming. Overall scores give each gradient equal weight.",
  "Only abundance-only forests were newly fitted here; other results reuse the verified identical-row trials.",
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
for (metric in c("Training", "Predictive", "Spearman")) {
  shown <- if (metric == "Training")
    first_training |> mutate(value = R2_training_equal_gradient) else
      if (metric == "Predictive") first_overall |> mutate(value = R2pred) else
        first_rank_overall |> mutate(value = Spearman)
  shown <- shown |>
    mutate(specification = factor(specification, levels = spec_order, labels = spec_labels),
           method = factor(method, levels = c("RF_without_species", "RF_with_species")),
           weighting = factor(weighting, levels = c("Unweighted", "Inverse variance")))
  print(ggplot(shown, aes(specification, model, fill = value)) +
    (if (metric == "Spearman") rank_style else style) +
    facet_grid(weighting ~ method, labeller = labeller(method = method_labels,
      weighting = c(Unweighted = "Unweighted", "Inverse variance" = "SE-weighted"))) +
    labs(y = if (metric == "Spearman") "Spearman rank correlation" else paste(metric, "R\u00b2"),
         fill = if (metric == "Spearman") "Spearman" else "R\u00b2",
         caption = if (metric == "Spearman") paste(first_caption, rank_caption, sep = "\n") else first_caption,
         tag = c(Training = "A", Predictive = "B", Spearman = "C")[[metric]]))
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
  labs(y = NULL, fill = "R\u00b2 or change", caption = split_caption, tag = "D"))
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
      sep = "\n"), tag = "E"))
shown <- split_changes |>
  mutate(split = factor(split, levels = rev(seq_len(n_splits))), value = R2_change)
print(ggplot(shown, aes(model, split, fill = value)) +
  style +
  labs(y = "Control-plot partition", fill = "Change in predictive R\u00b2",
    caption = paste(split_caption,
      "Positive values: independent abundance improved prediction in that partition.", sep = "\n"),
    tag = "F"))
dev.off()
message("Saved ", file.path(out_dir, "random_forest_abundance_checks.pdf"))
