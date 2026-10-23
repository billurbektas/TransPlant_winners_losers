#!/usr/bin/env Rscript

# Random forests: origin and destination climate, without a warming-difference term.
# Run: Rscript R/additional_checks/random_forest_origin_destination_climate.R
# Species-inclusive models predict only species represented in training.
# Every model is fitted with and without inverse-variance sampling weights.
# US_Arizona is excluded from all training data and all performance scores.
# All species: climate/trait predictions for every held-out observation.
# Known species: only held-out observations of species present in that fold's
# training data; unseen species are excluded. 

library(targets)
library(dplyr)
library(ggplot2)
library(ranger)

# 1. Data and model formulas ----------------------------------------------
# Temperature and log precipitation from BOTH sites, standardized.
# Trait_z is the existing standardized log trait. No explicit warming predictor.
# Both site temperatures still contain the information needed to derive their
# difference; removing W_z does not make that information unavailable.
# Keep the same transformations and complete-case rows as the mixed models.
traits <- tar_read(predictor_specifications)
datasets <- list(Climate = tar_read(resident_auc_data), LDMC = traits$LDMC$data,
                 LMA = traits$LMA$data, Height = traits$Height$data)
# Response, trait values and effect rows stay unchanged. Replace W_z with
# destination-site temperature and precipitation alongside origin-site climate.
# Use the same full resident cohort as the scaling reference for every trait subset.
site_climate <- tar_read(analysis_data) |>
  distinct(experiment_id, origin_T_JJA, origin_P_JJA,
           destination_T_JJA, destination_P_JJA)
stopifnot(!anyDuplicated(site_climate$experiment_id),
          all(is.finite(site_climate$origin_T_JJA)),
          all(is.finite(site_climate$origin_P_JJA) & site_climate$origin_P_JJA > 0),
          all(is.finite(site_climate$destination_T_JJA)),
          all(is.finite(site_climate$destination_P_JJA) & site_climate$destination_P_JJA > 0))
climate_predictors <- datasets$Climate |>
  select(effect_id, experiment_id, T_z, P_z) |>
  left_join(site_climate, by = "experiment_id", relationship = "many-to-one") |>
  mutate(origin_T_z = as.numeric(scale(origin_T_JJA)),
         origin_P_z = as.numeric(scale(log(origin_P_JJA))),
         destination_T_z = T_z, destination_P_z = P_z) |>
  select(-T_z, -P_z)
# Exclude this gradient before making folds or fitting full-data forests.
for (model in names(datasets)) {
  datasets[[model]] <- datasets[[model]] |>
    filter(as.character(gradient_id) != "US_Arizona") |>
    left_join(climate_predictors |> select(effect_id, origin_T_z, origin_P_z,
                                         destination_T_z, destination_P_z),
              by = "effect_id", relationship = "one-to-one")
}
formulas <- list(
  Climate = list(
    without_species = auc_effect ~ origin_T_z + origin_P_z + destination_T_z + destination_P_z,
    with_species = auc_effect ~ origin_T_z + origin_P_z + destination_T_z + destination_P_z + species_id),
  LDMC = list(
    without_species = auc_effect ~ origin_T_z + origin_P_z + destination_T_z + destination_P_z + Trait_z,
    with_species = auc_effect ~ origin_T_z + origin_P_z + destination_T_z + destination_P_z + Trait_z + species_id),
  LMA = list(
    without_species = auc_effect ~ origin_T_z + origin_P_z + destination_T_z + destination_P_z + Trait_z,
    with_species = auc_effect ~ origin_T_z + origin_P_z + destination_T_z + destination_P_z + Trait_z + species_id),
  Height = list(
    without_species = auc_effect ~ origin_T_z + origin_P_z + destination_T_z + destination_P_z + Trait_z,
    with_species = auc_effect ~ origin_T_z + origin_P_z + destination_T_z + destination_P_z + Trait_z + species_id)
)
base_datasets <- datasets

# 2. Separate results and checkpoints -------------------------------------
base_seed <- 20260914L
out_dir <- "output/additional_checks/random_forest_origin_destination_climate_without_US_Arizona"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(climate_predictors |> filter(effect_id %in% base_datasets$Climate$effect_id),
          file.path(out_dir, "origin_destination_climate_inputs.csv"), row.names = FALSE)

# 3. Fit each forest on the other groups; predict the omitted group --------
# Species categories are learned from training responses only.
# Inverse-variance weights (1/SE^2) change bootstrap sampling during training,
# not the weighting of observations in the performance scores.
all_predictions <- list()
all_counts <- list()
input_signatures <- list()
holdout <- "Gradient"
datasets <- base_datasets
fold_column <- "gradient_id"
groups <- sort(unique(as.character(datasets$Climate$gradient_id)))
fit_dir <- file.path(out_dir, "checkpoints", "gradient")
input_signature <- digest::digest(list(
  data = datasets, formulas = formulas, seed = base_seed,
  trees = 500L, mtry = "floor(sqrt(p))", minimum_node_size = 5L,
  replace = TRUE, sample_fraction = 1, weights = NULL, splitrule = "variance",
  unordered_factors = "order", threads = 3L, prediction_rule_version = 1L,
  ranger_version = as.character(packageVersion("ranger"))
))
weighted_signature <- digest::digest(list(input_signature,
  models = names(datasets), sampling_weights = "1 / auc_se^2"))
input_signatures[[holdout]] <- input_signature
for (weighting in c("Unweighted", "Inverse variance")) {
  cache_dir <- if (weighting == "Unweighted") fit_dir else file.path(fit_dir, "weighted")
  signature <- if (weighting == "Unweighted") input_signature else weighted_signature
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  for (g in seq_along(groups)) {
    group <- groups[[g]]
    checkpoint <- file.path(cache_dir, paste0(group, ".rds"))
    if (file.exists(checkpoint)) {
      saved <- readRDS(checkpoint)
      if (!identical(saved$input_signature, signature)) {
        stop("Incompatible saved forest checkpoint: ", checkpoint)
      }
    } else {
      saved <- list(fits = list(), predictions = data.frame(), counts = data.frame())
    }
    missing_models <- setdiff(names(datasets), saved$counts$model)
    if (length(missing_models)) {
      message("Fitting ", weighting, " forests: ", group, " / ", paste(missing_models, collapse = ", "))
      fold_fits <- saved$fits
      fold_predictions <- list(saved$predictions)
      fold_counts <- list(saved$counts)
      for (m in match(missing_models, names(datasets))) {
        model <- names(datasets)[[m]]
        data <- as.data.frame(datasets[[m]])
        data$species_id <- as.character(data$species_id)
        train <- data[as.character(data[[fold_column]]) != group, , drop = FALSE]
        test <- data[as.character(data[[fold_column]]) == group, , drop = FALSE]
        train$species_id <- factor(train$species_id)
        known <- test$species_id %in% levels(train$species_id)
        sampling_weights <- if (weighting == "Inverse variance") 1 / train$auc_se^2 else NULL
        fold_counts[[model]] <- data.frame(
          heldout_group = group, model = model,
          n_train = nrow(train), n_test = nrow(test),
          n_test_species = length(unique(test$species_id)),
          n_known_effects = sum(known), n_unseen_effects = sum(!known),
          n_known_species = length(unique(test$species_id[known])),
          n_unseen_species = length(unique(test$species_id[!known])))

        for (s in seq_along(formulas[[model]])) {
          forest_type <- names(formulas[[model]])[[s]]
          formula <- formulas[[model]][[s]]
          predictor_names <- all.vars(formula)[-1L]
          fit <- ranger::ranger(
            formula = formula, data = train[, all.vars(formula), drop = FALSE],
            num.trees = 500L,
            mtry = max(1L, as.integer(floor(sqrt(length(predictor_names))))),
            min.node.size = 5L, splitrule = "variance",
            replace = TRUE, sample.fraction = 1, case.weights = sampling_weights,
            respect.unordered.factors = "order", importance = "none",
            write.forest = TRUE, oob.error = FALSE, na.action = "na.fail",
            seed = base_seed + 100L * g + 10L * m + s, num.threads = 3L)
          key <- paste(model, forest_type, sep = "__")
          fold_fits[[key]] <- list(fit = fit, formula = formula,
            training_effect_ids = as.character(train$effect_id),
            training_species = levels(train$species_id))

          # Never send an unseen species to the species-inclusive forest.
          selected <- if (forest_type == "with_species") known else rep(TRUE, nrow(test))
          prediction_data <- test[selected, , drop = FALSE]
          if (forest_type == "with_species") {
            prediction_data$species_id <- factor(prediction_data$species_id,
                                                 levels = levels(train$species_id))
          }
          predicted <- numeric()
          if (nrow(prediction_data)) {
            predicted <- predict(fit,
              data = prediction_data[, predictor_names, drop = FALSE],
              num.threads = 3L)$predictions
          }
          fold_predictions[[key]] <- data.frame(
            heldout_group = rep(group, nrow(prediction_data)),
            model = rep(model, nrow(prediction_data)),
            method = rep(paste0("RF_", forest_type), nrow(prediction_data)),
            effect_id = as.character(prediction_data$effect_id),
            species_id = as.character(prediction_data$species_id),
            observed = prediction_data$auc_effect, predicted = predicted,
            known_species = known[selected])
        }
      }
      saved <- list(input_signature = signature,
        ranger_version = as.character(packageVersion("ranger")),
        fits = fold_fits, predictions = bind_rows(fold_predictions),
        counts = bind_rows(fold_counts))
      saveRDS(saved, paste0(checkpoint, ".tmp"))
      stopifnot(file.rename(paste0(checkpoint, ".tmp"), checkpoint))
    }
    saved$predictions$gradient_id <- as.character(base_datasets$Climate$gradient_id[
      match(saved$predictions$effect_id, base_datasets$Climate$effect_id)])
    all_predictions[[paste(holdout, weighting, group)]] <-
      mutate(saved$predictions, holdout = holdout, weighting = weighting)
    all_counts[[paste(holdout, weighting, group)]] <-
      mutate(saved$counts, holdout = holdout, weighting = weighting)
  }
}
predictions <- bind_rows(all_predictions)
counts <- bind_rows(all_counts)

# 4. Predictive R2 --------------------------------------------------------
# All species: only forests without species identity, evaluated on every row.
# Known species: exclude held-out observations of species absent from training.
# Score both forests on the remaining identical rows, whether or not the forest
# uses species name as a predictor.
evaluation_rows <- bind_rows(
  predictions |> filter(method == "RF_without_species") |> mutate(scope = "All species"),
  predictions |> filter(known_species) |> mutate(scope = "Known species")
)
# Within each held-out group: R2 = 1 - SSE/SST, using that group's observed mean.
fold_scores <- evaluation_rows |>
  group_by(holdout, weighting, model, method, scope, heldout_group) |>
  summarise(n_effects = n(), n_species = n_distinct(species_id),
    SSE = sum((predicted - observed)^2),
    SST = sum((observed - mean(observed))^2),
    MSE = mean((predicted - observed)^2),
    mean_observed = mean(observed), mean_observed_squared = mean(observed^2),
    .groups = "drop")
comparisons <- data.frame(
  method = c("RF_without_species", "RF_without_species", "RF_with_species"),
  scope = c("All species", "Known species", "Known species"))
by_fold <- merge(counts, comparisons, by = NULL) |>
  left_join(fold_scores,
    by = c("holdout", "weighting", "model", "method", "scope", "heldout_group")) |>
  mutate(n_effects = coalesce(n_effects, 0L), n_species = coalesce(n_species, 0L),
    R2pred = if_else(!is.na(SST) & SST > 0, 1 - SSE / SST, NA_real_),
    status = case_when(n_effects == 0 ~ "Empty subset",
                      SST <= 0 ~ "Zero observed variance", TRUE ~ "OK"))

# Overall R2 gives each geographic gradient equal weight.
# Average gradient MSEs and observed moments before calculating R2.
# Do not average the individual gradient R2s.
gradient_scores <- evaluation_rows |>
  group_by(holdout, weighting, model, method, scope, gradient_id) |>
  summarise(n_effects = n(), MSE = mean((predicted - observed)^2),
    mean_observed = mean(observed), mean_observed_squared = mean(observed^2),
    .groups = "drop")
overall <- by_fold |>
  distinct(holdout, weighting, model, method, scope) |>
  left_join(
    gradient_scores |>
      group_by(holdout, weighting, model, method, scope) |>
      summarise(n_gradients_evaluated = n(), n_effects = sum(n_effects),
        MSE = mean(MSE), mean_observed = mean(mean_observed),
        mean_observed_squared = mean(mean_observed_squared), .groups = "drop"),
    by = c("holdout", "weighting", "model", "method", "scope")) |>
  mutate(n_gradients_evaluated = coalesce(n_gradients_evaluated, 0L),
    n_effects = coalesce(n_effects, 0L),
    n_gradients_total = n_distinct(base_datasets$Climate$gradient_id),
    observed_variance = mean_observed_squared - mean_observed^2,
    R2pred = if_else(is.finite(observed_variance) & observed_variance > 0,
                     1 - MSE / observed_variance, NA_real_),
    aggregation = "Equal geographic-gradient weight",
    status = case_when(n_gradients_evaluated == 0 ~ "Empty subset",
                      !is.finite(R2pred) ~ "Zero observed variance",
                      n_gradients_evaluated < n_gradients_total ~ "Empty gradients omitted",
                      TRUE ~ "OK"))
write.csv(predictions, file.path(out_dir, "random_forest_predictions.csv"), row.names = FALSE)
write.csv(by_fold, file.path(out_dir, "random_forest_by_fold.csv"), row.names = FALSE)
write.csv(overall, file.path(out_dir, "random_forest_overall.csv"), row.names = FALSE)
write.csv(counts, file.path(out_dir, "random_forest_fold_counts.csv"), row.names = FALSE)

# 5. Full-data training R2 -------------------------------------------------
# Fit all retained observations once for each specification and weighting.
# Training R2 is in-sample fit, not held-out or out-of-bag prediction.
datasets <- base_datasets
input_signature <- input_signatures[["Gradient"]]
weighted_signature <- digest::digest(list(input_signature,
  models = names(datasets), sampling_weights = "1 / auc_se^2"))
fit_dir <- file.path(out_dir, "checkpoints", "gradient")
training_scores <- list()
for (weighting in c("Unweighted", "Inverse variance")) {
  cache_dir <- if (weighting == "Unweighted") fit_dir else file.path(fit_dir, "weighted")
  signature <- if (weighting == "Unweighted") input_signature else weighted_signature
  full_fit_dir <- file.path(cache_dir, "full_data")
  dir.create(full_fit_dir, recursive = TRUE, showWarnings = FALSE)
  for (m in seq_along(datasets)) {
    model <- names(datasets)[[m]]
    training_data <- as.data.frame(datasets[[m]])
    training_data$species_id <- factor(as.character(training_data$species_id))
    sampling_weights <- if (weighting == "Inverse variance") 1 / training_data$auc_se^2 else NULL
    for (s in seq_along(formulas[[model]])) {
      forest_type <- names(formulas[[model]])[[s]]
      formula <- formulas[[model]][[s]]
      predictor_names <- all.vars(formula)[-1L]
      key <- paste(model, forest_type, sep = "__")
      full_seed <- base_seed + 10000L + 10L * m + s
      full_signature <- digest::digest(list(signature, key, full_seed,
                                           evaluation = "full_data_training_v1"))
      checkpoint <- file.path(full_fit_dir, paste0(key, ".rds"))
      if (file.exists(checkpoint)) {
        full_saved <- readRDS(checkpoint)
        if (!identical(full_saved$input_signature, full_signature)) {
          stop("Incompatible full-data forest checkpoint: ", checkpoint)
        }
        fit <- full_saved$fit
      } else {
        message("Fitting full-data forest: ", key)
        fit <- ranger::ranger(
          formula = formula, data = training_data[, all.vars(formula), drop = FALSE],
          num.trees = 500L,
          mtry = max(1L, as.integer(floor(sqrt(length(predictor_names))))),
          min.node.size = 5L, splitrule = "variance",
          replace = TRUE, sample.fraction = 1, case.weights = sampling_weights,
          respect.unordered.factors = "order", importance = "none",
          write.forest = TRUE, oob.error = FALSE, na.action = "na.fail",
          seed = full_seed, num.threads = 3L)
        full_saved <- list(input_signature = full_signature, fit = fit, formula = formula,
          training_effect_ids = as.character(training_data$effect_id),
          ranger_version = as.character(packageVersion("ranger")))
        saveRDS(full_saved, paste0(checkpoint, ".tmp"))
        stopifnot(file.rename(paste0(checkpoint, ".tmp"), checkpoint))
      }
      fitted_values <- predict(fit,
        data = training_data[, predictor_names, drop = FALSE], num.threads = 3L)$predictions
      y <- training_data$auc_effect
      SSE <- sum((y - fitted_values)^2)
      SST <- sum((y - mean(y))^2)
      sizes <- table(as.character(training_data$gradient_id))
      weights <- 1 / (length(sizes) * as.numeric(sizes[as.character(training_data$gradient_id)]))
      weighted_mean <- sum(weights * y)
      weighted_mse <- sum(weights * (y - fitted_values)^2)
      weighted_variance <- sum(weights * (y - weighted_mean)^2)
      training_scores[[paste(weighting, key)]] <- data.frame(
        model = model, method = paste0("RF_", forest_type), weighting = weighting,
        evaluation = "In-sample: full training dataset", n_effects = length(y),
        n_species = nlevels(training_data$species_id), n_gradients = length(sizes),
        SSE = SSE, SST = SST,
        R2_training = if (SST > 0) 1 - SSE / SST else NA_real_,
        R2_training_equal_gradient = if (weighted_variance > 0)
          1 - weighted_mse / weighted_variance else NA_real_)
    }
  }
}
training_scores <- bind_rows(training_scores)
write.csv(training_scores, file.path(out_dir, "random_forest_training_R2.csv"), row.names = FALSE)

# 6. One PDF: training fit, overall predictions, then each held-out group ---
# Colours saturate below -1; the numbers show the actual R2.
labels <- c(
  "Unweighted:RF_without_species" = "Without\nspecies ID\nUnweighted",
  "Unweighted:RF_with_species" = "With\nspecies ID\nUnweighted",
  "Inverse variance:RF_without_species" = "Without\nspecies ID\nSE-weighted",
  "Inverse variance:RF_with_species" = "With\nspecies ID\nSE-weighted")
scope_labels <- c(
  "All species" = "All held-out observations",
  "Known species" = "Held-out observations:\nspecies seen in training")
model_caption <- paste(
  "Climate: Delta AUC ~ T_origin + P_origin + T_destination + P_destination; trait models add + Trait.",
  "Models with species ID also include + Species. Trait = LDMC, LMA or height, fitted separately.",
  "T = temperature; P = precipitation. No explicit warming-difference predictor.",
  "Leave-one-gradient-out; US_Arizona excluded. Overall R2 gives each geographic gradient equal weight.", sep = "\n")
cells <- list(
  geom_tile(colour = "white", linewidth = 0.7),
  geom_text(aes(label = ifelse(is.na(R2), "NA", sprintf("%.3f", R2)),
    colour = I(ifelse(!is.na(R2) & abs(R2) >= 0.65, "white", "black"))), size = 5.5),
  scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC",
    midpoint = 0, limits = c(-1, 1), oob = scales::squish, na.value = "grey80",
    breaks = c(-1, -0.5, 0, 0.5, 1), labels = c("<= -1", "-0.5", "0", "0.5", "1"),
    name = expression(R^2), guide = guide_colourbar(barwidth = grid::unit(10, "cm"))),
  scale_x_discrete(labels = labels, position = "top"),
  scale_y_discrete(limits = rev),
  labs(caption = model_caption),
  theme_minimal(base_size = 17),
  theme(panel.grid = element_blank(), axis.ticks = element_blank(),
    axis.title = element_text(size = 21), axis.text = element_text(size = 16),
    legend.position = "bottom", legend.text = element_text(size = 16),
    strip.text = element_text(size = 18, face = "bold"),
    plot.caption = element_text(size = 15, hjust = 0),
    plot.caption.position = "plot",
    plot.tag = element_text(size = 24, face = "bold"),
    plot.tag.position = "topleft",
    panel.spacing = grid::unit(1, "cm")))
pdf(file.path(out_dir, "random_forest_model_performance.pdf"), width = 16, height = 15)
print(ggplot(training_scores |> mutate(R2 = R2_training_equal_gradient,
    column = factor(paste(weighting, method, sep = ":"), levels = names(labels))),
    aes(column, model, fill = R2)) +
  cells + labs(x = NULL, y = expression(Training~R^2), tag = "A"))
print(ggplot(overall |> mutate(R2 = R2pred,
    column = factor(paste(weighting, method, sep = ":"), levels = names(labels))),
    aes(column, model, fill = R2)) + cells +
  facet_wrap(~ scope, labeller = labeller(scope = as_labeller(scope_labels))) +
  labs(x = NULL, y = expression(Predictive~R^2), tag = "B"))
figure_number <- 3L
for (scope in c("All species", "Known species")) {
  plot_rows <- by_fold |>
    filter(.data$scope == .env$scope) |>
    mutate(R2 = R2pred,
      column = factor(paste(weighting, method, sep = ":"), levels = names(labels)))
  print(ggplot(plot_rows, aes(column, heldout_group, fill = R2)) +
    cells + facet_wrap(~ model + scope, ncol = 2,
      labeller = labeller(scope = as_labeller(scope_labels))) +
    labs(x = NULL, y = "Held-out gradient",
         fill = expression(Predictive~R^2), tag = LETTERS[figure_number]))
  figure_number <- figure_number + 1L
}
dev.off()


