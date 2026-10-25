#!/usr/bin/env Rscript

# Hold out every observation of one species at a time.
# Test species recorded in at least two experiments; others stay in training.
# Run: Rscript R/additional_checks/random_forest_species_holdout.R

library(targets)
library(dplyr)
library(ggplot2)
library(ranger)

# 1. Data and visible model formulas --------------------------------------
# Delta AUC is not divided by duration. US_Arizona is excluded throughout.
# T_z, P_z, W_z and Trait_z retain their existing transformations.
# Gradient is a categorical predictor. Species identity is NOT a predictor:
# every test species is entirely absent from that fold's training data.
traits <- tar_read(predictor_specifications)
datasets <- list(Climate = tar_read(resident_auc_data), LDMC = traits$LDMC$data,
                 LMA = traits$LMA$data, Height = traits$Height$data)
for (model in names(datasets)) {
  datasets[[model]] <- datasets[[model]] |>
    filter(as.character(gradient_id) != "US_Arizona")
}
formulas <- list(
  Climate = auc_effect ~ T_z + P_z + W_z + duration + plot_area + gradient_id,
  LDMC = auc_effect ~ T_z + P_z + W_z + Trait_z + duration + plot_area + gradient_id,
  LMA = auc_effect ~ T_z + P_z + W_z + Trait_z + duration + plot_area + gradient_id,
  Height = auc_effect ~ T_z + P_z + W_z + Trait_z + duration + plot_area + gradient_id
)
base_seed <- 20260914L
out_dir <- "output/additional_checks/random_forest_species_holdout_without_US_Arizona"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Eligibility is checked within each model's complete-data subset.
species_roles <- list()
for (model in names(datasets)) {
  species_roles[[model]] <- datasets[[model]] |>
    group_by(species_id) |>
    summarise(n_effects = n(), n_experiments = n_distinct(experiment_id),
              n_gradients = n_distinct(gradient_id), .groups = "drop") |>
    mutate(model = model, held_out = n_experiments > 1L)
}
species_roles <- bind_rows(species_roles)
write.csv(species_roles, file.path(out_dir, "species_roles.csv"), row.names = FALSE)

# 2. Leave one species out ------------------------------------------------
# 500 trees; no tuning. SE weights affect bootstrap sampling, not R2 scoring.
# Keep species found in only one experiment in every training fold.
all_predictions <- list()
training_scores <- list()
for (weighting in c("Unweighted", "SE-weighted")) {
  for (m in seq_along(datasets)) {
    model <- names(datasets)[[m]]
    data <- as.data.frame(datasets[[model]])
    data$species_id <- as.character(data$species_id)
    data$gradient_id <- as.character(data$gradient_id)
    formula <- formulas[[model]]
    predictor_names <- all.vars(formula)[-1L]
    test_species <- species_roles |>
      filter(.data$model == .env$model, held_out) |>
      pull(species_id) |> as.character() |> sort()
    cache_dir <- file.path(out_dir, "checkpoints", weighting, model)
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    signature <- digest::digest(list(
      data, formula, test_species, weighting, base_seed,
      trees = 500L, mtry = "floor(sqrt(p))", node_size = 5L,
      splitrule = "variance", replace = TRUE, sample_fraction = 1,
      unordered_factors = "order", threads = 3L,
      ranger_version = as.character(packageVersion("ranger"))))

    # The final iteration fits all observations to give in-sample training R2.
    for (g in seq_len(length(test_species) + 1L)) {
      full_data <- g > length(test_species)
      species <- if (full_data) "Full training data" else test_species[[g]]
      checkpoint <- file.path(cache_dir,
        if (full_data) "full_data.rds" else paste0(digest::digest(species), ".rds"))
      train <- if (full_data) data else data[data$species_id != species, , drop = FALSE]
      test <- if (full_data) data else data[data$species_id == species, , drop = FALSE]
      train$gradient_id <- factor(train$gradient_id)
      test$gradient_id <- factor(test$gradient_id, levels = levels(train$gradient_id))
      if (anyNA(test$gradient_id)) {
        stop("Held-out species has a gradient absent from training: ", species)
      }

      if (file.exists(checkpoint)) {
        saved <- readRDS(checkpoint)
        if (!identical(saved$signature, signature)) {
          stop("Incompatible forest checkpoint: ", checkpoint)
        }
      } else {
        if (g == 1L || g %% 25L == 0L || full_data) {
          message(weighting, " / ", model, ": ", g, "/", length(test_species) + 1L)
        }
        fit <- ranger::ranger(
          formula = formula, data = train[, all.vars(formula), drop = FALSE],
          num.trees = 500L, mtry = as.integer(floor(sqrt(length(predictor_names)))),
          min.node.size = 5L, splitrule = "variance",
          replace = TRUE, sample.fraction = 1,
          case.weights = if (weighting == "SE-weighted") 1 / train$auc_se^2 else NULL,
          respect.unordered.factors = "order", importance = "none",
          write.forest = TRUE, oob.error = FALSE, na.action = "na.fail",
          seed = base_seed + 100L * g + 10L * m + 1L, num.threads = 3L)
        predicted <- predict(fit, data = test[, predictor_names, drop = FALSE],
                             num.threads = 3L)$predictions
        saved <- list(signature = signature, fit = fit, formula = formula,
          training_effect_ids = as.character(train$effect_id),
          training_species = unique(train$species_id),
          predictions = data.frame(model = model, weighting = weighting,
            effect_id = as.character(test$effect_id), species_id = test$species_id,
            experiment_id = as.character(test$experiment_id),
            gradient_id = as.character(test$gradient_id),
            observed = test$auc_effect, predicted = predicted))
        saveRDS(saved, paste0(checkpoint, ".tmp"))
        stopifnot(file.rename(paste0(checkpoint, ".tmp"), checkpoint))
      }
      if (full_data) {
        training_scores[[paste(weighting, model)]] <- saved$predictions
      } else {
        all_predictions[[paste(weighting, model, species)]] <- saved$predictions
      }
    }
  }
}
predictions <- bind_rows(all_predictions)
training_predictions <- bind_rows(training_scores)
write.csv(predictions, file.path(out_dir, "random_forest_predictions.csv"), row.names = FALSE)

# 3. R2 for each held-out species -----------------------------------------
# Compare squared prediction errors with variation around that species' mean.
# With just two observations this score can be very unstable. Zero observed
# variation makes R2 undefined; keep it NA, not zero.
by_species <- predictions |>
  group_by(weighting, model, species_id) |>
  summarise(n_effects = n(), n_experiments = n_distinct(experiment_id),
    n_gradients = n_distinct(gradient_id),
    SSE = sum((observed - predicted)^2),
    SST = sum((observed - mean(observed))^2), .groups = "drop") |>
  mutate(R2pred = if_else(SST > 0, 1 - SSE / SST, NA_real_))
write.csv(by_species, file.path(out_dir, "random_forest_by_species.csv"), row.names = FALSE)

# 4. Overall R2: continue giving geographic gradients equal weight ---------
# This is not the average of species R2s. First average squared errors within
# each gradient, then average across gradients. The denominator is the
# observed variance with the SAME equal-gradient weights and a grand mean.
# Recalculate it on the held-out observations available for each trait.
scoring_rows <- bind_rows(
  predictions |> mutate(evaluation = "Held-out species"),
  training_predictions |> mutate(evaluation = "Training"))
gradient_scores <- scoring_rows |>
  group_by(evaluation, weighting, model, gradient_id) |>
  summarise(n_effects = n(), MSE = mean((predicted - observed)^2),
    mean_observed = mean(observed), mean_observed_squared = mean(observed^2),
    .groups = "drop")
overall <- gradient_scores |>
  group_by(evaluation, weighting, model) |>
  summarise(n_gradients = n(), n_effects = sum(n_effects), MSE = mean(MSE),
    mean_observed = mean(mean_observed),
    mean_observed_squared = mean(mean_observed_squared), .groups = "drop") |>
  mutate(observed_variance = mean_observed_squared - mean_observed^2,
    R2 = if_else(observed_variance > 0, 1 - MSE / observed_variance, NA_real_))
write.csv(filter(overall, evaluation == "Held-out species"),
          file.path(out_dir, "random_forest_overall.csv"), row.names = FALSE)
write.csv(filter(overall, evaluation == "Training"),
          file.path(out_dir, "random_forest_training_R2.csv"), row.names = FALSE)

# 5. Compact figures instead of hundreds of species names -----------------
# Exact species results are in random_forest_by_species.csv.
# The last page counts species in R2 ranges, including undefined scores.
range_labels <- c("Below -1", "-1 to <0", "0 to <0.25", "0.25 to <0.5",
                  "0.5 to <0.75", "0.75 to 1", "Undefined")
species_ranges <- by_species |>
  mutate(R2_range = as.character(cut(R2pred,
    breaks = c(-Inf, -1, 0, 0.25, 0.5, 0.75, Inf), right = FALSE,
    labels = range_labels[1:6])),
    R2_range = factor(coalesce(R2_range, "Undefined"), levels = range_labels)) |>
  count(weighting, model, R2_range, .drop = FALSE) |>
  group_by(weighting, model) |>
  mutate(percent = 100 * n / sum(n)) |>
  ungroup()
write.csv(species_ranges, file.path(out_dir, "random_forest_species_R2_ranges.csv"),
          row.names = FALSE)

model_caption <- paste(
  "Delta AUC ~ T + P + W + Duration + Plot area + Gradient ID; trait models add + Trait.",
  "Trait = LDMC, LMA or height, fitted separately. No species identity or colonizer AUC predictor.",
  "Each test species is completely excluded from training. Species in only one experiment remain training-only.",
  "Overall R2 gives geographic gradients equal weight. SE-weighting changes forest training, not score weighting.",
  sep = "\n")
plot_theme <- theme_minimal(base_size = 19) +
  theme(panel.grid = element_blank(), axis.title = element_text(size = 22),
    axis.text = element_text(size = 19), strip.text = element_text(size = 21, face = "bold"),
    legend.position = "bottom", plot.caption = element_text(size = 15, hjust = 0),
    plot.caption.position = "plot", plot.tag = element_text(size = 26, face = "bold"),
    plot.tag.position = "topleft")
pdf(file.path(out_dir, "random_forest_model_performance.pdf"), width = 17, height = 12)
for (evaluation in c("Training", "Held-out species")) {
  print(ggplot(filter(overall, .data$evaluation == .env$evaluation),
    aes(factor(weighting, levels = c("Unweighted", "SE-weighted")), model, fill = R2)) +
    geom_tile(colour = "white", linewidth = 1) +
    geom_text(aes(label = sprintf("%.3f", R2),
      colour = I(ifelse(abs(R2) >= 0.65, "white", "black"))), size = 8) +
    scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC",
      midpoint = 0, limits = c(-1, 1), oob = scales::squish, name = expression(R^2)) +
    scale_x_discrete(position = "top") + scale_y_discrete(limits = rev) +
    labs(x = NULL, y = if (evaluation == "Training") expression(Training~R^2)
         else expression(Predictive~R^2), caption = model_caption,
         tag = if (evaluation == "Training") "A" else "B") + plot_theme)
}
print(ggplot(species_ranges, aes(R2_range, model, fill = percent)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.0f%%\nn=%d", percent, n),
    colour = I(ifelse(percent > 50, "white", "black"))), size = 6) +
  scale_fill_gradient(low = "white", high = "#2166AC", limits = c(0, 100),
                      name = "Species (%)") +
  scale_y_discrete(limits = rev) + facet_wrap(~ weighting, ncol = 1) +
  labs(x = expression(Held-out~species~R^2), y = NULL, tag = "C",
    caption = paste(model_caption,
      "Cells show the percentage and count of species. Species-level R2 can be unstable with few observations.",
      "Undefined means no observed variation. These percentages do not determine the overall R2 on page B.",
      sep = "\n")) + plot_theme)
dev.off()
