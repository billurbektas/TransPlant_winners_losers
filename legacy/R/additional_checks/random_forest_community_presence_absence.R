#!/usr/bin/env Rscript

# Presence-absence community composition and leave-one-gradient-out random forests.
# Run: Rscript R/additional_checks/random_forest_community_presence_absence.R
library(targets)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ranger)

# 1. Species presence-absence in all observed communities ----------------------
# Each row is one physical plot in one observed year, counted once even when
# the same controls are used by several experiments. Include all species in
# the reviewed community data, not just the resident-model species.
community <- tar_read(analysis_community)
taxonomy <- tar_read(community_taxonomy_species) |>
  distinct(name_original, accepted_name_reviewed)
experiments <- tar_read(experiments) |> filter(included)
community_membership <- list()
for (i in seq_len(nrow(experiments))) {
  experiment <- experiments[i, ]
  community_membership[[i]] <- community |>
    filter(region == experiment$region) |>
    mutate(role = case_when(
      origin_site == experiment$origin_site &
        destination_site == experiment$origin_site & treatment == "LocalControl" ~ "Origin control",
      origin_site == experiment$destination_site &
        destination_site == experiment$destination_site & treatment == "LocalControl" ~ "Destination control",
      origin_site == experiment$origin_site &
        destination_site == experiment$destination_site & treatment == "Warm" ~ "Warmed")) |>
    filter(!is.na(role)) |>
    distinct(plot_year_id, stable_plot_id, year, relative_year, role) |>
    mutate(experiment_id = as.character(experiment$experiment_id))
}
community_membership <- bind_rows(community_membership)
community_rows <- community |>
  filter(plot_year_id %in% community_membership$plot_year_id) |>
  left_join(taxonomy, by = c("species" = "name_original"), relationship = "many-to-one")
stopifnot(!anyNA(community_rows$accepted_name_reviewed))
# A species is present if any record under its accepted name has positive cover.
# Combine synonyms before recording 0/1. Do not normalise these binary rows or
# invent unsurveyed plot-years; absent species receive zero.
presences <- community_rows |>
  group_by(plot_year_id, accepted_name_reviewed) |>
  summarise(presence = as.integer(any(relative_cover > 0)), .groups = "drop") |>
  pivot_wider(names_from = accepted_name_reviewed, values_from = presence,
              values_fill = 0) |>
  arrange(plot_year_id)
presence_matrix <- as.matrix(presences[, -1])
rownames(presence_matrix) <- presences$plot_year_id
plot_metadata <- community_rows |>
  distinct(plot_year_id, stable_plot_id, region, year, relative_year, treatment)
stopifnot(!anyDuplicated(plot_metadata$plot_year_id),
          all(presence_matrix %in% c(0, 1)), all(rowSums(presence_matrix) > 0))

out_dir <- "output/additional_checks/random_forest_community_presence_absence"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(community_membership, file.path(out_dir, "community_membership.csv"), row.names = FALSE)

# 2. Binary Bray-Curtis (Sorensen) PCoA; axes reaching at least 50% --------------------
# Lingoes adds a constant to squared distances if needed to remove negative
# eigenvalues. The threshold refers to the corrected positive eigenvalues.
# Keep distance-preserving coordinates; do NOT standardise axes individually.
ordination_signature <- digest::digest(list(
  presence_matrix, distance = "binary_bray_sorensen", correction = "lingoes",
  vegan_version = as.character(packageVersion("vegan"))))
ordination_file <- file.path(out_dir, "community_pcoa.rds")
if (file.exists(ordination_file)) {
  saved_ordination <- readRDS(ordination_file)
  if (!identical(saved_ordination$signature, ordination_signature)) {
    stop("Community ordination checkpoint does not match the current data.")
  }
  ordination <- saved_ordination$ordination
} else {
  message("PCoA: ", nrow(presence_matrix), " plot-year communities, ",
          ncol(presence_matrix), " species.")
  community_distance <- vegan::vegdist(presence_matrix, method = "bray", binary = TRUE)
  ordination <- vegan::wcmdscale(community_distance, eig = TRUE, add = "lingoes")
  saveRDS(list(signature = ordination_signature, ordination = ordination,
              metadata = plot_metadata, species = colnames(presence_matrix)),
          ordination_file)
}
positive_eigenvalues <- ordination$eig[ordination$eig > 0]
axis_table <- data.frame(
  axis = seq_along(positive_eigenvalues),
  eigenvalue = positive_eigenvalues,
  variance_fraction = positive_eigenvalues / sum(positive_eigenvalues),
  cumulative_fraction = cumsum(positive_eigenvalues) / sum(positive_eigenvalues))
n_axes <- which(axis_table$cumulative_fraction >= 0.5)[1]
axis_table$retained <- axis_table$axis <= n_axes
# Lingoes shifts each of the n-1 centred eigenvalues by its constant.
# Report its share of corrected inertia so the 50% threshold is not mistaken
# for 50% of uncorrected ecological variation.
correction_fraction <- ordination$ac * (nrow(presence_matrix) - 1) / sum(ordination$eig)
coordinates <- ordination$points[, seq_len(n_axes), drop = FALSE]
colnames(coordinates) <- paste0("PCoA", seq_len(n_axes))
axis_names <- colnames(coordinates)
scores <- data.frame(plot_year_id = rownames(coordinates), coordinates,
                     check.names = FALSE)
write.csv(axis_table, file.path(out_dir, "pcoa_axes.csv"), row.names = FALSE)
write.csv(scores, file.path(out_dir, "community_pcoa_scores.csv"), row.names = FALSE)
write.csv(data.frame(n_communities = nrow(presence_matrix),
  n_species = ncol(presence_matrix), retained_axes = n_axes,
  cumulative_fraction = axis_table$cumulative_fraction[n_axes],
  lingoes_constant = ordination$ac,
  correction_fraction_of_inertia = correction_fraction),
  file.path(out_dir, "pcoa_summary.csv"), row.names = FALSE)
message("Retained ", n_axes, " axes: ",
        round(100 * axis_table$cumulative_fraction[n_axes], 2), "%.")

# 3. Origin-control centroids and first-year spread -----------------------
# First year = first recorded control year, not necessarily a pre-treatment year.
# All-year centroid = average of all observed control plot-years: each plot-year
# has equal weight. Shared controls can legitimately give identical centroids
# to several experiments; each experiment remains one reference in comparisons.
controls <- community_membership |>
  filter(role == "Origin control") |>
  left_join(scores, by = "plot_year_id", relationship = "many-to-one")
first_controls <- controls |> group_by(experiment_id) |>
  filter(year == min(year)) |> ungroup()
first_centroids <- first_controls |>
  group_by(experiment_id) |>
  summarise(across(all_of(axis_names), mean), .groups = "drop")
all_centroids <- controls |>
  group_by(experiment_id) |>
  summarise(across(all_of(axis_names), mean), .groups = "drop")
first_distance <- as.matrix(dist(as.matrix(first_centroids[, axis_names])))
all_distance <- as.matrix(dist(as.matrix(all_centroids[, axis_names])))
diag(first_distance) <- NA_real_
diag(all_distance) <- NA_real_
community_features <- first_centroids |>
  select(experiment_id) |>
  mutate(first_mean_distance = rowMeans(first_distance, na.rm = TRUE)) |>
  left_join(all_centroids |> select(experiment_id) |>
    mutate(all_year_mean_distance = rowMeans(all_distance, na.rm = TRUE)),
    by = "experiment_id")
# Multivariate SD = sqrt(sum of sample variances along the retained axes).
# Equivalently: sqrt(sum(squared distances to own centroid) / (n_plots - 1)).
# It is invariant to rotation; CV is unsuitable for centred coordinate means.
first_spread <- first_controls |>
  group_by(experiment_id) |>
  summarise(first_control_year = min(year), n_first_plots = n_distinct(stable_plot_id),
    first_spread_sd = if (n() > 1L)
      sqrt(sum(vapply(pick(all_of(axis_names)), var, numeric(1)))) else NA_real_,
    .groups = "drop")
community_features <- community_features |>
  left_join(first_spread, by = "experiment_id") |>
  left_join(controls |> group_by(experiment_id) |>
    summarise(n_control_plot_years = n(), n_control_years = n_distinct(year),
              .groups = "drop"), by = "experiment_id")
stopifnot(nrow(community_features) == nrow(experiments),
          all(is.finite(as.matrix(community_features[, c(
            "first_mean_distance", "first_spread_sd", "all_year_mean_distance")]))))
write.csv(first_centroids, file.path(out_dir, "first_year_control_centroids.csv"), row.names = FALSE)
write.csv(all_centroids, file.path(out_dir, "all_year_control_centroids.csv"), row.names = FALSE)
write.csv(community_features, file.path(out_dir, "community_predictors.csv"), row.names = FALSE)

# 4. Add three community summaries to the ORIGINAL climate/trait models ---
# Destination temperature, log destination precipitation and warming remain.
# Ordination and summaries above use all 42 included experiments. RF excludes
# US_Arizona as before. These are EXPLORATORY scores: held-out communities
# influence the global ordination and reference distances; all-year controls
# also share observations with the control contribution to Delta AUC.
traits <- tar_read(predictor_specifications)
datasets <- list(Climate = tar_read(resident_auc_data), LDMC = traits$LDMC$data,
                 LMA = traits$LMA$data, Height = traits$Height$data)
for (model in names(datasets)) {
  datasets[[model]] <- datasets[[model]] |>
    filter(as.character(gradient_id) != "US_Arizona") |>
    left_join(community_features |> select(experiment_id, first_mean_distance,
                first_spread_sd, all_year_mean_distance),
              by = "experiment_id", relationship = "many-to-one")
}
formulas <- list(
  Climate = list(
    without_species = auc_effect ~ T_z + P_z + W_z + first_mean_distance + first_spread_sd + all_year_mean_distance,
    with_species = auc_effect ~ T_z + P_z + W_z + first_mean_distance + first_spread_sd + all_year_mean_distance + species_id),
  LDMC = list(
    without_species = auc_effect ~ T_z + P_z + W_z + Trait_z + first_mean_distance + first_spread_sd + all_year_mean_distance,
    with_species = auc_effect ~ T_z + P_z + W_z + Trait_z + first_mean_distance + first_spread_sd + all_year_mean_distance + species_id),
  LMA = list(
    without_species = auc_effect ~ T_z + P_z + W_z + Trait_z + first_mean_distance + first_spread_sd + all_year_mean_distance,
    with_species = auc_effect ~ T_z + P_z + W_z + Trait_z + first_mean_distance + first_spread_sd + all_year_mean_distance + species_id),
  Height = list(
    without_species = auc_effect ~ T_z + P_z + W_z + Trait_z + first_mean_distance + first_spread_sd + all_year_mean_distance,
    with_species = auc_effect ~ T_z + P_z + W_z + Trait_z + first_mean_distance + first_spread_sd + all_year_mean_distance + species_id)
)
base_datasets <- datasets
base_seed <- 20260914L

# 5. Fit each forest on the other groups; predict the omitted group --------
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

# 6. Predictive R2 --------------------------------------------------------
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

# 7. Full-data training R2 -------------------------------------------------
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

# 8. One PDF: training fit, overall predictions, then each held-out group ---
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
  "Delta AUC ~ destination T + destination P + warming + three community summaries; trait models add + Trait.",
  "Models with species ID add + Species. Trait = LDMC, LMA or height, fitted separately.",
  "Community summaries: first-year mean centroid distance, first-year spread, all-year mean centroid distance.",
  paste0("Presence-absence / Sorensen PCoA: ", n_axes, " axes reaching ", round(100 * axis_table$cumulative_fraction[n_axes], 1), "% of corrected variance."),
  paste0("Caution: the Lingoes correction contributes ", round(100 * correction_fraction, 1), "% of total corrected inertia."),
  "EXPLORATORY: all communities define the ordination; all-year control observations also contribute to Delta AUC.",
  "RF leaves out geographic gradients; Arizona excluded from RF, retained in the all-community ordination.",
  sep = "\n")
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
# PCoA display uses two axes only; centroid calculations use ALL retained axes.
ordination_display <- data.frame(plot_year_id = rownames(ordination$points),
  Axis1 = ordination$points[, 1], Axis2 = ordination$points[, 2]) |>
  left_join(plot_metadata, by = "plot_year_id", relationship = "one-to-one")
print(ggplot(ordination_display, aes(Axis1, Axis2, colour = treatment)) +
  geom_point(alpha = 0.35, size = 1) + coord_equal() +
  scale_colour_manual(values = c(LocalControl = "#3E6B8A", Warm = "#B54A43"),
                      labels = c(LocalControl = "Control", Warm = "Warmed"), name = NULL) +
  labs(x = sprintf("PCoA 1 (%.1f%%)", 100 * axis_table$variance_fraction[1]),
       y = sprintf("PCoA 2 (%.1f%%)", 100 * axis_table$variance_fraction[2]),
       tag = "E", caption = model_caption) +
  theme_classic(base_size = 18) + theme(legend.position = "bottom",
    plot.caption = element_text(size = 12, hjust = 0), plot.tag = element_text(size = 24)))
print(ggplot(axis_table, aes(axis, 100 * cumulative_fraction)) +
  geom_line(linewidth = 1) + geom_hline(yintercept = 50, linetype = 2, colour = "#B54A43") +
  geom_vline(xintercept = n_axes, linetype = 2, colour = "#B54A43") +
  annotate("text", x = n_axes, y = 55, label = paste(n_axes, "axes retained"), hjust = 0) +
  labs(x = "Number of PCoA axes", y = "Cumulative corrected variance (%)", tag = "F",
       caption = paste("Presence-absence / Sorensen distances; Lingoes constant =", signif(ordination$ac, 5),
         "\nThe correction contributes", round(100 * correction_fraction, 1),
         "% of corrected inertia; this is not a count of independent ecological gradients.")) +
  theme_classic(base_size = 18) + theme(plot.tag = element_text(size = 24)))
dev.off()
