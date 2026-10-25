#!/usr/bin/env Rscript

# Random forests using AUC up to recorded year 2 and year 4.
# Run: Rscript R/additional_checks/random_forest_fixed_horizon.R
# Destination climate + warming, with traits and/or species as before.
# No community-PCoA predictors; US_Arizona is excluded from every fit.
library(targets)
library(dplyr)
library(ggplot2)
library(ranger)

# 1. Cached responses and unchanged predictor values -----------------------
traits <- tar_read(predictor_specifications)
original_datasets <- list(
  Climate = tar_read(resident_auc_data), LDMC = traits$LDMC$data,
  LMA = traits$LMA$data, Height = traits$Height$data)
horizon_effects <- tar_read(auc_contrasts) |>
  ungroup() |>
  filter(transform == "raw", contrast == "to_origin", horizon %in% c("2", "4")) |>
  transmute(horizon = as.character(horizon),
    experiment_id = as.character(experiment_id),
    source_species = as.character(species),
    horizon_delta_auc = as.numeric(estimate), horizon_se = as.numeric(SE))
# These are the existing trapezoidal AUCs truncated at relative year 2 or 4,
# NOT maximum-duration AUC divided by duration. Year zero was not observed.
# Resident membership and exclusions remain those of the current analysis.
# Species pools were assigned from the complete experiment, not reclassified
# separately at each horizon. Climate and trait transformations stay unchanged.
# Horizon-specific SEs, not maximum-duration SEs, supply inverse-variance weights.

# Record the actual observation window: a horizon is an end year, not a duration.
community <- tar_read(analysis_community)
experiments <- tar_read(experiments) |> filter(included)
observation_windows <- list()
for (i in seq_len(nrow(experiments))) {
  experiment <- experiments[i, ]
  observation_windows[[i]] <- community |>
    filter(region == experiment$region) |>
    mutate(ODT = case_when(
      origin_site == experiment$origin_site &
        destination_site == experiment$origin_site & treatment == "LocalControl" ~ "origin",
      origin_site == experiment$origin_site &
        destination_site == experiment$destination_site & treatment == "Warm" ~ "warmed")) |>
    filter(!is.na(ODT)) |>
    group_by(ODT, species) |>
    filter(n() >= 2L) |>
    group_by(ODT) |>
    summarise(first_year = min(relative_year), last_year = max(relative_year),
              .groups = "drop") |>
    mutate(experiment_id = as.character(experiment$experiment_id))
}
observation_windows <- bind_rows(observation_windows) |>
  tidyr::pivot_wider(names_from = ODT, values_from = c(first_year, last_year))

formulas <- list(
  Climate = list(
    without_species = auc_effect ~ T_z + P_z + W_z,
    with_species = auc_effect ~ T_z + P_z + W_z + species_id),
  LDMC = list(
    without_species = auc_effect ~ T_z + P_z + W_z + Trait_z,
    with_species = auc_effect ~ T_z + P_z + W_z + Trait_z + species_id),
  LMA = list(
    without_species = auc_effect ~ T_z + P_z + W_z + Trait_z,
    with_species = auc_effect ~ T_z + P_z + W_z + Trait_z + species_id),
  Height = list(
    without_species = auc_effect ~ T_z + P_z + W_z + Trait_z,
    with_species = auc_effect ~ T_z + P_z + W_z + Trait_z + species_id)
)

out_dir <- "output/additional_checks/random_forest_fixed_horizon_without_US_Arizona"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
base_seed <- 20260914L
all_horizon_predictions <- all_horizon_folds <- all_horizon_overall <- list()
all_horizon_training <- all_horizon_counts <- all_horizon_inputs <- list()
availability <- list()

pdf(file.path(out_dir, "random_forest_model_performance.pdf"), width = 16, height = 15)
for (horizon in c("2", "4")) {
  # 2. Replace only each response and its SE; retain eligible resident rows ----
  selected_effects <- horizon_effects |> filter(.data$horizon == .env$horizon)
  datasets <- original_datasets
  for (model in names(datasets)) {
    original <- datasets[[model]] |>
      filter(as.character(gradient_id) != "US_Arizona")
    datasets[[model]] <- original |>
      select(-horizon) |>
      inner_join(selected_effects, by = c("experiment_id", "source_species"),
                 relationship = "one-to-one") |>
      mutate(maximum_duration_effect_id = effect_id,
        effect_id = paste(horizon, experiment_id, species_id, sep = "::"),
        auc_effect = horizon_delta_auc, auc_se = horizon_se) |>
      select(-horizon_delta_auc, -horizon_se)
    stopifnot(all(is.finite(datasets[[model]]$auc_effect)),
              all(is.finite(datasets[[model]]$auc_se) & datasets[[model]]$auc_se > 0))
    availability[[paste(horizon, model)]] <- original |>
      group_by(experiment_id, gradient_id) |>
      summarise(n_maximum_duration_effects = n(),
        n_horizon_effects = sum(effect_id %in% datasets[[model]]$maximum_duration_effect_id),
        .groups = "drop") |>
      left_join(observation_windows, by = "experiment_id", relationship = "many-to-one") |>
      mutate(horizon = .env$horizon, model = .env$model, included = n_horizon_effects > 0,
        origin_integration_span = if_else(included, as.numeric(horizon) - first_year_origin, NA_real_),
        warmed_integration_span = if_else(included, as.numeric(horizon) - first_year_warmed, NA_real_))
  }
  base_datasets <- datasets
  all_horizon_inputs[[horizon]] <- bind_rows(datasets, .id = "model")
  page_offset <- if (horizon == "2") 0L else 4L

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
  fit_dir <- file.path(out_dir, "checkpoints", paste0("year_", horizon), "gradient")
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
  
  # 5. Full-data training R2 -------------------------------------------------
  # Fit all retained observations once for each specification and weighting.
  # Training R2 is in-sample fit, not held-out or out-of-bag prediction.
  datasets <- base_datasets
  input_signature <- input_signatures[["Gradient"]]
  weighted_signature <- digest::digest(list(input_signature,
    models = names(datasets), sampling_weights = "1 / auc_se^2"))
  fit_dir <- file.path(out_dir, "checkpoints", paste0("year_", horizon), "gradient")
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
    "Climate: Delta AUC ~ T + P + W. Trait models: Delta AUC ~ T + P + W + Trait.",
    "Models with species ID also include + Species. Trait = LDMC, LMA or height, fitted separately.",
    "T = destination temperature; P = destination precipitation; W = warming.",
    paste0("Response: raw warmed-minus-origin AUC up to recorded year ", horizon, "; SEs use the same horizon."),
    "Integration starts at recorded year 1 or 2, not year zero. Eligible experiments differ between horizons.",
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
  print(ggplot(training_scores |> mutate(R2 = R2_training_equal_gradient,
      column = factor(paste(weighting, method, sep = ":"), levels = names(labels))),
      aes(column, model, fill = R2)) +
    cells + labs(x = NULL, y = expression(Training~R^2), tag = LETTERS[page_offset + 1L]))
  print(ggplot(overall |> mutate(R2 = R2pred,
      column = factor(paste(weighting, method, sep = ":"), levels = names(labels))),
      aes(column, model, fill = R2)) + cells +
    facet_wrap(~ scope, labeller = labeller(scope = as_labeller(scope_labels))) +
    labs(x = NULL, y = expression(Predictive~R^2), tag = LETTERS[page_offset + 2L]))
  figure_number <- page_offset + 3L
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

  all_horizon_predictions[[horizon]] <- mutate(predictions, horizon = horizon)
  all_horizon_folds[[horizon]] <- mutate(by_fold, horizon = horizon)
  all_horizon_overall[[horizon]] <- mutate(overall, horizon = horizon)
  all_horizon_training[[horizon]] <- mutate(training_scores, horizon = horizon)
  all_horizon_counts[[horizon]] <- mutate(counts, horizon = horizon)
}
dev.off()

# 7. One set of tables: horizon identifies each separate analysis ------------
write.csv(bind_rows(all_horizon_inputs), file.path(out_dir, "model_inputs.csv"), row.names = FALSE)
write.csv(bind_rows(availability), file.path(out_dir, "horizon_availability.csv"), row.names = FALSE)
write.csv(bind_rows(all_horizon_predictions), file.path(out_dir, "random_forest_predictions.csv"), row.names = FALSE)
write.csv(bind_rows(all_horizon_folds), file.path(out_dir, "random_forest_by_fold.csv"), row.names = FALSE)
write.csv(bind_rows(all_horizon_overall), file.path(out_dir, "random_forest_overall.csv"), row.names = FALSE)
write.csv(bind_rows(all_horizon_training), file.path(out_dir, "random_forest_training_R2.csv"), row.names = FALSE)
write.csv(bind_rows(all_horizon_counts), file.path(out_dir, "random_forest_fold_counts.csv"), row.names = FALSE)
# A change in R2 can reflect the response window AND the different eligible
# experiments/gradients. These two runs alone cannot separate those explanations.
