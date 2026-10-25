#!/usr/bin/env Rscript

# Relative-abundance contrasts, followed by the same abundance/climate/trait
# forest comparisons. Year is continuous; only origin controls and warmed plots.
# Run: Rscript R/additional_checks/random_forest_year_group_contrasts_all_species.R

library(targets)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ranger)

out_dir <- "output/additional_checks/random_forest_year_group_contrasts_all_species"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
base_seed <- 20260914L

# 1. Observed plot-years and resident species -----------------------------
# Keep the preceding all-species trial's cohort, climate/trait transformations,
# and first-year origin abundance. US_Arizona remains excluded.
previous_inputs <- readRDS(paste0(
  "output/additional_checks/random_forest_abundance_checks_all_species/",
  "analysis_inputs.rds"))
datasets <- previous_inputs$datasets
cohort <- datasets$Climate
community <- tar_read(analysis_community)
experiments <- tar_read(experiments) |>
  filter(included, experiment_id %in% cohort$experiment_id)

plot_rows <- list()
year_rows <- list()
for (i in seq_len(nrow(experiments))) {
  experiment_row <- experiments[i, ]
  observations <- bind_rows(
    origin_controls = community |>
      filter(region == experiment_row$region,
        origin_site == experiment_row$origin_site,
        destination_site == experiment_row$origin_site,
        treatment == "LocalControl"),
    warmed = community |>
      filter(region == experiment_row$region,
        origin_site == experiment_row$origin_site,
        destination_site == experiment_row$destination_site,
        treatment == "Warm"),
    .id = "experimental_group")
  surveyed_plots <- observations |>
    distinct(experimental_group, plot_year_id, stable_plot_id, year,
             experimental_year = relative_year)
  focal_species <- cohort |>
    filter(experiment_id == experiment_row$experiment_id) |>
    select(effect_id, experiment_id, gradient_id, species_id, source_species)
  recorded_abundance <- observations |>
    group_by(experimental_group, plot_year_id, source_species = species) |>
    summarise(relative_abundance = sum(relative_cover), .groups = "drop")

  # A species absent from a surveyed plot has abundance zero. A missing survey
  # is NOT a zero observation and is not imputed. Keep the existing full-
  # community relative-cover denominator, not a resident-only denominator.
  plot_rows[[i]] <- focal_species |>
    cross_join(surveyed_plots) |>
    left_join(recorded_abundance,
      by = c("experimental_group", "plot_year_id", "source_species"),
      relationship = "many-to-one") |>
    mutate(relative_abundance = coalesce(relative_abundance, 0))
  origin_years <- sort(unique(surveyed_plots$experimental_year[
    surveyed_plots$experimental_group == "origin_controls"]))
  warmed_years <- sort(unique(surveyed_plots$experimental_year[
    surveyed_plots$experimental_group == "warmed"]))
  shared_years <- intersect(origin_years, warmed_years)
  year_rows[[i]] <- data.frame(experiment_id = experiment_row$experiment_id,
    origin_years = paste(origin_years, collapse = ", "),
    warmed_years = paste(warmed_years, collapse = ", "),
    shared_years = paste(shared_years, collapse = ", "),
    n_shared_years = length(shared_years),
    origin_min_year = min(origin_years), origin_max_year = max(origin_years),
    warmed_min_year = min(warmed_years), warmed_max_year = max(warmed_years))
}
plot_data <- bind_rows(plot_rows) |>
  mutate(experimental_year = as.numeric(experimental_year),
    experimental_group = factor(experimental_group,
      levels = c("origin_controls", "warmed")),
    species = factor(species_id), experiment = factor(experiment_id))
# emmeans' default for a continuous covariate: its arithmetic mean across ALL
# model-frame rows, not a separate mean within each experiment or species.
reference_year <- mean(plot_data$experimental_year)
contrast_years <- bind_rows(year_rows) |>
  mutate(reference_year = reference_year,
    outside_origin_year_range = reference_year < origin_min_year | reference_year > origin_max_year,
    outside_warmed_year_range = reference_year < warmed_min_year | reference_year > warmed_max_year)
stopifnot(!anyDuplicated(plot_data[c("effect_id", "experimental_group", "plot_year_id")]),
          all(is.finite(plot_data$relative_abundance)),
          all(is.finite(plot_data$experimental_year)))
write.csv(contrast_years, file.path(out_dir, "contrast_years.csv"), row.names = FALSE)

# 2. The requested four-way model and species x experiment contrasts ------
response_formula <- relative_abundance ~ experimental_year *
  experimental_group * species * experiment

# All interactions are present: each observed species x experiment x group
# has its own straight line. The following is the ordinary least-squares
# solution to that formula, without building columns for absent combinations.
# It retains ONE residual variance pooled over the whole model, not separate
# per-species error variances. No gradient or plot random effect is added.
# Repeated measurements and controls shared by experiments are consequently
# not independent in reality: these model-based SEs are exploratory weights,
# not a repeated-measures correction or a test of treatment significance.
cell_coefficients <- plot_data |>
  group_by(effect_id, experiment_id, species_id, experimental_group) |>
  summarise(n = n(), n_years = n_distinct(experimental_year),
    mean_year = mean(experimental_year), mean_abundance = mean(relative_abundance),
    Sxx = sum((experimental_year - mean(experimental_year))^2),
    Sxy = sum((experimental_year - mean(experimental_year)) *
                (relative_abundance - mean(relative_abundance))), .groups = "drop") |>
  mutate(year_slope = Sxy / Sxx,
         intercept = mean_abundance - year_slope * mean_year)
stopifnot(all(cell_coefficients$n_years >= 2L), all(cell_coefficients$Sxx > 0))
plot_data <- plot_data |>
  left_join(cell_coefficients |>
    select(effect_id, experimental_group, mean_year, mean_abundance, year_slope),
    by = c("effect_id", "experimental_group"), relationship = "many-to-one") |>
  mutate(fitted_abundance = mean_abundance + year_slope * (experimental_year - mean_year),
         residual = relative_abundance - fitted_abundance)
model_rank <- 2L * nrow(cell_coefficients)
residual_df <- nrow(plot_data) - model_rank
residual_sum_squares <- sum(plot_data$residual^2)
residual_variance <- residual_sum_squares / residual_df

# Same convention as:
#   emmeans(fit, ~ experimental_group | species * experiment)
#   contrast(means, list("warmed - origin" = c(-1, 1)))
# Year stays at the one overall model-frame mean calculated above, with no
# `at` or `cov.reduce` override. Some experiments may require extrapolation;
# contrast_years.csv records this without changing the chosen convention.
# The result is an abundance difference, NOT a slope difference or AUC.
contrast_predictions <- cell_coefficients |>
  left_join(contrast_years |> select(experiment_id, reference_year),
    by = "experiment_id", relationship = "many-to-one") |>
  mutate(fitted_mean = mean_abundance + year_slope * (reference_year - mean_year),
    mean_variance = residual_variance * (1 / n + (reference_year - mean_year)^2 / Sxx))
contrasts <- contrast_predictions |>
  select(effect_id, experiment_id, species_id, reference_year,
         experimental_group, fitted_mean, mean_variance) |>
  pivot_wider(names_from = experimental_group, values_from = c(fitted_mean, mean_variance)) |>
  mutate(abundance_change = fitted_mean_warmed - fitted_mean_origin_controls,
    contrast_se = sqrt(mean_variance_warmed + mean_variance_origin_controls),
    direction = case_when(abundance_change > 0 ~ "Winner",
                          abundance_change < 0 ~ "Loser", TRUE ~ "Zero"))
stopifnot(nrow(contrasts) == nrow(cohort), !anyDuplicated(contrasts$effect_id),
          all(is.finite(contrasts$abundance_change)),
          all(is.finite(contrasts$contrast_se) & contrasts$contrast_se > 0))
write.csv(contrasts, file.path(out_dir, "species_experiment_contrasts.csv"), row.names = FALSE)
write.csv(cell_coefficients, file.path(out_dir, "response_model_coefficients.csv"), row.names = FALSE)
saveRDS(list(formula = response_formula, coefficients = cell_coefficients,
  rank = model_rank, n = nrow(plot_data), residual_df = residual_df,
  residual_sum_squares = residual_sum_squares, residual_variance = residual_variance,
  reference_year = reference_year,
  year_convention = "emmeans default: mean experimental_year across the entire model frame",
  contrast_years = contrast_years, contrasts = contrasts),
  file.path(out_dir, "response_model.rds"))

cohort_counts <- list()
for (model in names(datasets)) {
  datasets[[model]] <- datasets[[model]] |>
    select(-auc_effect, -auc_se) |>
    left_join(contrasts |> select(effect_id, abundance_change, contrast_se, reference_year),
      by = "effect_id", relationship = "one-to-one")
  data <- datasets[[model]]
  cohort_counts[[model]] <- data.frame(model, n_effects = nrow(data),
    n_species = n_distinct(data$species_id), n_experiments = n_distinct(data$experiment_id),
    n_gradients = n_distinct(data$gradient_id))
}
write.csv(bind_rows(cohort_counts), file.path(out_dir, "cohort_counts.csv"), row.names = FALSE)
saveRDS(list(datasets = datasets, plot_data = plot_data, contrasts = contrasts,
             contrast_years = contrast_years), file.path(out_dir, "analysis_inputs.rds"))
message("Estimated ", nrow(contrasts), " species x experiment contrasts from ",
        nrow(plot_data), " plot-species-year observations at experimental year ",
        signif(reference_year, 6), " (emmeans default).")

# 3. Forest specifications ------------------------------------------------
# T_z/P_z: destination climate; W_z: warming. Trait_z is the existing log-trait
# standardisation. First-origin abundance is unchanged from the previous trial.
formulas <- list(
  Climate = list(
    "Abundance only" = list(
      without_species = abundance_change ~ first_origin_abundance,
      with_species = abundance_change ~ first_origin_abundance + species_id),
    "Climate / trait" = list(
      without_species = abundance_change ~ T_z + P_z + W_z,
      with_species = abundance_change ~ T_z + P_z + W_z + species_id),
    "Climate / trait + abundance" = list(
      without_species = abundance_change ~ T_z + P_z + W_z + first_origin_abundance,
      with_species = abundance_change ~ T_z + P_z + W_z + first_origin_abundance + species_id)),
  LDMC = list(
    "Abundance only" = list(
      without_species = abundance_change ~ first_origin_abundance,
      with_species = abundance_change ~ first_origin_abundance + species_id),
    "Climate / trait" = list(
      without_species = abundance_change ~ T_z + P_z + W_z + Trait_z,
      with_species = abundance_change ~ T_z + P_z + W_z + Trait_z + species_id),
    "Climate / trait + abundance" = list(
      without_species = abundance_change ~ T_z + P_z + W_z + Trait_z + first_origin_abundance,
      with_species = abundance_change ~ T_z + P_z + W_z + Trait_z + first_origin_abundance + species_id)),
  LMA = list(
    "Abundance only" = list(
      without_species = abundance_change ~ first_origin_abundance,
      with_species = abundance_change ~ first_origin_abundance + species_id),
    "Climate / trait" = list(
      without_species = abundance_change ~ T_z + P_z + W_z + Trait_z,
      with_species = abundance_change ~ T_z + P_z + W_z + Trait_z + species_id),
    "Climate / trait + abundance" = list(
      without_species = abundance_change ~ T_z + P_z + W_z + Trait_z + first_origin_abundance,
      with_species = abundance_change ~ T_z + P_z + W_z + Trait_z + first_origin_abundance + species_id)),
  Height = list(
    "Abundance only" = list(
      without_species = abundance_change ~ first_origin_abundance,
      with_species = abundance_change ~ first_origin_abundance + species_id),
    "Climate / trait" = list(
      without_species = abundance_change ~ T_z + P_z + W_z + Trait_z,
      with_species = abundance_change ~ T_z + P_z + W_z + Trait_z + species_id),
    "Climate / trait + abundance" = list(
      without_species = abundance_change ~ T_z + P_z + W_z + Trait_z + first_origin_abundance,
      with_species = abundance_change ~ T_z + P_z + W_z + Trait_z + first_origin_abundance + species_id)))

# 4. Leave one gradient out; keep every forest separate --------------------
# Species-ID models predict only species seen in that fold's training data.
# No-ID models predict all test species, and are also scored on the identical
# known-species subset for comparison with the species-ID forest.
groups <- sort(unique(as.character(datasets$Climate$gradient_id)))
input_signature <- digest::digest(list(datasets, formulas, response_formula,
  base_seed, trees = 500L, node_size = 5L, mtry = "floor(sqrt(p))",
  sampling = "bootstrap; case.weights = 1 / contrast_se^2 or NULL",
  unordered_factors = "order", ranger_version = as.character(packageVersion("ranger"))))
comparison_predictions <- list()
training_scores <- list()
for (weighting in c("Unweighted", "Inverse variance")) {
  checkpoint_dir <- file.path(out_dir, "checkpoints",
    if (weighting == "Unweighted") "unweighted" else "weighted")
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  for (g in seq_along(c(groups, "Full_training"))) {
    group <- c(groups, "Full_training")[[g]]
    full_training <- group == "Full_training"
    checkpoint <- file.path(checkpoint_dir, paste0(group, ".rds"))
    if (file.exists(checkpoint)) {
      saved <- readRDS(checkpoint)
      stopifnot(identical(saved$signature, input_signature))
    } else {
      message("Relative-abundance forests: ", weighting, " / ", group)
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
            # The response changed, so ALL forests are fitted anew. Never use
            # the old AUC forest checkpoints. Weights use the new contrast SE.
            fit <- ranger::ranger(
              formula, data = train[, all.vars(formula), drop = FALSE],
              num.trees = 500L, mtry = floor(sqrt(length(predictors))),
              min.node.size = 5L, splitrule = "variance", replace = TRUE,
              sample.fraction = 1,
              case.weights = if (weighting == "Inverse variance") 1 / train$contrast_se^2 else NULL,
              respect.unordered.factors = "order", importance = "none",
              write.forest = TRUE, oob.error = FALSE, na.action = "na.fail",
              seed = seed, num.threads = 3L)
            key <- paste(specification, model, method, sep = "__")
            fits[[key]] <- list(fit = fit, formula = formula,
              training_effect_ids = as.character(train$effect_id),
              training_species = levels(train$species_id))
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
              observed = test_rows$abundance_change, predicted = predicted)
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
      saved <- list(signature = input_signature, fits = fits,
                    predictions = bind_rows(predictions), training = bind_rows(training))
      saveRDS(saved, paste0(checkpoint, ".tmp"))
      stopifnot(file.rename(paste0(checkpoint, ".tmp"), checkpoint))
    }
    comparison_predictions[[paste(weighting, group)]] <- saved$predictions
    training_scores[[paste(weighting, group)]] <- saved$training
  }
}
comparison_predictions <- bind_rows(comparison_predictions)
training_scores <- bind_rows(training_scores)
saveRDS(comparison_predictions, file.path(out_dir, "forest_predictions.rds"))
write.csv(training_scores, file.path(out_dir, "training_R2.csv"), row.names = FALSE)

# 5. Predictive R2, species ranks and winner/loser scores ------------------
evaluation_predictions <- bind_rows(
  comparison_predictions |> filter(method == "RF_without_species") |>
    mutate(scope = "Every held-out species"),
  comparison_predictions |> filter(known_species) |>
    mutate(scope = "Species seen in training"))
# Within each held-out gradient, compare every species x experiment pair.
# No averaging across experiments: R2 = 1 - SSE / SST, with SST centred on
# the mean observed contrast across all evaluated pairs in that gradient.
fold_scores <- evaluation_predictions |>
  group_by(specification, weighting, model, method, scope, gradient_id) |>
  summarise(n_effects = n(), n_species = n_distinct(species_id),
    SSE = sum((predicted - observed)^2), SST = sum((observed - mean(observed))^2),
    MSE = mean((predicted - observed)^2), mean_y = mean(observed),
    mean_y2 = mean(observed^2), .groups = "drop") |>
  mutate(R2pred = if_else(SST > 0, 1 - SSE / SST, NA_real_)) |>
  complete(nesting(specification, weighting, model, method, scope), gradient_id = groups,
           fill = list(n_effects = 0L, n_species = 0L))
# Give each gradient equal weight. This is NOT mean(gradient R2).
overall <- fold_scores |>
  group_by(specification, weighting, model, method, scope) |>
  summarise(n_gradients = sum(n_effects > 0L), n_effects = sum(n_effects),
    MSE = mean(MSE[is.finite(MSE)]), mean_y = mean(mean_y[is.finite(mean_y)]),
    mean_y2 = mean(mean_y2[is.finite(mean_y2)]), .groups = "drop") |>
  mutate(R2pred = if_else(mean_y2 - mean_y^2 > 0, 1 - MSE / (mean_y2 - mean_y^2), NA_real_))
write.csv(fold_scores, file.path(out_dir, "effect_R2_by_gradient.csv"), row.names = FALSE)
write.csv(overall, file.path(out_dir, "equal_gradient_R2.csv"), row.names = FALSE)
fold_counts <- comparison_predictions |>
  filter(specification == "Abundance only", method == "RF_without_species",
         weighting == "Unweighted") |>
  group_by(model, gradient_id) |>
  summarise(n_effects = n(), n_species = n_distinct(species_id),
    n_known_effects = sum(known_species), n_known_species = n_distinct(species_id[known_species]),
    n_excluded_from_species_id = sum(!known_species), .groups = "drop")
write.csv(fold_counts, file.path(out_dir, "fold_counts.csv"), row.names = FALSE)

# Only rank and category scores use species averages within each gradient.
# R2 comes from the individual species x experiment pairs in fold_scores.
species_predictions <- evaluation_predictions |>
  group_by(specification, weighting, model, method, scope, gradient_id, species_id) |>
  summarise(n_effects = n(), observed = mean(observed), predicted = mean(predicted),
            .groups = "drop")
gradient_metrics <- species_predictions |>
  group_by(specification, weighting, model, method, scope, gradient_id) |>
  summarise(n_species = n(),
    Spearman = if (n() >= 2L && n_distinct(observed) > 1L && n_distinct(predicted) > 1L)
      cor(observed, predicted, method = "spearman") else NA_real_,
    observed_winners = sum(observed > 0), observed_losers = sum(observed < 0),
    observed_zero = sum(observed == 0),
    correct_winners = sum(observed > 0 & predicted > 0),
    correct_losers = sum(observed < 0 & predicted < 0), .groups = "drop") |>
  complete(nesting(specification, weighting, model, method, scope), gradient_id = groups,
           fill = list(n_species = 0L)) |>
  left_join(fold_scores |>
      select(specification, weighting, model, method, scope, gradient_id,
             n_effects, SSE, SST, R2pred),
    by = c("specification", "weighting", "model", "method", "scope", "gradient_id"),
    relationship = "one-to-one") |>
  mutate(winner_recall = if_else(observed_winners > 0, correct_winners / observed_winners, NA_real_),
    loser_recall = if_else(observed_losers > 0, correct_losers / observed_losers, NA_real_),
    balanced_accuracy = (winner_recall + loser_recall) / 2,
    winner_loser_skill = winner_recall + loser_recall - 1,
    view = case_when(
      scope == "Every held-out species" ~ "Without species ID: every held-out species",
      method == "RF_without_species" ~ "Without species ID: species seen in training",
      TRUE ~ "With species ID: species seen in training"),
    status = case_when(n_species == 0L ~ "Empty subset",
      !is.finite(R2pred) ~ "R2 undefined", !is.finite(Spearman) ~ "Spearman undefined",
      !is.finite(winner_loser_skill) ~ "Only one observed category", TRUE ~ "OK"))
# Categories use sign, NOT significance. A predicted zero detects neither
# winner nor loser. An observed zero is excluded from the category score.
write.csv(species_predictions, file.path(out_dir, "species_predictions.csv"), row.names = FALSE)
write.csv(gradient_metrics, file.path(out_dir, "predictive_metrics_by_gradient.csv"), row.names = FALSE)

# 6. Full-training R2, then the held-out density comparisons ---------------
# Show SE-weighted forests; both weighting options are available in the tables.
specification_order <- c("Abundance only", "Climate / trait", "Climate / trait + abundance")
model_colours <- c("Abundance only" = "#808080", "Climate / trait" = "#0072B2",
                   "Climate / trait + abundance" = "#D55E00")
model_labels <- c("Abundance only" = "Abundance only",
                  "Climate / trait" = "Climate (+ trait)",
                  "Climate / trait + abundance" = "Climate (+ trait) + abundance")
metric_views <- c("Without species ID: every held-out species",
                 "Without species ID: species seen in training",
                 "With species ID: species seen in training")
density_metrics <- gradient_metrics |>
  filter(weighting == "Inverse variance") |>
  select(view, specification, model, gradient_id, R2pred, Spearman, winner_loser_skill) |>
  pivot_longer(c(R2pred, Spearman, winner_loser_skill), names_to = "metric", values_to = "value") |>
  mutate(specification = factor(specification,
      levels = specification_order),
    model = factor(model, levels = c("Climate", "Height", "LDMC", "LMA")),
    metric = factor(metric, levels = c("R2pred", "Spearman", "winner_loser_skill"),
      labels = c("Predictive R\u00b2", "Spearman rank correlation", "Winner-loser\nprediction skill")))
# Add the earlier overall equal-gradient R2 beside the mean gradient R2.
# This uses all effect pairs; it does not average responses within species.
overall_annotations <- overall |>
  filter(weighting == "Inverse variance") |>
  left_join(gradient_metrics |> distinct(model, method, scope, view),
    by = c("model", "method", "scope"), relationship = "many-to-one") |>
  transmute(view, model = factor(model, levels = levels(density_metrics$model)),
    specification = factor(specification, levels = specification_order),
    metric = factor("Predictive R\u00b2", levels = levels(density_metrics$metric)),
    overall_R2 = R2pred)
references <- tibble(metric = factor(levels(density_metrics$metric),
  levels = levels(density_metrics$metric)), reference = 0)
pdf(file.path(out_dir, "predictive_metric_densities.pdf"), width = 20, height = 15)

# These are predictions on the same observations used to fit the full-data
# forests, NOT out-of-bag or held-out predictions. The existing scores give
# each geographic gradient equal weight and use individual effect rows.
training_table <- training_scores |>
  filter(weighting == "Inverse variance") |>
  mutate(specification = factor(specification, levels = specification_order),
    model = factor(model, levels = rev(c("Climate", "Height", "LDMC", "LMA"))),
    method = factor(method, levels = c("RF_without_species", "RF_with_species"),
      labels = c("Without species ID", "With species ID")))
print(ggplot(training_table, aes(specification, model)) +
  geom_tile(aes(fill = specification), alpha = 0.18, colour = "white", linewidth = 2) +
  geom_text(aes(label = sprintf("%.3f", R2_training_equal_gradient)), size = 10) +
  facet_wrap(vars(method), nrow = 1) +
  scale_fill_manual(values = model_colours, breaks = specification_order,
    labels = model_labels, guide = "none") +
  scale_x_discrete(position = "top", expand = expansion(add = 0),
    labels = c("Abundance only" = "Abundance only", "Climate / trait" = "Climate (+ trait)",
      "Climate / trait + abundance" = "Climate (+ trait)\n+ abundance")) +
  scale_y_discrete(expand = expansion(add = 0)) +
  labs(x = NULL, y = NULL, title = "Full-training R\u00b2",
    subtitle = "SE-weighted forests: relative-abundance contrasts", tag = "A",
    caption = paste(
      "Training R2: predict the same observations used to fit each full-data forest. These are not held-out or out-of-bag scores.",
      "R2 = 1 - weighted SSE / weighted SST; each geographic gradient has equal total weight. Scores use individual species x experiment effects.",
      "Trait rows use their own complete observations for all three model specifications. US_Arizona is excluded. Cell colour identifies the predictor set, not R2 magnitude.",
      "The following pages show held-out scores: R2 uses every species x experiment pair; rank and category scores use species averages within each gradient.",
      sep = "\n")) +
  theme_minimal(base_size = 25) +
  theme(panel.grid = element_blank(),
    panel.spacing = grid::unit(2, "lines"), strip.text = element_text(size = 27, face = "bold"),
    axis.text.x = element_text(size = 22), axis.text.y = element_text(size = 26),
    plot.title = element_text(size = 30, face = "bold"), plot.subtitle = element_text(size = 24),
    plot.tag = element_text(size = 30, face = "bold"), plot.tag.position = "topleft",
    plot.caption = element_text(size = 15, hjust = 0),
    plot.margin = margin(20, 25, 20, 25)))

for (page in seq_along(metric_views)) {
  page_metrics <- density_metrics |> filter(view == metric_views[[page]])
  density_counts <- page_metrics |>
    group_by(model, metric, specification) |>
    summarise(n_defined = sum(is.finite(value)), .groups = "drop") |>
    arrange(model, metric, specification) |>
    group_by(model, metric) |>
    summarise(label = paste0("n = ", paste(n_defined, collapse = " | ")), .groups = "drop")
  page_metrics <- page_metrics |> filter(is.finite(value))
  density_means <- page_metrics |>
    group_by(model, metric, specification) |>
    summarise(mean_score = mean(value), .groups = "drop") |>
    left_join(overall_annotations |> filter(view == metric_views[[page]]) |> select(-view),
      by = c("model", "metric", "specification"), relationship = "one-to-one") |>
    mutate(label = if_else(metric == "Predictive R\u00b2",
        sprintf("Mean = %.3f | Overall = %.3f", mean_score, overall_R2),
        sprintf("Mean = %.3f", mean_score)),
      label_line = 1.2 + (as.integer(specification) - 1L) * 1.4)
  density_ranges <- page_metrics |>
    group_by(metric) |>
    summarise(low = case_when(first(metric) == "Predictive R\u00b2" ~
      min(value) - 0.1 * diff(range(value)), TRUE ~ -1), high = 1, .groups = "drop") |>
    pivot_longer(c(low, high), values_to = "value") |>
    cross_join(tibble(model = factor(c("Climate", "Height", "LDMC", "LMA"),
      levels = levels(density_metrics$model))))
  print(ggplot(page_metrics, aes(value, fill = specification, colour = specification)) +
    geom_blank(data = density_ranges, aes(value, 0), inherit.aes = FALSE) +
    geom_density(data = filter(page_metrics, metric == "Predictive R\u00b2"),
      bounds = c(-Inf, 1), n = 4096, alpha = 0.25, linewidth = 0.9) +
    geom_density(data = filter(page_metrics, metric == "Spearman rank correlation"),
      bounds = c(-1, 1), n = 4096, alpha = 0.25, linewidth = 0.9) +
    geom_density(data = filter(page_metrics, metric == "Winner-loser\nprediction skill"),
      bounds = c(-1, 1), n = 4096, alpha = 0.25, linewidth = 0.9) +
    geom_rug(sides = "b", alpha = 0.5, linewidth = 0.4, show.legend = FALSE) +
    geom_vline(data = references, aes(xintercept = reference),
      linetype = "dashed", colour = "grey20", linewidth = 0.6) +
    geom_text(data = density_counts, aes(label = label), x = Inf, y = Inf,
      hjust = 1.1, vjust = 1.3, size = 5, inherit.aes = FALSE) +
    geom_text(data = density_means, aes(label = label, colour = specification, vjust = label_line),
      x = -Inf, y = Inf, hjust = -0.08, size = 4.8, inherit.aes = FALSE, show.legend = FALSE) +
    facet_wrap(vars(model, metric), ncol = 3, scales = "free",
      labeller = label_wrap_gen(width = 50, multi_line = FALSE)) +
    scale_fill_manual(values = model_colours, breaks = specification_order, labels = model_labels, name = NULL) +
    scale_colour_manual(values = model_colours, breaks = specification_order, labels = model_labels, name = NULL) +
    guides(fill = guide_legend(override.aes = list(alpha = 0.5))) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 5), expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 4), expand = expansion(mult = c(0, 0.2))) +
    labs(x = NULL, y = "Density", title = "SE-weighted forests: relative-abundance contrasts",
      subtitle = metric_views[[page]], tag = LETTERS[[page + 1L]],
      caption = paste(
        sprintf("Response: fitted warmed minus origin-control relative abundance at experimental year %.4f (overall model-frame mean; emmeans default for continuous year).", reference_year),
        "Within each held-out gradient, R2 compares every observed species x experiment contrast with its prediction, without averaging experiments.",
        "Mean = average of valid gradient scores. Overall = equal-gradient R2 calculated across all effect pairs, not the mean of gradient R2 values.",
        "Rank and category scores retain one mean contrast per species per gradient; their coloured annotations show mean gradient scores only.",
        "Climate = destination temperature + log precipitation + warming; trait rows retain matched complete observations. First-origin abundance is unchanged. US_Arizona excluded.",
        "Prediction skill = winner detection rate + loser detection rate - 1; winner/loser uses the contrast sign, not significance. Zero is the reference for all three metrics.",
        "Weights = inverse squared contrast SE, from the requested ordinary linear model; repeated-plot/shared-control dependence is not corrected in these exploratory SEs.",
        "Automatic density smoothing; ticks show gradient scores; n gives valid counts in legend order. Undefined scores are omitted separately; all negative R2 values are retained.",
        sep = "\n")) +
    theme_bw(base_size = 20) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "grey92"), panel.spacing = grid::unit(1, "lines"),
      strip.text = element_text(size = 17, face = "bold"),
      strip.background = element_rect(fill = "grey96", colour = "grey75"),
      axis.title.y = element_text(size = 22), axis.text = element_text(size = 16),
      plot.title = element_text(size = 24, face = "bold"), plot.subtitle = element_text(size = 20),
      plot.tag = element_text(size = 26, face = "bold"), plot.tag.position = "topleft",
      plot.caption = element_text(size = 13, hjust = 0), legend.position = "top",
      legend.text = element_text(size = 19), legend.key.width = grid::unit(1.1, "cm"),
      plot.margin = margin(15, 20, 15, 20)))
}
dev.off()
write.csv(density_metrics |>
  group_by(view, model, metric, specification) |>
  summarise(n_defined = sum(is.finite(value)), mean_score = mean(value[is.finite(value)]),
    .groups = "drop") |>
  left_join(overall_annotations, by = c("view", "model", "metric", "specification"),
    relationship = "one-to-one"),
  file.path(out_dir, "density_means.csv"), row.names = FALSE)
message("Saved ", file.path(out_dir, "predictive_metric_densities.pdf"))
