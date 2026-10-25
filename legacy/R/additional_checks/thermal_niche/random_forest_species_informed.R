#!/usr/bin/env Rscript

# Can a species' response in an unseen gradient be predicted from how it
# responded in the other gradients? Species-informed random forests, scored
# per species instead of per gradient.
# Run after random_forest_thermal_niche.R (it reads that script's inputs):
#   Rscript R/additional_checks/thermal_niche/random_forest_species_informed.R
#
# Hold-out: one whole gradient at a time, as in the other checks; nothing from
# the held-out gradient (its species or its sites) enters training or tuning.
# Training uses every niche-cohort effect of the other gradients.
# Scored species: those predicted with at least 4 other gradients in training,
# i.e. found in >= 5 gradients. Species with >= 3 other gradients (>= 4 in
# total) are scored as well and reported separately.
#
# Scores: each species' experiments are averaged within each gradient, giving
# one held-out point per species x gradient. Across its points, each species
# gets two R2 values, one Spearman rank correlation and one winner/loser
# prediction skill (winner recall + loser recall - 1; undefined when a species
# only wins or only loses):
#   within-species R2, around the species' own mean: does the model know in
#     which gradients the species does better or worse? Knowing the species'
#     typical response earns nothing here, as with Spearman.
#   absolute R2, around the overall mean of all scored points: also rewards
#     knowing whether the species is usually a winner or a loser.
#
# Models (response: the mean-year warmed - origin contrast):
#   Species mean                 the species' mean contrast in the training
#                                gradients; no forest.
#   Abundance + Niche            forest without species information.
#   Species ID + Abundance + Niche
#                                species ID as a factor; ranger orders the
#                                species by mean training response
#                                (respect.unordered.factors = "order"), so
#                                trees can split on groups of species and then
#                                on abundance or niche within them.
#   Species ID + Climate + Abundance + Niche
#   Species memory + Abundance + Niche
#                                species identity as a leakage-safe number: a
#                                training row gets its species' mean from the
#                                OTHER training gradients, a held-out row from
#                                all of them, both shrunk toward the training
#                                mean by n / (n + 5). Checks the factor encoding,
#                                whose ordering uses the rows' own responses.
# Forests: no weights; min.node.size and mtry tuned by an inner
# leave-one-gradient-out loop inside each training set, scored on held-out rows
# of species seen in the inner training data. Half-samples without replacement.
# Because each species' prediction comes from training data that excludes the
# gradient it is scored in, the species-mean reference predicts (S - y_g) /
# (n - 1) for gradient g: a decreasing function of the held-out value. Its
# within-species Spearman is -1 and its within-species R2 is 1 - (n / (n - 1))^2
# by construction; species-ID forests inherit part of this, because they learn
# a species' level the same way. Within-species scores therefore measure only
# context dependence; the absolute R2 is the one that credits species memory.

library(dplyr)
library(tidyr)
library(ggplot2)
library(ranger)
library(parallel)

input_path <- "output/additional_checks/thermal_niche/analysis_inputs.rds"
if (!file.exists(input_path)) stop(input_path, " is missing: run random_forest_thermal_niche.R first.")
out_dir <- "output/additional_checks/thermal_niche/species_informed"
dir.create(file.path(out_dir, "checkpoints"), recursive = TRUE, showWarnings = FALSE)
base_seed <- 20261025L
n_cores <- max(1L, detectCores() - 2L)
memory_shrinkage <- 5
min_other_gradients <- 4L

# 1. Data and scored species -------------------------------------------------
cohort <- readRDS(input_path) |>
  mutate(species_id = as.character(species_id))
gradients <- sort(unique(cohort$gradient_id))
species_gradients <- cohort |>
  group_by(species_id) |>
  summarise(n_gradients = n_distinct(gradient_id), n_experiments = n_distinct(experiment_id),
            n_effects = n(), .groups = "drop")
scored <- list(
  "At least 4 other gradients" = species_gradients$species_id[species_gradients$n_gradients >= min_other_gradients + 1L],
  "At least 3 other gradients" = species_gradients$species_id[species_gradients$n_gradients >= min_other_gradients])
write.csv(species_gradients |>
            mutate(scored_primary = species_id %in% scored[[1]], scored_secondary = species_id %in% scored[[2]]),
          file.path(out_dir, "species_gradient_counts.csv"), row.names = FALSE)
primary <- names(scored)[[1]]

# 2. Models ------------------------------------------------------------------
climate <- c("T_z", "P_z", "W_z")
abundance <- c("first_origin_abundance", "dest_abundance", "dest_present")
niche <- c("dMAT_dest", "dAP_dest_log", "niche_suitability_change")
forests <- list(
  "Abundance + Niche" = c(abundance, niche),
  "Species ID + Abundance + Niche" = c("species_id", abundance, niche),
  "Species ID + Climate + Abundance + Niche" = c("species_id", climate, abundance, niche),
  "Species memory + Abundance + Niche" = c("species_memory", abundance, niche))
model_levels <- c("Species mean", names(forests))

# Species mean of `from`, shrunk toward the mean of `from`, for the species of `to`.
species_mean <- function(from, to, k) {
  means <- from |> group_by(species_id) |> summarise(n = n(), total = sum(response), .groups = "drop")
  overall <- mean(from$response)
  matched <- means[match(to$species_id, means$species_id), ]
  ifelse(is.na(matched$n), overall, (matched$total + k * overall) / (matched$n + k))
}
# Training rows get their species' mean from the other training gradients only.
add_memory <- function(train, test, k = memory_shrinkage) {
  train$species_memory <- NA_real_
  for (g in unique(train$gradient_id)) {
    held <- train$gradient_id == g
    train$species_memory[held] <- species_mean(train[!held, ], train[held, ], k)
  }
  test$species_memory <- species_mean(train, test, k)
  list(train = train, test = test)
}
# Species ID is a factor over the training species; unseen species become NA
# and are never predicted by the species-ID forests.
as_species_factor <- function(train, test) {
  train$species_id <- factor(train$species_id)
  test$species_id <- factor(test$species_id, levels = levels(train$species_id))
  list(train = train, test = test)
}
prepare <- function(train, test, predictors) {
  if ("species_memory" %in% predictors) {
    prepared <- add_memory(train, test)
    train <- prepared$train
    test <- prepared$test
  }
  if ("species_id" %in% predictors) {
    prepared <- as_species_factor(train, test)
    train <- prepared$train
    test <- prepared$test
  }
  list(train = train, test = test, known = !is.na(match(test$species_id, unique(train$species_id))))
}

# 3. Nested leave-one-gradient-out forests ----------------------------------
node_grid <- c(5L, 25L, 75L, 200L, 500L)
mtry_grid <- function(p) sort(unique(pmin(p, c(1L, floor(sqrt(p)), ceiling(p / 2), p))))
inner_trees <- 200L
outer_trees <- 1000L

equal_gradient_R2 <- function(observed, predicted, gradient) {
  mse <- tapply((predicted - observed)^2, gradient, mean)
  mean_y <- tapply(observed, gradient, mean)
  mean_y2 <- tapply(observed^2, gradient, mean)
  1 - mean(mse) / (mean(mean_y2) - mean(mean_y)^2)
}
grow <- function(train, predictors, node, mtry, trees, seed, threads = 1L) {
  ranger(x = train[, predictors, drop = FALSE], y = train$response,
         num.trees = trees, mtry = mtry, min.node.size = node,
         sample.fraction = 0.5, replace = FALSE, splitrule = "variance",
         respect.unordered.factors = "order", importance = "none",
         seed = seed, num.threads = threads)
}
forecast <- function(fit, data, predictors, threads = 1L) {
  predict(fit, data = data[, predictors, drop = FALSE], num.threads = threads)$predictions
}
# Inner score: equal-gradient R2 on the inner held-out rows of species seen in
# the inner training data, the rows the outer scores are built from.
tune <- function(train, predictors, seed) {
  grid <- expand.grid(min_node_size = node_grid, mtry = mtry_grid(length(predictors)))
  inner <- sort(unique(train$gradient_id))
  splits <- lapply(inner, function(g) {
    held <- train$gradient_id == g
    prepare(train[!held, ], train[held, ], predictors)
  })
  grid$inner_R2 <- vapply(seq_len(nrow(grid)), function(k) {
    observed <- predicted <- gradient <- list()
    for (j in seq_along(inner)) {
      split <- splits[[j]]
      fit <- grow(split$train, predictors, grid$min_node_size[k], grid$mtry[k], inner_trees, seed + j)
      rows <- split$test[split$known, ]
      observed[[j]] <- rows$response
      predicted[[j]] <- forecast(fit, rows, predictors)
      gradient[[j]] <- rows$gradient_id
    }
    equal_gradient_R2(unlist(observed), unlist(predicted), unlist(gradient))
  }, numeric(1))
  grid
}

ranger_version <- as.character(packageVersion("ranger"))
fold_results <- list()
for (s in seq_along(forests)) {
  model <- names(forests)[[s]]
  predictors <- forests[[s]]
  checkpoint <- file.path(out_dir, "checkpoints", paste0(gsub("[^A-Za-z0-9]+", "_", model), ".rds"))
  signature <- digest::digest(list(
    data = cohort[c("effect_id", "gradient_id", "species_id", "response",
                    setdiff(predictors, "species_memory"))],
    predictors = predictors, node_grid = node_grid, mtry = "1, sqrt, p/2, p",
    inner_trees = inner_trees, outer_trees = outer_trees, shrinkage = memory_shrinkage,
    sample_fraction = 0.5, replace = FALSE, unordered = "order", seed = base_seed + 1000L * s,
    ranger_version = ranger_version))
  if (file.exists(checkpoint) && identical(readRDS(checkpoint)$signature, signature)) {
    fold_results[[model]] <- readRDS(checkpoint)
    next
  }
  message("Nested forests: ", model)
  folds <- mclapply(seq_along(gradients), function(g) {
    held <- cohort$gradient_id == gradients[[g]]
    seed <- base_seed + 1000L * s + 10L * g
    grid <- tune(cohort[!held, ], predictors, seed)
    # Ties go to the larger node size, then the smaller mtry (the smoother forest).
    best <- grid[order(-grid$inner_R2, -grid$min_node_size, grid$mtry), ][1, ]
    split <- prepare(cohort[!held, ], cohort[held, ], predictors)
    fit <- grow(split$train, predictors, best$min_node_size, best$mtry, outer_trees, seed)
    rows <- split$test[split$known, ]
    list(predictions = data.frame(effect_id = rows$effect_id, predicted = forecast(fit, rows, predictors)),
         grid = mutate(grid, heldout_gradient = gradients[[g]]),
         tuning = data.frame(heldout_gradient = gradients[[g]], min_node_size = best$min_node_size,
                             mtry = best$mtry, inner_R2 = best$inner_R2))
  }, mc.cores = n_cores, mc.preschedule = FALSE)
  failed <- vapply(folds, inherits, logical(1), "try-error")
  if (any(failed)) stop(model, ": ", folds[failed][[1]])
  saved <- list(signature = signature,
    predictions = bind_rows(lapply(folds, `[[`, "predictions")),
    grid = bind_rows(lapply(folds, `[[`, "grid")),
    tuning = bind_rows(lapply(folds, `[[`, "tuning")))
  saveRDS(saved, paste0(checkpoint, ".tmp"))
  stopifnot(file.rename(paste0(checkpoint, ".tmp"), checkpoint))
  fold_results[[model]] <- saved
}
stack_results <- function(part) {
  bind_rows(lapply(names(fold_results), function(model) {
    mutate(fold_results[[model]][[part]], model = model)
  }))
}

# The species-mean reference, from the training gradients of each fold.
reference <- bind_rows(lapply(gradients, function(g) {
  held <- cohort$gradient_id == g
  test <- cohort[held, ]
  seen <- test$species_id %in% cohort$species_id[!held]
  data.frame(effect_id = test$effect_id[seen], model = "Species mean",
             predicted = species_mean(cohort[!held, ], test[seen, ], k = 0))
}))
predictions <- bind_rows(stack_results("predictions"), reference) |>
  left_join(cohort |> select(effect_id, gradient_id, experiment_id, species_id, observed = response),
            by = "effect_id", relationship = "many-to-one")
saveRDS(predictions, file.path(out_dir, "predictions.rds"))
tuning <- stack_results("tuning")
write.csv(tuning, file.path(out_dir, "tuning.csv"), row.names = FALSE)
write.csv(stack_results("grid"), file.path(out_dir, "tuning_grid.csv"), row.names = FALSE)
for (species in scored[[2]]) {
  for (model in model_levels) {
    n_rows <- sum(predictions$model == model & predictions$species_id == species)
    stopifnot(n_rows == sum(cohort$species_id == species))
  }
}

# 4. Per-species scores -----------------------------------------------------
# One point per species x held-out gradient: experiments averaged within it.
points <- predictions |>
  filter(species_id %in% scored[[2]]) |>
  group_by(model, species_id, gradient_id) |>
  summarise(n_experiments = n(), observed = mean(observed), predicted = mean(predicted),
            .groups = "drop")
score_species <- function(points, overall_mean) {
  points |>
    group_by(model, species_id) |>
    summarise(n_gradients = n(),
      SSE = sum((predicted - observed)^2), SST = sum((observed - mean(observed))^2),
      MSE = mean((predicted - observed)^2), variance = mean((observed - mean(observed))^2),
      variance_absolute = mean((observed - overall_mean)^2),
      R2 = 1 - SSE / SST, R2_absolute = 1 - MSE / variance_absolute,
      Spearman = if (n_distinct(observed) > 1L && n_distinct(predicted) > 1L)
        cor(observed, predicted, method = "spearman") else NA_real_,
      winner_recall = if (any(observed > 0)) mean(predicted[observed > 0] > 0) else NA_real_,
      loser_recall = if (any(observed < 0)) mean(predicted[observed < 0] < 0) else NA_real_,
      .groups = "drop") |>
    mutate(winner_loser_skill = winner_recall + loser_recall - 1)
}
species_scores <- bind_rows(lapply(names(scored), function(cohort_name) {
  cohort_points <- filter(points, species_id %in% scored[[cohort_name]])
  overall_mean <- mean(cohort_points$observed[cohort_points$model == "Species mean"])
  score_species(cohort_points, overall_mean) |>
    mutate(scored_cohort = cohort_name, overall_mean = overall_mean)
})) |>
  mutate(model = factor(model, levels = model_levels))
write.csv(species_scores, file.path(out_dir, "species_scores.csv"), row.names = FALSE)

# Summaries. Equal-species R2: every species weighted equally, each around its
# own mean (1 - mean MSE / mean variance). Pooled R2: all species x gradient
# points around their grand mean, so it also rewards knowing which species
# win overall. No between-species score is reported: averaging a species'
# held-out predictions over all its gradients recombines every gradient (the
# mean of its leave-one-out species means is its overall mean), so it would
# not be a held-out test.
n_boot <- 2000L
set.seed(base_seed)
summarise_models <- function(cohort_name) {
  scores <- species_scores |> filter(scored_cohort == cohort_name)
  cohort_points <- points |> filter(species_id %in% scored[[cohort_name]])
  species <- sort(unique(scores$species_id))
  draws <- replicate(n_boot, sample(species, replace = TRUE))
  boot <- scores |>
    group_by(model) |>
    group_modify(~ {
      rows <- lapply(seq_len(n_boot), function(b) match(draws[, b], .x$species_id))
      data.frame(draw = seq_len(n_boot),
        R2 = vapply(rows, function(i) 1 - mean(.x$MSE[i]) / mean(.x$variance[i]), numeric(1)),
        R2_absolute = vapply(rows, function(i) 1 - mean(.x$MSE[i]) / mean(.x$variance_absolute[i]),
                             numeric(1)))
    }) |>
    ungroup()
  pooled <- cohort_points |>
    group_by(model) |>
    summarise(pooled_R2 = 1 - sum((predicted - observed)^2) / sum((observed - mean(observed))^2),
              .groups = "drop")
  summary <- scores |>
    group_by(model) |>
    summarise(n_species = n(), n_points = sum(n_gradients),
      equal_species_R2 = 1 - mean(MSE) / mean(variance), median_species_R2 = median(R2),
      share_species_R2_positive = mean(R2 > 0),
      equal_species_R2_absolute = 1 - mean(MSE) / mean(variance_absolute),
      median_species_R2_absolute = median(R2_absolute),
      mean_Spearman = mean(Spearman, na.rm = TRUE), n_Spearman_defined = sum(is.finite(Spearman)),
      mean_skill = mean(winner_loser_skill, na.rm = TRUE),
      n_skill_defined = sum(is.finite(winner_loser_skill)), .groups = "drop") |>
    left_join(boot |> group_by(model) |>
                summarise(R2_low = quantile(R2, 0.025), R2_high = quantile(R2, 0.975),
                          R2_absolute_low = quantile(R2_absolute, 0.025),
                          R2_absolute_high = quantile(R2_absolute, 0.975), .groups = "drop"),
              by = "model", relationship = "one-to-one") |>
    left_join(pooled, by = "model", relationship = "one-to-one") |>
    mutate(scored_cohort = cohort_name)
  list(summary = summary, boot = mutate(boot, scored_cohort = cohort_name))
}
summaries <- lapply(names(scored), summarise_models)
summary_table <- bind_rows(lapply(summaries, `[[`, "summary"))
boot <- bind_rows(lapply(summaries, `[[`, "boot"))

# Training R2 of forests refitted to all gradients with each model's most often
# tuned settings, scored the same way on the primary species (not held out).
modal_settings <- tuning |>
  count(model, min_node_size, mtry) |>
  group_by(model) |>
  slice_max(n, n = 1L, with_ties = FALSE) |>
  ungroup()
training <- list()
for (s in seq_along(forests)) {
  model <- names(forests)[[s]]
  predictors <- forests[[s]]
  setting <- modal_settings[modal_settings$model == model, ]
  data <- cohort
  if ("species_memory" %in% predictors) data <- add_memory(data, data)$train
  if ("species_id" %in% predictors) data$species_id <- factor(data$species_id)
  fit <- grow(data, predictors, setting$min_node_size, setting$mtry, outer_trees, base_seed + s,
              threads = n_cores)
  fitted_points <- data.frame(species_id = as.character(data$species_id), gradient_id = data$gradient_id,
                              observed = data$response, predicted = forecast(fit, data, predictors, n_cores)) |>
    filter(species_id %in% scored[[primary]]) |>
    group_by(species_id, gradient_id) |>
    summarise(observed = mean(observed), predicted = mean(predicted), .groups = "drop") |>
    group_by(species_id) |>
    summarise(MSE = mean((predicted - observed)^2), variance = mean((observed - mean(observed))^2),
              .groups = "drop")
  primary_mean <- species_scores$overall_mean[species_scores$scored_cohort == primary][1]
  absolute <- data.frame(species_id = as.character(data$species_id), gradient_id = data$gradient_id,
                         observed = data$response) |>
    filter(species_id %in% scored[[primary]]) |>
    group_by(species_id, gradient_id) |>
    summarise(observed = mean(observed), .groups = "drop") |>
    group_by(species_id) |>
    summarise(variance_absolute = mean((observed - primary_mean)^2), .groups = "drop")
  fitted_points <- left_join(fitted_points, absolute, by = "species_id", relationship = "one-to-one")
  training[[s]] <- data.frame(model, scored_cohort = primary,
    training_equal_species_R2 = 1 - mean(fitted_points$MSE) / mean(fitted_points$variance),
    training_equal_species_R2_absolute = 1 - mean(fitted_points$MSE) / mean(fitted_points$variance_absolute))
}
summary_table <- summary_table |>
  left_join(bind_rows(training) |> mutate(model = factor(model, levels = model_levels)),
            by = c("model", "scored_cohort"), relationship = "one-to-one")
write.csv(summary_table, file.path(out_dir, "species_based_R2.csv"), row.names = FALSE)

# Paired comparisons on the same species and bootstrap draws.
comparisons <- tribble(
  ~comparison,                          ~with,                                       ~without,
  "Species ID added",                   "Species ID + Abundance + Niche",            "Abundance + Niche",
  "Forest vs species mean",             "Species ID + Abundance + Niche",            "Species mean",
  "Climate added",                      "Species ID + Climate + Abundance + Niche",  "Species ID + Abundance + Niche",
  "Leakage-safe memory vs species ID",  "Species memory + Abundance + Niche",        "Species ID + Abundance + Niche")
comparison_table <- bind_rows(lapply(names(scored), function(cohort_name) {
  scores <- species_scores |> filter(scored_cohort == cohort_name)
  bind_rows(lapply(seq_len(nrow(comparisons)), function(i) {
    with <- scores |> filter(model == comparisons$with[i])
    without <- scores |> filter(model == comparisons$without[i])
    without <- without[match(with$species_id, without$species_id), ]
    boot_with <- boot |> filter(scored_cohort == cohort_name, model == comparisons$with[i])
    boot_without <- boot |> filter(scored_cohort == cohort_name, model == comparisons$without[i])
    paired <- match(boot_with$draw, boot_without$draw)
    change <- boot_with$R2 - boot_without$R2[paired]
    change_absolute <- boot_with$R2_absolute - boot_without$R2_absolute[paired]
    data.frame(scored_cohort = cohort_name, comparisons[i, ],
      R2_change = (1 - mean(with$MSE) / mean(with$variance)) - (1 - mean(without$MSE) / mean(without$variance)),
      R2_change_low = quantile(change, 0.025), R2_change_high = quantile(change, 0.975),
      R2_absolute_change = (1 - mean(with$MSE) / mean(with$variance_absolute)) -
        (1 - mean(without$MSE) / mean(without$variance_absolute)),
      R2_absolute_change_low = quantile(change_absolute, 0.025),
      R2_absolute_change_high = quantile(change_absolute, 0.975),
      n_species = nrow(with), n_species_lower_error = sum(with$MSE < without$MSE),
      Spearman_change = mean(with$Spearman - without$Spearman, na.rm = TRUE),
      n_species_higher_Spearman = sum(with$Spearman > without$Spearman, na.rm = TRUE),
      skill_change = mean(with$winner_loser_skill - without$winner_loser_skill, na.rm = TRUE),
      row.names = NULL)
  }))
}))
write.csv(comparison_table, file.path(out_dir, "comparisons.csv"), row.names = FALSE)

# 5. Figures -----------------------------------------------------------------
ink <- "#0b0b0b"
ink_2 <- "#52514e"
grid_colour <- "#e4e3df"
model_colours <- c("Species mean" = "#8a8984", "Abundance + Niche" = "#2a78d6",
  "Species ID + Abundance + Niche" = "#eb6834", "Species ID + Climate + Abundance + Niche" = "#1baf7a",
  "Species memory + Abundance + Niche" = "#6250d6")
theme_check <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = grid_colour, linewidth = 0.3),
        axis.text = element_text(colour = ink_2), axis.title = element_text(colour = ink_2),
        plot.title = element_text(colour = ink, face = "bold"),
        plot.subtitle = element_text(colour = ink_2), legend.position = "top",
        legend.title = element_blank(), strip.text = element_text(colour = ink, hjust = 0))
scale_labels <- c(within = "Within species (around each species' own mean)",
                  absolute = "Absolute (around the overall mean)")
primary_summary <- summary_table |> filter(scored_cohort == primary)
primary_scores <- species_scores |> filter(scored_cohort == primary)
cohort_text <- sprintf("%d species predicted with >= %d other gradients in training (%d species x gradient points)",
                       primary_summary$n_species[1], min_other_gradients, primary_summary$n_points[1])

r2_long <- bind_rows(
  primary_summary |> transmute(model, scale = scale_labels[["within"]], value = equal_species_R2,
                               low = R2_low, high = R2_high, median = median_species_R2),
  primary_summary |> transmute(model, scale = scale_labels[["absolute"]], value = equal_species_R2_absolute,
                               low = R2_absolute_low, high = R2_absolute_high, median = median_species_R2_absolute)) |>
  mutate(model = factor(model, levels = rev(model_levels)), scale = factor(scale, levels = scale_labels))
r2_plot <- ggplot(r2_long, aes(value, model, colour = model)) +
  geom_vline(xintercept = 0, colour = ink_2, linewidth = 0.4) +
  geom_linerange(aes(xmin = low, xmax = high), linewidth = 0.7) +
  geom_point(size = 2.6) +
  geom_text(aes(label = sprintf("%.3f (median species %.2f)", value, median)), colour = ink_2, size = 3,
            vjust = -1) +
  scale_colour_manual(values = model_colours, guide = "none") +
  facet_wrap(~ scale, scales = "free_x") +
  labs(x = "Equal-species R² across each species' held-out gradients (95% bootstrap over species)",
       y = NULL, title = "Predicting a species in an unseen gradient from its other gradients",
       subtitle = paste0("Mean-year contrast, leave one gradient out. ", cohort_text, ".")) +
  theme_check

comparison_long <- bind_rows(
  comparison_table |> transmute(scored_cohort, comparison, scale = scale_labels[["within"]],
    change = R2_change, low = R2_change_low, high = R2_change_high,
    label = sprintf("%+.3f; Spearman %+.2f; %d/%d species lower error", R2_change, Spearman_change,
                    n_species_lower_error, n_species)),
  comparison_table |> transmute(scored_cohort, comparison, scale = scale_labels[["absolute"]],
    change = R2_absolute_change, low = R2_absolute_change_low, high = R2_absolute_change_high,
    label = sprintf("%+.3f; skill %+.2f", R2_absolute_change, skill_change))) |>
  filter(scored_cohort == primary) |>
  mutate(comparison = factor(comparison, levels = rev(comparisons$comparison)),
         scale = factor(scale, levels = scale_labels))
comparison_plot <- ggplot(comparison_long, aes(change, comparison)) +
  geom_vline(xintercept = 0, colour = ink_2, linewidth = 0.4) +
  geom_linerange(aes(xmin = low, xmax = high), colour = ink_2, linewidth = 0.7) +
  geom_point(size = 2.6, colour = ink) +
  geom_text(aes(label = label), colour = ink_2, size = 2.8, vjust = -1.1) +
  facet_wrap(~ scale, scales = "free_x") +
  labs(x = "Change in equal-species R² (95% bootstrap over species)", y = NULL,
       title = "What species information adds", subtitle = cohort_text) +
  theme_check

gradients_plot <- primary_scores |>
  filter(model %in% c("Abundance + Niche", "Species ID + Abundance + Niche")) |>
  select(model, species_id, n_gradients, R2, R2_absolute) |>
  pivot_longer(c(R2, R2_absolute), names_to = "scale", values_to = "value") |>
  mutate(scale = factor(scale_labels[if_else(scale == "R2", "within", "absolute")], levels = scale_labels)) |>
  ggplot(aes(n_gradients, pmax(value, -3), colour = model)) +
  geom_hline(yintercept = 0, colour = ink_2, linewidth = 0.4) +
  geom_point(size = 2, alpha = 0.8, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5,
                                                                    jitter.height = 0, seed = 1)) +
  scale_colour_manual(values = model_colours) +
  facet_wrap(~ scale) +
  labs(x = "Gradients the species occurs in", y = "Species R² (values below −3 drawn at −3)",
       title = "Does more information about a species help?", subtitle = "One point per species") +
  theme_check + theme(panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.3))

cairo_pdf(file.path(out_dir, "species_informed.pdf"), width = 11, height = 6.5, onefile = TRUE)
print(r2_plot)
print(comparison_plot)
print(gradients_plot)
invisible(dev.off())

# Predicted vs observed: one held-out point per species x gradient and model,
# the 1:1 line, and each model's least-squares fit of observed on predicted
# (slope 1 = calibrated; flatter = predictions too extreme, steeper = too timid).
calibration_points <- points |>
  filter(species_id %in% scored[[primary]]) |>
  mutate(model = factor(model, levels = model_levels))
calibration_fits <- calibration_points |>
  group_by(model) |>
  summarise(intercept = coef(lm(observed ~ predicted))[[1]], slope = coef(lm(observed ~ predicted))[[2]],
            x_min = min(predicted), x_max = max(predicted), .groups = "drop") |>
  mutate(y_min = intercept + slope * x_min, y_max = intercept + slope * x_max) |>
  left_join(primary_summary |> select(model, equal_species_R2_absolute), by = "model",
            relationship = "one-to-one") |>
  mutate(model = factor(as.character(model), levels = model_levels))
model_legend <- setNames(sprintf("%s (slope %.2f)", calibration_fits$model, calibration_fits$slope),
                         as.character(calibration_fits$model))
axis_range <- range(c(calibration_points$observed, calibration_points$predicted))
calibration_plot <- ggplot(calibration_points, aes(predicted, observed, colour = model)) +
  geom_hline(yintercept = 0, colour = grid_colour, linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = grid_colour, linewidth = 0.4) +
  geom_abline(slope = 1, intercept = 0, colour = ink_2, linewidth = 0.5, linetype = "dashed") +
  geom_point(size = 1.3, alpha = 0.25) +
  geom_segment(data = calibration_fits, aes(x = x_min, xend = x_max, y = y_min, yend = y_max),
               linewidth = 1.1) +
  scale_colour_manual(values = model_colours, labels = model_legend) +
  guides(colour = guide_legend(ncol = 2, override.aes = list(alpha = 1, size = 2.5))) +
  coord_equal(xlim = axis_range, ylim = axis_range) +
  labs(x = "Held-out prediction (species × gradient mean)", y = "Observed contrast (species × gradient mean)",
       title = "Predicted vs observed, all models",
       subtitle = paste0(cohort_text, ".\nDashed: 1:1. Solid: least-squares fit of observed on predicted.")) +
  theme_check + theme(panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.3))
calibration_facets <- ggplot(calibration_points, aes(predicted, observed)) +
  geom_hline(yintercept = 0, colour = grid_colour, linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = grid_colour, linewidth = 0.4) +
  geom_abline(slope = 1, intercept = 0, colour = ink_2, linewidth = 0.5, linetype = "dashed") +
  geom_point(aes(colour = model), size = 1.4, alpha = 0.5) +
  geom_segment(data = calibration_fits, aes(x = x_min, xend = x_max, y = y_min, yend = y_max, colour = model),
               linewidth = 1.1) +
  geom_text(data = calibration_fits, aes(label = sprintf("slope %.2f\nabsolute R² %.2f", slope,
                                                          equal_species_R2_absolute)),
            x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, size = 3, colour = ink_2) +
  scale_colour_manual(values = model_colours, guide = "none") +
  coord_equal(xlim = axis_range, ylim = axis_range) +
  facet_wrap(~ model, nrow = 2) +
  labs(x = "Held-out prediction (species × gradient mean)", y = "Observed contrast (species × gradient mean)",
       title = "Predicted vs observed, by model", subtitle = paste0(cohort_text, ". Dashed: 1:1.")) +
  theme_check + theme(panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.3))
cairo_pdf(file.path(out_dir, "predicted_vs_observed.pdf"), width = 10, height = 8, onefile = TRUE)
print(calibration_plot)
print(calibration_facets)
invisible(dev.off())
write.csv(calibration_fits |> mutate(model = as.character(model)),
          file.path(out_dir, "predicted_vs_observed_fits.csv"), row.names = FALSE)

# Species-based densities, laid out as the earlier predictive_metric_densities.pdf:
# rows are comparisons, columns are metrics, densities run over species.
r2_floor <- -3
density_sets <- list(
  "Species information" = c("Species mean", "Abundance + Niche", "Species ID + Abundance + Niche"),
  "Climate and encoding" = c("Species ID + Abundance + Niche", "Species ID + Climate + Abundance + Niche",
                             "Species memory + Abundance + Niche"))
metric_labels <- c(R2 = "Species R², within species", R2_absolute = "Species R², absolute",
                   Spearman = "Spearman rank correlation", winner_loser_skill = "Winner-loser prediction skill")
r2_metrics <- metric_labels[c("R2", "R2_absolute")]
density_metrics <- bind_rows(lapply(names(density_sets), function(set_name) {
  primary_scores |>
    filter(model %in% density_sets[[set_name]]) |>
    mutate(set = set_name, position = match(as.character(model), density_sets[[set_name]]))
})) |>
  select(set, position, model, species_id, R2, R2_absolute, Spearman, winner_loser_skill) |>
  pivot_longer(c(R2, R2_absolute, Spearman, winner_loser_skill), names_to = "metric", values_to = "value") |>
  mutate(metric = factor(metric_labels[metric], levels = metric_labels),
         set = factor(set, levels = names(density_sets)))
density_counts <- density_metrics |>
  group_by(set, metric, position) |>
  summarise(n_defined = sum(is.finite(value)), n_floor = sum(value < r2_floor, na.rm = TRUE),
            .groups = "drop") |>
  arrange(set, metric, position) |>
  group_by(set, metric) |>
  summarise(label = paste0("n = ", paste(n_defined, collapse = " | "),
                           if (first(metric) %in% r2_metrics && any(n_floor > 0))
                             paste0("\nbelow −3: ", paste(n_floor, collapse = " | ")) else ""),
            .groups = "drop")
density_metrics <- density_metrics |>
  filter(is.finite(value)) |>
  mutate(value = if_else(metric %in% r2_metrics, pmax(value, r2_floor), value))
overall_lookup <- bind_rows(
  primary_summary |> transmute(model, metric = r2_metrics[["R2"]], overall = equal_species_R2,
                               median = median_species_R2),
  primary_summary |> transmute(model, metric = r2_metrics[["R2_absolute"]], overall = equal_species_R2_absolute,
                               median = median_species_R2_absolute)) |>
  mutate(metric = factor(metric, levels = metric_labels))
density_means <- density_metrics |>
  group_by(set, metric, model, position) |>
  summarise(mean_score = mean(value), .groups = "drop") |>
  left_join(overall_lookup, by = c("model", "metric"), relationship = "many-to-one") |>
  mutate(label = if_else(is.na(overall), sprintf("Mean = %.3f", mean_score),
                         sprintf("Median %.2f | Equal-species %.3f", median, overall)),
         label_line = 1.2 + (position - 1L) * 1.4)
write.csv(density_means |> mutate(across(c(set, metric, model), as.character)),
          file.path(out_dir, "density_means.csv"), row.names = FALSE)

cairo_pdf(file.path(out_dir, "predictive_metric_densities.pdf"), width = 24, height = 11, onefile = TRUE)
tiles <- primary_summary |>
  select(model, "Training, within" = training_equal_species_R2, "Held-out, within" = equal_species_R2,
         "Training, absolute" = training_equal_species_R2_absolute,
         "Held-out, absolute" = equal_species_R2_absolute) |>
  pivot_longer(-model, names_to = "score", values_to = "value") |>
  mutate(model = factor(model, levels = rev(model_levels)),
         score = factor(score, levels = c("Training, within", "Held-out, within",
                                          "Training, absolute", "Held-out, absolute")))
print(ggplot(tiles, aes(score, model)) +
  geom_tile(aes(fill = model), alpha = 0.18, colour = "white", linewidth = 2) +
  geom_text(aes(label = if_else(is.na(value), "—", sprintf("%.3f", value))), size = 10) +
  scale_fill_manual(values = model_colours, guide = "none") +
  scale_x_discrete(position = "top", expand = expansion(add = 0)) +
  scale_y_discrete(expand = expansion(add = 0)) +
  labs(x = NULL, y = NULL, title = "Equal-species R²: training and held-out", subtitle = cohort_text,
       tag = "A", caption = paste(
         "Training R2: forests refitted to all gradients (each model's most often tuned settings) predict the same rows; not held-out or out-of-bag.",
         "Held-out R2: leave one gradient out. Both average each species' experiments within each gradient, weight species equally, and centre each species",
         "on its own mean (within) or on the overall mean of all scored points (absolute). The species mean has no training fit. Cell colour identifies the model.",
         sep = "\n")) +
  theme_minimal(base_size = 25) +
  theme(panel.grid = element_blank(), axis.text = element_text(size = 22),
        plot.title = element_text(size = 30, face = "bold"), plot.subtitle = element_text(size = 22),
        plot.tag = element_text(size = 30, face = "bold"), plot.tag.position = "topleft",
        plot.caption = element_text(size = 15, hjust = 0), plot.margin = margin(20, 25, 20, 25)))
density_ranges <- density_metrics |>
  group_by(set, metric) |>
  summarise(low = if (first(metric) %in% r2_metrics) min(value) - 0.1 * diff(range(value)) else -1,
            high = 1, .groups = "drop") |>
  pivot_longer(c(low, high), values_to = "value")
print(ggplot(density_metrics, aes(value, fill = model, colour = model)) +
  geom_blank(data = density_ranges, aes(value, 0), inherit.aes = FALSE) +
  geom_density(data = filter(density_metrics, metric %in% r2_metrics),
               bounds = c(-Inf, 1), n = 4096, alpha = 0.25, linewidth = 0.9) +
  geom_density(data = filter(density_metrics, !metric %in% r2_metrics),
               bounds = c(-1, 1), n = 4096, alpha = 0.25, linewidth = 0.9) +
  geom_rug(sides = "b", alpha = 0.5, linewidth = 0.4, show.legend = FALSE) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey20", linewidth = 0.6) +
  geom_text(data = density_counts, aes(label = label), x = Inf, y = Inf,
            hjust = 1.05, vjust = 1.3, size = 5, inherit.aes = FALSE) +
  geom_text(data = density_means, aes(label = label, colour = model, vjust = label_line),
            x = -Inf, y = Inf, hjust = -0.05, size = 4.6, inherit.aes = FALSE, show.legend = FALSE) +
  facet_wrap(vars(set, metric), ncol = 4, scales = "free",
             labeller = label_wrap_gen(width = 50, multi_line = FALSE)) +
  scale_fill_manual(values = model_colours, name = NULL) +
  scale_colour_manual(values = model_colours, name = NULL) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.5), nrow = 1)) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 5), expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 4), expand = expansion(mult = c(0, 0.2))) +
  labs(x = NULL, y = "Density", title = "Species-based scores: predicting each species in its held-out gradients",
       subtitle = cohort_text, tag = "B",
       caption = paste(
         "Response: fitted warmed minus origin-control relative abundance at the mean model-frame year (emmeans default). Leave one gradient out; nothing from the held-out gradient is used.",
         "Each species' experiments are averaged within each held-out gradient; each species then gets two R2 values, one Spearman and one prediction skill across its gradients.",
         "Within-species R2 centres each species on its own mean (does the model know where the species does better or worse?); absolute R2 centres on the overall mean (does it also know the species' typical response?).",
         "Densities run over species. Species R2 values below -3 are drawn at -3 (counts in the corner). Equal-species = 1 - mean MSE / mean variance across species.",
         "Prediction skill = winner detection rate + loser detection rate - 1 across the species' gradients; undefined when the species only wins or only loses.",
         "Species mean = the species' mean in the training gradients: its within-species Spearman is -1 and within-species R2 1 - (n / (n - 1))^2 by construction.",
         "Forests: nested leave-one-gradient-out tuning, no weights; species ID ordered by mean training response. Species memory = shrunken species mean from other gradients.",
         sep = "\n")) +
  theme_bw(base_size = 20) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92"), panel.spacing = grid::unit(1, "lines"),
        strip.text = element_text(size = 16, face = "bold"),
        strip.background = element_rect(fill = "grey96", colour = "grey75"),
        axis.title.y = element_text(size = 22), axis.text = element_text(size = 15),
        plot.title = element_text(size = 24, face = "bold"), plot.subtitle = element_text(size = 20),
        plot.tag = element_text(size = 26, face = "bold"), plot.tag.position = "topleft",
        plot.caption = element_text(size = 13, hjust = 0), legend.position = "top",
        legend.text = element_text(size = 16), legend.key.width = grid::unit(1.1, "cm"),
        plot.margin = margin(15, 20, 15, 20)))
invisible(dev.off())

print(as.data.frame(summary_table), digits = 3)
print(as.data.frame(comparison_table), digits = 3)
