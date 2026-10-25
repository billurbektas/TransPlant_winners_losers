#!/usr/bin/env Rscript

# Does a species' first-year abundance at the origin and at the destination
# predict winners and losers in an unseen gradient?
# Run: Rscript R/additional_checks/origin_destination_abundance/random_forest_origin_destination_abundance.R
#
# Response: the warmed - origin relative-abundance contrasts of
# random_forest_year_group_contrasts_all_species.R (at the mean model-frame year)
# and random_forest_year4_and_max_contrasts_all_species.R (at year 4 and at each
# experiment's last year), read from their saved inputs; run those two first.
# US_Arizona is excluded, as there.
# Origin abundance: mean relative cover of the focal species over the origin-control
# plots in their first recorded year, surveyed absences as zero (identical to the
# all-species abundance check). These plots also enter the origin line of the contrast.
# Destination abundance: the same over the destination-control plots in their first
# recorded year. Those plots never enter the contrast; they show whether the species
# already grows in the warmer community. dest_present = destination abundance > 0.
#
# Forests carry no species ID and no weights: species ID lowered transfer in every
# earlier check, and contrast_se is constant within experiments. min.node.size and
# mtry are tuned by leave-one-gradient-out inside each outer training set, so the
# held-out gradient never informs its own tuning. Each forest is reported beside a
# linear model on the same predictors, the training mean, and an oracle that
# knows each held-out experiment's true mean.
#
# Independent-plot check: the 20 plot splits of
# random_forest_abundance_checks_all_species.R. Origin abundance comes from half A
# of the origin-control plots; each species x experiment origin line of the contrast
# is refitted from half B alone, the warmed line is unchanged, and the contrast is
# re-evaluated at the same year. Origin abundance and the response then share no
# plots. Each split reuses the outer folds' tuned settings.

library(targets)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ranger)
library(parallel)

out_dir <- "output/additional_checks/origin_destination_abundance"
dir.create(file.path(out_dir, "checkpoints"), recursive = TRUE, showWarnings = FALSE)
base_seed <- 20261024L
n_cores <- max(1L, detectCores() - 2L)
contrast_dir <- "output/additional_checks/random_forest_year_group_contrasts_all_species"
split_dir <- "output/additional_checks/random_forest_abundance_checks_all_species"

# 1. Effects, contrasts and first-year abundances ----------------------------
effects <- tar_read(resident_auc_data) |>
  filter(as.character(gradient_id) != "US_Arizona") |>
  mutate(across(c(experiment_id, gradient_id, species_id), as.character)) |>
  select(-auc_effect, -auc_se)
experiments <- tar_read(experiments) |>
  filter(included, experiment_id %in% effects$experiment_id)
community <- tar_read(analysis_community)

contrast_files <- c("Mean year" = "analysis_inputs.rds", "Year 4" = "analysis_inputs_year4.rds",
                    "Last year" = "analysis_inputs_max_year.rds")
contrast_names <- names(contrast_files)
contrast_inputs <- lapply(contrast_files, function(file) {
  path <- file.path(contrast_dir, file)
  if (!file.exists(path)) stop(path, " is missing: run the year-group contrast scripts first.")
  readRDS(path)
})
contrasts <- bind_rows(lapply(contrast_names, function(contrast) {
  contrast_inputs[[contrast]]$contrasts |>
    transmute(contrast, effect_id, reference_year = as.numeric(reference_year),
              fitted_mean_origin_controls, fitted_mean_warmed, response = abundance_change)
}))
for (contrast in contrast_names) {
  stopifnot(setequal(contrasts$effect_id[contrasts$contrast == contrast], effects$effect_id))
}

first_year_abundance <- function(experiment, site) {
  controls <- community |>
    filter(region == experiment$region, origin_site == site,
           destination_site == site, treatment == "LocalControl")
  controls <- controls |> filter(year == min(year))
  plots <- distinct(controls, plot_year_id)
  stopifnot(nrow(plots) > 0L)
  cover <- controls |>
    group_by(plot_year_id, source_species = species) |>
    summarise(relative_abundance = sum(relative_cover), .groups = "drop")
  effects |>
    filter(experiment_id == experiment$experiment_id) |>
    select(effect_id, source_species) |>
    cross_join(plots) |>
    left_join(cover, by = c("plot_year_id", "source_species"), relationship = "many-to-one") |>
    mutate(relative_abundance = coalesce(relative_abundance, 0)) |>
    group_by(effect_id) |>
    summarise(abundance = mean(relative_abundance), n_plots = n(),
              n_plots_present = sum(relative_abundance > 0), .groups = "drop") |>
    mutate(first_year = min(controls$year))
}
abundance <- list()
for (i in seq_len(nrow(experiments))) {
  experiment <- experiments[i, ]
  origin <- first_year_abundance(experiment, experiment$origin_site)
  destination <- first_year_abundance(experiment, experiment$destination_site)
  abundance[[i]] <- origin |>
    transmute(effect_id, first_origin_abundance = abundance, origin_first_year = first_year,
              n_origin_control_plots = n_plots) |>
    left_join(destination |> transmute(effect_id, dest_abundance = abundance,
                                       dest_first_year = first_year, n_dest_control_plots = n_plots,
                                       n_dest_plots_present = n_plots_present),
              by = "effect_id", relationship = "one-to-one")
}
abundance <- bind_rows(abundance) |> mutate(dest_present = as.numeric(dest_abundance > 0))
stopifnot(!anyDuplicated(abundance$effect_id), nrow(abundance) == nrow(effects))
cohort <- effects |>
  left_join(abundance, by = "effect_id", relationship = "one-to-one") |>
  as.data.frame()
gradients <- sort(unique(cohort$gradient_id))
write.csv(abundance, file.path(out_dir, "first_year_abundance.csv"), row.names = FALSE)
write.csv(cohort |> summarise(n_effects = n(), n_species = n_distinct(species_id),
    n_experiments = n_distinct(experiment_id), n_gradients = n_distinct(gradient_id),
    share_present_at_destination = mean(dest_present),
    share_absent_at_origin = mean(first_origin_abundance == 0)),
  file.path(out_dir, "cohort_counts.csv"), row.names = FALSE)
write.csv(contrasts |> group_by(contrast) |>
            summarise(min_reference_year = min(reference_year),
                      max_reference_year = max(reference_year), .groups = "drop"),
          file.path(out_dir, "contrast_years.csv"), row.names = FALSE)
saveRDS(list(cohort = cohort, contrasts = contrasts), file.path(out_dir, "analysis_inputs.rds"))
# One data frame per contrast; `response` is that contrast.
response_data <- lapply(setNames(contrast_names, contrast_names), function(contrast) {
  cohort |>
    left_join(contrasts |> filter(contrast == .env$contrast) |> select(effect_id, response),
              by = "effect_id", relationship = "one-to-one")
})

# 2. Specifications ---------------------------------------------------------
climate <- c("T_z", "P_z", "W_z")
origin <- "first_origin_abundance"
destination <- c("dest_abundance", "dest_present")
specifications <- list(
  "Climate" = climate,
  "Origin abundance" = origin,
  "Destination abundance" = destination,
  "Origin + destination abundance" = c(origin, destination),
  "Climate + origin abundance" = c(climate, origin),
  "Climate + destination abundance" = c(climate, destination),
  "Climate + origin + destination abundance" = c(climate, origin, destination))

# 3. Nested leave-one-gradient-out forests ----------------------------------
# Half-samples without replacement decorrelate trees and stop each tree from
# memorising repeated rows; tuning searches node sizes from the earlier checks'
# 5 up to leaves of several hundred effects.
node_grid <- c(5L, 25L, 75L, 200L, 500L)
mtry_grid <- function(p) sort(unique(pmin(p, c(1L, floor(sqrt(p)), ceiling(p / 2), p))))
inner_trees <- 200L
outer_trees <- 1000L
split_trees <- 500L
n_permutations <- 5L

# The scripts' equal-gradient R2: every held-out gradient gets equal weight.
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
         importance = "none", seed = seed, num.threads = threads)
}
forecast <- function(fit, data, predictors, threads = 1L) {
  predict(fit, data = data[, predictors, drop = FALSE], num.threads = threads)$predictions
}
# Score every grid cell by leave-one-gradient-out within `train` alone. All cells
# share one seed per inner fold, so they differ only by their settings.
tune <- function(train, predictors, seed) {
  grid <- expand.grid(min_node_size = node_grid, mtry = mtry_grid(length(predictors)))
  inner <- sort(unique(train$gradient_id))
  grid$inner_R2 <- vapply(seq_len(nrow(grid)), function(k) {
    predicted <- numeric(nrow(train))
    for (j in seq_along(inner)) {
      held <- train$gradient_id == inner[j]
      fit <- grow(train[!held, ], predictors, grid$min_node_size[k], grid$mtry[k],
                  inner_trees, seed + j)
      predicted[held] <- forecast(fit, train[held, ], predictors)
    }
    equal_gradient_R2(train$response, predicted, train$gradient_id)
  }, numeric(1))
  grid
}

ranger_version <- as.character(packageVersion("ranger"))
fold_results <- list()
for (k in seq_along(contrast_names)) {
  contrast <- contrast_names[[k]]
  data <- response_data[[contrast]]
  for (s in seq_along(specifications)) {
    specification <- names(specifications)[[s]]
    predictors <- specifications[[s]]
    key <- paste(contrast, specification, sep = " | ")
    checkpoint <- file.path(out_dir, "checkpoints", paste0(gsub("[^A-Za-z0-9]+", "_", key), ".rds"))
    signature <- digest::digest(list(
      data = data[c("effect_id", "gradient_id", "response", predictors)],
      predictors = predictors, node_grid = node_grid, mtry = "1, sqrt, p/2, p",
      inner_trees = inner_trees, outer_trees = outer_trees, permutations = n_permutations,
      sample_fraction = 0.5, replace = FALSE, seed = base_seed + 10000L * k + 1000L * s,
      ranger_version = ranger_version))
    if (file.exists(checkpoint) && identical(readRDS(checkpoint)$signature, signature)) {
      fold_results[[key]] <- readRDS(checkpoint)
      next
    }
    message("Nested forests: ", key)
    folds <- mclapply(seq_along(gradients), function(g) {
      held <- data$gradient_id == gradients[[g]]
      train <- data[!held, ]
      test <- data[held, ]
      seed <- base_seed + 10000L * k + 1000L * s + 10L * g
      grid <- tune(train, predictors, seed)
      # Ties go to the larger node size, then the smaller mtry (the smoother forest).
      best <- grid[order(-grid$inner_R2, -grid$min_node_size, grid$mtry), ][1, ]
      fit <- grow(train, predictors, best$min_node_size, best$mtry, outer_trees, seed)
      linear <- lm(reformulate(predictors, "response"), data = train)
      # Transfer importance: shuffle one predictor among the held-out gradient's rows.
      set.seed(seed + 7L)
      permuted <- list()
      for (predictor in predictors) {
        for (r in seq_len(n_permutations)) {
          shuffled <- test
          shuffled[[predictor]] <- sample(test[[predictor]])
          permuted[[paste(predictor, r)]] <- data.frame(predictor = predictor, permutation = r,
            effect_id = test$effect_id, predicted = forecast(fit, shuffled, predictors))
        }
      }
      list(predictions = data.frame(effect_id = test$effect_id,
             "Random forest" = forecast(fit, test, predictors),
             "Linear model" = unname(predict(linear, newdata = test)), check.names = FALSE),
           grid = mutate(grid, heldout_gradient = gradients[[g]]),
           tuning = data.frame(heldout_gradient = gradients[[g]], min_node_size = best$min_node_size,
                               mtry = best$mtry, inner_R2 = best$inner_R2),
           permuted = bind_rows(permuted))
    }, mc.cores = n_cores, mc.preschedule = FALSE)
    failed <- vapply(folds, inherits, logical(1), "try-error")
    if (any(failed)) stop(key, ": ", folds[failed][[1]])
    saved <- list(signature = signature,
      predictions = bind_rows(lapply(folds, `[[`, "predictions")),
      grid = bind_rows(lapply(folds, `[[`, "grid")),
      tuning = bind_rows(lapply(folds, `[[`, "tuning")),
      permuted = bind_rows(lapply(folds, `[[`, "permuted")))
    stopifnot(nrow(saved$predictions) == nrow(data),
              setequal(saved$predictions$effect_id, data$effect_id))
    saveRDS(saved, paste0(checkpoint, ".tmp"))
    stopifnot(file.rename(paste0(checkpoint, ".tmp"), checkpoint))
    fold_results[[key]] <- saved
  }
}
result_index <- expand.grid(specification = names(specifications), contrast = contrast_names,
                            stringsAsFactors = FALSE) |>
  mutate(key = paste(contrast, specification, sep = " | "))
stack_results <- function(part) {
  bind_rows(lapply(seq_len(nrow(result_index)), function(i) {
    mutate(fold_results[[result_index$key[i]]][[part]],
           contrast = result_index$contrast[i], specification = result_index$specification[i])
  }))
}

keys <- bind_rows(lapply(contrast_names, function(contrast) {
  response_data[[contrast]] |>
    transmute(contrast, effect_id, gradient_id, experiment_id, species_id, observed = response)
}))
predictions <- stack_results("predictions") |>
  pivot_longer(c(`Random forest`, `Linear model`), names_to = "method", values_to = "predicted")
# References: the training mean of the other gradients, and the oracle.
reference <- keys |>
  group_by(contrast) |>
  group_modify(~ {
    training_mean <- vapply(.x$gradient_id, function(g) mean(.x$observed[.x$gradient_id != g]), numeric(1))
    bind_rows(
      data.frame(effect_id = .x$effect_id, specification = "Training mean",
                 method = "Reference", predicted = training_mean),
      data.frame(effect_id = .x$effect_id, specification = "Oracle: held-out experiment mean",
                 method = "Reference", predicted = ave(.x$observed, .x$experiment_id)))
  }) |>
  ungroup()
predictions <- bind_rows(predictions, reference) |>
  left_join(keys, by = c("contrast", "effect_id"), relationship = "many-to-one")
saveRDS(predictions, file.path(out_dir, "predictions.rds"))
tuning <- stack_results("tuning")
write.csv(tuning, file.path(out_dir, "tuning.csv"), row.names = FALSE)
write.csv(stack_results("grid"), file.path(out_dir, "tuning_grid.csv"), row.names = FALSE)

# 4. Scores -----------------------------------------------------------------
# R2 on every species x experiment pair of a held-out gradient; Spearman and
# winner/loser skill on species means within that gradient (as in earlier checks).
by_gradient <- predictions |>
  group_by(contrast, specification, method, gradient_id) |>
  summarise(n_effects = n(), MSE = mean((predicted - observed)^2),
            mean_y = mean(observed), mean_y2 = mean(observed^2),
            R2pred = 1 - sum((predicted - observed)^2) / sum((observed - mean(observed))^2),
            .groups = "drop") |>
  left_join(predictions |>
      group_by(contrast, specification, method, gradient_id, species_id) |>
      summarise(observed = mean(observed), predicted = mean(predicted), .groups = "drop") |>
      group_by(contrast, specification, method, gradient_id) |>
      summarise(n_species = n(),
        Spearman = if (n() >= 2L && n_distinct(observed) > 1L && n_distinct(predicted) > 1L)
          cor(observed, predicted, method = "spearman") else NA_real_,
        winner_recall = if (any(observed > 0)) mean(predicted[observed > 0] > 0) else NA_real_,
        loser_recall = if (any(observed < 0)) mean(predicted[observed < 0] < 0) else NA_real_,
        .groups = "drop") |>
      mutate(winner_loser_skill = winner_recall + loser_recall - 1),
    by = c("contrast", "specification", "method", "gradient_id"), relationship = "one-to-one")
by_gradient <- by_gradient |>
  left_join(by_gradient |> filter(specification == "Climate") |>
              select(contrast, method, gradient_id, climate_MSE = MSE),
            by = c("contrast", "method", "gradient_id"), relationship = "many-to-one")
write.csv(by_gradient, file.path(out_dir, "metrics_by_gradient.csv"), row.names = FALSE)

# Whole-gradient bootstrap: resample the 19 gradients, recompute equal-gradient R2,
# and its change from Climate with the same method on the same draws.
n_boot <- 2000L
set.seed(base_seed)
draws <- replicate(n_boot, sample(gradients, replace = TRUE))
boot_R2 <- function(stats) {
  apply(draws, 2, function(draw) {
    s <- stats[match(draw, stats$gradient_id), ]
    1 - mean(s$MSE) / (mean(s$mean_y2) - mean(s$mean_y)^2)
  })
}
boot <- by_gradient |>
  group_by(contrast, specification, method) |>
  group_modify(~ data.frame(draw = seq_len(n_boot), R2 = boot_R2(.x))) |>
  ungroup()
boot <- boot |>
  left_join(boot |> filter(specification == "Climate") |>
              select(contrast, method, draw, climate_R2 = R2),
            by = c("contrast", "method", "draw"), relationship = "many-to-one")
overall <- by_gradient |>
  group_by(contrast, specification, method) |>
  summarise(n_gradients = n(), n_effects = sum(n_effects),
    R2pred = 1 - mean(MSE) / (mean(mean_y2) - mean(mean_y)^2),
    Spearman = mean(Spearman, na.rm = TRUE), n_Spearman_defined = sum(is.finite(Spearman)),
    winner_loser_skill = mean(winner_loser_skill, na.rm = TRUE),
    n_gradients_better_than_climate = sum(MSE < climate_MSE), .groups = "drop") |>
  left_join(boot |>
      group_by(contrast, specification, method) |>
      summarise(R2_low = quantile(R2, 0.025), R2_high = quantile(R2, 0.975),
                change_vs_climate_low = quantile(R2 - climate_R2, 0.025, na.rm = TRUE),
                change_vs_climate_high = quantile(R2 - climate_R2, 0.975, na.rm = TRUE), .groups = "drop"),
    by = c("contrast", "specification", "method"), relationship = "one-to-one")
overall <- overall |>
  left_join(overall |> filter(specification == "Climate") |>
              select(contrast, method, climate_R2 = R2pred),
            by = c("contrast", "method"), relationship = "many-to-one") |>
  mutate(change_vs_climate = R2pred - climate_R2) |>
  select(-climate_R2)
reference_rows <- overall$method == "Reference"
overall[reference_rows, c("n_gradients_better_than_climate", "change_vs_climate",
                          "change_vs_climate_low", "change_vs_climate_high")] <- NA
# Training R2 of forests refitted to all gradients with each specification's most
# frequently chosen settings; the gap to held-out R2 shows how much is memorised.
modal_settings <- tuning |>
  count(contrast, specification, min_node_size, mtry) |>
  group_by(contrast, specification) |>
  slice_max(n, n = 1L, with_ties = FALSE) |>
  ungroup()
full_fits <- list()
training <- list()
for (i in seq_len(nrow(result_index))) {
  contrast <- result_index$contrast[i]
  specification <- result_index$specification[i]
  setting <- modal_settings[modal_settings$contrast == contrast &
                              modal_settings$specification == specification, ]
  fit <- grow(response_data[[contrast]], specifications[[specification]], setting$min_node_size,
              setting$mtry, outer_trees, base_seed + i, threads = n_cores)
  full_fits[[result_index$key[i]]] <- fit
  fitted <- forecast(fit, response_data[[contrast]], specifications[[specification]], n_cores)
  training[[i]] <- data.frame(contrast, specification, method = "Random forest",
    training_R2 = equal_gradient_R2(response_data[[contrast]]$response, fitted,
                                    response_data[[contrast]]$gradient_id))
}
overall <- overall |>
  left_join(bind_rows(training), by = c("contrast", "specification", "method"),
            relationship = "one-to-one")
write.csv(overall, file.path(out_dir, "equal_gradient_R2.csv"), row.names = FALSE)

# Drop in held-out equal-gradient R2 when one predictor is shuffled.
heldout_R2 <- overall |>
  filter(method == "Random forest") |>
  select(contrast, specification, heldout_R2 = R2pred)
importance <- stack_results("permuted") |>
  left_join(keys, by = c("contrast", "effect_id"), relationship = "many-to-one") |>
  group_by(contrast, specification, predictor, permutation) |>
  summarise(R2 = equal_gradient_R2(observed, predicted, gradient_id), .groups = "drop") |>
  left_join(heldout_R2, by = c("contrast", "specification"), relationship = "many-to-one") |>
  group_by(contrast, specification, predictor, heldout_R2) |>
  summarise(R2_drop = mean(heldout_R2 - R2), R2_drop_min = min(heldout_R2 - R2),
            R2_drop_max = max(heldout_R2 - R2), .groups = "drop")
write.csv(importance, file.path(out_dir, "transfer_importance.csv"), row.names = FALSE)

# 5. Direction of the abundance effects --------------------------------------
# Partial dependence from the all-gradient forests: the mean prediction with one
# predictor set to each value and the others left as observed.
full_specification <- "Climate + origin + destination abundance"
pd_values <- function(x) {
  sort(unique(c(0, quantile(x[x > 0], seq(0, 0.975, length.out = 30L), names = FALSE))))
}
pd <- list()
for (contrast in contrast_names) {
  X <- response_data[[contrast]][, specifications[[full_specification]], drop = FALSE]
  for (variable in c(origin, "dest_abundance")) {
    values <- pd_values(cohort[[variable]])
    pd[[paste(contrast, variable)]] <- data.frame(contrast, specification = full_specification,
      variable, value = values, predicted = vapply(values, function(v) {
        X[[variable]] <- v
        if (variable == "dest_abundance") X$dest_present <- as.numeric(v > 0)
        mean(forecast(full_fits[[paste(contrast, full_specification, sep = " | ")]], X,
                      specifications[[full_specification]], n_cores))
      }, numeric(1)))
  }
}
pd <- bind_rows(pd)
write.csv(pd, file.path(out_dir, "partial_dependence.csv"), row.names = FALSE)
# Observed means: absent species, then deciles of the species present.
binned <- bind_rows(lapply(contrast_names, function(contrast) {
  bind_rows(lapply(c(origin, "dest_abundance"), function(variable) {
    x <- response_data[[contrast]][[variable]]
    response_data[[contrast]] |>
      mutate(bin = if_else(x > 0, ntile(if_else(x > 0, x, NA_real_), 10L), 0L)) |>
      group_by(bin) |>
      summarise(value = median(.data[[variable]]), observed = mean(response),
                se = sd(response) / sqrt(n()), n_effects = n(), .groups = "drop") |>
      mutate(contrast, variable)
  }))
}))
write.csv(binned, file.path(out_dir, "observed_by_abundance_decile.csv"), row.names = FALSE)

# Linear models on all gradients with standardised predictors; CR2 standard
# errors clustered by gradient, because experiments within a gradient share sites.
coefficient_specifications <- c(full_specification, "Origin + destination abundance")
coefficients <- list()
for (contrast in contrast_names) {
  for (specification in coefficient_specifications) {
    predictors <- specifications[[specification]]
    scaled <- response_data[[contrast]]
    scaled[predictors] <- lapply(scaled[predictors], function(x) (x - mean(x)) / sd(x))
    fit <- lm(reformulate(predictors, "response"), data = scaled)
    test <- clubSandwich::coef_test(fit, vcov = "CR2", cluster = scaled$gradient_id,
                                    test = "Satterthwaite")
    coefficients[[paste(contrast, specification)]] <- data.frame(contrast, specification,
      term = test$Coef, estimate = test$beta, SE = test$SE, df = test$df_Satt,
      p_value = test$p_Satt, row.names = NULL)
  }
}
coefficients <- bind_rows(coefficients)
write.csv(coefficients, file.path(out_dir, "linear_coefficients.csv"), row.names = FALSE)

# 6. Independent-plot check -------------------------------------------------
plot_data <- contrast_inputs[["Mean year"]]$plot_data |>
  filter(experimental_group == "origin_controls", effect_id %in% cohort$effect_id) |>
  select(effect_id, stable_plot_id, experimental_year, relative_abundance)
assignments <- read.csv(file.path(split_dir, "split_plot_assignments.csv"), stringsAsFactors = FALSE)
split_abundance <- read.csv(file.path(split_dir, "split_responses.csv"), stringsAsFactors = FALSE) |>
  select(split, effect_id, split_origin_abundance = first_origin_abundance)
stopifnot(all(plot_data$stable_plot_id %in% assignments$stable_plot_id),
          setequal(split_abundance$effect_id, cohort$effect_id))
# The origin line of each species x experiment cell, by least squares as in the
# contrast script, evaluated at each contrast's reference year.
origin_lines <- function(rows) {
  rows |>
    group_by(effect_id) |>
    summarise(n_years = n_distinct(experimental_year), mean_year = mean(experimental_year),
      mean_abundance = mean(relative_abundance),
      Sxx = sum((experimental_year - mean(experimental_year))^2),
      Sxy = sum((experimental_year - mean(experimental_year)) *
                  (relative_abundance - mean(relative_abundance))), .groups = "drop") |>
    mutate(year_slope = Sxy / Sxx)
}
refit_contrasts <- function(lines) {
  stopifnot(nrow(lines) == nrow(cohort), all(lines$n_years >= 2L), all(lines$Sxx > 0))
  contrasts |>
    left_join(lines, by = "effect_id", relationship = "many-to-one") |>
    mutate(refitted_origin = mean_abundance + year_slope * (reference_year - mean_year),
           refitted_response = fitted_mean_warmed - refitted_origin)
}
# With every plot, the refit must reproduce the saved contrasts.
check <- refit_contrasts(origin_lines(plot_data))
stopifnot(max(abs(check$refitted_origin - check$fitted_mean_origin_controls)) < 1e-10,
          max(abs(check$refitted_response - check$response)) < 1e-10)
splits <- sort(unique(assignments$split))
split_scores <- mclapply(splits, function(split) {
  B <- assignments$stable_plot_id[assignments$split == split & assignments$half == "B"]
  refitted <- refit_contrasts(origin_lines(filter(plot_data, stable_plot_id %in% B)))
  rows <- list()
  for (k in seq_along(contrast_names)) {
    contrast <- contrast_names[[k]]
    data <- cohort |>
      select(-first_origin_abundance) |>
      left_join(refitted |> filter(contrast == .env$contrast) |>
                  select(effect_id, response = refitted_response),
                by = "effect_id", relationship = "one-to-one") |>
      left_join(split_abundance |> filter(split == .env$split) |>
                  select(effect_id, first_origin_abundance = split_origin_abundance),
                by = "effect_id", relationship = "one-to-one")
    for (s in seq_along(specifications)) {
      specification <- names(specifications)[[s]]
      predictors <- specifications[[s]]
      predicted <- matrix(NA_real_, nrow(data), 2L,
                          dimnames = list(NULL, c("Random forest", "Linear model")))
      for (g in seq_along(gradients)) {
        held <- data$gradient_id == gradients[[g]]
        setting <- tuning[tuning$contrast == contrast & tuning$specification == specification &
                            tuning$heldout_gradient == gradients[[g]], ]
        fit <- grow(data[!held, ], predictors, setting$min_node_size, setting$mtry, split_trees,
                    base_seed + 100000L * split + 10000L * k + 1000L * s + 10L * g)
        predicted[held, "Random forest"] <- forecast(fit, data[held, ], predictors)
        predicted[held, "Linear model"] <- predict(lm(reformulate(predictors, "response"),
                                                      data = data[!held, ]), newdata = data[held, ])
      }
      rows[[paste(contrast, specification)]] <- data.frame(split, contrast, specification,
        method = colnames(predicted),
        R2pred = apply(predicted, 2, function(p) equal_gradient_R2(data$response, p, data$gradient_id)))
    }
    training_mean <- numeric(nrow(data))
    for (g in gradients) training_mean[data$gradient_id == g] <- mean(data$response[data$gradient_id != g])
    rows[[paste(contrast, "Training mean")]] <- data.frame(split, contrast,
      specification = "Training mean", method = "Reference",
      R2pred = equal_gradient_R2(data$response, training_mean, data$gradient_id))
  }
  bind_rows(rows)
}, mc.cores = n_cores, mc.preschedule = FALSE)
failed <- vapply(split_scores, inherits, logical(1), "try-error")
if (any(failed)) stop("Split check: ", split_scores[failed][[1]])
split_scores <- bind_rows(split_scores)
rownames(split_scores) <- NULL
# Paired changes on the same split responses.
paired <- split_scores |>
  filter(method != "Reference") |>
  select(split, contrast, method, specification, R2pred) |>
  pivot_wider(names_from = specification, values_from = R2pred) |>
  transmute(split, contrast, method,
    "Origin added to climate" = `Climate + origin abundance` - Climate,
    "Destination added to climate" = `Climate + destination abundance` - Climate,
    "Origin added to climate + destination" =
      `Climate + origin + destination abundance` - `Climate + destination abundance`,
    "Destination added to climate + origin" =
      `Climate + origin + destination abundance` - `Climate + origin abundance`) |>
  pivot_longer(-c(split, contrast, method), names_to = "comparison", values_to = "R2_change")
split_summary <- paired |>
  group_by(contrast, method, comparison) |>
  summarise(n_splits = n(), n_improved = sum(R2_change > 0), median_change = median(R2_change),
            min_change = min(R2_change), max_change = max(R2_change), .groups = "drop")
write.csv(split_scores, file.path(out_dir, "split_R2.csv"), row.names = FALSE)
write.csv(paired, file.path(out_dir, "split_paired_changes.csv"), row.names = FALSE)
write.csv(split_summary, file.path(out_dir, "split_summary.csv"), row.names = FALSE)

# 7. Figures ----------------------------------------------------------------
ink <- "#0b0b0b"
ink_2 <- "#52514e"
grid_colour <- "#e4e3df"
method_colours <- c("Random forest" = "#2a78d6", "Linear model" = "#eb6834", "Reference" = "#8a8984")
theme_check <- theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = grid_colour, linewidth = 0.3),
        axis.text = element_text(colour = ink_2), axis.title = element_text(colour = ink_2),
        plot.title = element_text(colour = ink, face = "bold"),
        plot.subtitle = element_text(colour = ink_2), legend.position = "top",
        legend.title = element_blank(), strip.text = element_text(colour = ink, hjust = 0))
contrast_labels <- contrasts |>
  group_by(contrast) |>
  summarise(label = if (min(reference_year) == max(reference_year))
    paste0(contrast[1], " (", format(round(min(reference_year), 2), nsmall = 2), ")") else
    paste0(contrast[1], " (", min(reference_year), "–", max(reference_year), ")"), .groups = "drop")
label_contrast <- function(x) {
  factor(contrast_labels$label[match(x, contrast_labels$contrast)],
         levels = contrast_labels$label[match(contrast_names, contrast_labels$contrast)])
}
specification_levels <- rev(c("Training mean", "Oracle: held-out experiment mean", names(specifications)))
oracle <- overall |> filter(specification == "Oracle: held-out experiment mean")
response_label <- "Warmed − origin relative abundance"

r2_plot <- overall |>
  mutate(specification = factor(specification, levels = specification_levels),
         method = factor(method, levels = names(method_colours)),
         contrast = label_contrast(contrast)) |>
  ggplot(aes(R2pred, specification, colour = method)) +
  geom_vline(xintercept = 0, colour = ink_2, linewidth = 0.4) +
  geom_vline(data = mutate(oracle, contrast = label_contrast(contrast)),
             aes(xintercept = R2pred), colour = ink_2, linewidth = 0.4, linetype = "dashed") +
  geom_linerange(aes(xmin = R2_low, xmax = R2_high), linewidth = 0.6,
                 position = position_dodge(width = 0.6)) +
  geom_point(size = 2, position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = method_colours) +
  facet_wrap(~ contrast, nrow = 1) +
  labs(x = "Equal-gradient predictive R², leave one gradient out", y = NULL,
       title = "Does first-year abundance at origin and destination predict the warmed − origin contrast?",
       subtitle = paste0(nrow(cohort), " effects, ", n_distinct(cohort$species_id), " species, ",
                         length(gradients), " gradients. Bars: 95% whole-gradient bootstrap. ",
                         "Dashed line: oracle ceiling for experiment-level predictors.")) +
  theme_check

importance_plot <- importance |>
  filter(specification == full_specification) |>
  mutate(contrast = label_contrast(contrast)) |>
  ggplot(aes(R2_drop, predictor)) +
  geom_vline(xintercept = 0, colour = ink_2, linewidth = 0.4) +
  geom_linerange(aes(xmin = R2_drop_min, xmax = R2_drop_max), colour = method_colours[["Random forest"]],
                 linewidth = 0.6) +
  geom_point(colour = method_colours[["Random forest"]], size = 2) +
  facet_wrap(~ contrast, nrow = 1) +
  labs(x = "Drop in held-out R² when the predictor is shuffled within the held-out gradient",
       y = NULL, title = "Which predictors carry the transfer?",
       subtitle = paste0(full_specification, "; mean and range over five shuffles, nested forests")) +
  theme_check

variable_labels <- c(first_origin_abundance = "First-year origin-control abundance",
                     dest_abundance = "First-year destination-control abundance")
pd_plot <- ggplot() +
  geom_hline(yintercept = 0, colour = ink_2, linewidth = 0.4) +
  geom_pointrange(data = mutate(binned, variable = factor(variable_labels[variable], levels = variable_labels),
                                contrast = label_contrast(contrast)),
                  aes(value, observed, ymin = observed - 2 * se, ymax = observed + 2 * se),
                  colour = method_colours[["Reference"]], size = 0.2, linewidth = 0.4) +
  geom_line(data = mutate(pd, variable = factor(variable_labels[variable], levels = variable_labels),
                          contrast = label_contrast(contrast)),
            aes(value, predicted), colour = method_colours[["Random forest"]], linewidth = 0.7) +
  scale_x_sqrt() +
  facet_grid(contrast ~ variable, scales = "free_x") +
  labs(x = "Relative abundance (square-root scale)", y = response_label,
       title = "Direction of the abundance effects",
       subtitle = paste0("Line: partial dependence of the all-gradient forest (", full_specification,
                         "). Grey: observed mean ± 2 SE for absent species and deciles of present species.")) +
  theme_check + theme(panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.3))

coefficient_plot <- coefficients |>
  filter(term != "(Intercept)") |>
  mutate(specification = factor(specification, levels = coefficient_specifications),
         contrast = label_contrast(contrast)) |>
  ggplot(aes(estimate, term)) +
  geom_vline(xintercept = 0, colour = ink_2, linewidth = 0.4) +
  geom_linerange(aes(xmin = estimate - qt(0.975, df) * SE, xmax = estimate + qt(0.975, df) * SE),
                 colour = method_colours[["Linear model"]], linewidth = 0.6) +
  geom_point(colour = method_colours[["Linear model"]], size = 2) +
  facet_grid(specification ~ contrast, scales = "free_y", space = "free_y") +
  labs(x = "Change in the contrast per SD of the predictor (95% CI, CR2 clustered by gradient)",
       y = NULL, title = "Linear abundance effects on all gradients") +
  theme_check

split_plot <- paired |>
  mutate(contrast = label_contrast(contrast)) |>
  ggplot(aes(R2_change, comparison, colour = method)) +
  geom_vline(xintercept = 0, colour = ink_2, linewidth = 0.4) +
  geom_point(size = 1.8, alpha = 0.7, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = method_colours) +
  facet_wrap(~ contrast, nrow = 1) +
  labs(x = "Change in held-out R² on the same split response", y = NULL,
       title = "Independent-plot check: origin abundance from half A, origin line of the contrast from half B",
       subtitle = "One point per split (20 splits). Destination plots never enter the contrast.") +
  theme_check

pdf(file.path(out_dir, "origin_destination_abundance.pdf"), width = 11, height = 7)
print(r2_plot)
print(importance_plot)
print(pd_plot)
print(coefficient_plot)
print(split_plot)
invisible(dev.off())

# 8. Per-gradient score densities, laid out as the earlier checks'
# predictive_metric_densities.pdf: rows are contrasts, columns are metrics, and
# each page overlays three specifications fitted with one method.
density_sets <- list(
  "Abundance alone" = c("Origin abundance", "Destination abundance", "Origin + destination abundance"),
  "Origin abundance added to climate" =
    c("Climate", "Climate + origin abundance", "Climate + origin + destination abundance"),
  "Destination abundance added to climate" =
    c("Climate", "Climate + destination abundance", "Climate + origin + destination abundance"))
density_caption <- paste(
  "Response: fitted warmed minus origin-control relative abundance (year-group contrast scripts) at the mean model-frame year, year 4, and each experiment's last year.",
  "Leave one gradient out; forests tuned by an inner leave-one-gradient-out loop, no species ID, no weights. Linear models use the same predictors.",
  "Within each held-out gradient, R2 compares every species x experiment contrast with its prediction. Rank and category scores use one mean contrast per species per gradient.",
  "Mean = average of valid gradient scores. Overall = equal-gradient R2 across all effect pairs, not the mean of gradient R2 values.",
  "Prediction skill = winner detection rate + loser detection rate - 1; winner/loser uses the contrast sign, not significance. Zero is the reference for all three metrics.",
  "Origin / destination abundance = first recorded year of origin / destination controls. Automatic density smoothing; ticks show gradient scores; n gives valid counts in legend order. US_Arizona excluded.",
  sep = "\n")
plot_density_pages <- function(by_gradient, overall, density_sets, caption, file) {
  series_colours <- c("#2a78d6", "#eb6834", "#1baf7a")
  metric_labels <- c(R2pred = "Predictive R²", Spearman = "Spearman rank correlation",
                     winner_loser_skill = "Winner-loser prediction skill")
  metrics <- by_gradient |>
    filter(method != "Reference") |>
    select(contrast, specification, method, gradient_id, R2pred, Spearman, winner_loser_skill) |>
    pivot_longer(c(R2pred, Spearman, winner_loser_skill), names_to = "metric", values_to = "value") |>
    mutate(metric = factor(metric_labels[metric], levels = metric_labels),
           contrast = label_contrast(contrast))
  theme_density <- theme_bw(base_size = 20) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "grey92"), panel.spacing = grid::unit(1, "lines"),
      strip.text = element_text(size = 17, face = "bold"),
      strip.background = element_rect(fill = "grey96", colour = "grey75"),
      axis.title.y = element_text(size = 22), axis.text = element_text(size = 16),
      plot.title = element_text(size = 24, face = "bold"), plot.subtitle = element_text(size = 20),
      plot.tag = element_text(size = 26, face = "bold"), plot.tag.position = "topleft",
      plot.caption = element_text(size = 13, hjust = 0), legend.position = "top",
      legend.text = element_text(size = 19), legend.key.width = grid::unit(1.1, "cm"),
      plot.margin = margin(15, 20, 15, 20))
  pdf(file, width = 20, height = 13)
  # Training R2 predicts the rows used to fit each all-gradient forest; it is
  # not held-out or out-of-bag. Held-out R2 is the leave-one-gradient-out score.
  tiles <- overall |>
    filter(method == "Random forest") |>
    mutate(contrast = label_contrast(contrast),
           specification = factor(specification, levels = rev(names(specifications))))
  print(ggplot(tiles, aes(contrast, specification)) +
    geom_tile(fill = "grey96", colour = "white", linewidth = 2) +
    geom_text(aes(label = sprintf("%.3f\nheld-out %.3f", training_R2, R2pred)), size = 6.5,
              lineheight = 0.9) +
    scale_x_discrete(position = "top", expand = expansion(add = 0)) +
    scale_y_discrete(expand = expansion(add = 0)) +
    labs(x = NULL, y = NULL, title = "Full-training R² of the random forests", tag = "A",
         subtitle = "Top: training R² of the all-gradient forest. Bottom: leave-one-gradient-out R².",
         caption = paste(
           "Training R2: predict the same observations used to fit each all-gradient forest (each specification's most often tuned settings). Not held-out or out-of-bag.",
           "R2 = 1 - SSE / SST with each gradient weighted equally, on individual species x experiment effects. US_Arizona is excluded.",
           sep = "\n")) +
    theme_minimal(base_size = 20) +
    theme(panel.grid = element_blank(), axis.text = element_text(size = 18),
          plot.title = element_text(size = 24, face = "bold"), plot.subtitle = element_text(size = 19),
          plot.tag = element_text(size = 26, face = "bold"), plot.tag.position = "topleft",
          plot.caption = element_text(size = 13, hjust = 0), plot.margin = margin(15, 20, 15, 20)))
  page_number <- 1L
  summaries <- list()
  for (method in c("Random forest", "Linear model")) {
    for (set_name in names(density_sets)) {
      set <- density_sets[[set_name]]
      page_number <- page_number + 1L
      colours <- setNames(series_colours[seq_along(set)], set)
      page <- metrics |>
        filter(method == .env$method, specification %in% set) |>
        mutate(specification = factor(specification, levels = set))
      counts <- page |>
        group_by(contrast, metric, specification) |>
        summarise(n_defined = sum(is.finite(value)), .groups = "drop") |>
        arrange(contrast, metric, specification) |>
        group_by(contrast, metric) |>
        summarise(label = paste0("n = ", paste(n_defined, collapse = " | ")), .groups = "drop")
      page <- page |> filter(is.finite(value))
      means <- page |>
        group_by(contrast, metric, specification) |>
        summarise(mean_score = mean(value), .groups = "drop") |>
        left_join(overall |>
            filter(method == .env$method, specification %in% set) |>
            transmute(contrast = label_contrast(contrast),
                      specification = factor(specification, levels = set),
                      metric = factor(metric_labels[["R2pred"]], levels = metric_labels),
                      overall_R2 = R2pred),
          by = c("contrast", "metric", "specification"), relationship = "one-to-one") |>
        mutate(label = if_else(is.na(overall_R2), sprintf("Mean = %.3f", mean_score),
                               sprintf("Mean = %.3f | Overall = %.3f", mean_score, overall_R2)),
               label_line = 1.2 + (as.integer(specification) - 1L) * 1.4)
      summaries[[paste(method, set_name)]] <- mutate(means, method, comparison = set_name)
      ranges <- page |>
        group_by(contrast, metric) |>
        summarise(low = if (first(metric) == metric_labels[["R2pred"]])
                    min(value) - 0.1 * diff(range(value)) else -1,
                  high = 1, .groups = "drop") |>
        pivot_longer(c(low, high), values_to = "value")
      print(ggplot(page, aes(value, fill = specification, colour = specification)) +
        geom_blank(data = ranges, aes(value, 0), inherit.aes = FALSE) +
        geom_density(data = filter(page, metric == metric_labels[["R2pred"]]),
                     bounds = c(-Inf, 1), n = 4096, alpha = 0.25, linewidth = 0.9) +
        geom_density(data = filter(page, metric != metric_labels[["R2pred"]]),
                     bounds = c(-1, 1), n = 4096, alpha = 0.25, linewidth = 0.9) +
        geom_rug(sides = "b", alpha = 0.5, linewidth = 0.4, show.legend = FALSE) +
        geom_vline(xintercept = 0, linetype = "dashed", colour = "grey20", linewidth = 0.6) +
        geom_text(data = counts, aes(label = label), x = Inf, y = Inf,
                  hjust = 1.1, vjust = 1.3, size = 5, inherit.aes = FALSE) +
        geom_text(data = means, aes(label = label, colour = specification, vjust = label_line),
                  x = -Inf, y = Inf, hjust = -0.08, size = 4.8, inherit.aes = FALSE, show.legend = FALSE) +
        facet_wrap(vars(contrast, metric), ncol = 3, scales = "free",
                   labeller = label_wrap_gen(width = 50, multi_line = FALSE)) +
        scale_fill_manual(values = colours, name = NULL) +
        scale_colour_manual(values = colours, name = NULL) +
        guides(fill = guide_legend(override.aes = list(alpha = 0.5))) +
        scale_x_continuous(breaks = scales::breaks_pretty(n = 5), expand = expansion(mult = c(0.02, 0.02))) +
        scale_y_continuous(breaks = scales::breaks_pretty(n = 4), expand = expansion(mult = c(0, 0.2))) +
        labs(x = NULL, y = "Density", title = paste0(method, "s: ", set_name),
             subtitle = "Warmed − origin relative-abundance contrasts, every held-out species",
             tag = LETTERS[[page_number]], caption = caption) +
        theme_density)
    }
  }
  invisible(dev.off())
  bind_rows(summaries)
}
density_means <- plot_density_pages(by_gradient, overall, density_sets, density_caption,
                                    file.path(out_dir, "predictive_metric_densities.pdf"))
write.csv(density_means |> mutate(across(c(contrast, metric, specification), as.character)),
          file.path(out_dir, "density_means.csv"), row.names = FALSE)

print(overall |>
  arrange(contrast, method, desc(R2pred)) |>
  select(contrast, specification, method, R2pred, R2_low, R2_high, change_vs_climate,
         n_gradients_better_than_climate, Spearman, winner_loser_skill, training_R2) |>
  as.data.frame(), digits = 3)
print(as.data.frame(split_summary), digits = 3)
