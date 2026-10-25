#!/usr/bin/env Rscript
# Run: Rscript R/additional_checks/investigations/community_and_prediction_investigations.R
# Species heatmaps, gradient prediction checks and controls-only CCA in one PDF.
library(targets)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(vegan)

out_dir <- "output/additional_checks/investigations"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Species occurrence, winner/loser counts and mean Delta AUC -------------
# Maximum-duration resident effects, including Arizona: no species model here.
effects <- tar_read(resident_auc_data)
classifications <- tar_read(winner_loser) |>
  ungroup() |>
  filter(horizon == "max", transform == "raw", contrast == "to_origin", pool == "resident") |>
  transmute(experiment_id = as.character(experiment_id),
    source_species = as.character(species), category = as.character(category))
effects <- effects |>
  left_join(classifications, by = c("experiment_id", "source_species"),
            relationship = "one-to-one")
stopifnot(!anyNA(effects$category))

# Parentheses count gradients, not plots. A species has one estimate per
# experiment; a gradient can contain several experiments.
species_order <- effects |>
  group_by(species_id) |>
  summarise(n_gradients = n_distinct(gradient_id), n_effects = n(), .groups = "drop") |>
  arrange(desc(n_gradients), desc(n_effects), species_id) |>
  mutate(rank = row_number(), label = paste0(species_id, " (", n_gradients, ")"))
species_cells <- effects |>
  group_by(species_id, gradient_id) |>
  summarise(n_effects = n(), winners = sum(category == "winner"),
    losers = sum(category == "loser"),
    mean_delta_auc = mean(auc_effect), .groups = "drop") |>
  mutate(winner_fraction = winners / n_effects, loser_fraction = losers / n_effects,
    annotation = paste0("W", winners, "/", n_effects, "\nL", losers, "/", n_effects)) |>
  left_join(species_order |> rename(n_total_effects = n_effects),
            by = "species_id", relationship = "many-to-one")
write.csv(species_order, file.path(out_dir, "species_occurrence_ranking.csv"), row.names = FALSE)
write.csv(species_cells, file.path(out_dir, "species_gradient_summary.csv"), row.names = FALSE)

# A blank cell means no resident estimate, not a measured zero effect.
# W/L use the original significant winner/loser labels (mvt-adjusted contrasts).
# Non-significant estimates remain in the denominator and in the mean colour.
gradients <- sort(unique(as.character(effects$gradient_id)))
colour_limit <- max(abs(species_cells$mean_delta_auc))
heatmap_pages <- list()
for (page in seq_len(ceiling(nrow(species_order) / 33))) {
  shown <- species_order |> slice(((page - 1L) * 33L + 1L):min(page * 33L, n()))
  plot_data <- expand_grid(species_id = shown$species_id, gradient_id = gradients) |>
    left_join(species_cells |> select(species_id, gradient_id, mean_delta_auc, annotation),
              by = c("species_id", "gradient_id"), relationship = "one-to-one") |>
    mutate(species_id = factor(species_id, levels = rev(shown$species_id)),
           gradient_id = factor(gradient_id, levels = gradients))
  heatmap_pages[[page]] <- ggplot(plot_data, aes(gradient_id, species_id, fill = mean_delta_auc)) +
    geom_tile(colour = "grey85", linewidth = 0.2) +
    geom_text(aes(label = annotation,
      colour = I(ifelse(!is.na(mean_delta_auc) & abs(mean_delta_auc) > 0.6 * colour_limit,
                        "white", "black"))), size = 2.8, na.rm = TRUE) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
      limits = c(-colour_limit, colour_limit), na.value = "grey92", name = "Mean Delta AUC",
      guide = guide_colourbar(barwidth = grid::unit(8,"cm"))) +
    scale_y_discrete(labels = setNames(shown$label, shown$species_id)) +
    labs(x = "Geographic gradient", y = NULL,
      title = paste0("Species ranked by gradient occurrence: ", min(shown$rank), "-", max(shown$rank)),
      caption = "Species (number of gradients). W = winners/estimates; L = losers/estimates. Colour = arithmetic mean Delta AUC.\nGrey = no resident estimate. Counts include non-significant effects; all 20 gradients are included.") +
    theme_minimal(base_size = 15) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 55, hjust = 1),
      axis.text.y = element_text(size = 12), legend.position = "bottom",
      plot.caption = element_text(hjust = 0, size = 12))
}

# 2. Why do Lautaret and Susalps appear easier to predict? ------------------
# Saved leave-gradient-out predictions; no forests are refitted.
focus <- c("FR_Lautaret", "DE_Susalps")

# 2.1 First test the premise across the existing RF trials -------------------
# Never pool prediction errors across different responses or holdout designs.
trial_folders <- c(
  Original = "random_forest_without_US_Arizona",
  Duration_normalized = "random_forest_duration_normalized_without_US_Arizona",
  Experiment_predictors = "random_forest_experiment_predictors_without_US_Arizona",
  Origin_climate = "random_forest_origin_climate_without_US_Arizona",
  Origin_destination_climate = "random_forest_origin_destination_climate_without_US_Arizona",
  Community_abundance = "random_forest_community_composition",
  Community_presence_absence = "random_forest_community_presence_absence",
  Fixed_horizon = "random_forest_fixed_horizon_without_US_Arizona")
trial_scores <- list()
for (trial in names(trial_folders)) {
  path <- file.path("output/additional_checks", trial_folders[[trial]], "random_forest_by_fold.csv")
  if (!file.exists(path)) next
  rows <- read.csv(path) |> filter(holdout == "Gradient")
  if (!"horizon" %in% names(rows)) rows$horizon <- "max"
  rows$horizon <- as.character(rows$horizon)
  trial_scores[[trial]] <- rows |>
    mutate(trial = if (trial == "Fixed_horizon") paste0(trial, "_", horizon) else trial)
}
trial_scores <- bind_rows(trial_scores) |>
  group_by(trial, model, method, scope, weighting) |>
  mutate(rank_within_model = min_rank(desc(R2pred)),
    n_gradients_scored = sum(is.finite(R2pred)),
    performance_percentile = 100 * (1 - (rank_within_model - 0.5) / n_gradients_scored)) |>
  ungroup()
trial_summary <- trial_scores |>
  group_by(trial, heldout_group, method, scope) |>
  summarise(n_model_scores = sum(is.finite(R2pred)),
    median_R2 = median(R2pred, na.rm = TRUE),
    fraction_positive = mean(R2pred > 0, na.rm = TRUE),
    median_percentile = median(performance_percentile, na.rm = TRUE), .groups = "drop")
write.csv(trial_scores, file.path(out_dir, "cross_trial_gradient_scores.csv"), row.names = FALSE)
write.csv(trial_summary, file.path(out_dir, "cross_trial_gradient_summary.csv"), row.names = FALSE)

# 2.2 Response spread, species overlap and climate support in original models
traits <- tar_read(predictor_specifications)
datasets <- list(Climate = tar_read(resident_auc_data), LDMC = traits$LDMC$data,
                 LMA = traits$LMA$data, Height = traits$Height$data)
for (model in names(datasets)) datasets[[model]] <- datasets[[model]] |>
  filter(as.character(gradient_id) != "US_Arizona")
original_predictions <- read.csv("output/additional_checks/random_forest_without_US_Arizona/random_forest_predictions.csv") |>
  filter(holdout == "Gradient")
diagnostics <- list()
support <- list()
for (model in names(datasets)) {
  data <- datasets[[model]]
  for (gradient in sort(unique(as.character(data$gradient_id)))) {
    train <- data |> filter(as.character(gradient_id) != gradient)
    test <- data |> filter(as.character(gradient_id) == gradient)
    known <- as.character(test$species_id) %in% as.character(train$species_id)
    frequencies <- table(as.character(train$species_id))
    training_occurrences <- as.numeric(frequencies[as.character(test$species_id)])
    training_occurrences[is.na(training_occurrences)] <- 0
    # Each experiment is one climate point; standardise using training experiments.
    train_climate <- train |> distinct(experiment_id, T_z, P_z, W_z)
    test_climate <- test |> distinct(experiment_id, T_z, P_z, W_z)
    train_matrix <- as.matrix(train_climate[, c("T_z", "P_z", "W_z")])
    test_matrix <- as.matrix(test_climate[, c("T_z", "P_z", "W_z")])
    centre <- colMeans(train_matrix); spread <- apply(train_matrix, 2, sd)
    train_matrix <- sweep(sweep(train_matrix, 2, centre), 2, spread, "/")
    test_matrix <- sweep(sweep(test_matrix, 2, centre), 2, spread, "/")
    nearest <- inside <- numeric(nrow(test_matrix))
    for (j in seq_len(nrow(test_matrix))) {
      nearest[j] <- min(sqrt(rowSums(sweep(train_matrix, 2, test_matrix[j, ])^2)))
      inside[j] <- all(test_matrix[j, ] >= apply(train_matrix, 2, min) &
                       test_matrix[j, ] <= apply(train_matrix, 2, max))
    }
    support[[paste(model, gradient)]] <- data.frame(model, gradient_id = gradient,
      n_effects = nrow(test), n_experiments = n_distinct(test$experiment_id),
      n_species = n_distinct(test$species_id), known_effect_fraction = mean(known),
      known_species_fraction = mean(unique(test$species_id) %in% train$species_id),
      mean_training_occurrences = mean(training_occurrences),
      nearest_climate_distance = mean(nearest), climate_range_fraction = mean(inside),
      observed_mean = mean(test$auc_effect), observed_SD = sd(test$auc_effect),
      median_SE = median(test$auc_se), mean_sampling_variance = mean(test$auc_se^2),
      median_duration = median(test$duration), median_plot_area = median(test$plot_area))
    for (method in c("RF_without_species", "RF_with_species")) {
      for (weighting in c("Unweighted", "Inverse variance")) {
        rows <- original_predictions |>
          filter(.data$model == .env$model, heldout_group == gradient,
                 .data$method == .env$method, .data$weighting == .env$weighting)
        for (scope in if (method == "RF_with_species") "Known species" else c("All species", "Known species")) {
          z <- if (scope == "Known species") rows |> filter(known_species) else rows
          y <- z$observed; yhat <- z$predicted
          SSE <- sum((y - yhat)^2); SST <- sum((y - mean(y))^2)
          # Compare with predictions that ignore all predictors: zero response,
          # or the training-set mean. These are diagnostics, not extra fitted models.
          zero_SSE <- sum(y^2); mean_SSE <- sum((y - mean(train$auc_effect))^2)
          largest <- if (length(y)) which.max(abs(y)) else integer()
          kept <- setdiff(seq_along(y), largest)
          trimmed_SST <- sum((y[kept] - mean(y[kept]))^2)
          diagnostics[[paste(model, gradient, method, weighting, scope)]] <- data.frame(
            model, gradient_id = gradient, method, weighting, scope, n = length(y),
            R2pred = if (SST > 0) 1 - SSE/SST else NA_real_, RMSE = sqrt(mean((y-yhat)^2)),
            response_SD = if (length(y)) sqrt(SST/length(y)) else NA_real_,
            prediction_SD = sd(yhat), observed_mean = mean(y), predicted_mean = mean(yhat),
            R2_zero = if (SST > 0) 1-zero_SSE/SST else NA_real_,
            R2_training_mean = if (SST > 0) 1-mean_SSE/SST else NA_real_,
            skill_vs_zero = if (zero_SSE > 0) 1-SSE/zero_SSE else NA_real_,
            skill_vs_training_mean = if (mean_SSE > 0) 1-SSE/mean_SSE else NA_real_,
            largest_effect_SST_fraction = if (SST > 0) max((y-mean(y))^2)/SST else NA_real_,
            R2_without_largest_effect = if (trimmed_SST > 0)
              1-sum((y[kept]-yhat[kept])^2)/trimmed_SST else NA_real_)
        }
      }
    }
  }
}
support <- bind_rows(support)
diagnostics <- bind_rows(diagnostics) |>
  left_join(support, by = c("model", "gradient_id"), suffix = c("", "_all_effects"),
            relationship = "many-to-one")
write.csv(support, file.path(out_dir, "gradient_support.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(out_dir, "gradient_prediction_diagnostics.csv"), row.names = FALSE)

# 2.3 Figures: relative rankings do not imply useful absolute prediction ------
performance_pages <- list()
trial_labels <- c(Original="Original",Duration_normalized="AUC / duration",
  Experiment_predictors="Duration + plot area",Origin_climate="Origin climate",
  Origin_destination_climate="Origin + destination climate",
  Community_abundance="Community: abundance",Community_presence_absence="Community: presence-absence",
  Fixed_horizon_2="AUC to year 2",Fixed_horizon_4="AUC to year 4")
for (scope in c("All species", "Known species")) {
  shown <- trial_summary |> filter(.data$scope == .env$scope)
  performance_pages[[length(performance_pages)+1L]] <- ggplot(shown,
      aes(trial, heldout_group, fill = median_percentile)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f", median_R2)), size = 3.3) +
    facet_wrap(~method, labeller = as_labeller(c(RF_without_species = "Without species identity",
                                                RF_with_species = "With species identity"))) +
    scale_fill_gradient(low = "white", high = "#70B7A8", limits = c(0,100), name = "Performance percentile") +
    scale_x_discrete(labels=trial_labels) +
    scale_y_discrete(labels = function(x) ifelse(x %in% focus, paste0(x," *"),x)) +
    labs(x = NULL, y = "Held-out gradient", title = paste("Across RF trials:",
      if (scope=="All species") "all held-out observations" else "held-out species represented in training"),
      caption = "Text: median R² across four predictor sets and two fitting weights. Colour: median within-model rank percentile (higher = better).\nRanks are relative to the other gradients in each trial; a high rank can still have negative R². * = focal gradients. No RF refitting.") +
    theme_minimal(base_size = 15) + theme(axis.text.x = element_text(angle=55,hjust=1),
      panel.grid=element_blank(), legend.position="bottom",plot.caption=element_text(hjust=0))
}
shown <- diagnostics |> filter(model=="Climate",method=="RF_without_species",
  weighting=="Unweighted",scope=="All species") |>
  mutate(focal=if_else(gradient_id %in% focus,gradient_id,"Other gradients"))
colours <- c(FR_Lautaret="#B54A43", DE_Susalps="#2166AC", "Other gradients"="grey60")
error_plot <- ggplot(shown,aes(response_SD,RMSE,colour=focal)) +
  geom_abline(slope=1,intercept=0,linetype=2) + geom_point(size=3) +
  geom_text_repel(aes(label=gradient_id),size=3.2,max.overlaps=Inf,seed=1) +
  scale_colour_manual(values=colours,name=NULL) +
  labs(x="Held-out response SD",y="Prediction RMSE",
    subtitle="Below dashed line: positive R²; above: negative R²") + theme_classic(base_size=17)
baseline_plot <- shown |> select(gradient_id,focal,R2pred,R2_zero,R2_training_mean) |>
  pivot_longer(starts_with("R2"),names_to="prediction",values_to="R2") |>
  ggplot(aes(R2,gradient_id,colour=prediction)) + geom_vline(xintercept=0,linetype=2) +
  geom_point(size=3) + labs(x="Predictive R²",y=NULL,colour=NULL) +
  scale_colour_manual(values=c(R2pred="#2166AC",R2_zero="grey55",R2_training_mean="#B54A43"),
    labels=c(R2pred="Climate forest",R2_zero="Always predict zero",R2_training_mean="Always predict training mean")) +
  theme_classic(base_size=16)
performance_pages[[length(performance_pages)+1L]] <- (error_plot + baseline_plot) +
  plot_annotation(title="Original climate forest: is apparent success mostly a near-zero response?",
    caption="Unweighted climate-only forests, all held-out observations. Both constant predictions use no held-out responses. R² uses the held-out mean only in its denominator.")

support_long <- shown |> select(gradient_id,focal,R2pred,known_effect_fraction,nearest_climate_distance,
  median_SE,median_duration,mean_training_occurrences,largest_effect_SST_fraction) |>
  pivot_longer(-c(gradient_id,focal,R2pred),names_to="measure",values_to="value")
measure_labels <- c(known_effect_fraction="Fraction of effects with known species",
  nearest_climate_distance="Nearest training climate distance",
  median_SE="Median response SE", median_duration="Median recorded experiment duration",
  mean_training_occurrences="Mean training occurrences per test species",
  largest_effect_SST_fraction="Largest share of response variation from one effect")
performance_pages[[length(performance_pages)+1L]] <- ggplot(support_long,aes(value,R2pred,colour=focal)) +
  geom_hline(yintercept=0,linetype=2) + geom_point(size=3) +
  geom_text_repel(data=filter(support_long,focal!="Other gradients"),aes(label=gradient_id),
    size=3.5,max.overlaps=Inf,seed=1) +
  facet_wrap(~measure,scales="free_x",ncol=3,labeller=as_labeller(measure_labels)) +
  scale_colour_manual(values=colours,name=NULL) + labs(x=NULL,y="Predictive R²",
    title="Possible explanations: species support, climate coverage and response information",
    caption="Descriptive comparisons across 19 gradients, not causal tests. Climate distance uses T, log(P), W scaled among training experiments.\nSampling SE measures estimation uncertainty, not biological variation. Removing one extreme effect is tabulated separately; no model was refitted.") +
  theme_classic(base_size=16) + theme(legend.position="bottom",plot.caption=element_text(hjust=0))

# 2.4 Which species actually benefit from adding species identity? -----------
# Paired errors on EXACTLY the same known-species observations. Training means
# below use other gradients only; they describe transferable response signs.
paired_rows <- list()
for (model in names(datasets)) {
  for (gradient in focus) {
    train_species <- datasets[[model]] |> filter(as.character(gradient_id) != gradient) |>
      group_by(species_id) |>
      summarise(training_n = n(), training_gradients = n_distinct(gradient_id),
        training_mean = mean(auc_effect),
        training_weighted_mean = weighted.mean(auc_effect, 1/auc_se^2), .groups="drop") |>
      mutate(species_id=as.character(species_id))
    z <- original_predictions |> filter(.data$model==.env$model, heldout_group==gradient, known_species) |>
      select(model,heldout_group,weighting,effect_id,species_id,observed,method,predicted) |>
      pivot_wider(names_from=method,values_from=predicted) |>
      left_join(train_species,by="species_id",relationship="many-to-one") |>
      mutate(error_without=(observed-RF_without_species)^2,
        error_with=(observed-RF_with_species)^2,
        error_reduction=error_without-error_with)
    paired_rows[[paste(model,gradient)]] <- z
  }
}
paired_rows <- bind_rows(paired_rows)
stopifnot(all(is.finite(paired_rows$RF_without_species)),all(is.finite(paired_rows$RF_with_species)))
paired_summary <- paired_rows |> group_by(model,heldout_group,weighting) |>
  summarise(n=n(), SST=sum((observed-mean(observed))^2),
    R2_without=1-sum(error_without)/SST, R2_with=1-sum(error_with)/SST,
    R2_gain=R2_with-R2_without, .groups="drop")
species_gains <- paired_rows |> group_by(model,heldout_group,weighting,species_id) |>
  summarise(n_effects=n(),observed_mean=mean(observed), predicted_with=mean(RF_with_species),
    predicted_without=mean(RF_without_species), error_reduction=sum(error_reduction),
    training_n=first(training_n),training_gradients=first(training_gradients),
    training_mean=first(training_mean),training_weighted_mean=first(training_weighted_mean),
    .groups="drop") |> arrange(model,heldout_group,weighting,desc(error_reduction))
write.csv(paired_summary,file.path(out_dir,"focal_species_identity_comparison.csv"),row.names=FALSE)
write.csv(species_gains,file.path(out_dir,"focal_species_error_contributions.csv"),row.names=FALSE)
write.csv(paired_rows,file.path(out_dir,"focal_paired_predictions.csv"),row.names=FALSE)
shown <- paired_rows |> filter(model=="Climate",weighting=="Inverse variance")
labels <- shown |> group_by(heldout_group) |> slice_max(abs(error_reduction),n=5,with_ties=FALSE) |> ungroup()
prediction_plot <- ggplot(shown,aes(observed,RF_with_species,colour=heldout_group)) +
  geom_abline(slope=1,intercept=0,linetype=2) + geom_hline(yintercept=0,colour="grey80") +
  geom_vline(xintercept=0,colour="grey80") + geom_point(size=2,alpha=.6) +
  geom_text_repel(data=labels,aes(label=species_id),size=3,max.overlaps=Inf,seed=1) +
  facet_wrap(~heldout_group,ncol=1) + scale_colour_manual(values=colours,guide="none") +
  labs(x="Observed Delta AUC",y="Predicted Delta AUC",subtitle="Species-inclusive climate forest; SE-weighted") +
  theme_classic(base_size=16)
gain_rows <- species_gains |> filter(model=="Climate",weighting=="Inverse variance") |>
  group_by(heldout_group) |> slice_max(abs(error_reduction),n=10,with_ties=FALSE) |> ungroup() |>
  mutate(species_panel=paste(species_id,heldout_group,sep="___"))
gain_rows$species_panel <- factor(gain_rows$species_panel,
  levels=gain_rows$species_panel[order(gain_rows$error_reduction)])
gain_plot <- ggplot(gain_rows,aes(error_reduction,species_panel,fill=heldout_group)) +
  geom_col() + geom_vline(xintercept=0,linetype=2) + facet_wrap(~heldout_group,ncol=1,scales="free_y") +
  scale_y_discrete(labels=setNames(gain_rows$species_id,as.character(gain_rows$species_panel))) +
  scale_fill_manual(values=colours,guide="none") +
  labs(x="Squared error without species - with species",y=NULL,
    subtitle="Ten most influential species per gradient") + theme_classic(base_size=14)
performance_pages[[length(performance_pages)+1L]] <- (prediction_plot + gain_plot) +
  plot_annotation(title="Does species identity improve predictions of the same held-out observations?",
    caption="Known species only; original climate forests with inverse-variance training weights. Right: positive bars mean species identity reduces error.\nEach species' contribution is summed across its held-out effects. These are descriptive contributions, not independent tests of species effects.",
    theme=theme(plot.title=element_text(size=20),plot.caption=element_text(size=14,hjust=0)))

# 2.5 Numerical findings, with limits on what they establish -----------------
findings <- c("LAUTARET AND SUSALPS: WHAT THE SAVED PREDICTIONS SHOW",
  "", "Higher rank does not necessarily mean useful prediction; compare R2 with zero.",
  "The statements below refer to the original maximum-duration RF trial.", "")
for (gradient in focus) {
  z <- diagnostics |> filter(gradient_id==gradient,method=="RF_with_species",
    weighting=="Inverse variance",scope=="Known species")
  c <- diagnostics |> filter(gradient_id==gradient,model=="Climate",
    method=="RF_without_species",weighting=="Unweighted",scope=="All species")
  s <- support |> filter(model=="Climate",gradient_id==gradient)
  influence <- z |> filter(model=="Climate")
  best <- species_gains |> filter(heldout_group==gradient,model=="Climate",weighting=="Inverse variance") |>
    slice_max(error_reduction,n=3,with_ties=FALSE)
  findings <- c(findings, gradient,
    sprintf("  Climate-only R2: %.3f. Always predicting zero gives R2: %.3f.",c$R2pred,c$R2_zero),
    sprintf("  Adding species (SE-weighted): R2 %.3f to %.3f across the four predictor sets.",min(z$R2pred),max(z$R2pred)),
    sprintf("  Training contains %.1f%% of test species; %.1f%% of test effects have known species.",
      100*s$known_species_fraction,100*s$known_effect_fraction),
    sprintf("  Nearest training-climate distance: %.3f (training-standardized T, log(P), W).",s$nearest_climate_distance),
    sprintf("  Removing the largest absolute response from SCORING changes species/climate R2 from %.3f to %.3f.",
      influence$R2pred,influence$R2_without_largest_effect),
    paste("  Biggest error reductions from species identity:",paste(best$species_id,collapse=", ")), "")
}
findings <- c(findings,
  "INTERPRETATION",
  "Lautaret: the weighted species models capture some repeat-species responses; climate alone does not.",
  "The signal is partly influenced by a large response, but remains positive when it is removed from scoring.",
  "Susalps: relatively favourable ranks mostly mean less-negative R2, not consistently successful prediction.",
  "Species overlap and proximity to training climates are possible contributors, not proof of a cause.",
  "Lautaret's advantage does not persist in the year-2/year-4 trials; those also change the eligible data.",
  "No extra forests or mixed models were fitted. Diagnostics use held-out data only for assessment.")
writeLines(findings,file.path(out_dir,"gradient_investigation_findings.txt"))
findings_page <- ggplot() + annotate("text",x=0,y=1,label=paste(findings,collapse="\n"),
  hjust=0,vjust=1,size=5.2,lineheight=1.15) + xlim(0,1) + ylim(0,1) + theme_void()
performance_pages <- c(list(findings_page),performance_pages)

# 3. Control community composition and climate ----------------------------
# 3.1 All included origin and destination controls; no warmed communities ----
# Presence-absence of reviewed taxa, including non-residents and Arizona.
# Retain species present in controls in at least six distinct gradients.
# One row per physical plot-year. Shared controls are not duplicated in the fit.
community <- tar_read(analysis_community)
experiments <- tar_read(experiments) |> filter(included)
taxonomy <- tar_read(community_taxonomy_species) |>
  distinct(name_original, accepted_name_reviewed)
memberships <- list()
for (i in seq_len(nrow(experiments))) {
  experiment <- experiments[i, ]
  for (role in c("Origin", "Destination")) {
    site <- if (role == "Origin") experiment$origin_site else experiment$destination_site
    memberships[[paste(i,role)]] <- community |>
      filter(region == experiment$region, origin_site == site,
             destination_site == site, treatment == "LocalControl") |>
      distinct(plot_year_id, stable_plot_id, region, destination_site, year) |>
      mutate(experiment_id = as.character(experiment$experiment_id), role = role)
  }
}
memberships <- bind_rows(memberships)
control_metadata <- memberships |>
  distinct(plot_year_id, stable_plot_id, region, destination_site, year)
control_species <- community |>
  filter(plot_year_id %in% control_metadata$plot_year_id, treatment == "LocalControl") |>
  left_join(taxonomy, by = c("species" = "name_original"), relationship = "many-to-one") |>
  group_by(plot_year_id, accepted_name_reviewed) |>
  summarise(cover = sum(relative_cover), .groups = "drop") |>
  mutate(presence = as.integer(cover > 0))
stopifnot(!anyNA(control_species$accepted_name_reviewed))
species_occurrence <- control_species |>
  filter(presence == 1L) |>
  left_join(control_metadata |> select(plot_year_id, region),
            by = "plot_year_id", relationship = "many-to-one") |>
  group_by(accepted_name_reviewed) |>
  summarise(n_gradients = n_distinct(region), .groups = "drop") |>
  mutate(retained = n_gradients >= 6L) |>
  arrange(desc(n_gradients), accepted_name_reviewed)
write.csv(species_occurrence, file.path(out_dir, "controls_cca_species_selection.csv"),
          row.names = FALSE)
community_table <- control_species |>
  semi_join(filter(species_occurrence, retained), by = "accepted_name_reviewed") |>
  select(plot_year_id, accepted_name_reviewed, presence) |>
  pivot_wider(names_from = accepted_name_reviewed, values_from = presence, values_fill = 0) |>
  arrange(plot_year_id)
community_matrix <- as.matrix(community_table[, -1])
rownames(community_matrix) <- community_table$plot_year_id
# CCA cannot use a community with no retained species.
community_matrix <- community_matrix[rowSums(community_matrix) > 0, , drop = FALSE]
control_metadata <- control_metadata |>
  mutate(included_in_cca = plot_year_id %in% rownames(community_matrix))
write.csv(control_metadata, file.path(out_dir, "controls_cca_community_selection.csv"),
          row.names = FALSE)
control_coverage <- memberships |>
  distinct(region, role, plot_year_id) |>
  left_join(control_metadata |> select(plot_year_id, included_in_cca),
            by = "plot_year_id", relationship = "many-to-one") |>
  group_by(region, role) |>
  summarise(n_plot_years = n(), n_retained = sum(included_in_cca), .groups = "drop")
write.csv(control_coverage, file.path(out_dir, "controls_cca_gradient_coverage.csv"),
          row.names = FALSE)

# 3.2 The same site's climate for either control role -----------------------
# Summer means over the common cached 2008-2021 period, not treatment warming,
# not observed weather during each sampling year, and no extrapolated years.
climate <- tar_read(site_year_climate) |>
  group_by(region, site_id) |>
  summarise(temperature = mean(T_JJA), precipitation = mean(P_JJA),
    climate_first_year = min(year), climate_last_year = max(year), .groups = "drop")
environment <- control_metadata |>
  left_join(climate, by = c("region", "destination_site" = "site_id"), relationship = "many-to-one") |>
  mutate(log_precipitation = log(precipitation))
environment <- environment[match(rownames(community_matrix), environment$plot_year_id), ]
stopifnot(!anyDuplicated(environment$plot_year_id),
  identical(rownames(community_matrix), environment$plot_year_id),
  all(is.finite(environment$temperature)), all(is.finite(environment$log_precipitation)))

# 3.3 Literal model: species presence-absence ~ site climate -----------------
cca_fit <- vegan::cca(community_matrix ~ temperature + log_precipitation,
                     data = environment)
stopifnot(cca_fit$CCA$rank == 2L)
saveRDS(cca_fit, file.path(out_dir, "controls_cca.rds"))
write.csv(environment, file.path(out_dir, "controls_cca_inputs.csv"), row.names = FALSE)
write.csv(memberships, file.path(out_dir, "controls_cca_membership.csv"), row.names = FALSE)
axis_inertia <- data.frame(axis = names(cca_fit$CCA$eig), eigenvalue = cca_fit$CCA$eig,
  percent_total_inertia = 100 * cca_fit$CCA$eig / cca_fit$tot.chi,
  percent_constrained_inertia = 100 * cca_fit$CCA$eig / sum(cca_fit$CCA$eig))
write.csv(axis_inertia, file.path(out_dir, "controls_cca_axis_inertia.csv"), row.names = FALSE)

# Community (weighted-average) scores, not sites predicted solely from climate.
# Mean coordinates give each observed plot-year equal weight within a role.
# If a control serves several experiments, count it once per gradient and role.
site_scores <- as.data.frame(scores(cca_fit, display = "sites", choices = 1:2, scaling = 1))
site_scores$plot_year_id <- rownames(site_scores)
centroids <- memberships |> distinct(region, role, plot_year_id) |>
  semi_join(site_scores, by = "plot_year_id") |>
  left_join(site_scores, by = "plot_year_id", relationship = "many-to-one") |>
  group_by(region, role) |>
  summarise(CCA1 = mean(CCA1), CCA2 = mean(CCA2), n_plot_years = n(), .groups = "drop")
links <- centroids |> select(region, role, CCA1, CCA2) |>
  pivot_wider(names_from = role, values_from = c(CCA1, CCA2)) |>
  filter(complete.cases(pick(everything())))
# Draw an arrow only where both control roles still have a centroid.
centroid_labels <- centroids |>
  arrange(region, desc(role == "Destination")) |>
  distinct(region, .keep_all = TRUE)
write.csv(centroids, file.path(out_dir, "controls_cca_gradient_centroids.csv"), row.names = FALSE)

# 3.4 Species contributions, not merely the most extreme plotted scores -----
# Standard species coordinates are weighted by the species' column mass.
# Contribution of species s to axis k = mass_s * score_sk^2 (sums to one).
# Inertia contribution = that proportion * the axis eigenvalue.
species_coordinates <- scores(cca_fit, display = "species", choices = 1:2, scaling = 0)
species_mass <- weights(cca_fit, display = "species")[rownames(species_coordinates)]
axis_contribution <- sweep(species_coordinates^2, 1, species_mass, "*")
stopifnot(max(abs(colSums(axis_contribution) - 1)) < 1e-8)
species_inertia <- sweep(axis_contribution, 2, cca_fit$CCA$eig[1:2], "*")
contributions <- data.frame(species = rownames(species_coordinates),
  CCA1_inertia = species_inertia[,1], CCA2_inertia = species_inertia[,2],
  CCA1_contribution_percent = 100 * axis_contribution[,1],
  CCA2_contribution_percent = 100 * axis_contribution[,2],
  combined_percent = 100 * rowSums(species_inertia) / sum(cca_fit$CCA$eig[1:2])) |>
  arrange(desc(combined_percent))
stopifnot(abs(sum(contributions$combined_percent)-100) < 1e-8)
write.csv(contributions, file.path(out_dir, "controls_cca_species_contributions.csv"), row.names = FALSE)
top_species <- head(contributions, 20)
bar_data <- top_species |> select(species, CCA1_inertia, CCA2_inertia) |>
  pivot_longer(-species, names_to = "axis", values_to = "inertia") |>
  mutate(percent = 100 * inertia / sum(cca_fit$CCA$eig[1:2]),
         species = factor(species, levels = rev(top_species$species)))
gradient_colours <- setNames(hcl.colors(n_distinct(centroids$region), "Dark 3"), sort(unique(centroids$region)))
ordination_plot <- ggplot(centroids, aes(CCA1, CCA2, colour = region)) +
  geom_hline(yintercept=0,colour="grey85") + geom_vline(xintercept=0,colour="grey85") +
  geom_segment(data=links, aes(x=CCA1_Origin,y=CCA2_Origin,xend=CCA1_Destination,yend=CCA2_Destination),
    inherit.aes=FALSE, colour=gradient_colours[links$region],
    arrow=grid::arrow(length=grid::unit(0.16,"cm")),linewidth=0.8) +
  geom_point(aes(shape=role),size=3.5) +
  geom_text_repel(data=centroid_labels,aes(label=region),
    size=3.5,max.overlaps=Inf,seed=1,min.segment.length=0) +
  scale_colour_manual(values=gradient_colours,guide="none") +
  scale_shape_manual(values=c(Origin=16,Destination=17),name=NULL) + coord_equal() +
  labs(x=sprintf("CCA1: %.2f%% of total inertia",axis_inertia$percent_total_inertia[1]),
       y=sprintf("CCA2: %.2f%% of total inertia",axis_inertia$percent_total_inertia[2]),tag="A") +
  theme_classic(base_size=17) + theme(legend.position="bottom")
contribution_plot <- ggplot(bar_data,aes(percent,species,fill=axis)) + geom_col() +
  scale_fill_manual(values=c(CCA1_inertia="#2166AC",CCA2_inertia="#B54A43"),
                    labels=c("CCA1","CCA2"),name=NULL) +
  labs(x="Contribution to CCA1 + CCA2 inertia (%)",y=NULL,tag="B") +
  theme_classic(base_size=15) + theme(legend.position="bottom")
cca_page <- (ordination_plot + contribution_plot) + plot_layout(widths=c(1.3,1)) +
  plot_annotation(title="Control communities ~ temperature + log(precipitation)",
    subtitle=sprintf("Climate accounts for %.2f%% of total community inertia; %s control plot-years, %s species.",
      100*cca_fit$CCA$tot.chi/cca_fit$tot.chi,nrow(community_matrix),ncol(community_matrix)),
    caption=paste0("Arrows: origin to destination control centroid within each gradient; points are means of community scores, not fitted climate scores.\nPresence-absence; species recorded in controls in at least six gradients; unique physical plot-years. ",
      sum(!control_metadata$included_in_cca), " plot-years omitted with no retained species.\nNo centroid after filtering: ",
      paste(with(filter(control_coverage, n_retained == 0), paste(region, role)),
            collapse = "; "),
      ".\nSite climate = mean summer CHELSA climate over 2008-2021. Bars: 20 largest eigenvalue-weighted species contributions to the two constrained axes.\nContributions are not percentages of each species' own variation. Descriptive CCA, no permutation tests."),
    theme=theme(plot.title=element_text(size=22),plot.subtitle=element_text(size=17),
                plot.caption=element_text(size=13,hjust=0)))

# 4. One PDF, with CCA and species contributions together on the last page ---
pdf(file.path(out_dir,"investigations.pdf"),width=22,height=14)
for (page in heatmap_pages) print(page)
for (page in performance_pages) print(page)
print(cca_page)
dev.off()
message("Saved ",file.path(out_dir,"investigations.pdf"))
