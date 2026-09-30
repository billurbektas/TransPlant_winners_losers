# MODEL ANALYSIS: READ-ONLY SUMMARY
# =================================
#
# This file brings the model analysis together in normal script order so that
# the data choices, formulas, random effects, comparisons, and classifications
# can be reviewed in one place. The workflow itself is run with targets, which
# is useful for reproducibility but harder to read from beginning to end. Nothing
# in this file runs because all code is inside `if (FALSE)`.

if (FALSE) {

  # TERMS USED BELOW
  # ----------------
  # Delta AUC / auc_effect
  #   Difference in maximum-duration, raw relative-cover area under the time
  #   curve: warmed treatment minus the origin control.
  #
  # T / T_z
  #   Destination-site June-August temperature / its standardized value.
  #   Standardization uses the mean and SD across distinct experiments, so an
  #   experiment with more species does not receive extra weight.
  #
  # P / P_z
  #   Natural log of destination-site June-August precipitation / its
  #   standardized value across distinct experiments.
  #
  # W / W_z
  #   June-August temperature difference representing experimental warming /
  #   its standardized value across distinct experiments.
  #
  # TP_z
  #   T_z * P_z. In the random effects this represents differences among species
  #   or families in the temperature-by-precipitation response.
  #
  # sampling_variance_matrix
  #   A diagonal matrix containing auc_se^2 in effect-row order. The models
  #   assume independent sampling errors between Delta AUC estimates.
  #
  # GEN
  #   An unstructured covariance for the species or family T_z, P_z, and TP_z
  #   random coefficients: their three variances and three pairwise covariances
  #   are estimated. The intercept and W_z random effects are separate blocks,
  #   so their covariance with the T/P/TP block is fixed at zero.
  #
  # CR1
  #   Small-sample-adjusted, cluster-robust fixed-effect uncertainty, clustered
  #   by experiment. It changes standard errors, confidence intervals, tests,
  #   and denominator degrees of freedom, but not coefficient estimates.
  #
  # ML / REML
  #   Maximum likelihood / restricted maximum likelihood. ML is used to compare
  #   the climate-only and full fixed-effect formulas. REML is used for final
  #   coefficients, random effects, predictions, and heterogeneity.
  #
  # PEV1_z, PEV2_z, PEV3_z
  #   Standardized coordinates on the first three phylogenetic eigenvectors.
  #   They describe phylogenetic position; they are not measured plant traits.


  # 1. PREPARING THE RESIDENT-SPECIES RESPONSE DATA
  # =================================================
  # We use maximum-duration, raw relative-cover AUC contrasts for resident
  # species. One row represents one accepted species in one experiment. The
  # current validated data contain 1,604 effects from 526 species and 42
  # experiments.
  #
  # These six deep-lineage taxa are currently retained while relative cover and
  # Delta AUC are calculated, then excluded from the resident models. They are
  # incompatible with the focal functional comparison and strongly influenced
  # the early phylogenetic axes. Juniperus communis is the only gymnosperm in the
  # analysed vegetation; the other five are non-flowering vascular lineages.
  resident_model_exclusions <- c(
    "Gymnocarpium dryopteris",
    "Athyrium spinulosum",
    "Botrychium lunaria",
    "Equisetum arvense",
    "Selaginella selaginoides",
    "Juniperus communis"
  )
  #
  # A possible future change is to remove these taxa before relative-cover and
  # AUC calculation. That would change the response itself, so it is deliberately
  # not done here. It requires a separate scientific decision and a full rebuild.
  #
  # Climate variables are destination summer temperature (T), log destination
  # summer precipitation (P), and the summer warming difference (W). They are
  # standardized across distinct experiments before being copied to species rows.
  resident_auc_data <- winner_loser_results |>
    dplyr::filter(
      transform == "raw",
      contrast == "to_origin",
      pool == "resident",
      horizon == "max"
    ) |>
    dplyr::transmute(
      experiment_id,
      source_species = species,
      auc_effect = estimate,
      auc_se = SE
    ) |>
    dplyr::left_join(experiment_and_species_metadata) |>
    dplyr::mutate(
      T = destination_T_JJA,
      P = log(destination_P_JJA),
      W = diff_T_JJA,
      effect_id = paste("max", experiment_id, species_id, sep = "::")
    ) |>
    dplyr::filter(!species_id %in% resident_model_exclusions)
  #
  # Pipeline code:
  #   config/analysis.R                  resident_model_exclusions()
  #   config/workflow.R                  resident_analysis_spec()
  #   R/responses/response_interfaces.R  build_resident_auc_data_from_objects()
  #   R/responses/auc_model_helpers.R    scale_resident_auc_data()


  # 2. SAMPLING VARIANCES
  # ======================
  # The winner/loser analysis has already estimated each warmed-minus-origin
  # Delta AUC and its standard error. The meta-regressions use that information
  # directly: sampling variance is auc_se^2, in exactly the same order as the
  # effect rows, with zero off-diagonals.
  resident_sampling_variance_matrix <- build_sampling_variance_matrix(
    resident_auc_data,
    se_column = "auc_se",
    id_column = "effect_id"
  )
  #
  # This means that sampling errors for different species effects are assumed
  # independent, including effects from the same experiment. Dependence within
  # experiments is addressed separately for fixed-effect inference by CR1
  # clustering in Section 3. The plot-level AUC model is not refitted here and
  # no off-diagonal sampling covariance is reconstructed.
  #
  # The original winner/loser classification contrasts retain their existing
  # `adjust = "mvt"` multiplicity adjustment. That adjustment concerns the
  # classification contrasts and is unrelated to this sampling-variance matrix.
  #
  # Pipeline code:
  #   R/responses/auc_model_helpers.R  build_sampling_variance_matrix()
  #   R/responses/auc_model_helpers.R  validate_sampling_variance_matrix()
  #   R/responses/auc_model_helpers.R  subset_sampling_variance_matrix()


  # 3. RESIDENT CLIMATE MODEL
  # =========================
  # Question: how does the mean warming response in AUC vary with destination
  # temperature, destination precipitation, and warming magnitude, while
  # allowing species to differ in their overall and climate-dependent responses?
  resident_climate_formula <- ~ T_z * P_z * W_z

  resident_species_random_effects <- list(
    ~ 1 | species_id,
    ~ 0 + W_z | species_id,
    ~ 0 + T_z + P_z + TP_z | species_id
  )
  #
  # The first block is the among-species intercept variance. The second is the
  # among-species W_z-slope variance. In the third, `struct = "GEN"` estimates
  # the full 3 x 3 covariance of the T_z, P_z, and TP_z species coefficients.
  resident_climate_fit <- metafor::rma.mv(
    yi = auc_effect,
    V = resident_sampling_variance_matrix,
    mods = resident_climate_formula,
    random = resident_species_random_effects,
    struct = "GEN",
    data = resident_auc_data,
    method = "REML",
    sparse = TRUE
  )
  #
  # The same formula is also fitted with ML when it is the climate-only baseline
  # for a predictor comparison. Fixed-effect uncertainty for the final REML fit
  # is additionally calculated with experiment-clustered CR1:
  resident_climate_CR1 <- metafor::robust(
    resident_climate_fit,
    cluster = resident_climate_fit$data$experiment_id,
    adjust = TRUE,
    clubSandwich = FALSE
  )
  #
  # Pipeline code:
  #   config/workflow.R              resident_analysis_spec()
  #   R/models/meta_regression.R     fit_meta_regression()
  #   R/models/meta_regression.R     cluster_robust_inference()
  #   R/models/meta_regression.R     random_effect_diagnostics()
  #   R/targets/definitions/models.R resident climate model definitions


  # 4. SPECIES CONTEXT CLASSIFICATION
  # ==================================
  # Question: within climates represented by the data for a species, is its
  # fitted warming response always positive, always negative, or climate-dependent?
  #
  # A species needs observations from at least four experiments and a full-rank
  # [1, T_z, P_z, T_z*P_z] design. Predictions are restricted to the convex hull
  # of the climates where that species was observed, with warming fixed at the
  # median across experiments.
  species_support <- species_context_support(resident_auc_data)
  eligible_species <- species_support |>
    dplyr::filter(eligible) |>
    dplyr::pull(species_id)

  conditional_species_coefficients <-
    extract_conditional_species_coefficients(resident_climate_fit)
  median_warming_z <- resident_auc_data |>
    dplyr::distinct(experiment_id, W_z) |>
    dplyr::pull(W_z) |>
    stats::median()

  one_species_grid <- species_observed_hull_grid(
    resident_auc_data,
    species_name = eligible_species[[1]],
    W_z_value = median_warming_z,
    points = 81L
  )
  one_species_prediction <- predict_species_on_observed_support(
    resident_climate_fit,
    resident_auc_data,
    resident_climate_formula,
    conditional_species_coefficients,
    species_name = eligible_species[[1]],
    W_z_value = median_warming_z,
    points = 81L
  )
  #
  # Predictions combine the fixed climate coefficients with that species'
  # conditional intercept, W, T, P, and TP coefficients. Prediction uncertainty
  # includes fixed-effect uncertainty, random-coefficient prediction error, and
  # their covariance. Classification uses fitted extrema, not p-values or CIs:
  species_class <- dplyr::case_when(
    minimum_fitted_delta_auc > 0 ~ "Consistent winner",
    maximum_fitted_delta_auc < 0 ~ "Consistent loser",
    TRUE ~ "Context-dependent"
  )
  species_context <- build_species_context_bundle(
    resident_climate_fit,
    resident_auc_data,
    resident_climate_formula,
    points = 81L
  )
  #
  # Pipeline code:
  #   R/models/species_context.R  species_context_support()
  #   R/models/species_context.R  species_observed_hull_grid()
  #   R/models/species_context.R  predict_species_on_observed_support()
  #   R/models/species_context.R  classify_supported_extrema()
  #   R/models/species_context.R  build_species_context_bundle()


  # 5. FAMILY-LEVEL RECONSTRUCTION AND MODEL
  # =========================================
  # Question: do plant families show consistent overall and climate-dependent
  # warming responses?
  #
  # Family cover is reconstructed from plot-level cover for the retained
  # resident taxa; family effects are not averages of species effects. Within
  # each experiment, `auc ~ ODT * family_id` provides the warmed-minus-origin
  # estimate and SE for every family. That model is fitted once because it is
  # needed to estimate the contrasts. Its off-diagonal `vcov()` is not used.
  # The reviewed taxonomy includes Viola palustris in Violaceae.
  family_effects <- build_family_auc_effects(
    family_plot_auc,
    family_membership,
    resident_auc_data
  )
  family_sampling_variance_matrix <- build_sampling_variance_matrix(
    family_effects$data,
    se_column = "family_auc_se"
  )

  family_random_effects <- list(
    ~ 1 | family_id,
    ~ 0 + W_z | family_id,
    ~ 0 + T_z + P_z + TP_z | family_id
  )
  family_fit <- metafor::rma.mv(
    yi = family_auc_effect,
    V = family_sampling_variance_matrix,
    mods = ~ T_z * P_z * W_z,
    random = family_random_effects,
    struct = "GEN",
    data = family_effects$data,
    method = "REML",
    sparse = TRUE
  )
  #
  # The pipeline first tries nlminb and may retry with bobyqa (`maxfun = 100000`)
  # without changing the model. Family support and classification then use the
  # same experiment-count, climate-hull, median-warming, and fitted-extrema rules
  # as the species analysis.
  #
  # Pipeline code:
  #   R/models/family.R  family cover, AUC effects, model, and classification
  #   R/models/family.R  build_family_auc_effects()
  #   R/models/family.R  fit_family_aggregated_model()


  # 6. LDMC, LMA, AND VEGETATIVE-HEIGHT MODELS
  # ============================================
  # Question: does LDMC, LMA, or vegetative height improve the explanation of
  # warming responses and their climate dependence beyond climate alone?
  #
  # Each trait is analysed separately. Only finite, positive species medians are
  # retained; missing values are not imputed. Medians are log transformed and
  # standardized across unique species. Each trait therefore has its own
  # complete-case effects, matching diagonal sampling variances, and matching
  # climate-only model.
  trait_climate_only_formula <- ~ T_z * P_z * W_z
  trait_full_formula <- ~ Trait_z * T_z * P_z * W_z
  #
  # The full formula includes all main effects, all two- and three-way
  # interactions, and Trait:T:P:W. The species random effects are unchanged and
  # no random trait slope is added. There is no intermediate additive model.
  # Climate-only ML is compared directly with full ML by LRT, AIC, and AICc;
  # climate-only and full models are also fitted with REML.
  LDMC_models <- fit_predictor_model_pair(
    LDMC_specification,
    climate_formula = trait_climate_only_formula
  )
  LMA_models <- fit_predictor_model_pair(
    LMA_specification,
    climate_formula = trait_climate_only_formula
  )
  height_models <- fit_predictor_model_pair(
    height_specification,
    climate_formula = trait_climate_only_formula
  )
  #
  # Pipeline code:
  #   R/models/predictors.R      prepare_trait_data()
  #   R/models/predictors.R      trait_predictor_specification()
  #   R/models/predictors.R      fit_predictor_model_pair()
  #   R/models/meta_regression.R fit_meta_regression()


  # 7. GRAMINOID AND WOODINESS MODELS
  # ==================================
  # Question: do graminoids differ from non-graminoids, or woody species from
  # non-woody species, in warming response and climate dependence?
  #
  # All resident effects are retained because both classifications are complete
  # and unique by species. Non-graminoid and non-woody are the reference levels.
  climate_only_formula <- ~ T_z * P_z * W_z
  graminoid_full_formula <- ~ GraminoidStatus * T_z * P_z * W_z
  woodiness_full_formula <- ~ Woodiness * T_z * P_z * W_z
  #
  # The full formulas include the group, T, P, and W main effects and all their
  # interactions. The species random effects remain unchanged. Climate-only ML
  # is compared directly with full ML; both roles also have REML fits.
  graminoid_models <- fit_predictor_model_pair(
    graminoid_specification,
    climate_formula = climate_only_formula,
    full_climate_ML = resident_climate_ML,
    full_climate_REML = resident_climate_fit
  )
  woodiness_models <- fit_predictor_model_pair(
    woodiness_specification,
    climate_formula = climate_only_formula,
    full_climate_ML = resident_climate_ML,
    full_climate_REML = resident_climate_fit
  )
  #
  # Pipeline code:
  #   R/models/predictors.R  prepare_functional_predictors()
  #   R/models/predictors.R  functional_predictor_specification()
  #   R/models/predictors.R  fit_predictor_model_pair()


  # 8. PHYLOGENETIC PEV MODEL
  # ==========================
  # Question: does phylogenetic position improve the explanation of mean or
  # climate-dependent warming responses beyond climate alone?
  #
  # The first three phylogenetic eigenvectors are standardized across unique
  # retained species before expansion to effect rows. Axis directions are fixed:
  # PEV1 has Poa trivialis positive, PEV2 has Erigeron formosissimus negative,
  # and PEV3 has Euphrasia minima positive.
  phylogeny_full_formula <-
    ~ T_z * P_z * W_z +
      PEV1_z + PEV2_z + PEV3_z +
      PEV1_z:W_z + PEV2_z:W_z + PEV3_z:W_z +
      PEV1_z:T_z + PEV2_z:T_z + PEV3_z:T_z +
      PEV1_z:P_z + PEV2_z:P_z + PEV3_z:P_z +
      PEV1_z:TP_z + PEV2_z:TP_z + PEV3_z:TP_z
  #
  # Each PEV is crossed with W, T, P, and T:P. The formula does not add further
  # PEV-by-W-by-T or PEV-by-W-by-P terms. The random effects are unchanged.
  # Climate-only ML is compared directly with full ML; both are also fitted by
  # REML. There is no intermediate additive PEV model.
  phylogeny_models <- fit_predictor_model_pair(
    phylogeny_specification,
    climate_formula = climate_only_formula,
    full_climate_ML = resident_climate_ML,
    full_climate_REML = resident_climate_fit
  )
  #
  # Pipeline code:
  #   R/models/phylogeny.R   phylogeny construction, validation, and PEVs
  #   R/models/predictors.R  standardize_unique_species_predictors()
  #   R/models/predictors.R  build_phylogeny_predictor_specification()
  #   R/models/predictors.R  fit_predictor_model_pair()


  # 9. HETEROGENEITY RESOLVED BY EACH PREDICTOR
  # =============================================
  # Question: how much smaller is the fitted among-species variation in the full
  # predictor model than in its climate-only model?
  #
  # The two REML fits use exactly the same effects and sampling variances.
  # Climate averaging uses distinct experiments rather than species-effect rows.
  baseline_fit <- one_predictor_pair$fits$climate_REML$fit
  full_fit <- one_predictor_pair$fits$full_REML$fit
  experiment_climates <- one_predictor_pair$data |>
    dplyr::distinct(experiment_id, T_z, P_z, W_z) |>
    dplyr::mutate(TP_z = T_z * P_z)

  baseline_H <- response_heterogeneity_metrics(
    baseline_fit,
    one_predictor_pair$data
  )
  full_H <- response_heterogeneity_metrics(
    full_fit,
    one_predictor_pair$data
  )
  percent_resolved <- heterogeneity_resolution(baseline_H, full_H)
  #
  # H_overall is the species intercept variance. H_W averages the W-slope
  # variance over experiment warming values. H_context averages x'Gx for the
  # correlated T, P, and T:P species coefficients. H_total is their sum.
  # A negative percentage means that a fitted heterogeneity component increased
  # in the full model; it does not contradict a significant fixed-effect LRT.
  direct_ML_comparison <- one_predictor_pair$comparisons
  delta_AIC_full_minus_climate <-
    one_predictor_pair$statistics$AIC[
      one_predictor_pair$statistics$model == "full predictor-climate"
    ] -
    one_predictor_pair$statistics$AIC[
      one_predictor_pair$statistics$model == "climate only"
    ]
  #
  # Pipeline code:
  #   R/evaluation/heterogeneity.R  response_heterogeneity_metrics()
  #   R/evaluation/heterogeneity.R  heterogeneity_resolution()
  #   R/evaluation/heterogeneity.R  heterogeneity_from_predictor_models()
  #   R/evaluation/heterogeneity.R  assemble_heterogeneity_synthesis()


  # 10. CROSS-VALIDATION AND MODEL DIAGNOSTICS
  # ===========================================
  # Question: does a predictor improve prediction for a complete geographic
  # gradient that was not used to fit the model?
  #
  # One gradient is held out at a time. The remaining effects and their diagonal
  # sampling variances are used to refit the climate-only and full REML models.
  # Trait models use their own complete-case effects; complete-data predictors
  # use all resident effects. Training and test effect IDs must be disjoint.
  one_partition <- model_performance_partition(
    one_predictor_specification$data,
    one_predictor_specification$sampling_variance_matrix,
    heldout_gradient = one_gradient
  )
  climate_fold_fit <- fit_model_performance_fold(
    one_partition,
    climate_only_formula,
    label = "held-out-gradient climate-only REML"
  )
  full_fold_fit <- fit_model_performance_fold(
    one_partition,
    one_predictor_specification$full_formula,
    label = "held-out-gradient full predictor REML",
    predictor_columns = one_predictor_specification$predictor_columns
  )
  climate_predictions <- predict_model_performance_fold(
    climate_fold_fit,
    one_partition,
    climate_only_formula,
    label = "climate-only"
  )
  full_predictions <- predict_model_performance_fold(
    full_fold_fit,
    one_partition,
    one_predictor_specification$full_formula,
    label = "full predictor"
  )
  #
  # General predictions use fixed effects and include new-species random
  # variation. For species present in training, a second prediction conditions
  # on their fitted species effects. Predictions are also marked according to
  # whether the held-out climate is inside the training climate support.
  prediction_error <- predicted_delta_auc - observed_delta_auc
  RMSE <- sqrt(mean(prediction_error^2))
  MAE <- mean(abs(prediction_error))
  bias <- mean(prediction_error)
  predictive_R2 <- 1 - mean(prediction_error^2) /
    stats::var(observed_delta_auc)
  coverage_95 <- mean(
    observed_delta_auc >= prediction_interval_low &
      observed_delta_auc <= prediction_interval_high
  )
  gaussian_NLPD <- model_performance_gaussian_score(
    observed_delta_auc,
    predicted_delta_auc,
    predictive_covariance
  )[["negative_log_predictive_density"]]
  #
  # The pipeline aggregates these quantities once with equal weight per gradient
  # and once by pooling effects. Predictor skill is calculated against the
  # matching climate-only fold, for example:
  RMSE_percent_improvement <- 100 * (1 - full_RMSE / climate_RMSE)
  incremental_predictive_R2 <- 1 - full_MSE / climate_MSE
  log_score_improvement <- climate_NLPD - full_NLPD
  #
  # In-sample marginal R2 is fixed-effect variance divided by fixed + species
  # random + mean sampling variance. Conditional R2 has fixed + species random
  # variance in the numerator. Fit checks cover design rank, convergence, finite
  # estimates, random effects, prediction covariance, warnings, and profiles.
  # Cross-validation models are evaluation-only and never replace final models.
  #
  # Pipeline code:
  #   R/evaluation/model_performance.R  model_performance_partition()
  #   R/evaluation/model_performance.R  fit_model_performance_fold()
  #   R/evaluation/model_performance.R  predict_model_performance_fold()
  #   R/evaluation/model_performance.R  model_performance_gaussian_score()
  #   R/evaluation/model_performance.R  aggregate_model_performance_statistics()
  #   R/evaluation/model_performance.R  summarise_model_predictive_skill()
  #   R/evaluation/model_performance.R  run_model_performance_gradient_fold()


  # NEXT STEPS AFTER THIS REVIEW
  # ============================
  # The corrected sampling variances and climate-only-versus-full comparisons
  # require a clean certification rebuild of models, classifications,
  # heterogeneity, cross-validation, reports, figures, release checks, and the
  # temporary historical comparisons. The archive remains in place until that
  # rebuild passes and its removal receives explicit approval.
}
