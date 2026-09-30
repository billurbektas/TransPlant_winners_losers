# MODEL ANALYSIS: READ-ONLY SUMMARY
# =================================
#
# This file brings the model analysis together in normal script order so that
# the data choices, formulas, random effects, comparisons, and classifications
# can be reviewed in one place. The workflow itself is run with targets, which
# is useful for reproducibility but harder to read from beginning to end. Here,
# I provide a summary for review!

if (FALSE) {

  # TERMS USED BELOW -- I will improve these for the final version of the code!
  # --------------------------------------------------------------------------
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
  # V_auc
  #   Sampling covariance matrix for Delta AUC estimates.
  #
  # GEN
  #   An unstructured covariance: variances and pairwise covariances are freely
  #   estimated for the T_z, P_z, and TP_z random coefficients.
  #
  # CR1
  #   Small-sample-adjusted, cluster-robust fixed-effect uncertainty, clustered
  #   by experiment. It changes standard errors, confidence intervals, tests,
  #   and denominator degrees of freedom, but not coefficient estimates.
  #
  # ML / REML
  #   Maximum likelihood / restricted maximum likelihood. ML is used to compare
  #   nested fixed-effect formulas. REML is used for final estimates of model
  #   coefficients, variance, predictions, and heterogeneity.
  #
  # M0 / M1 / M2 / M3
  #   These are current internal labels that will be replaced in the next naming
  #   pass. M0 is climate only WITH species random slopes and intercepts. 
  #   M1 adds a trait, functional group, or phylogenetic axis WITHOUT climate interactions.
  #   M2 contains the full trait- or functional-group-by-climate interactions.
  #   M3 contains the specified phylogenetic interactions.
  #
  # PEV1_z, PEV2_z, PEV3_z
  #   Standardized coordinates on the first three phylogenetic eigenvectors.
  #   They describe phylogenetic position; they are not measured plant traits.


  # 1. PREPARING THE RESIDENT-SPECIES RESPONSE DATA
  # ===========================================================================
  # We use maximum-duration, raw relative-cover AUC contrasts for resident
  # species. One row represents one accepted species in one experiment.
  #
  # The current validated data contain 1,604 effects from 526 species and 42
  # experiments. 
  #
  # Six taxa remain in the response calculations but are excluded when the
  # resident models and phylogeny are prepared:
  # Taken out because they are outliers in our functional grouping downstream,
  # also being deep lineages they were dominating PEVs too much. 
  # TODO: need to ommit these from the start!   
  resident_model_exclusions <- c(
    "Gymnocarpium dryopteris",
    "Athyrium spinulosum",
    "Botrychium lunaria",
    "Equisetum arvense",
    "Selaginella selaginoides",
    "Juniperus communis" # TODO check other gymnosperms
  )
  #
  # The climate variables joined to each effect are destination summer
  # temperature (T), log destination summer precipitation (P), and the summer
  # warming difference (W). They are standardized using the distinct experiments
  # before being copied to species rows.
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


  # 2. RECONSTRUCTING SAMPLING COVARIANCE
  # ===========================================================================
  # Delta AUC estimates from different species in the same experiment are not
  # statistically independent. To recover their shared sampling uncertainty, a
  # plot-level model is fitted separately within each experiment:
  covariance_reconstruction_formula <- auc ~ ODT * species
  #
  # ODT (experimental groups) is ordered as destination control, origin control, and warmed. Estimated
  # marginal means are calculated for each species and the contrast c(0, -1, 1)
  # gives warmed minus origin. Its `vcov()` matrix contains both the sampling
  # variances and the covariances among species contrasts in that experiment.
  plot_auc_model <- stats::lm(auc ~ ODT * species, data = one_experiment)
  species_means <- emmeans::emmeans(plot_auc_model, ~ ODT | species)
  warmed_minus_origin <- emmeans::contrast(
    species_means,
    method = list(warmed_minus_origin = c(0, -1, 1)),
    adjust = "none"
  )
  one_experiment_covariance <- stats::vcov(warmed_minus_origin)

  #
  # Pipeline code:
  #   R/responses/auc_model_helpers.R  reconstruct_auc_covariance()
  #   R/responses/auc_model_helpers.R  validate_auc_covariance()
  #   R/responses/auc_model_helpers.R  select_auc_covariance()
  #   R/responses/auc_model_helpers.R  subset_auc_v()


  # 3. RESIDENT CLIMATE MODEL
  # ===========================================================================
  # Question: how does the mean warming response in AUC vary with destination
  # temperature, destination precipitation, and warming magnitude, while
  # allowing species to differ in their overall (intercept) 
  # and climate-dependent responses (slopes)?
  #
  # The response is Delta AUC and the fixed-effects formula is:
  resident_climate_formula <- ~ T_z * P_z * W_z
  #
  # Species can differ in their overall response, their response to warming
  # magnitude, and their response across destination temperature and
  # precipitation (including their two-way interaction TP_z):
  resident_species_random_effects <- list(
    ~ 1 | species_id,
    ~ 0 + W_z | species_id,
    ~ 0 + T_z + P_z + TP_z | species_id
  )
  #
  # The intercept and W terms are separate variance components. The T, P, and TP
  # coefficients use `struct = "GEN"`, allowing their variances and covariances
  # to be estimated. Covariance between these three random-effect blocks is zero.
  #
  # The final model is fitted with REML and the full sampling covariance V_auc.
  # An otherwise identical ML fit is retained for later nested model comparisons.
  resident_climate_fit <- metafor::rma.mv(
    yi = auc_effect,
    V = V_auc,
    mods = ~ T_z * P_z * W_z,
    random = resident_species_random_effects,
    struct = "GEN",
    data = resident_auc_data,
    method = "REML",
    sparse = TRUE
  )
  #
  # Fixed-effect uncertainty is shown both from the fitted model and with a CR1
  # adjustment clustered by experiment:
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
  # ===========================================================================
  # Question: within climates represented by the data for a species, is its
  # fitted warming response always positive, always negative, or dependent on
  # climate context?
  #
  # A species is classified only when it occurs in at least four experiments and
  # its [1, T_z, P_z, T_z*P_z] design has rank four. Predictions are made over a
  # fine T/P grid inside that species' observed destination-climate hull, with
  # W_z fixed at the median warming among experiments. Observed climate points
  # are included in addition to the grid.
  #
  # Predictions combine the fixed climate coefficients with the species'
  # conditional intercept, W, T, P, and TP coefficients. Their SEs include fixed
  # uncertainty, random-coefficient prediction error, and fixed/random covariance.
  #
  # The classification uses fitted extrema, not p-values or confidence intervals:
  species_class <- dplyr::case_when(
    minimum_fitted_delta_auc > 0 ~ "Consistent winner",
    maximum_fitted_delta_auc < 0 ~ "Consistent loser",
    TRUE ~ "Context-dependent"
  )
  #
  # Species without enough support remain in the support table but are not forced
  # into a class. The outputs contain the support assessment, conditional
  # coefficients, prediction grid and uncertainty, fitted extrema and their
  # climate locations, and the final class.
  #
  # Pipeline code:
  #   R/models/species_context.R  species_context_support()
  #   R/models/species_context.R  species_observed_hull_grid()
  #   R/models/species_context.R  predict_species_on_observed_support()
  #   R/models/species_context.R  build_species_context_bundle()
  #   R/models/species_context.R  classify_supported_extrema()


  # 5. FAMILY-LEVEL RECONSTRUCTION AND MODEL
  # ===========================================================================
  # Question: do plant families show consistent overall and climate-dependent
  # warming responses?
  #
  # Family cover is reconstructed directly from raw plot-level cover for the
  # retained resident taxa. Family effects are therefore not averages of species
  # effect sizes. Within each experiment, `auc ~ ODT * family_id` and the same
  # c(0, -1, 1) contrast provide family Delta AUC effects and their sampling
  # covariance. The reviewed taxonomy includes Viola palustris in Violaceae.
  #
  # The fixed formula is the same climate formula as for species:
  family_climate_formula <- ~ T_z * P_z * W_z
  #
  # Families have the corresponding random-effects structure:
  family_random_effects <- list(
    ~ 1 | family_id,
    ~ 0 + W_z | family_id,
    ~ 0 + T_z + P_z + TP_z | family_id
  )
  #
  # The family model is fitted with REML and the reconstructed family covariance.
  # It first uses nlminb and can retry with bobyqa (`maxfun = 100000`) without
  # changing the formula or random structure.
  family_fit <- metafor::rma.mv(
    yi = family_auc_effect,
    V = family_auc_covariance,
    mods = ~ T_z * P_z * W_z,
    random = family_random_effects,
    struct = "GEN",
    data = family_auc_data,
    method = "REML",
    sparse = TRUE
  )
  #
  # Family support and classification follow the species rules: at least four
  # experiments, rank four for [1, T, P, T*P], predictions inside the observed
  # T/P hull at median warming, and classes based on the signs of fitted extrema.
  #
  # Pipeline code:
  #   R/models/family.R  family cover and AUC reconstruction
  #   R/models/family.R  fit_family_aggregated_model()
  #   R/models/family.R  family_support_table() and classification functions


  # 6. LDMC, LMA, AND VEGETATIVE-HEIGHT MODELS
  # ===========================================================================
  # Question: does LDMC, LMA, or vegetative height explain either the mean
  # warming response or the way that response changes with climate?
  #
  # Each trait is analyzed separately. Only species with a finite, positive
  # species-level median for that trait are retained; missing traits are not
  # imputed. Trait medians are log transformed and standardized across unique
  # retained species before being copied to their effect rows. Each trait
  # therefore has its own complete-case data, matching V subset, and matching
  # climate-only model.
  #
  # Three fixed-effects formulas are compared:
  trait_climate_only_formula <- ~ T_z * P_z * W_z
  trait_additive_formula <- ~ Trait_z + T_z * P_z * W_z
  trait_full_formula <- ~ Trait_z * T_z * P_z * W_z
  #
  # The full formula contains the trait and climate main effects; all six two-way
  # interactions among Trait, T, P, and W; all four three-way interactions; and
  # the Trait:T:P:W four-way interaction. The species random effects remain
  # exactly as in the resident climate model; no random trait slope is added.
  #
  # The climate-only, additive-trait, and full-interaction formulas are fitted by
  # ML for two nested comparisons:
  #   climate only versus additive trait;
  #   additive trait versus trait-by-climate interactions.
  # The climate-only and full models are also fitted by REML for final estimates,
  # heterogeneity, predictions, and experiment-clustered CR1 uncertainty.
  LDMC_models <- fit_predictor_sequence(
    LDMC_complete_rows_and_V,
    climate_formula = ~ T_z * P_z * W_z
  )
  LMA_models <- fit_predictor_sequence(
    LMA_complete_rows_and_V,
    climate_formula = ~ T_z * P_z * W_z
  )
  height_models <- fit_predictor_sequence(
    height_complete_rows_and_V,
    climate_formula = ~ T_z * P_z * W_z
  )
  #
  # A comparison is valid only when both models use identical effect rows and
  # the identical covariance matrix. The outputs include trait scaling and
  # missingness, likelihood-ratio tests, AIC/AICc, the full REML and CR1 tables,
  # random-effect diagnostics, heterogeneity, and held-out performance.
  #
  # Pipeline code:
  #   R/models/predictors.R      prepare_trait_data()
  #   R/models/predictors.R      predictor_formulas()
  #   R/models/predictors.R      trait_predictor_specification()
  #   R/models/predictors.R      fit_predictor_sequence()
  #   R/models/meta_regression.R fit_meta_regression()


  # 7. GRAMINOID AND WOODINESS MODELS
  # ===========================================================================
  # Question: do graminoids differ from non-graminoids, or woody species from
  # non-woody species, in their mean warming response or climate dependence?
  #
  # All resident effects are retained because both classifications are complete
  # and unique by species. Non-graminoid and non-woody are the reference levels.
  # These binary variables are not standardized.
  graminoid_additive_formula <- ~ GraminoidStatus + T_z * P_z * W_z
  graminoid_full_formula <- ~ GraminoidStatus * T_z * P_z * W_z
  woodiness_additive_formula <- ~ Woodiness + T_z * P_z * W_z
  woodiness_full_formula <- ~ Woodiness * T_z * P_z * W_z
  #
  # For either group G, the full formula contains the G, T, P, and W main effects;
  # all their two- and three-way interactions; and G:T:P:W. The species random
  # effects are unchanged and there is no random slope for group membership.
  #
  # ML compares climate only with the additive group effect, then the additive
  # model with all group-by-climate interactions. Full REML fits provide final
  # estimates, heterogeneity, predictions, and experiment-clustered CR1
  # uncertainty. Both analyses use the complete-data climate baselines and V_auc.
  graminoid_models <- fit_predictor_sequence(
    graminoid_specification,
    climate_formula = ~ T_z * P_z * W_z,
    full_M0_ML = full_climate_ML,
    full_M0_REML = full_climate_REML
  )
  woodiness_models <- fit_predictor_sequence(
    woodiness_specification,
    climate_formula = ~ T_z * P_z * W_z,
    full_M0_ML = full_climate_ML,
    full_M0_REML = full_climate_REML
  )
  #
  # Pipeline code:
  #   R/models/predictors.R  prepare_functional_predictors()
  #   R/models/predictors.R  functional_predictor_specification()
  #   R/models/predictors.R  predictor_formulas()
  #   R/models/predictors.R  fit_predictor_sequence()


  # 8. PHYLOGENETIC PEV MODEL
  # ===========================================================================
  # Question: does phylogenetic position explain the mean or climate-dependent
  # warming response after accounting for destination climate and warming?
  #
  # The resident phylogeny provides the first three phylogenetic eigenvectors.
  # Each axis is standardized across unique retained species before expansion to
  # effect rows. Axis directions are fixed for stable interpretation: PEV1 has
  # Poa trivialis positive, PEV2 has Erigeron formosissimus negative, and PEV3
  # has Euphrasia minima positive.
  phylogeny_additive_formula <-
    ~ T_z * P_z * W_z + PEV1_z + PEV2_z + PEV3_z
  phylogeny_full_formula <-
    ~ T_z * P_z * W_z +
      PEV1_z + PEV2_z + PEV3_z +
      PEV1_z:W_z + PEV2_z:W_z + PEV3_z:W_z +
      PEV1_z:T_z + PEV2_z:T_z + PEV3_z:T_z +
      PEV1_z:P_z + PEV2_z:P_z + PEV3_z:P_z +
      PEV1_z:TP_z + PEV2_z:TP_z + PEV3_z:TP_z
  #
  # The full formula contains the climate terms, three PEV main effects, and each
  # PEV crossed separately with W, T, P, and TP. It does not add PEV-by-W-by-T or
  # PEV-by-W-by-P terms beyond that list. The species random effects are unchanged.
  #
  # ML compares climate only with additive PEV effects, then additive PEV effects
  # with the specified PEV-by-context terms. The full formula is fitted by REML
  # for final estimates, heterogeneity, predictions, and experiment-clustered
  # CR1 uncertainty. It uses the complete resident cohort and full V_auc.
  phylogeny_models <- fit_predictor_sequence(
    phylogeny_specification,
    climate_formula = ~ T_z * P_z * W_z,
    full_M0_ML = full_climate_ML,
    full_M0_REML = full_climate_REML
  )
  #
  # Pipeline code:
  #   R/models/phylogeny.R   phylogeny construction, validation, and PEVs
  #   R/models/predictors.R  standardize_unique_species_predictors()
  #   R/models/predictors.R  build_phylogeny_predictor_specification()
  #   R/models/predictors.R  predictor_formulas()
  #   R/models/predictors.R  fit_predictor_sequence()


  # 9. HOW MUCH DOES EACH PREDICTOR REDUCE SPECIES HETEROGENEITY?
  # ===========================================================================
  # Question: how much smaller is the fitted among-species variation after adding
  # a trait, functional group, or phylogenetic position to the climate model?
  #
  # Each full predictor REML model is compared with its climate-only REML model
  # fitted to exactly the same rows and covariance matrix. The family model is
  # kept separate because this comparison concerns information about species.
  # Climate averaging uses distinct experiments rather than effect rows.
  #
  # The fitted species variation is divided into:
  #   H_overall: species random-intercept variance;
  #   H_W:       average variation associated with species warming slopes;
  #   H_context: average variation associated with the correlated T, P, and TP
  #              species coefficients;
  #   H_total:   the sum of those three components.
  H_overall <- species_intercept_variance
  H_W <- mean(experiment_climates$W_z^2 * species_W_slope_variance)
  H_context <- mean(rowSums(
    (X_context %*% G_context) * X_context
  ))
  H_total <- H_overall + H_W + H_context
  percent_resolved <- 100 * (baseline_H - predictor_H) / baseline_H
  #
  # Negative percent resolved is possible and means that the fitted variation
  # increased after adding the predictor. The summary also retains the ML nested
  # tests and Delta AIC/AICc. CR1 is not used to recalculate random-effect
  # variation because it concerns fixed-effect uncertainty.
  #
  # Pipeline code:
  #   R/evaluation/heterogeneity.R  response_heterogeneity_metrics()
  #   R/evaluation/heterogeneity.R  heterogeneity_resolution()
  #   R/evaluation/heterogeneity.R  heterogeneity_from_predictor_sequence()
  #   R/evaluation/heterogeneity.R  assemble_heterogeneity_synthesis()


  # 10. CROSS-VALIDATION AND MODEL DIAGNOSTICS
  # ===========================================================================
  # Question: does a predictor improve prediction for a complete geographic
  # gradient that was not used to fit the model?
  #
  # One gradient is held out at a time. Trait models use their matched
  # trait-complete rows and climate baselines; complete-data predictors use the
  # full resident climate baseline. Every fold refits the same full formula and
  # species random structure with REML using the training covariance submatrix.
  # Train/test sampling covariance must be effectively zero when whole gradients
  # are separated.
  one_gradient_fold <- run_model_performance_gradient_fold(
    heldout_gradient = one_gradient,
    predictor_specifications = all_predictor_specifications,
    full_data = resident_auc_data,
    full_V = V_auc,
    climate_formula = ~ T_z * P_z * W_z
  )
  #
  # Predictions are summarized both for a new species and, when a species was
  # seen in training, conditional on its fitted species effects. Results are also
  # separated according to whether the held-out climate lies inside the training
  # climate support.
  #
  # Full predictor models are compared with their matching climate baselines
  # using predictive R2, RMSE, MAE, bias, Gaussian negative log predictive
  # density, and 95% interval coverage. Results are summarized with equal weight
  # per gradient and, separately, equal weight per effect.
  #
  # The in-sample summaries use:
  #   marginal R2 = fixed-effect variance / total response variance;
  #   conditional R2 = (fixed + species-random variance) / total response variance.
  # Total response variance also includes mean sampling variance.
  #
  # Fit checks cover design rank, optimizer status, finite coefficients and
  # variance components, the expected random-effect structure, predictive
  # covariance, warnings, and variance profiles. Cross-validation fits are used
  # only for evaluation and never replace the final models.
  #
  # Pipeline code:
  #   R/evaluation/model_performance.R  model_performance_fit_r2()
  #   R/evaluation/model_performance.R  model_performance_partition()
  #   R/evaluation/model_performance.R  fit_model_performance_fold()
  #   R/evaluation/model_performance.R  run_model_performance_gradient_fold()
  #   R/evaluation/model_performance.R  metric and comparison functions
  #   R/models/meta_regression.R        random_effect_diagnostics()
  #   R/models/meta_regression.R        cluster_robust_inference()


  # AFTER THIS SUMMARY IS REVIEWED
  # ===========================================================================
  # The next pass will centralize all formulas and model settings in
  # `analysis_model_specifications()`, make the pipeline use that registry, and
  # simplify cryptic wording and old labels throughout the workflow. Temporary
  # historical certification code will be separated from the current models,
  # followed by a clean certification of tests, models, reports, and figures.
  #
  # The archive and old validation machinery will remain in place until the
  # revised workflow passes certification and receives explicit approval for
  # their removal.
}
