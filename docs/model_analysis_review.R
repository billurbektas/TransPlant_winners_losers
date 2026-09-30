# MODEL ANALYSIS: READ-ONLY SUMMARY
# =================================
#
# This file brings the model analysis together in normal script order so that
# the data choices, formulas, random effects, comparisons, and classifications
# can be reviewed in one place. The workflow itself is run with targets, which
# is useful for reproducibility but harder to read from beginning to end.

if (FALSE) {

  # TERMS USED BELOW
  # ----------------
  # Delta AUC / auc_effect
  #   Difference in raw relative-cover area under the time curve: warmed
  #   treatment minus the origin control. Plot-level AUC is calculated by the
  #   trapezoid rule over all available relative years for that experiment. It
  #   is not divided by study duration; "maximum duration" means that each
  #   experiment contributes its complete available time series.
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
  #   Destination minus origin June-August temperature, representing the
  #   magnitude of experimental warming / its standardized value across
  #   distinct experiments. W_z = 0 is the mean warming magnitude in these
  #   experiments, not an unwarmed treatment.
  #
  # TP_z
  #   T_z * P_z. In the random effects this represents differences among species
  #   or families in the temperature-by-precipitation response.
  #
  # sampling_variance_matrix
  #   A diagonal matrix containing auc_se^2 in effect-row order. The models
  #   assume independent sampling errors between Delta AUC estimates.
  #
  # ODT
  #   The three plot groups used in the plot-level AUC model: destination
  #   controls, origin controls, and warmed transplants. The response analysed
  #   below is the warmed-minus-origin contrast.
  #
  # GEN
  #   An unstructured covariance for the species or family T_z, P_z, and TP_z
  #   random coefficients: their three variances and three pairwise covariances
  #   are estimated. The intercept and W_z random effects are separate blocks,
  #   so their covariance with the T/P/TP block is fixed at zero.
  #
  # CR1
  #   Fixed-effect standard errors and tests adjusted for effects coming from
  #   the same experiment. It does not change the fitted coefficients or random
  #   effects.
  #
  # ML / REML
  #   Maximum likelihood / restricted maximum likelihood. ML is used to compare
  #   the climate-only and full fixed-effect formulas. REML is used for final
  #   coefficients, predictions, and heterogeneity. This rule is used for all
  #   six predictor analyses and is not repeated below. Fixed-effect uncertainty
  #   is reported with experiment-clustered CR1.
  #
  # LRT / AIC / AICc
  #   Likelihood-ratio test / Akaike information criterion / small-sample-
  #   corrected AIC. The LRT jointly tests all predictor terms added to the
  #   climate-only model. Lower AIC or AICc indicates the better-supported model.
  #
  # PEV1_z, PEV2_z, PEV3_z
  #   Standardized coordinates on the first three positive axes from a principal
  #   coordinates analysis of patristic distances in the species tree. They
  #   describe phylogenetic position; they are not measured plant traits.


  # 1. PREPARING THE RESIDENT-SPECIES RESPONSE DATA
  # =================================================
  # We use maximum-duration, raw relative-cover AUC contrasts for resident
  # species. Here, "resident" means that the species occurred in the origin
  # control of that experiment after the observation-count filter. One row
  # represents one accepted species in one experiment. The final dataset
  # contains 1,604 effects from 526 species and 42 experiments.
  #
  # These six deep-lineage taxa are currently retained while relative cover and
  # Delta AUC are calculated, then excluded from the resident models. They are
  # incompatible with the focal functional comparison and strongly influenced
  # the early phylogenetic axes. Juniperus communis is the only gymnosperm in the
  # analysed vegetation that's why removed at the end.
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
  # not done here. It requires a full rebuild.
  #
  # The response table is restricted to raw, warmed-minus-origin, resident,
  # maximum-duration contrasts. Reviewed accepted names and experiment metadata
  # are then joined, the six taxa above are removed, and completeness, positive
  # SEs, and one unique effect ID per species x experiment are required.
  resident_auc_unscaled <- build_resident_auc_data_from_objects(
    winner_loser_results,
    analysis_data,
    climate_context_experiments,
    horizons = "max"
  )

  # Climate variables are destination summer temperature (T), log destination
  # summer precipitation (P), and destination-minus-origin summer temperature
  # (W). P must be positive before it is logged. T, P, and W are standardized
  # using the mean and sample SD across the 42 distinct experiments, not across
  # the 1,604 species-effect rows.
  resident_scaled <- scale_resident_auc_data(resident_auc_unscaled)
  resident_auc_data <- resident_scaled$data |>
    dplyr::mutate(TP_z = T_z * P_z)
  #
  # Pipeline code:
  #   config/analysis.R                  resident_model_exclusions()
  #   config/workflow.R                  resident_analysis_spec()
  #   R/responses/response_interfaces.R  build_resident_auc_data_from_objects()
  #   R/responses/auc_model_helpers.R    scale_resident_auc_data()


  # 2. SAMPLING VARIANCES
  # ======================
  # The winner/loser analysis has already estimated each warmed-minus-origin
  # Delta AUC and its standard error. Within each experiment, plot AUC is fitted
  # as `auc ~ ODT * species`; `emmeans` then calculates, for every species,
  #
  #   auc_effect = mean AUC_warmed - mean AUC_origin
  #
  # and the corresponding OLS contrast SE. The meta-regressions use that
  # information directly: sampling variance is auc_se^2, in exactly the same
  # order as the effect rows, with zero off-diagonals.
  resident_sampling_variance_matrix <- build_sampling_variance_matrix(
    resident_auc_data,
    se_column = "auc_se",
    id_column = "effect_id"
  )
  #
  # This assumes independent sampling errors between species effects, including
  # effects from the same experiment. The separate CR1 calculation accounts for
  # experiment grouping when reporting fixed-effect uncertainty; it does not add
  # off-diagonal values to this matrix. No off-diagonal sampling
  # covariance is reconstructed.
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

  # All three predictors are standardized, and T:P is stored as TP_z for the
  # corresponding species random coefficient.

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
  # Fixed-effect uncertainty is calculated separately by experiment:
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
  # Question: at the median warming magnitude, and within destination climates
  # represented by the data for a species, is its fitted Delta AUC always
  # positive, always negative, or dependent on destination climate?
  #
  # First identify species for which climate dependence can actually be
  # estimated. A species needs observations from at least four experiments and
  # those experiments must provide four independent columns: 1, T, P, and T:P.
  # Thus, four experiments are necessary but are not sufficient if their climate
  # combinations are redundant.
  species_support <- resident_auc_data |>
    dplyr::distinct(species_id, experiment_id, T_z, P_z) |>
    dplyr::group_by(species_id) |>
    dplyr::group_modify(function(species_climates, species_name) {
      climate_design <- cbind(
        intercept = 1,
        T_z = species_climates$T_z,
        P_z = species_climates$P_z,
        TP_z = species_climates$T_z * species_climates$P_z
      )
      tibble::tibble(
        n_experiments = dplyr::n_distinct(species_climates$experiment_id),
        context_design_rank = qr(climate_design)$rank,
        eligible = n_experiments >= 4L & context_design_rank == 4L
      )
    }) |>
    dplyr::ungroup()
  eligible_species <- species_support |>
    dplyr::filter(eligible) |>
    dplyr::pull(species_id)

  # The model estimates five deviations for each species: an overall deviation
  # and deviations associated with W, T, P, and T:P. In plain language, these
  # describe how that species' fitted response differs from the average response
  # across all species. Species represented by less information are pulled more
  # strongly towards the overall response rather than being estimated entirely
  # independently.
  conditional_species_coefficients <-
    extract_conditional_species_coefficients(resident_climate_fit)

  # Classification is made at one representative warming magnitude: the median
  # standardized temperature difference between warmed and origin plots across
  # the distinct experiments.
  median_warming_z <- resident_auc_data |>
    dplyr::distinct(experiment_id, W_z) |>
    dplyr::pull(W_z) |>
    stats::median()

  # For each eligible species, make an 81 x 81 candidate grid between the
  # observed minimum and maximum T and P. Keep only combinations inside the
  # convex hull of that species' observed experiment climates, and add the
  # observed climate points themselves. We therefore do not classify a species
  # by extrapolating to climate combinations where it was never studied.
  one_species <- eligible_species[[1]]
  observed_species_climates <- resident_auc_data |>
    dplyr::filter(species_id == one_species) |>
    dplyr::distinct(T_z, P_z)
  candidate_climates <- tidyr::expand_grid(
    T_z = seq(
      min(observed_species_climates$T_z),
      max(observed_species_climates$T_z),
      length.out = 81L
    ),
    P_z = seq(
      min(observed_species_climates$P_z),
      max(observed_species_climates$P_z),
      length.out = 81L
    )
  )
  observed_climate_hull <- geometry::convhulln(
    as.matrix(observed_species_climates[c("T_z", "P_z")]),
    options = "Qt"
  )
  one_species_grid <- candidate_climates |>
    dplyr::filter(geometry::inhulln(
      observed_climate_hull,
      as.matrix(dplyr::select(candidate_climates, T_z, P_z))
    )) |>
    dplyr::bind_rows(observed_species_climates) |>
    dplyr::distinct(T_z, P_z) |>
    dplyr::mutate(
      W_z = median_warming_z,
      TP_z = T_z * P_z
    )

  # At every retained T-P combination, the fitted response is calculated as
  # the fixed climate response plus that species' conditional deviation:
  #
  # fixed(T,P,W) = beta_0 + beta_T*T + beta_P*P + beta_W*W
  #                + beta_TP*T*P + beta_TW*T*W + beta_PW*P*W
  #                + beta_TPW*T*P*W
  #
  # species deviation(T,P,W) = b_0,species + b_W,species*W
  #                              + b_T,species*T + b_P,species*P
  #                              + b_TP,species*T*P
  #
  # fitted Delta AUC = fixed(T,P,W) + species deviation(T,P,W)
  #
  # This calculation is repeated at every grid point for every eligible species.
  # Prediction SEs and confidence intervals are saved for describing uncertainty,
  # but they do not determine the winner/loser/context-dependent class.
  one_species_prediction <- predict_species_on_observed_support(
    resident_climate_fit, resident_auc_data, resident_climate_formula,
    conditional_species_coefficients, one_species, median_warming_z,
    points = 81L
  )

  # Repeat that same calculation for every eligible species before finding the
  # extrema used for classification.
  all_species_predictions <- dplyr::bind_rows(lapply(
    eligible_species,
    function(species_name) {
      predict_species_on_observed_support(
        resident_climate_fit, resident_auc_data, resident_climate_formula,
        conditional_species_coefficients, species_name, median_warming_z,
        points = 81L
      )
    }
  ))

  # Finally take the minimum and maximum fitted Delta AUC within each species'
  # observed climate hull. The rule uses the fitted surface itself, not a
  # p-value and not whether a confidence interval overlaps zero.
  species_extrema <- all_species_predictions |>
    dplyr::group_by(species_id) |>
    dplyr::summarise(
      minimum_fitted_delta_auc = min(estimate),
      maximum_fitted_delta_auc = max(estimate),
      .groups = "drop"
    )
  species_classification <- species_extrema |>
    dplyr::mutate(
      fitted_class = dplyr::case_when(
        minimum_fitted_delta_auc > 0 ~ "Consistent winner",
        maximum_fitted_delta_auc < 0 ~ "Consistent loser",
        TRUE ~ "Context-dependent"
      )
    )
  # In words:
  #   Consistent winner: fitted Delta AUC is positive everywhere examined.
  #   Consistent loser:  fitted Delta AUC is negative everywhere examined.
  #   Context-dependent: the fitted surface spans or touches zero within the hull.
  # The pipeline repeats the grid prediction for all eligible species:
  species_context <- build_species_context_bundle(
    resident_climate_fit, resident_auc_data, resident_climate_formula,
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
  # Family effects are reconstructed from the community data; they are not
  # averages of the species-level Delta AUC estimates. For each experiment:
  #   1. assign every retained resident species to its reviewed family;
  #   2. sum their relative cover within family x plot x year;
  #   3. enter zero when a surveyed plot-year contains no member of that family;
  #   4. complete missing plot x year x family combinations within treatment,
  #      replacing them by the mean for that treatment x year x family, and use
  #      zero if no such mean is available;
  #   5. integrate family relative cover through time for each plot.
  # The reviewed taxonomy includes Viola palustris in Violaceae.
  #
  # The resulting plot-level family AUC is fitted separately in each experiment:
  #
  #   auc ~ ODT * family_id
  #
  # For every family, the estimate is mean AUC_warmed - mean AUC_origin and its
  # sampling variance is the squared contrast SE. There is one such contrast per
  # family; its p-value is not used in the later family model or classification.
  # The model is fitted once because it is needed to estimate these contrasts.
  # Its off-diagonal `vcov()` is not used.
  family_effects <- build_family_auc_effects(
    family_plot_auc_bundle$plot_auc,
    family_resident_membership,
    resident_model_data
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
  # The fixed and random formulas have the same meaning as in the resident
  # climate model, with family replacing species as the grouping level.
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
  # Trait records come from the reviewed TRY export. Experimental records and
  # duplicates are removed, quantitative records require TRY ErrorRisk < 4, and
  # tropical, arid, and missing climate-zone records are excluded. Units are
  # checked before summarising observations by species.
  # LMA is calculated as 1 / SLA for each observation before taking the species
  # median. LDMC and vegetative height use their species medians directly.
  #
  # Each trait is then analysed separately. Only finite, positive species
  # medians are retained; missing values are not imputed. Medians are log
  # transformed and standardized across the unique trait-complete species before
  # being copied to their effect rows. Each trait therefore has its own
  # complete-case effects, matching diagonal sampling variances, and climate-only
  # model fitted to exactly those same rows.
  trait_climate_only_formula <- ~ T_z * P_z * W_z
  trait_full_formula <- ~ Trait_z * T_z * P_z * W_z

  ldmc_model_pair <- fit_predictor_model_pair(
    ldmc_predictor_specification,
    climate_formula = trait_climate_only_formula
  )
  lma_model_pair <- fit_predictor_model_pair(
    lma_predictor_specification,
    climate_formula = trait_climate_only_formula
  )
  height_model_pair <- fit_predictor_model_pair(
    height_predictor_specification,
    climate_formula = trait_climate_only_formula
  )
  #
  # The climate-only and full models are compared by LRT, AIC, and AICc. The LRT
  # is a joint test of the trait main effect and every trait-climate interaction:
  one_trait_models <- ldmc_model_pair
  one_trait_LR <- 2 * (
    as.numeric(stats::logLik(one_trait_models$fits$full_ML$fit)) -
      as.numeric(stats::logLik(one_trait_models$fits$climate_ML$fit))
  )
  one_trait_df <-
    attr(stats::logLik(one_trait_models$fits$full_ML$fit), "df") -
    attr(stats::logLik(one_trait_models$fits$climate_ML$fit), "df")
  one_trait_p <- stats::pchisq(
    one_trait_LR, df = one_trait_df, lower.tail = FALSE
  )
  delta_AIC_full_minus_climate <-
    stats::AIC(one_trait_models$fits$full_ML$fit) -
    stats::AIC(one_trait_models$fits$climate_ML$fit)
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

  climate_only_formula <- ~ T_z * P_z * W_z
  graminoid_full_formula <- ~ GraminoidStatus * T_z * P_z * W_z
  woodiness_full_formula <- ~ Woodiness * T_z * P_z * W_z
  #
  # The full formulas include the group, T, P, and W main effects and all their
  # interactions. The species random effects remain unchanged.
  graminoid_model_pair <- fit_predictor_model_pair(
    graminoid_predictor_specification,
    climate_formula = climate_only_formula,
    full_climate_ML = resident_climate_model_ML,
    full_climate_REML = resident_climate_model
  )
  woodiness_model_pair <- fit_predictor_model_pair(
    woodiness_predictor_specification,
    climate_formula = climate_only_formula,
    full_climate_ML = resident_climate_model_ML,
    full_climate_REML = resident_climate_model
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
  # The species tree is converted into evolutionary distances between species.
  # A principal coordinates analysis condenses those distances into numerical
  # axes describing phylogenetic position. We retain and standardize the first
  # three axes. Delta AUC is not used to construct them. Their signs are fixed so
  # that rerunning the same analysis cannot arbitrarily reverse an axis.
  phylogeny_full_formula <-
    ~ T_z * P_z * W_z +
      PEV1_z + PEV2_z + PEV3_z +
      PEV1_z:W_z + PEV2_z:W_z + PEV3_z:W_z +
      PEV1_z:T_z + PEV2_z:T_z + PEV3_z:T_z +
      PEV1_z:P_z + PEV2_z:P_z + PEV3_z:P_z +
      PEV1_z:TP_z + PEV2_z:TP_z + PEV3_z:TP_z
  #
  # Each PEV can affect mean Delta AUC and can modify its relationship with W,
  # T, P, and T:P. No higher-order PEV interactions are added. Phylogeny enters
  # through these fixed predictors; it is not added as another random covariance
  # matrix.
  phylogeny_model_pair <- fit_predictor_model_pair(
    phylogeny_predictor_specification,
    climate_formula = climate_only_formula,
    full_climate_ML = resident_climate_model_ML,
    full_climate_REML = resident_climate_model
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
  # The calculation compares the two REML fits. They use exactly the same effect
  # estimates and sampling variances, so the only difference is whether the
  # predictor and its climate interactions are present as fixed effects.
  # The calculation below is repeated separately for LDMC, LMA, height,
  # graminoid status, woodiness, and phylogenetic position. Here one model pair
  # is used to show the calculation.
  one_predictor_pair <- ldmc_model_pair
  baseline_fit <- one_predictor_pair$fits$climate_REML$fit
  full_fit <- one_predictor_pair$fits$full_REML$fit

  # Each experiment enters the climate average once.
  experiment_climates <- one_predictor_pair$data |>
    dplyr::distinct(experiment_id, T_z, P_z, W_z) |>
    dplyr::mutate(TP_z = T_z * P_z)

  # The fitted random effects contain three independent blocks:
  #   1. a species intercept variance, sigma_intercept^2;
  #   2. a species warming-slope variance, sigma_W^2;
  #   3. a 3 x 3 covariance matrix G for species deviations in T, P, and T:P.
  #
  # Convert these parameters into variation in fitted Delta AUC at the observed
  # experiment climates. The calculation below is the complete calculation used
  # for each fit.
  calculate_response_heterogeneity <- function(fit, climates) {
    climate_response_covariance <- as.matrix(fit$H)
    experiment_climate_terms <-
      as.matrix(climates[c("T_z", "P_z", "TP_z")])

    H_overall <- as.numeric(fit$sigma2[[1]])
    H_warming <- mean(climates$W_z^2 * as.numeric(fit$tau2[[1]]))
    H_context <- mean(rowSums(
      (experiment_climate_terms %*% climate_response_covariance) *
        experiment_climate_terms
    ))

    c(
      H_overall = H_overall,
      H_W = H_warming,
      H_context = H_context,
      H_total = H_overall + H_warming + H_context
    )
  }

  # The three results have a direct interpretation:
  #
  #   H_overall  remaining differences among species at mean T, P, and W;
  #   H_W        differences caused by species having different W responses;
  #   H_context  differences caused by species having different T, P, and T:P
  #              responses;
  #   H_total    the sum of these three components.
  #
  # H_W and H_context must be evaluated at actual climate values. For example,
  # different species warming slopes create little separation when W is small
  # and more separation when W is large. The same logic applies to the species
  # T, P, and T:P responses. We therefore calculate the variation expressed at
  # each experiment's T, P, and W, then average across experiments.
  baseline_H <- calculate_response_heterogeneity(
    baseline_fit, experiment_climates
  )
  full_H <- calculate_response_heterogeneity(
    full_fit, experiment_climates
  )

  # For each component separately, the reported percentage is
  #
  #   100 * (climate-only heterogeneity - full-model heterogeneity)
  #         / climate-only heterogeneity.
  #
  # For example, 20% means the fitted component is 20% smaller after adding the
  # predictor; zero means no change; -10% means it is 10% larger. Percentages
  # can be unstable when the climate-only component is close to zero.
  percent_resolved <- 100 * (baseline_H - full_H) / baseline_H

  # These two results answer different questions. The model comparison asks
  # whether the predictor improves the average fitted response. The heterogeneity
  # percentages ask whether species become more similar around that average. The
  # model estimates the average relationship and the spread among species at the
  # same time; it does not first remove the predictor effect and then calculate a
  # fixed amount of leftover variation. The average fit can therefore improve
  # while the newly estimated spread among species stays similar or becomes
  # slightly larger. A negative percentage means that the fitted species spread
  # was larger in the full model, not that the predictor made species more variable
  # or that its average relationship was unsupported.
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


  # 10. MODEL R2 AND LEAVE-GRADIENT-OUT RMSE: APPENDIX FIGURE S7
  # =================================================================
  # Question: how accurately does each model predict a complete geographic
  # gradient that was not used to fit it?
  #
  # The 20 gradients are left out one at a time. The model is fitted to the other
  # 19 gradients and used to predict Delta AUC in the omitted gradient. The
  # prediction uses the fitted average response; it does not borrow a
  # species-specific response from the omitted gradient. For each predictor, the
  # full model and its climate-only comparison use exactly the same effect rows.
  one_predictor_specification <- ldmc_predictor_specification
  gradient_ids <- sort(unique(as.character(
    one_predictor_specification$data$gradient_id
  )))

  fold_results <- dplyr::bind_rows(lapply(
    gradient_ids,
    function(heldout_gradient) {
      partition <- model_performance_partition(
        one_predictor_specification$data,
        one_predictor_specification$sampling_variance_matrix,
        heldout_gradient = heldout_gradient
      )

      climate_fit <- fit_model_performance_fold(
        partition,
        climate_only_formula,
        label = "held-out-gradient climate-only REML"
      )
      full_fit <- fit_model_performance_fold(
        partition,
        one_predictor_specification$full_formula,
        label = "held-out-gradient full predictor REML",
        predictor_columns = one_predictor_specification$predictor_columns
      )

      climate_test_design <- auc_design_matrix(
        climate_only_formula,
        partition$test,
        names(stats::coef(climate_fit$fit))
      )
      full_test_design <- auc_design_matrix(
        one_predictor_specification$full_formula,
        partition$test,
        names(stats::coef(full_fit$fit))
      )
      observed <- partition$test$auc_effect
      climate_predicted <- as.numeric(
        climate_test_design %*% stats::coef(climate_fit$fit)
      )
      full_predicted <- as.numeric(
        full_test_design %*% stats::coef(full_fit$fit)
      )

      tibble::tibble(
        heldout_gradient = heldout_gradient,
        mean_observed = mean(observed),
        mean_observed_squared = mean(observed^2),
        climate_MSE = mean((climate_predicted - observed)^2),
        full_MSE = mean((full_predicted - observed)^2)
      )
    }
  ))

  # RMSE is calculated for each model by first giving every gradient equal
  # influence. We average the 20 gradient-level mean squared errors and then take
  # the square root. A gradient containing many species therefore does not
  # dominate the result merely because it contributes more effect rows.
  climate_RMSE <- sqrt(mean(fold_results$climate_MSE))
  full_RMSE <- sqrt(mean(fold_results$full_MSE))

  # Predictive R2 compares the same held-out prediction error with the variation
  # in observed Delta AUC across the held-out gradients. It can be negative when
  # prediction error is larger than that observed variation.
  heldout_mean <- mean(fold_results$mean_observed)
  heldout_variance <- mean(fold_results$mean_observed_squared) -
    heldout_mean^2
  climate_predictive_R2 <- 1 -
    mean(fold_results$climate_MSE) / heldout_variance
  full_predictive_R2 <- 1 -
    mean(fold_results$full_MSE) / heldout_variance

  # Appendix Figure S7 gives RMSE and predictive R2 a 95% uncertainty interval.
  # Complete gradients are resampled 5,000 times, preserving the study structure.
  set.seed(20260914L)
  performance_bootstrap <- replicate(5000L, {
    sampled_gradients <- sample.int(
      nrow(fold_results), nrow(fold_results), replace = TRUE
    )
    sampled <- fold_results[sampled_gradients, ]
    sampled_mean <- mean(sampled$mean_observed)
    sampled_variance <- mean(sampled$mean_observed_squared) - sampled_mean^2
    c(
      climate_RMSE = sqrt(mean(sampled$climate_MSE)),
      full_RMSE = sqrt(mean(sampled$full_MSE)),
      climate_predictive_R2 = 1 -
        mean(sampled$climate_MSE) / sampled_variance,
      full_predictive_R2 = 1 -
        mean(sampled$full_MSE) / sampled_variance
    )
  })
  performance_CI <- t(apply(
    performance_bootstrap,
    1L,
    stats::quantile,
    probs = c(0.025, 0.975), names = FALSE
  ))

  # Figure S7 also reports two summaries of the final fit. Marginal R2 is the
  # share of response variation attributed to the fixed climate and predictor
  # terms. Conditional R2 additionally includes the fitted differences among
  # species. Both include sampling variation in their denominator.
  final_fit <- one_predictor_pair$fits$full_REML$fit
  fit_weights <- model_performance_weights(
    final_fit$data, aggregation = "equal_gradient"
  )
  fixed_response <- as.numeric(final_fit$X %*% stats::coef(final_fit))
  fixed_mean <- sum(fit_weights * fixed_response)
  fixed_variance <- sum(fit_weights * (fixed_response - fixed_mean)^2)
  species_variance <- sum(
    fit_weights * model_performance_random_variance(final_fit, final_fit$data)
  )
  sampling_variance <- sum(
    fit_weights * as.numeric(Matrix::diag(final_fit$V))
  )
  total_variance <- fixed_variance + species_variance + sampling_variance
  marginal_R2 <- fixed_variance / total_variance
  conditional_R2 <-
    (fixed_variance + species_variance) / total_variance
  
  #
  # Pipeline code:
  #   R/evaluation/model_performance.R  model_performance_fit_r2()
  #   R/evaluation/model_performance.R  model_performance_partition()
  #   R/evaluation/model_performance.R  fit_model_performance_fold()
  #   R/evaluation/model_performance.R  aggregate_model_performance_statistics()
  #   R/evaluation/model_performance.R  summarise_model_cross_validation()
  }