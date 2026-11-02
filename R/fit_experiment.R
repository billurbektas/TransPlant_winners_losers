# One mixed model per experiment (PLAN.md Step A), shared by step 02
# (residents, reference origin_control) and step 07 (colonizers, reference
# destination_control). The model and the two emmeans calls are PLAN.md's code
# block. The residual variance is one per species (dispformula = ~ species,
# decided 30 Oct): with one variance per experiment, as lmer fits, every
# species of an experiment got the same SE, so V would weight experiments and
# not species. Added here: warnings and messages are captured for the step
# report; a fit that fails to converge or has a non-positive-definite Hessian
# stops. Degrees of freedom of every contrast: control plots + warmed plots - 2
# of the experiment (decided 30 Oct), because treatment is applied to whole
# plots, so the plot is the replicate; the yearly surveys of a plot are repeated
# measures. The plot x species random intercept is one variance shared by all
# species (kept, decided 30 Oct).
fit_experiment <- function(d, control_level) {
  d <- droplevels(d)
  df_plots <- dplyr::n_distinct(d$plot_id[d$treatment == control_level]) +
    dplyr::n_distinct(d$plot_id[d$treatment == "warmed"]) - 2
  notes <- character()
  capture <- function(expr) withCallingHandlers(expr,
    warning = function(w) { notes <<- c(notes, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) { notes <<- c(notes, trimws(conditionMessage(m))); invokeRestart("muffleMessage") })

  # an experiment with one species (four colonizer pools, step 07a) gets the same
  # model without the species terms, which a one-level factor cannot carry
  one_species <- nlevels(factor(d$species)) == 1
  m <- capture(glmmTMB::glmmTMB(
    if (one_species) relative_cover ~ treatment * year + (1 | plot_id) else relative_cover ~ treatment * year * species + (1 | plot_id:species),
    dispformula = if (one_species) ~ 1 else ~ species,
    data = d, REML = TRUE
  ))
  if (m$fit$convergence != 0 || !isTRUE(m$sdr$pdHess))
    stop("glmmTMB did not converge (code ", m$fit$convergence, ", pdHess ", m$sdr$pdHess, "): ",
         m$fit$message, call. = FALSE)

  # with one species the grid has no species column; it is added back by name
  add_species <- function(x) if (one_species) cbind(species = factor(unique(d$species)), x) else x
  em <- capture(emmeans::emmeans(m, if (one_species) ~ treatment | year else ~ treatment | species * year, df = df_plots))
  # year-specific contrast, used in Step B
  by_year <- capture(emmeans::contrast(em, method = "revpairwise")) |>
    as.data.frame() |> add_species() |>
    dplyr::filter(contrast == paste("warmed -", control_level))
  # contrast averaged over years, used as strength from Step C onward
  em_mean <- capture(emmeans::emmeans(m, if (one_species) ~ treatment else ~ treatment | species, df = df_plots))
  strength <- capture(emmeans::contrast(em_mean, method = "revpairwise")) |>
    as.data.frame() |> add_species() |>
    dplyr::filter(contrast == paste("warmed -", control_level))

  # residual SD per species (predict type "disp" is the SD for the Gaussian family)
  one_per_species <- d[!duplicated(d$species), ]
  sigma_species <- stats::setNames(stats::predict(m, newdata = one_per_species, type = "disp"),
                                   as.character(one_per_species$species))

  list(model = m, by_year = by_year, strength = strength, warnings = unique(notes),
       convergence = m$fit$convergence, pdHess = m$sdr$pdHess, df_plots = df_plots,
       sd_plot_species = attr(glmmTMB::VarCorr(m)$cond[[if (one_species) "plot_id" else "plot_id:species"]], "stddev")[[1]],
       sigma_species = sigma_species, one_species = one_species)
}
