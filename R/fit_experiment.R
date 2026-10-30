# One mixed model per experiment (PLAN.md Step A), shared by step 02
# (residents, reference origin_control) and step 07 (colonizers, reference
# destination_control). The model and the two emmeans calls are PLAN.md's code
# block. The residual variance is one per species (dispformula = ~ species,
# decided 30 Oct): with one variance per experiment, as lmer fits, every
# species of an experiment got the same SE, so V would weight experiments and
# not species. Added here: warnings and messages are captured for the step
# report; a fit that fails to converge or has a non-positive-definite Hessian
# stops; the contrasts use asymptotic degrees of freedom in every experiment
# (decided 30 Oct).
fit_experiment <- function(d, control_level) {
  d <- droplevels(d)
  notes <- character()
  capture <- function(expr) withCallingHandlers(expr,
    warning = function(w) { notes <<- c(notes, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) { notes <<- c(notes, trimws(conditionMessage(m))); invokeRestart("muffleMessage") })

  m <- capture(glmmTMB::glmmTMB(
    relative_cover ~ treatment * year * species + (1 | plot_id:species),
    dispformula = ~ species,
    data = d, REML = TRUE
  ))
  if (m$fit$convergence != 0 || !isTRUE(m$sdr$pdHess))
    stop("glmmTMB did not converge (code ", m$fit$convergence, ", pdHess ", m$sdr$pdHess, "): ",
         m$fit$message, call. = FALSE)

  em <- capture(emmeans::emmeans(m, ~ treatment | species * year, df = Inf))
  # year-specific contrast, used in Step B
  by_year <- capture(emmeans::contrast(em, method = "revpairwise")) |>
    as.data.frame() |>
    dplyr::filter(contrast == paste("warmed -", control_level))
  # contrast averaged over years, used as strength from Step C onward
  em_mean <- capture(emmeans::emmeans(m, ~ treatment | species, df = Inf))
  strength <- capture(emmeans::contrast(em_mean, method = "revpairwise")) |>
    as.data.frame() |>
    dplyr::filter(contrast == paste("warmed -", control_level))

  list(model = m, by_year = by_year, strength = strength, warnings = unique(notes),
       convergence = m$fit$convergence, pdHess = m$sdr$pdHess,
       sd_plot_species = attr(glmmTMB::VarCorr(m)$cond[["plot_id:species"]], "stddev")[[1]])
}
