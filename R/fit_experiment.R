# One mixed model per experiment (PLAN.md Step A), shared by step 02
# (residents, reference origin_control) and step 07 (colonizers, reference
# destination_control). The model and the two emmeans calls are PLAN.md's code
# block; added here: warnings and messages of lmer and emmeans are captured for
# the step report, a fit that fails to converge stops, and the contrasts use
# asymptotic degrees of freedom in every experiment (decided 30 Oct: emmeans
# would otherwise use Kenward-Roger below 3,000 rows only, so SEs would be
# computed two ways).
fit_experiment <- function(d, control_level) {
  d <- droplevels(d)
  notes <- character()
  capture <- function(expr) withCallingHandlers(expr,
    warning = function(w) { notes <<- c(notes, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) { notes <<- c(notes, trimws(conditionMessage(m))); invokeRestart("muffleMessage") })

  m <- capture(lme4::lmer(
    relative_cover ~ treatment * year * species + (1 | plot_id:species),
    data = d
  ))
  conv <- m@optinfo$conv$lme4$messages
  if (any(grepl("failed to converge", conv)))
    stop("lmer did not converge: ", paste(conv, collapse = "; "), call. = FALSE)

  em <- capture(emmeans::emmeans(m, ~ treatment | species * year, lmer.df = "asymptotic"))
  # year-specific contrast, used in Step B
  by_year <- capture(emmeans::contrast(em, method = "revpairwise")) |>
    as.data.frame() |>
    dplyr::filter(contrast == paste("warmed -", control_level))
  # contrast averaged over years, used as strength from Step C onward
  em_mean <- capture(emmeans::emmeans(m, ~ treatment | species, lmer.df = "asymptotic"))
  strength <- capture(emmeans::contrast(em_mean, method = "revpairwise")) |>
    as.data.frame() |>
    dplyr::filter(contrast == paste("warmed -", control_level))

  list(model = m, by_year = by_year, strength = strength, warnings = unique(notes),
       singular = lme4::isSingular(m))
}
