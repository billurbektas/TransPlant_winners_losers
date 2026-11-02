# Per-gradient prediction metrics (PLAN.md Step D, "Route 2"); used by steps 06
# and 07. As PLAN.md writes it, plus one rule from MILESTONES step 06: Spearman is
# NA where every prediction in the gradient is the same (cor() would return NA
# with a warning).
metrics_gradient <- function(observed, predicted, V) {
  tibble::tibble(
    r2       = 1 - mean((predicted - observed)^2) / V,
    spearman = if (length(observed) >= 5 && length(unique(predicted)) > 1)
                 cor(observed, predicted, method = "spearman") else NA_real_,
    skill    = {
      ow <- observed > 0; pw <- predicted > 0
      if (all(ow) || all(!ow)) NA_real_ else
        mean(pw[ow]) + mean(!pw[!ow]) - 1
    }
  )
}
