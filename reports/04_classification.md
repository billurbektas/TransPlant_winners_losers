# Step 04: Consistency across experiments (Step C, residents)

For the 74 resident species found in at least three gradients (523 effects in 39 experiments), one `rma.mv` estimated a pooled losing-winning strength per species and a species-specific among-experiment variance, on the strength from year 3 (step 03b). Each species' prediction interval for a new experiment, which carries both its among-experiment variance and the uncertainty of its pooled strength, places it as a consistent winner (interval above zero), a consistent loser (below zero) or context-dependent (spanning zero).

## Model as fitted (for reviewers)

```r
# the model: one fixed mean and one among-experiment variance per species
# rma.mv(strength, strength_se^2, mods = ~ 0 + species, random = ~ species | effect_id, struct = "DIAG")
# fitted as one rma.mv per species, which maximises the same REML likelihood (no shared parameter):
m_s <- rma.mv(yi = strength, V = strength_se^2, data = d, random = ~1 | effect_id, method = "REML")   # d = one species' effects
p <- predict(m_s)
# lower, upper = p$pi.lb, p$pi.ub = mu_s -/+ 1.96 * sqrt(tau2_s + se(mu_s)^2)
```

## Counts

| rule | n |
| --- | --- |
| effects (step 03b) | 1,004 |
| species in >= 3 gradients | 74 |
| their effects | 523 |
| experiments they occur in | 39 |
| species with tau2_s estimated at about 0 (< 1e-8) | 14 |

## Classes

| class | species |
| --- | --- |
| Consistent loser | 6 |
| Context-dependent | 66 |
| Consistent winner | 2 |

Consistent winners and losers:

| species | class | n_experiments | mu | se_mu | tau | lower | upper |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Helictochloa versicolor | Consistent loser | 5 | -0.0343 | 0.0059 | 0.0036 | -0.0478 | -0.0209 |
| Carex vaginata | Consistent loser | 3 | -0.0274 | 0.0084 | 0.0000 | -0.0437 | -0.0110 |
| Omalotheca supina | Consistent loser | 3 | -0.0121 | 0.0054 | 0.0000 | -0.0227 | -0.0014 |
| Veronica alpina | Consistent loser | 6 | -0.0084 | 0.0042 | 0.0000 | -0.0167 | -0.0001 |
| Soldanella alpina | Consistent loser | 5 | -0.0030 | 0.0007 | 0.0000 | -0.0045 | -0.0016 |
| Trifolium repens | Consistent loser | 9 | -0.0023 | 0.0007 | 0.0000 | -0.0037 | -0.0009 |
| Poa pratensis | Consistent winner | 12 | 0.0091 | 0.0032 | 0.0000 | 0.0028 | 0.0153 |
| Festuca ovina | Consistent winner | 4 | 0.0171 | 0.0067 | 0.0000 | 0.0040 | 0.0303 |

Species with tau2_s at about zero (their interval is the uncertainty of the mean alone): Campanula rotundifolia, Carex norvegica, Carex vaginata, Cerastium fontanum, Euphrasia minima, Festuca ovina, Myosotis alpestris, Omalotheca supina, Parnassia palustris, Poa pratensis, Soldanella alpina, Trifolium repens, Veronica alpina, Veronica serpyllifolia.

## Fig 2

![Fig 2](../figures/main/fig2.png)

Fig 2 b) and c): Festuca rubra and Bistorta vivipara, the highest and lowest pooled strength among species in at least 11 gradients; one row per experiment (country, origin to destination elevation), mean relative cover over years 3 and later in origin controls (hollow) and warmed plots (filled).
Fig 2 a): point = pooled strength, filled by class; line = 95% interval including the variation among experiments (a line that does not cross zero marks a consistent loser or winner), shaded by the number of gradients and labelled gradients (experiments). 18 context-dependent species with a pooled strength between -0.005 and +0.005 are not shown.

## Warnings captured

- The `label.size` argument of `geom_label()` is deprecated as of ggplot2 3.5.0.
ℹ Please use the `linewidth` argument instead.

## Deviations from PLAN.md

None: PLAN.md Step C was revised on 30 Oct (Bills) to fixed species means with predict() intervals, and Fig 2 was redesigned on 1 Nov (Bills: species overview with two species panels; PLAN.md updated). Decision approved with the plan: se_mu added to species_summary.
