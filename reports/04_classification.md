# Step 04: Consistency across experiments (Step C, residents)

For the 137 resident species found in at least four experiments (998 effects in 43 experiments), one `rma.mv` estimated a pooled losing-winning strength per species and a species-specific among-experiment variance, on the strength from year 3 (step 03b). Each species' prediction interval for a new experiment, which carries both its among-experiment variance and the uncertainty of its pooled strength, places it as a consistent winner (interval above zero), a consistent loser (below zero) or context-dependent (spanning zero).

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
| effects (step 03b) | 1,655 |
| species with >= 4 experiments | 137 |
| their effects | 998 |
| experiments they occur in | 43 |
| species with tau2_s estimated at about 0 (< 1e-8) | 40 |

## Classes

| class | species |
| --- | --- |
| Consistent loser | 14 |
| Context-dependent | 122 |
| Consistent winner | 1 |

Consistent winners and losers:

| species | class | n_experiments | mu | se_mu | tau | lower | upper |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Festuca violacea | Consistent loser | 4 | -0.0666 | 0.0164 | 0.0272 | -0.1288 | -0.0043 |
| Potentilla nivea | Consistent loser | 4 | -0.0375 | 0.0094 | 0.0138 | -0.0702 | -0.0048 |
| Helictochloa versicolor | Consistent loser | 5 | -0.0343 | 0.0058 | 0.0033 | -0.0474 | -0.0212 |
| Agrostis stolonifera | Consistent loser | 5 | -0.0300 | 0.0060 | 0.0000 | -0.0418 | -0.0182 |
| Alpagrostis alpina | Consistent loser | 4 | -0.0173 | 0.0042 | 0.0026 | -0.0269 | -0.0076 |
| Carex bigelowii | Consistent loser | 6 | -0.0088 | 0.0035 | 0.0000 | -0.0157 | -0.0018 |
| Salix reticulata | Consistent loser | 4 | -0.0082 | 0.0023 | 0.0000 | -0.0127 | -0.0037 |
| Veronica alpina | Consistent loser | 6 | -0.0078 | 0.0025 | 0.0000 | -0.0128 | -0.0029 |
| Rhinanthus minor | Consistent loser | 7 | -0.0078 | 0.0028 | 0.0018 | -0.0143 | -0.0013 |
| Omalotheca supina | Consistent loser | 5 | -0.0077 | 0.0029 | 0.0020 | -0.0146 | -0.0008 |
| Soldanella alpina | Consistent loser | 5 | -0.0028 | 0.0007 | 0.0000 | -0.0041 | -0.0015 |
| Trifolium badium | Consistent loser | 4 | -0.0022 | 0.0004 | 0.0000 | -0.0029 | -0.0014 |
| Silene flos-cuculi | Consistent loser | 5 | -0.0016 | 0.0005 | 0.0000 | -0.0026 | -0.0006 |
| Primula veris | Consistent loser | 6 | -0.0011 | 0.0003 | 0.0000 | -0.0017 | -0.0004 |
| Helianthemum nummularium | Consistent winner | 5 | 0.0024 | 0.0007 | 0.0000 | 0.0010 | 0.0038 |

Species with tau2_s at about zero (their interval is the uncertainty of the mean alone): Agrostis stolonifera, Alchemilla conjuncta, Alopecurus pratensis, Androsace septentrionalis, Bartsia alpina, Campanula rotundifolia, Carex bigelowii, Carex myosuroides, Carex norvegica, Centaurea jacea, Cerastium fontanum, Dactylorhiza viridis, Empetrum nigrum, Euphrasia minima, Festuca halleri, Fragaria virginiana, Gentiana nivalis, Gentianella campestris, Helianthemum nummularium, Hieracium alpinum, Kalmia procumbens, Leucanthemum ircutianum, Luzula spicata, Oreojuncus trifidus, Phleum alpinum, Pinguicula vulgaris, Primula veris, Prunella vulgaris, Ranunculus brotherusii, Rumex acetosa, Salix reticulata, Silene dioica, Silene flos-cuculi, Soldanella alpina, Thesium alpinum, Trifolium badium, Vaccinium uliginosum, Veronica alpina, Vicia cracca, Viola riviniana.

## Fig 2

![Fig 2](../figures/main/fig2.png)

Fig 2 b) and c): Festuca rubra and Bistorta vivipara, the highest and lowest pooled strength among species in at least 12 gradients; one row per experiment (country, origin to destination elevation), mean relative cover over years 3 and later in origin controls (hollow) and warmed plots (filled).
Fig 2 a): point = pooled strength, filled by class; line = 95% interval including the variation among experiments (a line that does not cross zero marks a consistent loser or winner), shaded by the number of gradients and labelled gradients (experiments). 63 context-dependent species with a pooled strength between -0.005 and +0.005 are not shown.

## Warnings captured

- Refitted with the bobyqa optimizer after nlminb did not declare convergence (tau2_s at 0): Bartsia alpina.
- The `label.size` argument of `geom_label()` is deprecated as of ggplot2 3.5.0.
ℹ Please use the `linewidth` argument instead.

## Deviations from PLAN.md

None: PLAN.md Step C was revised on 30 Oct (Bills) to fixed species means with predict() intervals, and Fig 2 was redesigned on 1 Nov (Bills: species overview with two species panels; PLAN.md updated). Decision approved with the plan: se_mu added to species_summary.
