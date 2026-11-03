# Step 05: The ladder, fitted (Step D, route 1, residents)

The model ladder was fitted on the strength from year 3 (step 03b), on one set of rows (2 Nov): the 512 effects of the 71 residents found in at least 3 gradients that carry every property (39 experiments, 20 gradients). Null, Start, Climate, Climate + start, Species-dependent and one property model per property (LDMC, LMA, height, graminoid, woody, phylogeny as PEV1 to PEV3, niche as T_opt and P_opt) with the species slopes kept, 12 REML fits, all on the same rows. Each model reports τ² (species × site), σ²_experiment (shared site shift) and the two H_species pieces as raw values. The T, P, TP slope block is diagonal in every slope fit (Bills, 1 Nov). Fixed effects carry model-based (REML Wald) standard errors (Bills, 2 Nov; CR1 before). Fits ran in parallel on 7 cores, 11.8 fit-minutes in all; this run took 1 min.

## Slope structure

The 1 Nov GEN timing run (Species-dependent on the 1,655 rows of the 28 Oct pools, started from the gate fit's own optimum) did not converge, so the T, P, TP slope block is diagonal in every slope fit; the case is in MILESTONES step 05 notes. It is not repeated on the new rows (Bills, 2 Nov). The experiments' T and P are correlated at r = -0.49 on the ladder rows.

## Model as fitted (for reviewers)

```r
re_base   <- list(~ 1 | species, ~ 1 | experiment_id, ~ 1 | effect_id)
re_slopes <- c(re_base, list(~ 0 + W | species, ~ 0 + T + P + TP | species))
rma.mv(yi = strength, V = strength_se^2, mods = <formula>, random = <re_base or re_slopes>,
       struct = c("GEN", "GDIAG"), data = <rowset>, method = "REML", sparse = TRUE)
# struct only acts on re_slopes: W slope; T, P, TP slopes diagonal. No starting values: metafor's default for every fit
# Null ~ 1; Start ~A0_origin + A0_dest + A0_origin:W + A0_dest:W; Climate ~W * T * P
# Climate + start and Species-dependent ~W * T * P + A0_origin + A0_dest + A0_origin:W + A0_dest:W
# property models: Climate + start + X + X:W + X:T + X:P + X:TP per variable, with re_slopes
coef(summary(m))   # model-based (REML Wald) intervals for Table S5 and Fig S11
```

## Counts per rowset

| rowset | effects | species | experiments |
| --- | --- | --- | --- |
| all | 512 | 71 | 39 |
| ldmc | 512 | 71 | 39 |
| lma | 512 | 71 | 39 |
| height | 512 | 71 | 39 |
| growth_form | 512 | 71 | 39 |
| phylogeny | 512 | 71 | 39 |
| niche | 512 | 71 | 39 |

Each matches the expectation checked in the script. Every property rowset equals the ladder rows, so Null to Species-dependent are fitted once, on `all` (marked `reused` elsewhere).

Step 02's standardisation, kept (Bills, 2 Nov), read on the ladder rows:

| variable | across | mean | sd |
| --- | --- | --- | --- |
| W | the 39 experiments | 0.07 | 1.00 |
| T | the 39 experiments | 0.02 | 1.03 |
| P | the 39 experiments | 0.03 | 0.92 |
| A0_origin | the 512 effects | 0.13 | 1.03 |
| A0_dest | the 512 effects | 0.15 | 1.10 |
| ldmc_z | the 71 species | 0.13 | 0.84 |
| lma_z | the 71 species | -0.08 | 0.84 |
| height_z | the 71 species | -0.15 | 0.86 |
| PEV1 | the 71 species | 0.22 | 1.06 |
| PEV2 | the 71 species | -0.12 | 1.13 |
| PEV3 | the 71 species | 0.02 | 1.04 |
| T_opt | the 71 species | -0.02 | 0.86 |
| P_opt | the 71 species | -0.16 | 0.70 |

## Heterogeneity per model (×10⁻⁴)

| rowset | model | fitted_on | effects | species | τ² | σ²_experiment | H intercept | H slope | struct | minutes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| all | Null | all | 512 | 71 | 4.802 | 0.201 | 1.165 | 0.000 | none (no slopes) | 0.2 |
| all | Start | all | 512 | 71 | 3.360 | 0.107 | 0.728 | 0.000 | none (no slopes) | 0.2 |
| all | Climate | all | 512 | 71 | 4.902 | 0.000 | 1.154 | 0.000 | none (no slopes) | 0.4 |
| all | Climate + start | all | 512 | 71 | 3.296 | 0.000 | 0.742 | 0.000 | none (no slopes) | 0.4 |
| all | Species-dependent | all | 512 | 71 | 2.842 | 0.000 | 0.584 | 0.692 | DIAG | 1.3 |
| ldmc | Null | all | 512 | 71 | 4.802 | 0.201 | 1.165 | 0.000 | none (no slopes) | reused |
| ldmc | Start | all | 512 | 71 | 3.360 | 0.107 | 0.728 | 0.000 | none (no slopes) | reused |
| ldmc | Climate | all | 512 | 71 | 4.902 | 0.000 | 1.154 | 0.000 | none (no slopes) | reused |
| ldmc | Climate + start | all | 512 | 71 | 3.296 | 0.000 | 0.742 | 0.000 | none (no slopes) | reused |
| ldmc | Species-dependent | all | 512 | 71 | 2.842 | 0.000 | 0.584 | 0.692 | DIAG | reused |
| ldmc | LDMC | ldmc | 512 | 71 | 2.834 | 0.000 | 0.643 | 0.749 | DIAG | 1.5 |
| lma | Null | all | 512 | 71 | 4.802 | 0.201 | 1.165 | 0.000 | none (no slopes) | reused |
| lma | Start | all | 512 | 71 | 3.360 | 0.107 | 0.728 | 0.000 | none (no slopes) | reused |
| lma | Climate | all | 512 | 71 | 4.902 | 0.000 | 1.154 | 0.000 | none (no slopes) | reused |
| lma | Climate + start | all | 512 | 71 | 3.296 | 0.000 | 0.742 | 0.000 | none (no slopes) | reused |
| lma | Species-dependent | all | 512 | 71 | 2.842 | 0.000 | 0.584 | 0.692 | DIAG | reused |
| lma | LMA | lma | 512 | 71 | 2.804 | 0.000 | 0.499 | 0.730 | DIAG | 1.4 |
| height | Null | all | 512 | 71 | 4.802 | 0.201 | 1.165 | 0.000 | none (no slopes) | reused |
| height | Start | all | 512 | 71 | 3.360 | 0.107 | 0.728 | 0.000 | none (no slopes) | reused |
| height | Climate | all | 512 | 71 | 4.902 | 0.000 | 1.154 | 0.000 | none (no slopes) | reused |
| height | Climate + start | all | 512 | 71 | 3.296 | 0.000 | 0.742 | 0.000 | none (no slopes) | reused |
| height | Species-dependent | all | 512 | 71 | 2.842 | 0.000 | 0.584 | 0.692 | DIAG | reused |
| height | Height | height | 512 | 71 | 2.930 | 0.000 | 0.564 | 0.540 | DIAG | 1.0 |
| growth_form | Null | all | 512 | 71 | 4.802 | 0.201 | 1.165 | 0.000 | none (no slopes) | reused |
| growth_form | Start | all | 512 | 71 | 3.360 | 0.107 | 0.728 | 0.000 | none (no slopes) | reused |
| growth_form | Climate | all | 512 | 71 | 4.902 | 0.000 | 1.154 | 0.000 | none (no slopes) | reused |
| growth_form | Climate + start | all | 512 | 71 | 3.296 | 0.000 | 0.742 | 0.000 | none (no slopes) | reused |
| growth_form | Species-dependent | all | 512 | 71 | 2.842 | 0.000 | 0.584 | 0.692 | DIAG | reused |
| growth_form | Graminoid | growth_form | 512 | 71 | 2.915 | 0.000 | 0.622 | 0.665 | DIAG | 1.5 |
| growth_form | Woody | growth_form | 512 | 71 | 2.913 | 0.000 | 0.607 | 0.683 | DIAG | 1.4 |
| phylogeny | Null | all | 512 | 71 | 4.802 | 0.201 | 1.165 | 0.000 | none (no slopes) | reused |
| phylogeny | Start | all | 512 | 71 | 3.360 | 0.107 | 0.728 | 0.000 | none (no slopes) | reused |
| phylogeny | Climate | all | 512 | 71 | 4.902 | 0.000 | 1.154 | 0.000 | none (no slopes) | reused |
| phylogeny | Climate + start | all | 512 | 71 | 3.296 | 0.000 | 0.742 | 0.000 | none (no slopes) | reused |
| phylogeny | Species-dependent | all | 512 | 71 | 2.842 | 0.000 | 0.584 | 0.692 | DIAG | reused |
| phylogeny | Phylogeny | phylogeny | 512 | 71 | 2.925 | 0.000 | 0.616 | 0.756 | DIAG | 1.5 |
| niche | Null | all | 512 | 71 | 4.802 | 0.201 | 1.165 | 0.000 | none (no slopes) | reused |
| niche | Start | all | 512 | 71 | 3.360 | 0.107 | 0.728 | 0.000 | none (no slopes) | reused |
| niche | Climate | all | 512 | 71 | 4.902 | 0.000 | 1.154 | 0.000 | none (no slopes) | reused |
| niche | Climate + start | all | 512 | 71 | 3.296 | 0.000 | 0.742 | 0.000 | none (no slopes) | reused |
| niche | Species-dependent | all | 512 | 71 | 2.842 | 0.000 | 0.584 | 0.692 | DIAG | reused |
| niche | Niche | niche | 512 | 71 | 2.860 | 0.000 | 0.575 | 0.750 | DIAG | 1.1 |

All 8 slope fits converged with the diagonal T, P, TP block, each from metafor's default start (Bills, 1 Nov). Across the 8 fits the T slope variance ranges 0.000 to 0.274 ×10⁻⁴; the TP slope variance 0.030 to 0.193 ×10⁻⁴; only the P slope variance is at the boundary (at most 9.1e-11).

## Checks

τ² at Climate against Null (climate's main effects cannot reach the species × site piece; stop beyond 5%):

| rowset | τ² Null | τ² Climate | change |
| --- | --- | --- | --- |
| all | 4.802 | 4.902 | +2.1% |

Each property model against Species-dependent on the same rows, as a share of that model's τ² + σ²_experiment (a property reaches only H_species; stop beyond 5%, Bills 1 Nov):

| property | rowset | τ² + σ²_experiment (×10⁻⁴) | τ² move | σ²_experiment move |
| --- | --- | --- | --- | --- |
| LDMC | ldmc | 2.842 | -0.27% | -0.00% |
| LMA | lma | 2.842 | -1.34% | -0.00% |
| Height | height | 2.842 | +3.11% | +0.00% |
| Graminoid | growth_form | 2.842 | +2.57% | -0.00% |
| Woody | growth_form | 2.842 | +2.48% | -0.00% |
| Phylogeny | phylogeny | 2.842 | +2.94% | -0.00% |
| Niche | niche | 2.842 | +0.64% | -0.00% |

σ²_species at Null on all rows is 1.165 ×10⁻⁴, off the boundary; the Step C spread of μ_s (step 04, 74 species) is 5.452 ×10⁻⁴. Random structure, property terms (X, X:W, X:T, X:P, X:TP) and A0 interactions (A0:W only) are as PLAN.md writes them. Subspecies niche check (26 Oct): neither Solidago virgaurea nor Ranunculus acris is a resident at CH_Lavey under the original pools, so there is nothing to refit.

## Effective weights (all rows)

| weights | max/min within experiment (median) | share, low-A0_origin half | share, top 10% of effects |
| --- | --- | --- | --- |
| raw 1/V | 9.9 | 0.95 | 0.84 |
| Null | 2.3 | 0.58 | 0.13 |
| Start | 2.8 | 0.59 | 0.14 |
| Climate | 2.4 | 0.58 | 0.13 |
| Climate + start | 2.8 | 0.60 | 0.14 |
| Species-dependent | 2.9 | 0.60 | 0.16 |

The low-A0_origin half is the effects at or below the median A0_origin. Species-dependent weights include the slope variance at each effect's climate.

## Figures

![Fig S11, coefficients](../figures/appendix/figS11.png)

Fig 6 (the ladder) is drawn by step 06; the fitted surfaces are not drawn (Bills, 2 Nov).

## Warnings captured

- all / Null: Model does not contain an '~ inner | outer' term, so 'struct' argument is disregaded.
- all / Start: Model does not contain an '~ inner | outer' term, so 'struct' argument is disregaded.
- all / Climate: Model does not contain an '~ inner | outer' term, so 'struct' argument is disregaded.
- all / Climate + start: Model does not contain an '~ inner | outer' term, so 'struct' argument is disregaded.
- `geom_errorbarh()` was deprecated in ggplot2 4.0.0.
ℹ Please use the `orientation` argument of `geom_errorbar()` instead.

## Deviations from PLAN.md

- H_species tests `m$withG` where PLAN.md tests `is.null(m$tau2)`: rma.mv returns tau2 = 0, not NULL, for a model without slopes. Same quantity.
- The coefficient figure is appendix Fig S11, written to `figures/appendix/figS11.png`.
- Decided 1 and 2 Nov (Bills) and in PLAN.md, not deviations: the T, P, TP slope block diagonal; no ML refits and no AICc column; no profile CIs; the 5% property check measured against τ² + σ²_experiment; Fig 4 with twelve columns, drawn by step 06; the ladder on residents in at least 3 gradients that carry every property, one set of rows for every model; step 02's standardisation kept.
