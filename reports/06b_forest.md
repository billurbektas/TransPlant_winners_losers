# Step 06b: Sensitivity, random forest (held out)

Every level of the ladder except Null was refitted as a random forest (ranger 0.18.0) on the 512 ladder rows of step 05, with the variables of the step 06 formulas as predictors and species as a factor in the Species-dependent and the seven property + species forests. Each forest was fitted 16 times, once without each scored gradient, and asked to predict it, with 3000 trees, no case weights and one parameter set for every level and fold: minimum node size 20 and mtry the square root of the number of predictors, chosen from 10 settings as the one with the highest median out-of-bag R² over the 18 levels (Bills, 2 Nov). Only held-out predictions are read. They are scored as in step 06: R² against the metafor Null on the same gradient, Spearman ρ within the gradient and loser/winner skill, and set beside step 06's held-out values per gradient and as medians. The forest Null is the mean of the training effects.

## Counts

| rule | n |
| --- | --- |
| strength from year 3 (step 03b) | 1004 |
| ladder rows (step 05) | 512 |
| species | 71 |
| experiments | 39 |
| gradients training | 20 |
| gradients held out and scored (at least 5 effects) | 16 |
| forest levels | 18 |
| forests (levels x scored gradients) | 288 |

## Tuning: one parameter set for all forests

Each setting was fitted once per level on all 512 rows, and scored by its out-of-bag R² (each effect predicted by the trees that did not draw it). The setting with the highest median over the 18 levels (node20_sqrt) is used for every level and fold. Out-of-bag R² uses none of the held-out scores, so these stay unbiased; it measures prediction within known gradients (the other effects of a gradient stay in training), so its values say nothing about transfer: the forests with species score highest out of bag and lowest held out. In-sample R² would always pick the most flexible forest.

| setting | minimum node size | mtry | out-of-bag R², median over levels | lowest level | highest level | chosen |
| --- | --- | --- | --- | --- | --- | --- |
| node5_sqrt | 5 | √p | 0.079 | -0.090 | 0.195 |  |
| node5_half | 5 | p / 2 | 0.029 | -0.110 | 0.165 |  |
| node10_sqrt | 10 | √p | 0.088 | -0.088 | 0.200 |  |
| node10_half | 10 | p / 2 | 0.044 | -0.108 | 0.192 |  |
| node20_sqrt | 20 | √p | 0.096 | -0.085 | 0.213 | yes |
| node20_half | 20 | p / 2 | 0.065 | -0.104 | 0.220 |  |
| node40_sqrt | 40 | √p | 0.095 | -0.074 | 0.211 |  |
| node40_half | 40 | p / 2 | 0.070 | -0.095 | 0.221 |  |
| node80_sqrt | 80 | √p | 0.089 | -0.052 | 0.192 |  |
| node80_half | 80 | p / 2 | 0.066 | -0.077 | 0.205 |  |

## No case weights

PLAN.md's case weights 1 / SE² were dropped (Bills, 2 Nov). Precision is tied to effect size here (Spearman ρ between |strength| and SE 0.65): the most precise effects are the near-zero ones. On the 512 rows the 1 / SE² weights span a factor of 56,834, their effective number of rows, (Σw)² / Σw², is 20.5, the 5 most precise effects carry 40% of the weight and the single most precise one (Crocus vernus at DE_TransAlps) 16%; ranger draws bootstrap rows in proportion to case weights, so the weighted forests predicted about zero everywhere (first run, 2 Nov). The meta-regressions do not weight by 1 / SE² either: rma.mv, given the sampling variances only, weights by the inverse of SE² plus the random-effect variances it estimates in the same fit, which is far more even.

Spread of the forest predictions: within a held-out gradient their standard deviation is a median 42% of the observed one, and 24 of 288 forests predict one sign for every effect of the gradient (skill exactly 0 there).

## Forest beside meta-regression, held out, median over gradients

| level | metric | metafor | forest | median difference | forest higher |
| --- | --- | --- | --- | --- | --- |
| Null | R² against Null | 0.000 | 0.004 | 0.004 | 9 of 16 |
| Null | Spearman ρ |  |  |  |  |
| Null | loser/winner skill | 0.000 | 0.000 | 0.000 | 0 of 16 |
| Initial abundances | R² against Null | 0.172 | 0.051 | -0.042 | 6 of 16 |
| Initial abundances | Spearman ρ | 0.371 | 0.319 | -0.008 | 7 of 16 |
| Initial abundances | loser/winner skill | 0.147 | 0.196 | -0.004 | 8 of 16 |
| Climate | R² against Null | -0.002 | 0.006 | 0.003 | 9 of 16 |
| Climate | Spearman ρ |  |  |  |  |
| Climate | loser/winner skill |  |  |  |  |
| Climate + initial abundances | R² against Null | 0.158 | 0.056 | -0.034 | 7 of 16 |
| Climate + initial abundances | Spearman ρ | 0.353 | 0.329 | 0.008 | 9 of 16 |
| Climate + initial abundances | loser/winner skill | 0.197 | 0.188 | 0.062 | 11 of 16 |
| Species-dependent | R² against Null | 0.118 | -0.016 | -0.069 | 6 of 16 |
| Species-dependent | Spearman ρ | 0.382 | 0.387 | 0.012 | 8 of 16 |
| Species-dependent | loser/winner skill | 0.239 | 0.148 | -0.025 | 5 of 16 |
| LDMC | R² against Null | 0.140 | 0.039 | 0.020 | 10 of 16 |
| LDMC | Spearman ρ | 0.334 | 0.289 | 0.001 | 8 of 16 |
| LDMC | loser/winner skill | 0.061 | 0.172 | 0.107 | 9 of 16 |
| LMA | R² against Null | 0.158 | 0.093 | -0.024 | 7 of 16 |
| LMA | Spearman ρ | 0.303 | 0.267 | -0.033 | 7 of 16 |
| LMA | loser/winner skill | 0.093 | 0.176 | 0.021 | 8 of 16 |
| Height | R² against Null | 0.129 | 0.148 | -0.001 | 8 of 16 |
| Height | Spearman ρ | 0.332 | 0.412 | 0.040 | 9 of 16 |
| Height | loser/winner skill | 0.214 | 0.170 | -0.009 | 8 of 16 |
| Graminoid | R² against Null | 0.138 | 0.080 | 0.009 | 9 of 16 |
| Graminoid | Spearman ρ | 0.313 | 0.298 | -0.037 | 7 of 16 |
| Graminoid | loser/winner skill | 0.125 | 0.179 | 0.013 | 9 of 16 |
| Woody | R² against Null | 0.168 | 0.066 | -0.024 | 7 of 16 |
| Woody | Spearman ρ | 0.387 | 0.312 | -0.021 | 8 of 16 |
| Woody | loser/winner skill | 0.222 | 0.236 | 0.004 | 10 of 16 |
| Phylogeny | R² against Null | 0.058 | 0.046 | 0.052 | 10 of 16 |
| Phylogeny | Spearman ρ | 0.297 | 0.299 | 0.018 | 9 of 16 |
| Phylogeny | loser/winner skill | 0.046 | 0.162 | 0.026 | 8 of 16 |
| Niche | R² against Null | 0.184 | 0.121 | 0.011 | 9 of 16 |
| Niche | Spearman ρ | 0.383 | 0.355 | -0.007 | 7 of 16 |
| Niche | loser/winner skill | 0.235 | 0.239 | 0.019 | 9 of 16 |
| LDMC + species | R² against Null | 0.114 | -0.028 | -0.092 | 7 of 16 |
| LDMC + species | Spearman ρ | 0.345 | 0.347 | 0.067 | 9 of 16 |
| LDMC + species | loser/winner skill | 0.212 | 0.181 | -0.009 | 6 of 16 |
| LMA + species | R² against Null | 0.138 | -0.020 | -0.095 | 6 of 16 |
| LMA + species | Spearman ρ | 0.359 | 0.343 | 0.015 | 8 of 16 |
| LMA + species | loser/winner skill | 0.149 | 0.170 | 0.015 | 9 of 16 |
| Height + species | R² against Null | 0.091 | 0.025 | -0.006 | 8 of 16 |
| Height + species | Spearman ρ | 0.421 | 0.388 | 0.040 | 9 of 16 |
| Height + species | loser/winner skill | 0.239 | 0.127 | -0.065 | 7 of 16 |
| Graminoid + species | R² against Null | 0.104 | 0.006 | -0.045 | 7 of 16 |
| Graminoid + species | Spearman ρ | 0.378 | 0.367 | 0.044 | 9 of 16 |
| Graminoid + species | loser/winner skill | 0.192 | 0.126 | -0.018 | 8 of 16 |
| Woody + species | R² against Null | 0.106 | -0.002 | -0.070 | 6 of 16 |
| Woody + species | Spearman ρ | 0.392 | 0.396 | 0.005 | 8 of 16 |
| Woody + species | loser/winner skill | 0.239 | 0.177 | -0.041 | 5 of 16 |
| Phylogeny + species | R² against Null | 0.058 | -0.002 | -0.034 | 7 of 16 |
| Phylogeny + species | Spearman ρ | 0.356 | 0.307 | 0.016 | 9 of 16 |
| Phylogeny + species | loser/winner skill | 0.126 | 0.175 | 0.032 | 8 of 16 |
| Niche + species | R² against Null | 0.122 | 0.034 | -0.005 | 8 of 16 |
| Niche + species | Spearman ρ | 0.382 | 0.284 | -0.021 | 8 of 16 |
| Niche + species | loser/winner skill | 0.215 | 0.160 | -0.071 | 7 of 16 |

R² against Null, median held out: the forest is above the meta-regression at 5 of 18 levels.

## Seed check

Every forest refitted with a second seed: the largest change in a per-gradient value is 0.400 and the largest change in a median 0.046. Largest change in a median per metric: R² 0.012, Spearman ρ 0.036, skill 0.046; skill moves most because a prediction near zero can change sign between seeds.

## Figures

![Fig S11, forest beside meta-regression](../figures/diagnostics/06b_figS11.png)

## Warnings captured

None.

## Deviations from PLAN.md

- No case weights, instead of 1 / SE² (Bills, 2 Nov; PLAN.md edited).
- Held-out values only, no fit values (Bills, 2 Nov; PLAN.md edited).
- No partial-dependence plots (Bills, 2 Nov; PLAN.md edited).
- Forests also for the seven property + species levels (Bills, 2 Nov).
- Fig S11 written to `figures/diagnostics/` so that the report can link it, as steps 05 and 06.
