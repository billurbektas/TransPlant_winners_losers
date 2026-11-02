# Step 06b: Sensitivity, random forest (held out)

Every level of the ladder except Null was refitted as a random forest (ranger 0.18.0) on the 512 ladder rows of step 05, with the variables of the step 06 formulas as predictors and species as a factor in the Species-dependent and the seven property + species forests. Each forest was fitted 16 times, once without each scored gradient, and asked to predict it, with fixed settings (1000 trees, minimum node size 40, mtry the square root of the number of predictors, case weights 1 / SE²). Only held-out predictions are read. They are scored as in step 06: R² against the metafor Null on the same gradient, Spearman ρ within the gradient and loser/winner skill, and set beside step 06's held-out values per gradient and as medians. The forest Null is the 1 / SE²-weighted mean of the training effects.

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

## Case weights

The case weights 1 / SE² (PLAN.md) span a factor of 56,834 on the 512 rows. Their effective number of rows, (Σw)² / Σw², is 20.5; the 5 most precise effects carry 40% of the weight and the single most precise one (Crocus vernus at DE_TransAlps) 16%. ranger draws each tree's bootstrap rows in proportion to these weights, so the trees are built mostly from a few effects. The meta-regressions weight by 1 / (SE² + the random-effect variances), which is far more even.

The forest predictions are strongly shrunk towards the weighted mean: within a held-out gradient their standard deviation is a median 16% of the observed one, and 94 of 288 forests predict one sign for every effect of the gradient (skill exactly 0 there). Low forest skill therefore reflects these settings (1 / SE² weights, minimum node size 40) rather than a test of shape.

## Forest beside meta-regression, held out, median over gradients

| level | metric | metafor | forest | median difference | forest higher |
| --- | --- | --- | --- | --- | --- |
| Null | R² against Null | 0.000 | -0.002 | -0.002 | 8 of 16 |
| Null | Spearman ρ |  |  |  |  |
| Null | loser/winner skill | 0.000 | 0.000 | 0.000 | 0 of 16 |
| Initial abundances | R² against Null | 0.172 | 0.115 | -0.039 | 7 of 16 |
| Initial abundances | Spearman ρ | 0.371 | 0.371 | -0.015 | 7 of 16 |
| Initial abundances | loser/winner skill | 0.147 | 0.071 | -0.104 | 6 of 16 |
| Climate | R² against Null | -0.002 | -0.003 | 0.006 | 9 of 16 |
| Climate | Spearman ρ |  |  |  |  |
| Climate | loser/winner skill |  |  |  |  |
| Climate + initial abundances | R² against Null | 0.158 | 0.117 | -0.001 | 8 of 16 |
| Climate + initial abundances | Spearman ρ | 0.353 | 0.340 | -0.045 | 7 of 16 |
| Climate + initial abundances | loser/winner skill | 0.197 | 0.017 | -0.123 | 5 of 16 |
| Species-dependent | R² against Null | 0.118 | 0.097 | 0.001 | 8 of 16 |
| Species-dependent | Spearman ρ | 0.382 | 0.327 | 0.016 | 9 of 16 |
| Species-dependent | loser/winner skill | 0.239 | -0.004 | -0.207 | 3 of 16 |
| LDMC | R² against Null | 0.140 | 0.098 | -0.021 | 7 of 16 |
| LDMC | Spearman ρ | 0.334 | 0.373 | 0.019 | 8 of 16 |
| LDMC | loser/winner skill | 0.061 | 0.000 | -0.123 | 6 of 16 |
| LMA | R² against Null | 0.158 | 0.098 | -0.018 | 7 of 16 |
| LMA | Spearman ρ | 0.303 | 0.339 | -0.010 | 6 of 16 |
| LMA | loser/winner skill | 0.093 | 0.000 | -0.123 | 6 of 16 |
| Height | R² against Null | 0.129 | 0.100 | -0.012 | 8 of 16 |
| Height | Spearman ρ | 0.332 | 0.400 | 0.013 | 8 of 16 |
| Height | loser/winner skill | 0.214 | 0.000 | -0.224 | 4 of 16 |
| Graminoid | R² against Null | 0.138 | 0.105 | 0.010 | 8 of 16 |
| Graminoid | Spearman ρ | 0.313 | 0.317 | -0.040 | 7 of 16 |
| Graminoid | loser/winner skill | 0.125 | 0.000 | -0.153 | 4 of 16 |
| Woody | R² against Null | 0.168 | 0.108 | -0.004 | 8 of 16 |
| Woody | Spearman ρ | 0.387 | 0.303 | -0.051 | 5 of 16 |
| Woody | loser/winner skill | 0.222 | 0.000 | -0.157 | 3 of 16 |
| Phylogeny | R² against Null | 0.058 | 0.089 | -0.035 | 7 of 16 |
| Phylogeny | Spearman ρ | 0.297 | 0.316 | 0.018 | 8 of 16 |
| Phylogeny | loser/winner skill | 0.046 | 0.000 | -0.082 | 5 of 16 |
| Niche | R² against Null | 0.184 | 0.100 | -0.013 | 8 of 16 |
| Niche | Spearman ρ | 0.383 | 0.350 | -0.036 | 5 of 16 |
| Niche | loser/winner skill | 0.235 | 0.000 | -0.181 | 4 of 16 |
| LDMC + species | R² against Null | 0.114 | 0.092 | -0.024 | 7 of 16 |
| LDMC + species | Spearman ρ | 0.345 | 0.327 | -0.002 | 8 of 16 |
| LDMC + species | loser/winner skill | 0.212 | 0.000 | -0.212 | 7 of 16 |
| LMA + species | R² against Null | 0.138 | 0.089 | -0.028 | 7 of 16 |
| LMA + species | Spearman ρ | 0.359 | 0.348 | 0.033 | 9 of 16 |
| LMA + species | loser/winner skill | 0.149 | 0.000 | -0.175 | 4 of 16 |
| Height + species | R² against Null | 0.091 | 0.095 | 0.009 | 8 of 16 |
| Height + species | Spearman ρ | 0.421 | 0.348 | 0.051 | 9 of 16 |
| Height + species | loser/winner skill | 0.239 | 0.003 | -0.204 | 4 of 16 |
| Graminoid + species | R² against Null | 0.104 | 0.092 | 0.020 | 9 of 16 |
| Graminoid + species | Spearman ρ | 0.378 | 0.338 | 0.022 | 10 of 16 |
| Graminoid + species | loser/winner skill | 0.192 | 0.000 | -0.217 | 4 of 16 |
| Woody + species | R² against Null | 0.106 | 0.095 | 0.002 | 8 of 16 |
| Woody + species | Spearman ρ | 0.392 | 0.320 | -0.018 | 8 of 16 |
| Woody + species | loser/winner skill | 0.239 | 0.000 | -0.207 | 4 of 16 |
| Phylogeny + species | R² against Null | 0.058 | 0.094 | 0.032 | 9 of 16 |
| Phylogeny + species | Spearman ρ | 0.356 | 0.306 | 0.025 | 9 of 16 |
| Phylogeny + species | loser/winner skill | 0.126 | 0.000 | -0.141 | 6 of 16 |
| Niche + species | R² against Null | 0.122 | 0.090 | 0.029 | 9 of 16 |
| Niche + species | Spearman ρ | 0.382 | 0.323 | -0.033 | 5 of 16 |
| Niche + species | loser/winner skill | 0.215 | 0.000 | -0.152 | 7 of 16 |

R² against Null, median held out: the forest is above the meta-regression at 8 of 18 levels.

## Seed check

Every forest refitted with a second seed: the largest change in a per-gradient value is 0.277 and the largest change in a median 0.031.

## Figures

![Fig S11, forest beside meta-regression](../figures/diagnostics/06b_figS11.png)

## Warnings captured

None.

## Deviations from PLAN.md

- Held-out values only, no fit values (Bills, 2 Nov; PLAN.md edited).
- No partial-dependence plots (Bills, 2 Nov; PLAN.md edited).
- Forests also for the seven property + species levels (Bills, 2 Nov).
- Fig S11 written to `figures/diagnostics/` so that the report can link it, as steps 05 and 06.
