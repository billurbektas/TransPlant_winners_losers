# Step 12: Final review

A critical final pass over the whole repository (Bills, 3 Nov): every script of the prep pipeline and the analysis read line by line; everything run that does not fit a model or download data; the documentation checked against the code; every open item re-classified; and statistical weak points reported without changing any choice. `legacy/` was retired. What was fixed changes no result; what would change a result is listed for Bills under *Decisions*.

## Runs

| check | outcome |
| --- | --- |
| prep pipeline rebuilt from an empty store in a clean clone without `legacy/` | all 75 targets ran; every file in `data/derived/prep/` and every P report byte-identical to the working copy |
| all 15 analysis steps run from the cached fits in the clean clone | every step finished; every tracked output byte-identical, apart from figures whose two legends swapped order between runs (Figs 2, 3, 4, S10; fixed here) |
| two consecutive full runs in the clean clone, then the figure steps (04, 07b, 10) three more times | identical outputs; Figs 2, 3, 4, S3 and S10 byte-identical to the versions reviewed |
| cached ladder fits (05, 07c) against today's rows | 12 of 12 and 12 of 12 identical in strength and variance |
| cached prediction folds (06, 07d) against today's rows | 192 of 192 and 96 of 96 identical |
| package versions against renv.lock | every package the code loads is in the lock at the locked version (R 4.6.1); renv itself is not activated (no project library) |

## Fixed in this step (no result changes)

- **Nondeterministic figures.** Figs 2, 3, 4 and S10 swapped the order of their two legends between runs (ggplot2 orders them arbitrarily); the order is now fixed to the one in the reviewed figures. Label placement in Figs 4 and S3 now ends on an iteration count, not on ggrepel's half-second clock.
- **Table S5.** 30 resident rows were the five floor and ceiling fits repeated under each property's row set (all row sets are the ladder rows since 2 Nov); dropped, so Table S4 has one row per model and pool. The slope structure is labelled the same way for both pools.
- **Reports.** The 07a fit table showed 0 species for every experiment (`nlevels` on a character column). Steps 05 and 06 said the P and TP slope variances sit at the boundary; only P does (TP is of the size of T). Steps 06 and 07d dropped metafor's 'struct disregarded' warning from their reports, against CLAUDE.md; it is now listed with its count. Old figure numbers (Fig 4, S4, S7, S9, S11, S2) replaced by the final ones in headers, comments and report text. Step 03's report said step 10 draws Fig 1b; step 07a does.
- **README.** The forest sentence compared medians over 19 levels including the Null (the forest's Null is the training mean) and disagreed with the 06b report; it now uses the 06b definition over the 18 model levels. A sentence now states that the colonizer Null is a weak reference (see Statistics). The Fig S13 caption, the ranger snippet and the repository tree were corrected.
- **Prep text.** P2 said 'ErrorRisk >= 4' (the pipeline uses 3, written now from the counts table); P4 called experiments 'destination sites'; the `clean_try()` default threshold and the `gbif_niche.R` header count now match the pipeline. Rerunning the pipeline changed only the reports.
- **Dead code.** `R/yardstick.R` (the shared R² yardstick, replaced 2 Nov) removed.
- **`legacy/` retired.** Nothing read from it: the clean clone ran without it. Every comment citing `legacy/<file>:<lines>` now points to commit 7fed5b1, the last commit holding it; the old pipeline at its original paths is the tag `targets-pipeline-2026-10`. The four provenance files `data/README.md` and `data/source_manifest.csv` refer to moved to `prep/provenance/`. The `.gitignore` lines and the README tree entry removed. 67 git-ignored legacy figures remain on disk (see Decisions).

## Decisions for Bills (would change a result)

Bills' answers (3 Nov) were carried out in step 13: 1 and 4 kept as they are; 2, 3, 5, 6, 7 and 8 fixed; 10: renv activated, the data README and manifest corrected; and the colonizer destination abundance put on the residents' log scale. The items below describe the state the review found; step 13's notes give the effect of each fix.

1. **Colonizer Step C with measured absences (kept, Bills 3 Nov).** Effects with strength 0 and SE 0 enter `rma.mv` with V = 0. The likelihood then has a spike at τ² = 0 (table below), an artefact of declaring a value exactly known; the reported fits are the proper interior optimum, the one any small SE floor gives too (checked 3 Nov). What remains is that a precisely measured absence outweighs a noisy presence: 6 of 28 modelled colonizers have two or more absences, and 3 species (Euphrasia officinalis, Pilosella piloselloides, Rumex acetosa) have pooled means below 1e-6 although they arrived somewhere; Fig 5 counts them as colonizers (28; 25 with a tolerance), which matches 'present once is present'.

| species | effects with SE 0 | reported tau2 | REML log-lik at reported tau2 | REML log-lik at the spike, tau2 = 1e-12 |
| --- | --- | --- | --- | --- |
| Lotus corniculatus | 2 | 0.000313 | 5.05 | 9.91 |
| Phleum pratense | 2 | 0.00237 | 5.98 | 13.27 |

2. **CH_Calanda2 colonizer years.** Step 01 dropped a year when one of the two compared groups was not surveyed, which made sense for the destination-control contrast but not for the warmed-only colonizer model (2 Nov): the colonizers of CH_Calanda2::Cal::Nes lost year 4. Fixed in step 13: the colonizer rows there now have years 1, 2, 3, 4, as the residents' (1, 2, 3, 4).
3. **TRY names (fixed in step 13).** Records reached a species through TRY's accepted name even when their own name had been resolved and rejected (558 records), moving some records between species (step 15 finding 3). Now TRY's name is used only when the own name was never resolved or resolves to the same species: 1,665 kept records come in this way, 80 of them with an own name that resolves to the same species (`prep/R/traits.R` `attach_try_names()`).
4. **Climate years.** Experiment climate averages every calendar year from the first to the last survey, also years without a survey: 17 of 42 analysed experiments have such gaps (CN_Damxung::HIGH::LOW, CN_Heibei::3400::3200, CN_Heibei::3800::3200, CN_Heibei::3800::3400, ...). PLAN says 'averaged over the experiment's years', which reads either way.
5. **A clone could not run P6 (fixed in step 13).** `data/derived/prep/gbif_niche.rds`, a pipeline input, was git-ignored; it is now tracked (the csv agrees only to csv precision).
6. **Warnings in parallel workers were lost (fixed in step 13).** `parallel::mclapply` workers muffled warnings into their own copy of `warns`; the forest workers now return their warnings with their results.
7. **R² trimming in medians (fixed in step 13).** Taking out an R² below −100 moves the colonizer Woody held-out median from -0.55 to -0.28; 07e now takes the same gradient out of the forest's median too, so both are over the same gradients.
8. **Caches keyed by file existence (fixed in step 13).** Every cache now stores a key of its inputs (`R/cache.R`) and is redone when an upstream change alters them; every fit was redone once under the keys and reproduced every result.
9. **Smaller points.** Fig 2's panels b and c show plain mean cover, which differs in sign from the model strength for a few effects; 30 morphotaxon codes (`*_NID_*`) are labelled genus-level instead of unidentified (same treatment, different label); the IN_Kashmir exclusion is a region name, not a stated rule; Step C reports name the first optimizer's failure cause without reading the error.
10. **Environment and files (step 13).** renv is activated, its project library linked from the installed packages (versions as in the lock; ragg declared in `_dependencies.R` so figures keep their device); 61 locked packages that no code uses are not installed, nine of them used only by the old pipeline (clubSandwich, crew, geometry, ggalluvial, here, mclogit, tarchetypes, testthat, TNRS): a `renv::snapshot()` would drop them. Corrected `data/README.md` and `data/source_manifest.csv` are ready to copy in (a permission rule protects `data/`). CLAUDE.md still states rules about `legacy/` and `figures/diagnostics/`. Remove the leftover ignored files with `rm -rf legacy`.

## Open items

| item | status |
| --- | --- |
| Response scale (PLAN, open since 30 Oct) | still open: in Step C the top 10% of effects carry 84% of the weight 1/SE²; the ladder's effective weights are even (step 05) |
| GBIF download citation (deadline 28 Apr 2027) | open: DOI in README; goes into methods and data statement |
| [P1] reference in the README motivation | open: needs the network reference |
| Step 15, upstream-correction ledger | open (planned): finding 3 still live (decision 3 above) and finding 10 (Kashmir rule) still live; the others are marked fixed or decided in the ledger, not re-checked one by one here |
| Steps 10 and 11 | awaiting review |
| GEN slope check (additional_checks/05_gen_check.R) | obsolete: DIAG taken ladder-wide on 1 Nov |
| Fig S8 surfaces, Fig 5 per-gradient, results-dependent Fig 6, family rerun | closed: dropped 2 Nov |
| Colonizer A0_dest untransformed; SE 0 for absences; colonizer W-slope variance at the boundary | open (step 11); SE 0 is now decision 1 |

## Statistics: weak points not discussed before (reported, not changed)

- **Shared control plots across experiments.** 16 of 42 experiments share their origin controls with another experiment of the gradient; 62 species × origin-site groups hold 160 of the 512 ladder effects, whose strengths subtract the same control cover. Their sampling errors are correlated, V is diagonal and the ladder has no gradient term, so the fixed-effect intervals (Figs 7, S12, S16) are too narrow. Held-out prediction is not affected (whole gradients are held out).
- **The origin abundance shares data with the response.** A0_origin comes from the origin controls whose later cover is subtracted in strength: its correlation with the origin-control cover from year 3 is 0.71, and with strength -0.24. Part of the strong negative A0_origin effect is regression to the mean, on top of 'cannot lose more than you have'.
- **Standardisation includes the held-out gradient.** W, T and P are z-scored on all 42 experiments; leaving one gradient out would move their means by at most 0.15, 0.27 and 0.21 SD. Small, but the held-out scores use information from the held-out gradient.
- **Fig 4 classes are signs of noisy estimates.** 412 of 1,004 resident effects have |strength| below one SE; the loser and winner slices partly sort noise.
- **Consistent classes rest on τ² = 0.** 7 of 8 consistent residents have τ²_s at its boundary (14 of 74 species overall). With a t quantile (df = experiments − 1) in place of 1.96, 4 of the 8 consistent species become context-dependent.

| window | effects | median SE | mean |strength| |
| --- | --- | --- | --- |
| 1 | 209 | 0.011 | 0.028 |
| 2 | 402 | 0.0086 | 0.029 |
| 3 | 36 | 0.035 | 0.038 |
| 4 or more | 357 | 0.01 | 0.02 |

- **Unequal averaging windows** (table above): experiments contribute different numbers of years from year 3, so precision and the size of strength differ with the window.
- **Held-out medians are fragile.** Residents, median (interquartile range; gradients above 0): initial abundances 0.17 (-0.01 to 0.26; above 0 in 12 of 16); Climate + start 0.16 (0.04 to 0.25; above 0 in 12 of 16); Species-dependent 0.12 (-0.04 to 0.25; above 0 in 12 of 16). Colonizers on 8 gradients: destination abundance 0.22 (0.15 to 0.31; above 0 in 7 of 8).
- **The colonizer Null is a biased reference.** The forest's plain training mean reaches a median R² of 0.19 against the metafor Null on the colonizer gradients: the measured absences (strength and SE 0) pull the weighted intercept towards zero, so every colonizer R² against the Null is measured from a low baseline.
- **Colonizer climate scale.** The colonizer ladder carries W, T, P on the residents' 42-experiment scale; on its 110 rows T has mean 0.16 and SD 0.72; the Woody colonizer model rests on one woody species.

## Known and not reopened

Gaussian model on bounded cover (23 Oct); compositional covariance within an experiment (24 Oct); diagonal T, P, TP slopes (1 Nov); no uncertainty on τ_s (25 Oct); forest tuning by out-of-bag R² on all rows (2 Nov); multiple testing in Fig 7 (raised 2 Nov); colonizer A0_dest scale and the W-slope boundary (step 11).

## Warnings captured

- There are outcomes with non-positive sampling variances.
- 'V' appears to be not positive definite.

## Deviations from PLAN.md

None: this step changes no analysis choice.
