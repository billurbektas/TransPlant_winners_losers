# TransPlant warming winners and losers

This repository analyzes plant-community responses to downward transplantation along elevation gradients. A turf moved from a colder origin to a warmer destination is compared with an unmoved control at its origin. Species are classified from relative-cover trajectories by pool, contrast direction, and statistical significance. Significant destination-control colonizers and significant warmed-only outsiders are combined into the biological `colonizer` response. The pipeline asks whether responses persist through time, vary with destination climate, and are associated with functional groups or plant traits.

The reproducible cleaning, climate, winner/loser, temporal, and descriptive pipeline is defined in [`_targets.R`](_targets.R); its implementation is split by responsibility under [`R/pipeline/`](R/pipeline/). Modeling functions and standalone trial scripts are preserved, but model fitting is intentionally excluded from ordinary `targets::tar_make()` runs.

## Current checkpoint: where the project is left

The canonical data, TRY summaries, monthly CHELSA summaries, winner/loser classifications, temporal summaries, climate assignments, and descriptive figures are retained and ready for continued work. Generated model objects, model tables, exploratory model outputs, and superseded model archives were deliberately removed because no final inferential structure has been approved. They remain reproducible from the preserved scripts and `output/tables/analysis_dataset.csv`.

The only retained modeling figure is `plot/climate_response_data_model_probability_prototype.pdf`. It is a visual checkpoint for the pool-specific fixed-effects candidate, not an endorsement of that model as the final analysis. All descriptive, temporal, alluvial, legacy, and presentation figures remain available.

The legacy category-model trials remain review artifacts rather than final
inference. The separate resident continuous-AUC analysis below now provides the
approved dependence structure for resident responses; nonresident occurrence
and magnitude inference remains deferred pending a reviewed opportunity-pool
definition.

## Resident continuous-AUC hierarchical analysis

The standalone resident analysis in `R/13_frequentist_hierarchical_auc.R`
models the continuous raw warmed-minus-origin AUC contrast and its sampling
covariance. It is deliberately outside `_targets.R` and does not alter the
winner/loser/colonizer classifications. Maximum-duration AUC is primary and
the common first-three-year horizon is a sensitivity analysis. The model
prioritizes independent species random intercepts and warming slopes, uses
experiment-clustered robust inference, runs the prespecified fallback ladder,
and fits LDMC, LMA, and vegetative height separately. CR2 is attempted first;
for the realized crossed species/experiment hierarchy, the analysis records
that limitation and uses the small-sample-adjusted CR1 fallback.

After restoring the lockfile, run:

```r
source("R/13_frequentist_hierarchical_auc.R")
```

Outputs, including the audit, model ledger, predictions, diagnostics, figures,
and interpretive reports, are written to
`output/frequentist_hierarchical_auc/`. Set
`TRANSPLANT_AUC_SKIP_LOEO=true` for a faster smoke run without the
leave-one-experiment-out refits. The self-contained pre-model audit is written
before covariance reconstruction to `audit/data_audit.html`; the final visual
summary is `visual_report.html`. Set
`TRANSPLANT_AUC_SKIP_VISUAL_REPORT=true` to skip only the final visual render.
Both HTML files can be regenerated from saved outputs without model refits:

```r
source("R/14_render_frequentist_hierarchical_auc_report.R")
```

Colonizer inference is intentionally deferred.

## Current input and the adapter boundary

The current vegetation input is `data/TransPlantNetwork_060824.RData`. It remains the source until the revised metadata arrives. It contains:

| Object | Current form | Role |
|---|---:|---|
| `dat` | 61,430 long-format rows | Species cover and relative cover by region, site pair, treatment, plot, and year |
| `cleanedspecies` | 4,703 rows | Site-specific mapping from local species codes to matched botanical names |
| 22 regional list objects | `meta`, `community`, `cover`, `taxa`, and sometimes `trait` | Coordinates, elevations, experiment dates, plot sizes, and source tables |

The species data currently used are therefore the species observations in `dat`, translated from `SpeciesName` local codes with the region-and-destination-specific `cleanedspecies` table. Existing botanical corrections are applied before the canonical species list is sent to TNRS. The incoming data can keep different metadata internals as long as a replacement adapter produces these three tables:

- `sites`: one row per `region × site_id`, with longitude, latitude, elevation, establishment/survey years, plot size, and country.
- `community`: long-format observations with region, calendar and relative year, origin, destination, treatment, stable plot ID, plot-year ID, species, cover, and relative cover.
- `experiments`: one row per warmed origin-to-destination comparison, including control availability, analysis years, plot size, and exclusions.

No analysis module reads the regional objects directly after this adapter. A future input structure should only require a new function parallel to `read_current_transplant_input()` and the canonical constructors in [`community.R`](R/pipeline/community.R).

## Experiment and control structure

- **Warmed:** `origin_site != destination_site` and `treatment == "Warm"`.
- **Origin control:** `LocalControl` with origin and destination both equal to the source site.
- **Destination control:** `LocalControl` with origin and destination both equal to the destination site.
- **Stable plot ID:** region, origin, destination, treatment, and destination plot/turf. It deliberately excludes year.
- **Plot-year ID:** stable plot ID plus calendar year.

The adapter validates unique site metadata and both controls for every warmed pair. The current snapshot has 59 sites and 44 raw warmed pairs. Two are explicitly excluded:

- `IN_Kashmir`: insufficient duration.
- `FR_AlpeHuez`: warmed rows have missing `Year`.

This leaves 42 analysis comparisons. The nine ambiguous site/code keys in `SE_Abisko` are recorded in [`species_mapping_decisions.csv`](config/species_mapping_decisions.csv). Eight are accepted-name choices and `Equ pal/pra` is excluded as an unresolved aggregate. The complete candidate rows are written to `output/tables/species_mapping_review.csv` before analysis.

## Setup and execution

Use the explicit refresh entry point once, and whenever dependency updates are intentionally requested:

```r
source("R/00_refresh_dependencies.R")
```

It fetches the current tidyTRY version, archives the current `get_chelsa.R`, records both Git SHAs in [`external_sources.csv`](config/external_sources.csv), installs the required packages, and updates `renv.lock`. This is the only entry point that updates GitHub dependencies. The SHAs recorded when this revision was prepared are tidyTRY `3352199f…` and get_chelsa `3d410b8`.

Ordinary analysis uses the lockfile and does not update either repository:

```r
renv::restore()
targets::tar_make()
```

or:

```r
source("R/00_run_pipeline.R")
```

CHELSA branches use two local workers by default. Change the limit explicitly before starting R:

```sh
TRANSPLANT_CHELSA_WORKERS=4 Rscript R/00_run_pipeline.R
```

Useful staged runs are:

```r
targets::tar_make(names = c(sites_file, experiments_file, canonical_validation_file))
targets::tar_make(names = try_summary_file)
targets::tar_make(names = experiment_climate_file)
targets::tar_make(names = winner_loser_file)
```

TNRS, TRY cleaning, and CHELSA require network access. The 8.1-GB raw TRY file, monthly CHELSA cache, `_targets/` store, local `renv` library, and reconstructable model/intermediate objects are ignored by Git.

## TRY cleaning

The raw export is `data/try_47163.txt`. The project does not use `default_trait_map()`. [`try_trait_map.csv`](config/try_trait_map.csv) is the reviewed dictionary and is intentionally limited to:

- LDMC: two supplied leaf-dry-matter aliases.
- SLA: seven supplied SLA/1-LMA aliases covering petiole included, excluded, undefined, and rachis/midrib-excluded forms.
- Vegetative height: `Plant height vegetative` and `Plant vegetative height`.

The sequence in [`try_traits.R`](R/pipeline/try_traits.R) is:

1. Take canonical vegetation names after botanical corrections; preserve original and normalized matching fields.
2. TNRS-resolve the vegetation list and save the complete response.
3. Filter the raw export in chunks with original plus accepted names.
4. Remove experimental TRY datasets and TRY duplicates.
5. Clean quantitative records with `ErrorRisk < 4` through tidyTRY.
6. Apply the project dictionary and explicitly discard every trait except LDMC, SLA, and vegetative height.
7. Extract TRY observation coordinates.
8. Independently TNRS-resolve unique TRY-side `SpeciesName` and `AccSpeciesName` values.
9. Stop if a name is unresolved or has TNRS score below 0.90 unless it is accepted/excluded in [`taxonomy_overrides.csv`](config/taxonomy_overrides.csv).
10. Assign Köppen-Geiger zones and remove `Af`, `Am`, `Aw`, `BWh`, `BWk`, `BSh`, `BSk`, and missing assignments.
11. Confirm one standardized unit per trait and positive finite SLA, then calculate `LMA = 1 / SLA`.
12. Calculate species medians without imputation, retaining observation count, dataset count/list, unit, range, and climate-zone coverage.

TraitID 335, flowering phenology, Floraveg indicators, and dispersal downloads are deliberately out of scope. Filter losses are written to `try_filter_audit.csv`.

## Resumable CHELSA extraction

CHELSA uses the 59 TransPlant site coordinates, never TRY coordinates. The pipeline uses monthly CHELSA V2.1 `tas` and bias-corrected `pr`. This is an explicit data-availability decision: the daily bias-corrected `pr` archive has no files for 2020–2023, whereas the monthly archive extends through 2021. The project does not substitute or splice the separate uncorrected daily `prec` product. Monthly rasters are extracted for 2008–2021, the complete portion of the archive intersecting current experiment years. The URL and unit logic derive from `get_chelsa`; [`chelsa.R`](R/pipeline/chelsa.R) adds the interruption-safe behavior needed here.

Four experiments extend beyond monthly CHELSA coverage: `FR_Lautaret::G::L` requires 2022 and the three `US_Colorado` comparisons require 2022–2023. Their climate summaries use the available years only and are explicitly marked `climate_summary_status = "partial"` and `climate_complete = FALSE`; expected/available counts and missing years are retained in both climate tables. Lautaret therefore uses 2017–2021 (5/6 required years), and Colorado uses 2018–2021 (4/6 required years). The other 38 experiments have complete monthly coverage. This decision can be revisited if CHELSA extends the monthly archive or a different complete climate product is selected.

Before any raster extraction, the pipeline validates unique site keys, missing values, longitude/latitude bounds, exact coordinate duplicates, reuse across experiments, and origin/destination separation. Exact coordinate sharing is retained because it can represent the same physical site under different regional metadata; the reviewed groups are saved in `chelsa_duplicate_coordinates.csv` and all multi-experiment reuse in `chelsa_coordinate_reuse_across_experiments.csv`. Identical origin and destination coordinates within one warmed experiment stop the pipeline.

Work is split into year-month-variable branches. Each remote monthly raster is opened once and all site points are extracted together. Every month-variable result is checkpointed at:

```text
data/cache/chelsa_monthly/YYYY/MM/{variable}_YYYY-MM.rds
```

Each checkpoint is written to a temporary file, validated for schema, 59 unique sites, year, month, variable, unit, and finite values, and then atomically renamed. Restarting skips valid files. Invalid files are quarantined with an `.invalid-<timestamp>` suffix and re-extracted. Downloads use bounded exponential-backoff retries. If a required 2008–2021 raster still fails, the branch writes a missing-URL report under `data/cache/missing_monthly/` and stops.

For every site and year the pipeline calculates:

| Metric | Temperature | Precipitation |
|---|---|---|
| Annual | Days-weighted mean of monthly `tas` | Sum of monthly `pr` accumulations |
| JJA | Days-weighted mean of June–August monthly `tas` | Sum of June–August monthly `pr` accumulations |
| July | July monthly `tas` mean | July monthly `pr` accumulation |

For an experiment, the intended interval contains every calendar year from the first through last analyzed post-transplant year, including years without a vegetation survey. Available yearly values are averaged for origin and destination sites, then destination-minus-origin differences are calculated. If any required year is unavailable, the available-year average is retained but marked `partial`, with the expected/available counts and missing years recorded. Primary predictors are `diff_T_JJA`, `destination_T_JJA`, and `destination_P_JJA`; annual and July variables remain secondary outputs.

## Winner/loser calculation

[`winner_loser.R`](R/pipeline/winner_loser.R) retains the relative-cover AUC approach. Trapezoidal AUC is calculated per plot for raw and square-root relative cover at every supported annual horizon from year 2 through year 9, plus maximum duration. Within each experiment/horizon/transform, `auc ~ ODT * species` is fitted and warmed-minus-origin contrasts use multivariate-*t* adjustment. Each model group receives a stable key-derived random seed so randomized multivariate-*t* integration is invariant to the number or ordering of other horizons. Significance is `p < 0.05`.

**Provisional downstream-horizon decision:** the repeated-species, functional-group, trait, genus, destination-climate, and explanatory-model analyses currently use raw relative-cover results at `horizon == "max"`. Thus, every experiment contributes its complete available post-establishment survey interval rather than a common fixed horizon. The present 42 experiments span 3 years (9 experiments), 4 years (14), 5 years (2), 6 years (6), 7 years (8), or 9 years (3). Including `year_range` as a model covariate adjusts for duration statistically but does not make these response windows equivalent. This choice is being retained for now and may later be changed to a common horizon or a horizon-stratified analysis. The explicitly temporal figures and tests remain based on their stated annual horizons rather than `max`.

Before pool assignment, taxa explicitly excluded at the vegetation TNRS review gate (including `_NID_` morphotaxa) are removed, and the legacy eligibility rule is retained: a species×treatment group must contain at least two vegetation observations (`n >= 2`, equivalently `n > 1`). In the current data all 5,265 treatment×species groups already satisfy this rule. Increasing the threshold to `n >= 3` would remove 709 groups (13.5%); because this count combines plots and years, that stricter threshold is reserved for sensitivity analysis rather than used as the primary filter. Matching the original scripts, each species is assigned one fixed pool per experiment from observations across the complete analyzed duration, after the observation filter. The same pool is then joined to every annual horizon and the maximum-duration contrast. Validation stops if pool identity changes across horizons or if a colonizer/outsider is assigned a resident response category. Categories are:

Only vegetation records identified to species level or lower are admitted to the analysis community. One-word genus names, explicit `Genus sp./spp.` names, family-only names, and unresolved names are omitted before species pools and AUCs are calculated. The current resolution review omits 64 genus-only and four family-only names; the complete decisions are saved in `vegetation_taxonomic_resolution.csv` and the 68-name omission list in `omitted_non_species_taxa.csv`.

| Category | Definition |
|---|---|
| `loser` | Present in origin controls; adjusted `p < 0.05` and negative warmed-minus-origin contrast |
| `non_significant_loser` | Present in origin controls; adjusted `p >= 0.05` and negative contrast |
| `no_change` | Present in origin controls with a contrast exactly equal to zero; there are currently no such records |
| `non_significant_winner` | Present in origin controls; adjusted `p >= 0.05` and positive contrast |
| `winner` | Present in origin controls; adjusted `p < 0.05` and positive contrast |
| `colonizer` | Adjusted `p < 0.05` and positive contrast, and either (a) absent from origin controls but present in destination controls and warmed plots or (b) absent from both controls but present in warmed plots |
| `non_significant_colonizer` | Absent from origin controls, present in destination controls and warmed plots, and adjusted `p >= 0.05` |
| `non_significant_outsider` | Absent from both controls, present in warmed plots, and adjusted `p >= 0.05` |

Significant outsiders are therefore treated as colonizers. Non-significant outsiders remain explicitly recorded in the calculation-level `category` column and can always be recovered using `pool == "outsider"`. For downstream analyses, the separate `analysis_category` column merges `non_significant_outsider` into `non_significant_colonizer`; all other labels remain unchanged. This merge is therefore reversible without recalculating AUCs or contrasts. There are no extinction or colonizing-loser categories. A significant negative colonizer/outsider is treated as an invalid state and stops validation; accumulated origin AUC is zero for these pools, so such a result is not expected.

## Questions and tests

### Temporal consistency

The continuous test regresses strength at years 4 and 6 on first-two-year strength with crossed species and experiment random intercepts. Annotated transition matrices show exact persistence on the diagonal and every significance or response-direction change away from it. Pool changes are forbidden by validation rather than displayed as transitions. All temporal summaries use `analysis_category`, so destination-control colonizers and warmed-only outsiders are displayed together as `colonizer` or `non_significant_colonizer`; their underlying pool and calculation category remain recoverable from `winner_loser_results.csv`. The supplementary alluvial figure retains the First 2, First 4, and First 6 stacked bars and their original vertical proportions. Missing follow-up still occupies its required layout space so those proportions are not distorted, but that space, its label, and every ribbon entering it are fully transparent. Consequently, visible paths may end rather than being sent to a displayed unavailable state. Only loser↔winner response reversals are annotated, including transitions involving either significant or non-significant states. Each annotation is formatted as `count (source-conditional percentage)`. Bars contain counts only, with category identity supplied by the retained legend. Exact sequences, pairwise transition counts/proportions, Cohen's kappa, experiment-cluster bootstrap intervals, and Holm-adjusted *p*-values remain saved as tables.

### Climate consistency

Destination climates use five reviewed groups: cold/dry (`P < 550 mm`, `T < 12.8669 °C`), intermediate/dry (`P < 550`, `12.8669 ≤ T < 16.9427`), warm/dry (`P < 550`, `T ≥ 16.9427`), cold/intermediate (`550 ≤ P < 768`, `T < 16.9427`), and cold/wet (`P ≥ 768`, `T < 16.9427`). Their current experiment counts are 15, 9, 11, 2, and 5. A future warm, non-dry experiment is not forced into a misleading group: it is written to `climate_context_review.csv`, and only context-dependent targets stop pending review.

The retained signed species-consistency figure includes species occurring in more than three experiments, with significant winners/colonizers above zero, significant losers below, and all non-significant subcategories centered in gray. The climate-density panels use raw significant response counts: each experiment/category kernel and count-sized point is weighted by the number of species×experiment responses, not by a within-experiment proportion. All 42 experiments remain visible as faint background points, and the nonzero count points are filled by JJA warming. These densities describe where observations are concentrated; they are not response probabilities. Formal response probabilities instead come from the pool-specific models described below. The annotated warming-gradient map is descriptive and divides each observed destination JJA climate range into three equal-width intervals: temperature cutoffs are 12.8669 and 16.9427 °C, and precipitation cutoffs are 462.75 and 768.00 mm. The reviewed display assignments retain Cold/Intermediate, assign the NO Gudmedalen Rambera-to-Arhelleren comparison to Cold/Intermediate, the higher-precipitation NO Skjellingahaugen comparison to Cold/Wet, and the intermediate-precipitation DE Grainau and DE Susalps comparisons to Warm/Dry. The resulting occupied groups contain 14 Cold/Dry, 9 Intermediate/Dry, 11 Warm/Dry, 3 Cold/Intermediate, and 5 Cold/Wet experiments; the Warm/Wet cell remains explicitly unoccupied. These figure-specific assignments do not replace the five reviewed groups used by the context analyses and trait panels.

### Functional groups and traits

TNRS family and genus are used to derive two independent, binary functional statuses: `graminoid_status` (`graminoid`/`non_graminoid`) and `woody_status` (`woody`/`non_woody`). Poaceae, Cyperaceae, and Juncaceae define graminoids; a reviewed woody-genus list defines woody taxa. These statuses are tested in separate functional-information models and are not collapsed into one mutually exclusive forb/graminoid variable. The old `functional_group` field remains only as a compatibility field for preserved descriptive scripts. [`functional_group_overrides.csv`](config/functional_group_overrides.csv) is the explicit correction interface. The trait figure retains LDMC, LMA, and vegetative-height facets. Within every response category the five destination-climate groups are shown as side-by-side boxplots, with box and point colors identifying climate group. Trait values are log-transformed.

Six non-angiosperm taxa remain in the upstream community and winner/loser calculations but must be excluded before every downstream ecological model and phylogenetic analysis: *Gymnocarpium dryopteris*, *Athyrium spinulosum*, *Botrychium lunaria*, *Equisetum arvense*, *Selaginella selaginoides*, and *Juniperus communis*. The authoritative exclusion list and filter are in [`analysis.R`](config/analysis.R). Thus the preserved winner/loser table retains their calculated contrasts, while future downstream resident model runs use 526 rather than 532 species.

### Modeling checkpoint: trials completed, no final model selected

Modeling is now manual and exploratory. The complete model and plotting functions remain in [`climate_response_model.R`](R/pipeline/climate_response_model.R), [`biological_context_review.R`](R/pipeline/biological_context_review.R), and [`hierarchical_trait_trial.R`](R/pipeline/hierarchical_trait_trial.R). Standalone entry points remain under `R/`, but generated model objects, model tables, diagnostics, and review figures are intentionally not retained.

The upstream maximum-duration raw-cover snapshot contains 1,625 resident species×experiment records (132 losers, 69 winners, and 1,424 non-significant residents) and 614 nonresident records (34 significant colonizers and 580 non-significant nonresidents). The six declared downstream exclusions remain in those winner/loser calculations. The updated resident modeling boundary removes their 21 records, leaving 1,604 effects from 526 species across the same 42 experiments. Destination-control colonizers and warmed-only outsiders are combined for downstream colonizer analysis while their original calculation pool remains recoverable in the analysis dataset. A species can belong to different pools in different experiments but exactly one pool within a species×experiment.

The following structures have been tried:

- A fixed four-category multinomial combining residents, colonizers, and all non-significant responses converged, but it was retired because colonizers and residents do not share a biologically valid denominator.
- A fixed pool-specific candidate—resident multinomial plus nonresident colonizer binomial—converged and generated the retained probability prototype. Its ordinary Wald inference does not account for repeated species or observations sharing experiments, so it is not yet the approved final model.
- Separate loser, winner, and colonizer binomial sensitivity models with a species random intercept fitted for loser and winner; the colonizer species variance collapsed to zero and required a fixed-effects fallback.
- A simpler fixed climate-category model was explored for descriptive comparison, but it was not selected for continued inference.
- Family, functional-group, and one-trait-at-a-time mixed trials with species and experiment intercepts were only partly estimable. Resident multinomials timed out or failed, several binomials reached zero experiment variance, and sparse family/colonizer cells produced separation.
- The paper-inspired LDMC hierarchy used separate species intercept, species temperature, precipitation, and warming slopes plus experiment and regional-gradient intercepts. The resident PQL multinomial exceeded its feasibility limit; loser/winner binomials were singular and the colonizer fit had poor support and extreme variance estimates.
- Removing experiment and gradient effects allowed clean full-fixed loser and winner binomials with a species intercept. Full species random-slope versions remained singular; temperature-only slopes were numerically supported for both resident outcomes, precipitation-only for winners, and warming slopes reached boundaries. The sparse colonizer fits remained provisional.
- A final reduced resident `mblogit` trial with only separate species intercept, temperature, precipitation, and warming slopes produced inner singular convergence, near-zero random variances, and no usable fixed-effect covariance matrix.

None of these trials is the selected final inferential analysis. The next step is to choose one defensible base structure for dependence among repeated species and shared experiments, then use that same framework for family, functional-group, and trait extensions. Pool separation must remain explicit. The maximum-duration response window should be reviewed before final model fitting.

Run modeling only when explicitly requested. The focused pool-specific candidate and its prototype can be reconstructed with:

```sh
Rscript R/04_climate_response_model.R
```

The other preserved review scripts similarly recreate their own disposable outputs. Ordinary `targets::tar_make()` does not fit or save explanatory models.

## Outputs

Milestones are saved outside the targets cache so they remain inspectable:

| Path | Contents |
|---|---|
| `output/tables/input_inventory.csv` | Current input objects and sizes |
| `output/tables/sites.csv` | Canonical 59-site table |
| `output/tables/vegetation_taxonomic_resolution.csv` | Current vegetation species list with inferred identification level and include/omit decision |
| `output/tables/included_vegetation_species.csv` | Final species-level vegetation list admitted to the analysis |
| `output/tables/omitted_non_species_taxa.csv` | Genus-only and family-only taxa omitted before analysis |
| `output/intermediate/community.rds` | Canonical cleaned community table |
| `output/tables/experiments.csv` | Raw/included comparisons and exclusions |
| `output/tables/canonical_validation.csv` | Snapshot and pairing checks |
| `output/tables/species_mapping_review.csv` | Nine ambiguous local-code reviews |
| `output/tables/tnrs_vegetation_full.csv` | Complete vegetation TNRS response |
| `output/tables/tnrs_vegetation_review.csv` | Vegetation taxonomy gates/decisions |
| `output/tables/tnrs_try_full.csv` | Independent TRY-side TNRS response |
| `output/tables/tnrs_try_review.csv` | TRY taxonomy gates/decisions |
| `output/tables/try_filter_audit.csv` | Observation/species/dataset losses by stage |
| `output/tables/try_range_validation.csv` | Observed focal-trait ranges and plausible-range gates |
| `output/intermediate/try_traits_individual.rds` | Filtered observation-level focal traits |
| `output/tables/try_species_traits.csv` | Species medians and provenance summaries |
| `output/tables/site_year_climate.csv` | Annual/JJA/July site-year metrics |
| `output/tables/experiment_climate.csv` | Experiment-period origin/destination/difference climate |
| `output/tables/experiment_climate_coverage.csv` | Expected, available, and missing climate years plus the complete-coverage flag for every included experiment |
| `output/tables/chelsa_coordinate_validation.csv` | Pre-download coordinate checks and counts |
| `output/tables/chelsa_duplicate_coordinates.csv` | Exact coordinate groups assigned to multiple site IDs |
| `output/tables/chelsa_coordinate_reuse_across_experiments.csv` | Coordinate locations used by multiple included experiments |
| `output/intermediate/auc_by_plot.rds` | Plot-level AUC by horizon and transform |
| `output/tables/winner_loser_results.csv` | Contrasts, pools, categories, uncertainty |
| `output/tables/temporal_transitions.csv` | Flow counts and proportions |
| `output/tables/temporal_persistence.csv` | Matched flow totals and same-state persistence by horizon |
| `output/tables/temporal_state_counts.csv` | Left- and right-hand flow-box counts by state and horizon |
| `output/tables/temporal_linked_flows.csv` | Exact three-horizon category sequences and their counts |
| `output/tables/temporal_category_persistence.csv` | Same-category counts and proportions within each starting category |
| `output/tables/temporal_overall_persistence.csv` | Overall pairwise and all-three-horizon category persistence |
| `output/tables/temporal_transition_matrix.csv` | Complete starting-to-later category matrices with counts, row percentages, and transition types |
| `output/tables/temporal_category_composition.csv` | Counts, denominators, percentages, distinct species, and distinct experiments at every annual horizon from years 2 through 9; `no_change` is retained here for auditability but omitted from the temporal figure |
| `output/tables/temporal_kappa.csv` | Kappa estimates, bootstrap intervals, adjusted tests |
| `output/tables/functional_groups.csv` | Independent graminoid/non-graminoid and woody/non-woody assignments, review bases, and a deprecated compatibility group |
| `output/tables/analysis_dataset.csv` | Complete joined model/figure input |
| `output/tables/climate_context_experiments.csv` | Experiment assignments to the five reviewed destination JJA climate groups and their cutoffs |
| `output/tables/climate_context_review.csv` | Future warm, non-dry experiments requiring an explicit context decision |
| `output/tables/winner_loser_annotated.csv` | Maximum-duration winner/loser results annotated with accepted taxonomy, traits, functional groups, and experiment covariates; does not require CHELSA |
| `output/manifests/figures.csv` | Figure question, source, target, dimensions, and status |

Generated model directories and model-result tables are intentionally absent from this checkpoint. The preserved standalone modeling scripts recreate them when a trial is explicitly requested.

## Figures

The consolidated collection contains 15 PDFs directly under `plot/`. Six are generated by `_targets.R`:

- `plot/transient_winning_losing_strength.pdf`: First-2 versus later-horizon strength, with Pearson correlations and p-values.
- `plot/temporal_category_flows.pdf`: First 2 → First 4 → First 6 response flows.
- `plot/temporal_transition_matrices.pdf`: annual response composition and pairwise transition matrices.
- `plot/destination_climate_warming_gradient.pdf`: annotated experiments in destination JJA climate space.
- `plot/repeated_species_consistency.pdf`: significant-response consistency for species represented in more than three experiments. White integers inside orange, green, and purple segments give the corresponding significant loser, winner, or colonizer experiment count; gray non-significant segments are not labelled.
- `plot/compare_rel.abund_vs.sqrt.rel.abund.pdf`: raw- versus square-root-relative-cover sensitivity.

Six retained PDFs are produced by standalone scripts:

- `plot/climate_response_data_model_probability_prototype.pdf` from `R/04_climate_response_model.R`; this is the only retained modeling figure and is a provisional visual checkpoint.
- `plot/trait_categories_with_overall.pdf` from `R/07_climate_faceted_descriptive_review.R`.
- `plot/response_to_climate_alluvial.pdf`, `plot/response_by_climate_percentages.pdf`, `plot/response_by_climate_percentages_response_axis.pdf`, and `plot/response_to_family_alluvial.pdf` from `R/10_response_grouping_flows_gradients_review.R`. In the climate alluvial and both companion bar charts, percentages are response-conditional: the five climate values sum to 100% separately for the 132 losers, 69 winners, and 34 colonizers. The first bar chart places climate on the x-axis and colors responses; the alternative places response on the x-axis and colors climates. Both display the same underlying species×experiment counts. The family alluvial retains the resident/nonresident pool summary and places it immediately beside the response column before continuing through functional group to family. Response→functional-group labels retain their displayed-flow percentages. Functional-group→family labels use each complete response total and are shown for paths representing at least 5%; no additional response counts or percentages are printed to the right of the families.

The remaining root PDFs—`plot/eg.sp.1.pdf`, `plot/eg.sp.2.pdf`, and `plot/winners_losers.pdf`—are retained legacy artifacts. The eight existing presentation files remain unchanged under `plot/figures_presentation/`; `focal4_trait_phenology_warming.png` is included there even though its generating source is unavailable.

`output/manifests/figures.csv` is an explicit allowlist of these 15 root PDFs and the preserved presentation files. Modeling figures outside the retained prototype are excluded. Their generating functions and manual scripts remain available, but ordinary pipeline runs do not recreate them.

## Validation

Run local, non-network fixtures with:

```sh
Rscript tests/testthat.R
```

The suite covers the current 59-site/44-raw/42-included snapshot, stable plot IDs, controls, the nine mapping decisions, TRY alias coverage and LMA reciprocity, coordinate and monthly-checkpoint gates, days-weighted annual/JJA/July calculations, explicitly flagged partial experiment averages, fixed experiment-wide pools, outsider logic, flow conservation, the five destination-climate boundaries and review gate, exact resident/nonresident denominators, experiment-weighted scaling, model formulas and historical convergence decisions, deterministic probability calculations, denominator-selected review families, distinct-species trait scaling, and zero-event diagnostics. Pipeline targets additionally validate complete monthly coverage for 2008–2021, taxonomy thresholds, TRY units/climate filters, winner/loser key uniqueness, pool/category compatibility, and legacy AUC agreement within numerical tolerance.

The long-running TNRS, 8.1-GB TRY scan, full CHELSA extraction, and explanatory model refits are intentionally excluded from routine local execution. Model unit tests may use local fixtures, but `targets::tar_make()` does not fit explanatory models or write model artifacts.
