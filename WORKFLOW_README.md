# Canonical targets workflow

The pipeline is organized into three layers: configuration, reusable functions,
and target definitions. `_targets.R` connects those layers and returns the
complete dependency graph; it does not directly fit models itself.

## 1. Starting the pipeline

From the repository root, run:

```bash
Rscript scripts/run_pipeline.R
```

The entry point is [`scripts/run_pipeline.R`](scripts/run_pipeline.R). It
eventually calls:

```r
targets::tar_make()
```

`tar_make()` automatically finds and evaluates [`_targets.R`](_targets.R).

You can also build one target and all its prerequisites:

```bash
Rscript scripts/run_pipeline.R phylogeny_model_sequence_526
```

This builds the phylogenetic model branch, but not unrelated downstream reports
or figures.

## 2. Loading the targets framework

```r
library(targets)
```

This provides functions such as:

- `tar_target()` to define a computational node.
- `tar_option_set()` to set pipeline-wide options.
- `tar_make()` to execute outdated targets.
- `tar_read()` to retrieve completed objects.
- `tar_visnetwork()` to display the dependency graph.

At this point, no scientific analysis has run.

## 3. Loading configuration

```r
source("config/analysis.R")
source("config/canonical_analysis.R")
```

[`config/analysis.R`](config/analysis.R) defines general input and runtime
configuration:

- raw TransPlant input path;
- TRY input path;
- CHELSA cache path;
- output locations;
- climate years and variables;
- random seed;
- the six downstream exclusions.

[`config/canonical_analysis.R`](config/canonical_analysis.R) defines the settled
canonical analysis:

- expected 1,604 effects, 526 species, and 42 experiments;
- `T_z * P_z * W_z` climate formula;
- anchored species random-effects structure;
- classification rules;
- trait scaling;
- PEV axes and orientation;
- reference and canonical output paths.

These files define functions returning configuration lists. They do not load
data or fit models.

## 4. Loading all analytical functions

```r
active_pipeline_files <- canonical_pipeline_source_files()
```

This returns an explicit list of the six foundational upstream modules:

```text
R/functions/upstream_utils.R
R/functions/upstream_community.R
R/functions/upstream_climate.R
R/functions/upstream_traits.R
R/functions/upstream_winner_loser.R
R/functions/upstream_temporal.R
```

The existence check causes an immediate, clear failure if one is missing:

```r
if (any(!file.exists(active_pipeline_files))) {
  stop("One or more canonical upstream pipeline modules are missing.")
}
```

Then they are loaded:

```r
purrr::walk(active_pipeline_files, source)
```

The remaining files under `R/functions/` are found automatically:

```r
canonical_function_files <- setdiff(
  list.files("R/functions", pattern = "[.]R$", full.names = TRUE),
  active_pipeline_files
)
purrr::walk(canonical_function_files, source)
```

`setdiff()` prevents the six upstream modules from being sourced twice.

Important modules include:

- [`model_core.R`](R/functions/model_core.R): `M_both`, CR1, variance profiles,
  and fit summaries.
- [`species_context.R`](R/functions/species_context.R): species min/max
  classification and callouts.
- [`family.R`](R/functions/family.R): family aggregation, covariance, model,
  and classification.
- [`predictor_models.R`](R/functions/predictor_models.R): trait, functional,
  and PEV model sequences.
- [`phylogeny.R`](R/functions/phylogeny.R): tree reconstruction, PCoA,
  orientation, and equivalence.
- [`heterogeneity.R`](R/functions/heterogeneity.R): heterogeneity components
  and synthesis.
- [`reference_validation.R`](R/functions/reference_validation.R): numerical
  regression and release checks.
- [`figures.R`](R/functions/figures.R): canonical downstream figures.

These are primarily function definitions. The actual calls happen in targets.

## 5. Loading the target factories

```r
source("R/targets/data_targets.R")
source("R/targets/analysis_targets.R")
source("R/targets/publication_targets.R")
```

Each file defines one function that returns a list of targets:

| Factory | Targets | Responsibility |
|---|---:|---|
| `canonical_data_targets()` | 114 | Community, taxonomy, TRY, CHELSA, AUC, winner/loser, temporal data |
| `canonical_analysis_targets()` | 234 | 526-species models, phylogeny, heterogeneity, validation, exports |
| `canonical_publication_targets()` | 90 | Figures, reports, manifests, final release gate |

The files are:

- [`data_targets.R`](R/targets/data_targets.R)
- [`analysis_targets.R`](R/targets/analysis_targets.R)
- [`publication_targets.R`](R/targets/publication_targets.R)

Together they currently define 438 targets. These counts are obtained by
sourcing `_targets.R` and evaluating the lengths of the three target-factory
lists; update this table whenever a factory changes.

### Final main-manuscript figure map

The current numbered figures and their final file targets are:

| Figure | Output stem | PDF target | PNG target |
|---:|---|---|---|
| 1 | `fig01_winners_losers` | `figure01_winners_losers_pdf` | `figure01_winners_losers_png` |
| 2 | `fig02_repeated_species_consistency` | `figure02_repeated_species_consistency_pdf` | `figure02_repeated_species_consistency_png` |
| 3 | `fig03_temporal_category_flows` | `figure03_temporal_category_flows_file` | `figure03_temporal_category_flows_png` |
| 4 | `fig04_species_minmax_classification_526` | `figure04_species_minmax_classification_pdf` | `figure04_species_minmax_classification_png` |
| 5 | `fig05_climate_only_fixed_effect_prediction_surfaces_526` | `figure05_climate_only_fixed_effect_pdf` | `figure05_climate_only_fixed_effect_png` |
| 6 | `fig06_predictor_synthesis_support_and_heterogeneity_526` | `figure06_predictor_synthesis_pdf` | `figure06_predictor_synthesis_png` |
| 7 | `fig07_predictor_response_atlas_526` | `figure07_predictor_response_atlas_pdf` | `figure07_predictor_response_atlas_png` |

The files are written under `figures/manuscript/`. Figure-associated CSVs are
written only to the relevant `output/data/`, `output/tables/`, or
`output/audits/` directory. Approved appendices use descriptive `appendix_`
stems until their manuscript numbering is settled.

Figures 5–7 have standalone figure-only renderers:

- [`build_climate_only_fixed_effect_surfaces_526.R`](scripts/build_climate_only_fixed_effect_surfaces_526.R)
- [`build_predictor_synthesis_support_and_heterogeneity_526.R`](scripts/build_predictor_synthesis_support_and_heterogeneity_526.R)
- [`build_predictor_response_atlas_526.R`](scripts/build_predictor_response_atlas_526.R)

They read existing canonical exports and do not refit models or invoke
`targets`. In contrast, requesting a figure target through
`scripts/run_pipeline.R` lets `targets` rebuild any prerequisite it determines
is outdated.

## 6. Quarto availability

```r
if (!nzchar(Sys.which("quarto"))) {
  bundled_quarto <- "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto"
  if (file.exists(bundled_quarto)) Sys.setenv(QUARTO_PATH = bundled_quarto)
}
```

The pipeline first looks for a normal Quarto installation. If it cannot find
one, it tries the copy bundled with RStudio on macOS.

This affects report rendering only, not scientific calculations.

## 7. Parallel execution

```r
workers <- as.integer(Sys.getenv("TRANSPLANT_CHELSA_WORKERS", "2"))
```

The default is two local workers. You can change it when launching:

```bash
TRANSPLANT_CHELSA_WORKERS=4 Rscript scripts/run_pipeline.R
```

If `crew` is installed:

```r
controller <- crew::crew_controller_local(workers = workers)
```

`targets` can then run independent branches concurrently. This is particularly
useful for dynamically branched CHELSA jobs:

```text
chelsa_job_list
    └── chelsa_job [one branch per job]
          └── chelsa_monthly_file [parallel downloads/extractions]
```

The runner also allows parallel processing to be disabled:

```bash
TRANSPLANT_USE_CREW=false Rscript scripts/run_pipeline.R
```

## 8. Pipeline-wide options

```r
tar_option_set(
  packages = c(...),
  controller = controller,
  error = "stop",
  seed = 20260824L
)
```

This means:

- listed packages are available when target commands run;
- `crew` controls execution when enabled;
- the pipeline stops on an error;
- stochastic operations receive reproducible target-specific seeds.

The package list declares the runtime environment; it does not execute those
packages' analyses immediately.

## 9. Constructing the complete graph

The last expression is:

```r
c(
  canonical_data_targets(),
  canonical_analysis_targets(),
  canonical_publication_targets()
)
```

This returns the 438 target definitions to `targets`.

The order inside `c()` is mainly for organization. Actual execution order is
inferred from target names referenced inside each command.

For example:

```r
tar_target(
  community,
  canonical_community(current_input$dat, species_mapping, site_metadata)
)
```

From that expression, `targets` infers:

```text
current_input ─┐
species_mapping ├──> community
site_metadata ─┘
```

If `species_mapping` changes, `community` becomes outdated, followed by any
downstream target depending on `community`.

## 10. Main scientific dependency flow

A simplified version of the active graph is:

```text
Raw TransPlant input
    └── current_input
          ├── species_mapping
          ├── community
          ├── sites
          └── experiments
                 │
                 ├── taxonomy ──> analysis_community
                 │                    │
                 │                    └── plot AUC ──> contrasts ──> winner_loser
                 │                                                   └── temporal analysis
                 │
TRY export ──> filtering/taxonomy/range audit ──> species traits
                 │
CHELSA jobs ──> site-year climate ──> experiment climate
                 │
                 └─────────────┬─────────────────────────────────────┘
                               │
                         analysis_data
                               │
                  downstream 526-species boundary
                               │
                  downstream data + sampling V_auc
                               │
             ┌─────────────────┼────────────────────┐
             │                 │                    │
          M_both          family model       predictor branches
             │                                      │
       species context                 LDMC / LMA / Height
                                       Graminoid / Woodiness
                                       phylogeny / PEV1–PEV3
                                                  │
                                  heterogeneity synthesis
                                                  │
                               canonical validation assembly
                                                  │
                             numerical + figure release checks
                                                  │
                                  reports and manifests
```

## 11. The downstream boundary

A deliberate CSV round trip connects the upstream graph to the canonical
526-species analysis:

```text
winner_loser_file
analysis_data_file
climate_context_experiments_file
        │
        └── canonical_*_input
                └── downstream_auc_unscaled_526
                        └── downstream_scaled_526
                                └── downstream_auc_data
```

This preserves the exact numeric representation used by the anchored Step 22
analysis while making every file read visible and hash-tracked.

The covariance is separately reconstructed:

```text
downstream_auc_data + auc_inputs$auc
    └── downstream_sampling_covariance_526
          └── downstream_V_auc
```

These feed the actual canonical model:

```text
downstream_model_data_526
downstream_V_auc
canonical_downstream_settings
        └── mboth_526_refit
```

## 12. Parallel scientific branches

After the downstream data are established, branches can run independently.

For example:

```text
mboth_526_refit
    ├── variance profiles
    ├── CR1 results
    └── species_context_bundle_526
          └── species_classification_526
```

Family:

```text
community + experiments + resident membership
    └── family plot AUC
          └── family covariance
                └── family_model_bundle_526
                      └── family_classification_526
```

Traits and functional groups:

```text
TRY summaries ──> LDMC / LMA / Height specifications ──> model sequences

functional assignments
    ├── Graminoid specification ──> graminoid_model_sequence_526
    └── Woodiness specification ──> woodiness_model_sequence_526
```

Phylogeny:

```text
placement decisions + 526 species
    └── resident_phylogeny_rebuild
          └── resident_pcoa_candidate
                └── resident_pev_scores_oriented
                      └── standardized_pev_species_526
                            └── phylogeny_auc_data_526
                                  └── phylogeny_model_sequence_526
```

Thus, the active PEV model receives newly reconstructed PCoA values. Frozen PEV
values are confined to a separately named historical validation branch.

## 13. Reference objects versus scientific computation

`canonical_transition_reference.rds` and the frozen Step 22 fixture are loaded
into targets, but their role is comparison and migration validation.

For example:

```text
canonical reconstructed result ─┐
                                ├── equivalence check
frozen reference result ────────┘
```

The canonical `M_both` scientific object is `mboth_526_refit`. The target named
`mboth_526_model` is a historical model signature used for comparison, despite
its slightly ambiguous name.

Similarly:

```text
candidate_step22_results + step22_reference_results
    └── step22_reference_validation
```

The fixture is consumed by the validation target, not used to calculate the
candidate result.

## 14. Reports do not fit models

For example, [`09_mboth_526.qmd`](reports/09_mboth_526.qmd) reads completed
targets:

```r
bundle <- targets::tar_read(mboth_526_refit_bundle)
```

The report does not refit `M_both`.

Likewise, [`14_phylogeny.qmd`](reports/14_phylogeny.qmd) reads tree, PCoA, model,
and validation targets. Changing report prose should therefore not invalidate
the model branch.

## 15. What happens on subsequent runs

`targets` stores completed R objects and metadata in `_targets/`. On another
run, it checks:

- input-file content hashes;
- target commands;
- referenced function definitions;
- upstream target values;
- output files declared with `format = "file"`.

Only outdated targets and their affected descendants run again.

Useful commands from the repository root are:

```r
targets::tar_outdated()
targets::tar_visnetwork()
targets::tar_read(mboth_526_refit_bundle)
targets::tar_make(names = "release_validation")
```

The practical summary is:

```text
scripts/run_pipeline.R
    └── tar_make()
          └── _targets.R
                ├── configuration
                ├── reusable R/functions
                ├── upstream/data target definitions
                ├── downstream analysis target definitions
                └── publication/release target definitions
```
