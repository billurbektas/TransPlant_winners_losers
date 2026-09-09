# TransPlant warming winners and losers

This repository has one active scientific entry point: [`_targets.R`](_targets.R).
It builds the community, taxonomy, climate, trait, AUC, temporal, and resident
downstream branches, then validates the assembled downstream result against the
approved Step 22 numerical reference.

The current downstream analysis contains 1,604 resident effects from 526
species in 42 experiments. The six downstream exclusions remain in upstream
community, AUC, and winner/loser calculations:

- *Gymnocarpium dryopteris*
- *Athyrium spinulosum*
- *Botrychium lunaria*
- *Equisetum arvense*
- *Selaginella selaginoides*
- *Juniperus communis*

## Numerical authority and rollback

The immutable numerical authority is:

```text
output/step_by_step_modeling/downstream_526/downstream_526_reference_results.rds
SHA-256 8d5d7544a6e6bcb757827263763ffae95115fa10b52b187669eaa5627af6931f
```

A byte-identical regression fixture is stored at
`tests/fixtures/step22_reference_results.rds`. The pipeline never writes to the
immutable reference path.

Branch-level transition interfaces are condensed in the 180-KB
`tests/fixtures/canonical_transition_reference.rds`. This checksummed fixture
contains only validated tables, covariance matrices, compact model signatures,
and the Step 19 phylogeny/PCoA interface. It allows historical output trees to
be archived without making them hidden active dependencies.

The complete pre-refactor scientific source is recoverable from:

```text
branch checkpoint/pre-refactor-source-2026-09-09
tag    pre-refactor-source-2026-09-09
commit cc9ee54
```

Active refactoring is performed on `refactor/canonical-targets`. The source
checkpoint intentionally excludes raw data, caches, `_targets/`, generated
models and reports, bulk figures, whole-session snapshots, and
`repository_contents.txt`.

## Canonical scientific steps

| Step | Scientific stage | Report |
|---:|---|---|
| 01 | Inputs, configuration, hashes, and scope | `01_inputs_and_scope.qmd` |
| 02 | Community cleaning and adapter validation | `02_community_cleaning.qmd` |
| 03 | Canonical TNRS taxonomy and reviewed decisions | `03_taxonomy.qmd` |
| 04 | CHELSA coordinate audit and experiment climate | `04_climate.qmd` |
| 05 | TRY filtering and species trait summaries | `05_trait_preparation.qmd` |
| 06 | Plot AUC, effects, sampling covariance, validation | `06_auc_validation.qmd` |
| 07 | Upstream winner/loser and temporal analysis | `07_temporal_winner_loser.qmd` |
| 08 | Downstream exclusions and 526-species scope | `08_downstream_scope_526.qmd` |
| 09 | Anchored `M_both`, profiles, and CR1 | `09_mboth_526.qmd` |
| 10 | Species min/max context classification | `10_species_context.qmd` |
| 11 | Plot-cover family reconstruction and model | `11_family_resolution.qmd` |
| 12 | LDMC, LMA, and height information | `12_trait_information.qmd` |
| 13 | Separate Graminoid and Woodiness information | `13_functional_information.qmd` |
| 14 | Phylogeny, response-independent PCoA, PEV model | `14_phylogeny.qmd` |
| 15 | Response-scale heterogeneity synthesis | `15_heterogeneity_synthesis.qmd` |
| 16 | Numerical reference and release validation | `16_reference_and_release.qmd` |

Reports live in [`reports/`](reports/) and are rendered by `tarchetypes` to
`reports/rendered/`. Reports read narrowly scoped targets; they do not refit
models or search output directories.

## Target architecture

The root target file is deliberately small. It sources:

- active upstream targets from [`R/targets/data_targets.R`](R/targets/data_targets.R);
- downstream scientific and validation targets from
  [`R/targets/analysis_targets.R`](R/targets/analysis_targets.R);
- report, support-mask, and provenance targets from
  [`R/targets/publication_targets.R`](R/targets/publication_targets.R);
- pure analytical functions under [`R/functions/`](R/functions/).

The active pipeline-module list is explicit in
`canonical_pipeline_source_files()`. Obsolete review/model modules are not
loaded merely because they remain beside canonical upstream code during the
validation phase, and archived files are never reached through a recursive
source glob.

The major dependency branches are:

```text
community + taxonomy + TRY + CHELSA
                 |
           AUC effects + V_auc
             /             \
 temporal/upstream      downstream scope 526
                         /   |   |   |   \
                    M_both family traits status phylogeny
                         \   |   |   |   /
                       heterogeneity synthesis
                                  |
                  late Step 22 validation assembly
```

`candidate_step22_results` exists only for late regression validation. Ordinary
reports and figures depend on their smallest scientific target, so a plotting
change does not invalidate models and one predictor-model change does not
invalidate unrelated branches.

The settled formulas, exclusion list, expected counts, random structure,
classification rule, output interfaces, and PEV orientation are centralized in
[`config/canonical_analysis.R`](config/canonical_analysis.R). Analytical
functions accept objects and return objects. Only explicit file targets export
results.

## Settled downstream choices

- Response: maximum-duration raw relative-cover AUC, warmed minus origin.
- Sampling covariance: the validated `V_auc` in exact effect order.
- Base fixed effects: `T_z * P_z * W_z`.
- Species random effects: independent intercept, independent `W_z` slope, and
  correlated `T_z`/`P_z`/`TP_z` block (`GEN`).
- Robust inference: experiment-clustered, small-sample-adjusted CR1.
- Context class: minimum fitted response above zero is a consistent winner;
  maximum below zero is a consistent loser; otherwise context-dependent.
- Family analysis: family AUC and covariance are reconstructed from plot cover;
  species effects are not averaged.
- Continuous traits: positive species medians are log transformed and
  standardized across unique trait-complete species before AUC-row expansion.
- Functional information: Graminoid status and Woodiness are separate binary
  predictors.
- Phylogeny: response-independent PEV1–PEV3 only; ecological modeling scales
  them after unique-species PCoA.
- Heterogeneity: overall, warming, context, and total components retain the
  Step 22 definitions.

One downstream taxonomy exception is explicit in
[`config/downstream_species_overrides.csv`](config/downstream_species_overrides.csv):
the TNRS family lookup omits *Viola palustris*. The shared ledger augments the
active functional-group target and supplies its validated Violaceae,
non-graminoid, non-woody assignments to both the family and functional
branches. No Step 20 output is used to fill this row.

The Step 10 callout examples are selected before rendering with explicit rules:
*Potentilla erecta* (highest winner minimum), *Plantago lanceolata* (second
lowest loser maximum under the stable retained rule), and *Elymus repens* plus
*Silene acaulis* (two widest context-dependent spans).

## Phylogeny reproducibility

The lockfile pins `V.PhyloMaker2` 0.1.0 from
`jinyizju/V.PhyloMaker2` at commit
`7af3fb5152f691af2e4ec9d5e2e467d1b50505e9`.

The provenance gate records:

| Artifact | SHA-256 |
|---|---|
| `GBOTB.extended.LCVP`, R serialization v3 | `5a62592d963cfadea4138754a17f9d8154afe57cffc5666a0d7f138926159604` |
| `nodes.info.1.LCVP`, R serialization v3 | `a4716697bc63294920b6dd9535468522b33dab27571bfff198e440c7fc4d1604` |
| current final tree RDS | `a095615895ca5b4a99f7e4515a827de97f79079a4424d6b1f486cc591650b8db` |
| current final Newick | `4b14760c246c5bb7faec58e4c4298d9c539d7822df9511fb624b6180595b1891` |
| current `species_pev_scores.csv` | `9aed55a22a03d10436a64bb378cc61da3387540213ee440aa9eba3555cebf85a` |

The active phylogeny branch rebuilds the tree with the pinned package and the
526-row reviewed placement ledger in
[`config/phylogeny_placement_decisions.csv`](config/phylogeny_placement_decisions.csv).
It does not adopt the historical Step 19 tree RDS as its candidate result.
The historical RDS, Newick, and PEV table remain read-only regression
authorities until topology, branch lengths, documented shared placements,
eigenstructure, and oriented PEV1–PEV3 coordinates all pass.

PCoA signs use a response-independent convention: *Poa trivialis* is positive
on PEV1, *Erigeron formosissimus* is negative on PEV2, and *Euphrasia minima*
is positive on PEV3. Release validation checks species/order mapping, topology,
branch lengths, zero-distance placements, PEV1–PEV3 coordinates, and positive
variation percentages. Numerical equivalence is authoritative even when file
serialization metadata changes a checksum.

## Figures

Figure status is separate from numerical status:

- `CURRENT_CANONICAL_ANALYSIS_FIGURE`: generated from the current anchored
  analysis but open to aesthetic, mask, labeling, or panel refinement.
- `FINAL_MANUSCRIPT`: assigned only after explicit user approval.
- `KEEP_FOR_MANUSCRIPT`, `KEEP_HISTORICAL`, `REGENERATABLE_AUDIT`,
  `INTERMEDIATE_DISCARDABLE`, and `UNKNOWN_DO_NOT_TOUCH` protect other roles.

The complete inventory is
`output/manifests/figure_provenance.csv`. Eleven explicitly retained figures
are copied with checksum verification from `plot/` to `figures/retained/`;
their originals remain in place. The retained allowlist and expected hashes are
in [`config/retained_figures.csv`](config/retained_figures.csv).

The Woodiness surface uses a group-specific observed-climate support-mask
target. This target derives support separately for woody and non-woody rows and
cannot alter or refit the anchored Woodiness model. Regenerated scientific
figures are not subjected to pixel or file-hash regression tests.

## Installation and execution

Restore the exact R dependencies:

```r
renv::restore()
```

Quarto is a system requirement for reports. The pipeline uses `quarto` on
`PATH`; on macOS it also recognizes the executable bundled with RStudio. The R
packages `quarto` and `tarchetypes` are locked in `renv.lock`.

Run the complete active pipeline:

```r
targets::tar_make()
```

or:

```sh
Rscript scripts/run_pipeline.R
```

Useful scoped gates are:

```r
targets::tar_make(names = c(
  step22_reference_validation_file,
  phylogeny_reproduction_validation_file
))

targets::tar_make(names = c(
  published_report_10_species_context,
  published_report_13_functional_information,
  published_report_16_reference_and_release
))

targets::tar_make(names = release_validation)
```

Dependency refresh is intentionally separate from analysis:

```sh
Rscript scripts/refresh_dependencies.R
```

CHELSA and first-time TNRS/TRY construction require network access. Large raw
sources stay in their existing locations and are catalogued in
[`data/source_manifest.csv`](data/source_manifest.csv); they are not copied for
cosmetic restructuring. The CHELSA monthly cache is resumable and remains
ignored.

## Validation

Run local tests with:

```r
testthat::test_dir("tests/testthat")
```

The Step 22 comparator uses exact checks for identifiers, structures, counts,
classes, convergence flags, and PASS/FAIL fields. Prespecified tolerances are:

| Quantity | Absolute | Relative |
|---|---:|---:|
| deterministic preprocessing and scaling | `1e-12` | `1e-10` |
| fixed effects, CR1, likelihoods, IC, heterogeneity | `1e-8` | `1e-6` |
| optimizer-derived covariance and profiles | `1e-6` | `1e-4` |

Tolerance changes require explicit review. Figure pixels and regenerated
figure hashes are excluded from numerical regression.

## Outputs, provenance, and archive policy

Canonical exports are organized under:

```text
output/data/
output/models/
output/tables/
output/audits/
output/manifests/
```

`output/manifests/canonical_outputs.csv` is produced from the explicit
canonical file targets and records each exported artifact and checksum.
`output/manifests/active_scientific_source_inventory.csv` records the path,
role, size, and SHA-256 of every active source, configuration, report, script,
test, and compact regression fixture used by this lineage.
`output/manifests/canonical_run_provenance.csv` separately records the release
time, Git branch/commit and dirty state, R/platform details, key package
versions, and checksums of the lockfile, target script, Step 22 authority, and
compact transition fixture. It is downstream of the numerical release gate and
does not alter the immutable 12-field Step 22 schema.

The active graph does not read superseded Step output directories. Validated
Step 19–23 transition interfaces are held in the compact, immutable fixture
described above. Old files are not moved until their replacement targets pass
independently and the combined numerical, phylogeny, Step 22, and retained-
figure release gate passes.

Historical material is staged only under [`archive/`](archive/) with a manifest
recording original path, destination, category, checksum, producer,
dependencies, reason, and date. Archived R files are excluded from active
source globs, targets, reports, and tests.

`scripts/archive_validated_history.R` is dry-run-only by default and refuses to
proceed unless `release_validation` is current and fully passing. The explicit
`--apply` mode moves eligible source and generated history, verifies every
post-move SHA-256, and writes the archive ledger. It never moves `plot/`, raw
data, caches, `_targets/`, session snapshots, or the immutable Step 22 reference.

Nothing is deleted during the initial refactor. Deletion of snapshots, caches,
duplicate figures, generated outputs, or archived material requires an explicit
itemized approval with checksums and regeneration instructions. Unclassified
figures remain `UNKNOWN_DO_NOT_TOUCH`.
