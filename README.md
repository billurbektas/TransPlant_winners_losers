# TransPlant winners–losers analysis

The analysis for the TransPlant "winners and losers" manuscript, rebuilt as a
sequence of numbered R scripts. What the analysis is lives in
[`PLAN.md`](PLAN.md); the steps, their checks and their review status in
[`MILESTONES.md`](MILESTONES.md); how the work is done in
[`CLAUDE.md`](CLAUDE.md). Each step writes a review report under
[`reports/`](reports/), collected by `reports/REPORT.qmd`.

## Layout

- `prep/`: a `targets` pipeline for data preparation only, from the raw network
  SQLite, TRY, CHELSA, our GBIF download and the reviewed ledgers in
  `prep/config/` to the input tables under `data/derived/prep/`.
- `analysis/`: one script per step (`01_relative_cover.R` to `06b_forest.R`, the
  colonizers `07a` to `07f`, and `10_figures.R`), run in order; each reads only
  `data/derived/`.
  `analysis/additional_checks/` holds side checks that feed no step.
- `R/`: functions shared by two or more analysis scripts.
- `data/derived/NN_<name>/`: each step's outputs; `reports/` the review reports.
- `figures/main/fig1.png` to `fig6.png`: Sankey, residents' classification,
  colonizers' classification, climate space, traits (pathways), the residents'
  ladder. `figures/appendix/`: Figs S1 to S7. `tables/`: Tables S1 to S5.
  `figures/diagnostics/`: each step's review figures.
- `legacy/`: the old targets pipeline, kept read-only for reference.

## Running

```sh
Rscript -e 'renv::restore()'
Rscript -e 'targets::tar_make(script = "prep/_targets.R", store = "prep/_targets")'
LANG=en_US.UTF-8 Rscript analysis/01_relative_cover.R   # then 02, 03, 03b, 04, 05, 06, 06b, 07a-07f, 10
```

Run from the repository root. A UTF-8 locale keeps the figure labels (τ², ×10⁻⁴)
intact. Fitted models and fold results are cached under `data/derived/`, so a
rerun skips what exists.

## Data and citation

- Community data: the TransPlant network database
  (`data/transplant_network_clean.sqlite`, from
  TransPlantNetwork/transplant_network_data_cleaning).
- Species occurrences for the climatic niche: GBIF.org (28 October 2026) GBIF
  Occurrence Download https://doi.org/10.15468/dl.7p28y9 (download key
  0011587-261023074918325). This DOI must be cited in the methods and the data
  availability statement; GBIF deletes downloads that are not cited after
  28 April 2027.
- Climate: CHELSA V2.1. Traits: the TRY database export.
