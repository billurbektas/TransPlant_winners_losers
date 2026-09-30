# TransPlant winners–losers analysis

This repository contains the reproducible TransPlant community-response and
resident-model workflow. Active code is organized by scientific domain under
`R/`, the executable graph is defined in `_targets.R`, and the analysis is
reported in seven stage-level Quarto documents.

The authoritative architecture, target-stage table, naming rules, commands,
failure behavior, and two-phase history-retirement procedure are documented in
[`WORKFLOW_README.md`](WORKFLOW_README.md).

Quick start:

```sh
Rscript -e 'renv::restore()'
Rscript scripts/run_pipeline.R
```

Figure renderers are available through one dispatcher:

```sh
Rscript scripts/render_figure.R --list
```

The refactor is currently in Phase 1. Historical fixtures and `archive/` are
retained only for certification and will not be deleted until the rebuilt
workflow is cleanly validated and explicitly accepted.
