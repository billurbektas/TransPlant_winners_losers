# TransPlant winners–losers analysis

This repository contains the reproducible TransPlant community-response and
resident-model workflow. Start with [`analysis_workflow.qmd`](analysis_workflow.qmd):
it follows the science from the original data through response calculation,
explicit model specifications, inference, figures, and tables. Its main
analytical assignments are genuine targets through `tarchetypes::tar_tangle()`.

The authoritative architecture, target-stage table, naming rules, commands,
failure behavior, and two-phase history-retirement procedure are documented in
[`WORKFLOW_README.md`](WORKFLOW_README.md).

Quick start:

```sh
Rscript -e 'renv::restore()'
Rscript scripts/run_pipeline.R
```

Those commands are for the approved clean rebuild. The current review stage is
inspection-only; no scientific target should be run yet.

Figure renderers are available through one dispatcher:

```sh
Rscript scripts/render_figure.R --list
```

Tracked historical sources are preserved under the annotated recovery tag
`pre-reviewable-workflow-with-history-2026-09-30`. The archive, fixtures, and
old reports have been retired from the working tree. Current generated results
and target stores were preserved in the verified ZIP recorded in
`docs/pre_rebuild_output_snapshot.csv`. No clean rebuild has been run yet.
