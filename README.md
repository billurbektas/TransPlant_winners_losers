# TransPlant winners–losers analysis

This repository contains the reproducible TransPlant community-response and
resident-model workflow. Start with [`analysis_workflow.qmd`](analysis_workflow.qmd):
it follows the science from the original data through response calculation,
explicit model specifications, inference, figures, and tables. Its main
analytical assignments are genuine targets through `tarchetypes::tar_tangle()`.

The target-stage table, naming rules, commands, failure behavior, and backup
records are documented in
[`WORKFLOW_README.md`](WORKFLOW_README.md).

Quick start:

```sh
Rscript -e 'renv::restore()'
Rscript scripts/run_pipeline.R
```

Do not run the pipeline yet: the old outputs have been backed up and removed,
and the clean rebuild awaits approval. For that rebuild, set
`TRANSPLANT_CHELSA_OFFLINE=true` to require the 336 validated monthly files
without downloading any CHELSA data.

Figure renderers are available through one dispatcher:

```sh
Rscript scripts/render_figure.R --list
```

Tracked historical sources are preserved under the annotated recovery tag
`pre-reviewable-workflow-with-history-2026-09-30`. The archive, fixtures, and
old reports have been retired from the working tree. Previous generated results
and target stores are recoverable from the verified ZIP recorded in
[`docs/pre_rebuild_backup_receipts.csv`](docs/pre_rebuild_backup_receipts.csv).
No clean rebuild has been run yet.
