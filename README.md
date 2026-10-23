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
TRANSPLANT_CHELSA_OFFLINE=true Rscript scripts/run_pipeline.R
```

The clean rebuild completed on 2 October 2026 using the validated monthly
CHELSA cache. The final release checks and test suite passed with one explicitly
documented, unresolved model-profile interval; see
[`docs/rebuild_known_issues.md`](docs/rebuild_known_issues.md). The rebuilt
scientific outputs still await review and acceptance before a new release
reference is frozen. `TRANSPLANT_CHELSA_OFFLINE=true` requires the 336 monthly
checkpoints and prevents new CHELSA downloads.

Figure renderers are available through one dispatcher:

```sh
Rscript scripts/render_figure.R --list
```

Tracked historical sources are preserved under the annotated recovery tag
`pre-reviewable-workflow-with-history-2026-09-30`. The archive, fixtures, and
old reports have been retired from the working tree. Previous generated results
and target stores are recoverable from the verified ZIP recorded in
[`docs/pre_rebuild_backup_receipts.csv`](docs/pre_rebuild_backup_receipts.csv).
