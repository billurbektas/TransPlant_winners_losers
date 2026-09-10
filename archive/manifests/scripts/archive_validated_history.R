#!/usr/bin/env Rscript

# Gated archive migration for the canonical targets refactor.
#
# Default mode is a read-only plan. Use --apply only after `release_validation`
# is current and fully passing. This script moves files; it never deletes data.

args <- commandArgs(trailingOnly = TRUE)
do_apply <- identical(args, "--apply")
if (length(args) && !do_apply) {
  stop("Usage: Rscript scripts/archive_validated_history.R [--apply]",
       call. = FALSE)
}

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
expected_branch <- "refactor/canonical-targets"
branch <- system2("git", c("branch", "--show-current"), stdout = TRUE)
if (!identical(branch, expected_branch)) {
  stop("Archive migration requires branch ", expected_branch, ".", call. = FALSE)
}

reference <- file.path(
  root, "output", "step_by_step_modeling", "downstream_526",
  "downstream_526_reference_results.rds"
)
reference_hash <- if (file.exists(reference)) {
  strsplit(
    system2("shasum", c("-a", "256", reference), stdout = TRUE)[[1]],
    "  ", fixed = TRUE
  )[[1]][[1]]
} else {
  NA_character_
}
expected_reference_hash <-
  "8d5d7544a6e6bcb757827263763ffae95115fa10b52b187669eaa5627af6931f"
if (!identical(reference_hash, expected_reference_hash)) {
  stop("The immutable Step 22 reference is missing or has changed.", call. = FALSE)
}

outdated <- targets::tar_outdated(callr_function = NULL)
if ("release_validation" %in% outdated) {
  stop(
    "The release_validation target is stale; rebuild it before archiving.",
    call. = FALSE
  )
}

release <- tryCatch(
  targets::tar_read_raw("release_validation"),
  error = function(error) {
    stop(
      "Build the current release_validation target before archiving: ",
      conditionMessage(error), call. = FALSE
    )
  }
)
if (!is.data.frame(release) || !"pass" %in% names(release) ||
    any(is.na(release$pass) | !release$pass)) {
  stop("The current release validation is not fully passing.", call. = FALSE)
}

source_groups <- list(
  pretargets_workflow_history = c(
    "R/00_refresh_dependencies.R", "R/00_run_pipeline.R", "R/00_workflow.R",
    "R/01_import.R", "R/02_organize.R", "R/03_calculate_winners_losers.R",
    "R/figures_presentation.R", "R/functions.R", "R/theme_ggplot.R"
  ),
  legacy_modeling_v1 = c(
    sprintf("R/%02d", 4:14), "R/functional_groups.R", "R/reporting",
    "R/pipeline/biological_context_review.R",
    "R/pipeline/climate_response_model.R",
    "R/pipeline/additional_analysis.R",
    "R/pipeline/figures.R",
    "R/pipeline/frequentist_hierarchical_auc.R",
    "R/pipeline/frequentist_hierarchical_auc_reporting.R",
    "R/pipeline/hierarchical_trait_trial.R",
    "tests/testthat/test-biological-context-review.R",
    "tests/testthat/test-climate-response-model.R",
    "tests/testthat/test-frequentist-hierarchical-auc.R",
    "tests/testthat/test-hierarchical-ldmc-reductions.R",
    "tests/testthat/test-hierarchical-trait-trial.R"
  ),
  stepwise_model_development = sprintf("R/%02d", 15:29),
  canonical_source_snapshots = sprintf("R/%02d", 30:35)
)

resolve_prefixes <- function(prefixes) {
  all_files <- list.files(file.path(root, "R"), recursive = TRUE,
                          full.names = FALSE, all.files = TRUE)
  all_files <- file.path("R", all_files)
  unique(unlist(lapply(prefixes, function(prefix) {
    exact <- file.path(root, prefix)
    if (file.exists(exact) && !dir.exists(exact)) return(prefix)
    all_files[startsWith(all_files, paste0(prefix, "_")) |
                startsWith(all_files, paste0(prefix, "/"))]
  }), use.names = FALSE))
}

source_rows <- do.call(rbind, lapply(names(source_groups), function(category) {
  paths <- resolve_prefixes(source_groups[[category]])
  if (!length(paths)) return(NULL)
  destination <- file.path("archive", category, paths)
  data.frame(
    original_path = paths,
    destination_path = destination,
    category = category,
    producer = "historical scientific source",
    dependencies = "pre-refactor source checkpoint cc9ee54",
    reason = "inactive after canonical function/target extraction",
    stringsAsFactors = FALSE
  )
}))

checkpoint_tag <- "pre-refactor-source-2026-09-09"
checkpoint_files <- system2(
  "git", c("ls-tree", "-r", "--name-only", checkpoint_tag), stdout = TRUE
)
missing_from_checkpoint <- setdiff(source_rows$original_path, checkpoint_files)
if (length(missing_from_checkpoint)) {
  stop(
    "Historical source lacks checkpoint coverage: ",
    paste(missing_from_checkpoint, collapse = ", "), call. = FALSE
  )
}

output_groups <- c(
  "output/frequentist_hierarchical_auc" =
    "archive/legacy_modeling_v1/output/frequentist_hierarchical_auc",
  "output/step_by_step_modeling/auc_intercept_species" =
    "archive/stepwise_model_development/output/auc_intercept_species",
  "output/step_by_step_modeling/species_context" =
    "archive/stepwise_model_development/output/species_context",
  "output/step_by_step_modeling/species_summary_plots" =
    "archive/stepwise_model_development/output/species_summary_plots",
  "output/step_by_step_modeling/species_both" =
    "archive/stepwise_model_development/output/species_both",
  "output/step_by_step_modeling/species_both_traits" =
    "archive/stepwise_model_development/output/species_both_traits",
  "output/step_by_step_modeling/species_information_heterogeneity" =
    "archive/stepwise_model_development/output/species_information_heterogeneity",
  "output/step_by_step_modeling/family_aggregated_auc" =
    "archive/stepwise_model_development/output/family_aggregated_auc",
  "output/step_by_step_modeling/trait_data_audit" =
    "archive/canonical_source_snapshots/output/trait_data_audit",
  "output/step_by_step_modeling/resident_phylogeny_pcoa" =
    "archive/canonical_source_snapshots/output/resident_phylogeny_pcoa",
  "output/step_by_step_modeling/downstream_526" =
    "archive/canonical_source_snapshots/output/downstream_526"
)

output_rows <- do.call(rbind, lapply(names(output_groups), function(source_dir) {
  absolute_source <- file.path(root, source_dir)
  if (!dir.exists(absolute_source)) return(NULL)
  files <- list.files(absolute_source, recursive = TRUE, full.names = FALSE,
                      all.files = TRUE, no.. = TRUE)
  if (identical(source_dir, "output/step_by_step_modeling/downstream_526")) {
    files <- setdiff(files, "downstream_526_reference_results.rds")
  }
  if (!length(files)) return(NULL)
  category <- if (grepl("frequentist", source_dir, fixed = TRUE)) {
    "legacy_modeling_v1"
  } else if (grepl("trait_data_audit|resident_phylogeny_pcoa|downstream_526",
                   source_dir)) {
    "canonical_source_snapshots"
  } else {
    "stepwise_model_development"
  }
  data.frame(
    original_path = file.path(source_dir, files),
    destination_path = file.path(output_groups[[source_dir]], files),
    category = category,
    producer = "historical generated analysis output",
    dependencies = "condensed in canonical_transition_reference.rds where required",
    reason = "superseded by canonical targets exports after numerical equivalence",
    stringsAsFactors = FALSE
  )
}))

moves <- rbind(source_rows, output_rows)
moves <- moves[file.exists(file.path(root, moves$original_path)), , drop = FALSE]
if (!nrow(moves)) {
  message("No eligible files remain at their original paths.")
  quit(status = 0)
}
if (any(grepl("(^|/)plot/|data/cache|_targets|Environment_Save_", moves$original_path))) {
  stop("The move plan contains a prohibited data/cache/figure target.", call. = FALSE)
}
if (reference %in% file.path(root, moves$original_path)) {
  stop("The immutable Step 22 reference entered the move plan.", call. = FALSE)
}
if (any(file.exists(file.path(root, moves$destination_path)))) {
  stop("At least one archive destination already exists; refusing overwrite.",
       call. = FALSE)
}

hash_file <- function(path) {
  strsplit(
    system2("shasum", c("-a", "256", path), stdout = TRUE)[[1]],
    "  ", fixed = TRUE
  )[[1]][[1]]
}
moves$sha256 <- vapply(file.path(root, moves$original_path), hash_file, character(1))
moves$archive_date <- as.character(Sys.Date())
moves$status <- if (do_apply) "ARCHIVED_HASH_VERIFIED" else "PLANNED"
moves <- moves[, c(
  "original_path", "destination_path", "category", "sha256", "producer",
  "dependencies", "reason", "archive_date", "status"
)]

cat(nrow(moves), "files eligible for archive migration.\n")
print(table(moves$category))
cat("Total bytes:", sum(file.info(file.path(root, moves$original_path))$size), "\n")
if (!do_apply) {
  cat("Dry run only. Re-run with --apply after reviewing the release gate.\n")
  quit(status = 0)
}

manifest_path <- file.path(root, "archive", "manifests", "archive_moves.csv")
existing <- utils::read.csv(manifest_path, stringsAsFactors = FALSE,
                            check.names = FALSE)
write_manifest_atomic <- function(rows) {
  temporary <- paste0(manifest_path, ".tmp")
  utils::write.csv(rows, temporary, row.names = FALSE, na = "")
  if (!file.rename(temporary, manifest_path)) {
    unlink(temporary)
    stop("Could not install the archive manifest atomically.", call. = FALSE)
  }
}

for (index in seq_len(nrow(moves))) {
  source <- file.path(root, moves$original_path[[index]])
  destination <- file.path(root, moves$destination_path[[index]])
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (!file.rename(source, destination)) {
    stop("Archive move failed for ", moves$original_path[[index]], call. = FALSE)
  }
  if (!identical(hash_file(destination), moves$sha256[[index]])) {
    stop("Post-move checksum failed for ", moves$destination_path[[index]],
         call. = FALSE)
  }
  existing <- rbind(existing, moves[index, , drop = FALSE])
  write_manifest_atomic(existing)
}

cat("Archived", nrow(moves), "files with verified SHA-256 checksums.\n")
