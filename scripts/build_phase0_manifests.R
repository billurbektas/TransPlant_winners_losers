source("R/functions/provenance.R")

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
manifest_dir <- file.path(root, "output", "manifests")
dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)

figure_manifest <- build_figure_provenance(root)
utils::write.csv(
  figure_manifest,
  file.path(manifest_dir, "figure_provenance.csv"),
  row.names = FALSE,
  na = ""
)

reference_manifest <- verify_reference_artifacts(root)
utils::write.csv(
  reference_manifest,
  file.path(manifest_dir, "reference_artifact_verification.csv"),
  row.names = FALSE,
  na = ""
)

package_manifest <- phylogeny_package_provenance()
utils::write.csv(
  package_manifest,
  file.path(manifest_dir, "phylogeny_package_provenance.csv"),
  row.names = FALSE,
  na = ""
)

if (!all(reference_manifest$hash_matches %in% TRUE)) {
  stop("One or more immutable reference artifacts failed checksum verification.", call. = FALSE)
}

retained_rows <- figure_manifest$proposed_status == "KEEP_FOR_MANUSCRIPT"
if (any(!figure_manifest$expected_hash_matches[retained_rows])) {
  stop("One or more explicitly retained figures failed checksum verification.", call. = FALSE)
}

message("Wrote Phase 0 manifests for ", nrow(figure_manifest), " figures and ", nrow(reference_manifest), " reference artifacts.")
