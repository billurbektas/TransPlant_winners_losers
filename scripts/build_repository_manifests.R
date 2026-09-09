#!/usr/bin/env Rscript

# Rebuilds source/output inventories. It does not move, delete, or modify any
# scientific artifact outside output/manifests.

source("R/functions/provenance.R")
source("config/canonical_analysis.R")

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
manifest_dir <- file.path(root, "output", "manifests")
dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)

relative_files <- function(directory, pattern = NULL) {
  path <- file.path(root, directory)
  if (!dir.exists(path)) return(character())
  files <- list.files(path, recursive = TRUE, full.names = TRUE, all.files = TRUE)
  files <- files[file.exists(files) & !dir.exists(files)]
  if (!is.null(pattern)) files <- files[grepl(pattern, files)]
  relative_project_path(files, root)
}

hash_inventory <- function(paths) {
  full <- file.path(root, paths)
  info <- file.info(full)
  data.frame(
    path = paths,
    size_bytes = as.numeric(info$size),
    modified_utc = format(info$mtime, tz = "UTC", usetz = TRUE),
    sha256 = vapply(full, sha256_file, character(1)),
    stringsAsFactors = FALSE
  )
}

tracked <- system2("git", "ls-files", stdout = TRUE)
status_lines <- system2(
  "git", c("status", "--porcelain=v1", "--untracked-files=all"),
  stdout = TRUE
)
status_paths <- if (length(status_lines)) substring(status_lines, 4L) else character()
status_codes <- if (length(status_lines)) substring(status_lines, 1L, 2L) else character()
git_state <- function(path) {
  index <- match(path, status_paths)
  if (!is.na(index)) {
    if (identical(status_codes[[index]], "??")) return("untracked")
    return(paste0("tracked_", gsub(" ", "_", status_codes[[index]], fixed = TRUE)))
  }
  if (path %in% tracked) "tracked_clean" else "untracked"
}

root_source <- c(
  "_targets.R", ".gitignore", "README.md", "renv.lock",
  "TransPlant_winners_losers.Rproj", "external/get_chelsa.R"
)
source_paths <- sort(unique(c(
  root_source[file.exists(file.path(root, root_source))],
  relative_files("R", "[.][Rr]$"),
  relative_files("config", "[.](R|csv)$"),
  relative_files("tests", "[.](R|rds|csv|md)$"),
  relative_files("reports", "[.]qmd$"),
  relative_files("scripts", "[.]R$"),
  relative_files("data", "(README[.]md|source_manifest[.]csv)$"),
  relative_files("archive", "(README[.]md|archive_moves[.]csv)$")
)))
source_manifest <- hash_inventory(source_paths)
source_manifest$git_state <- vapply(source_paths, git_state, character(1))
active_paths <- relative_project_path(
  canonical_scientific_source_files(root), root
)
source_manifest$active_source <- source_manifest$path %in% active_paths
source_manifest$source_role <- dplyr::case_when(
  source_manifest$active_source ~ "ACTIVE_CANONICAL_SOURCE",
  grepl("^archive/", source_manifest$path) ~ "ARCHIVED_SOURCE_OR_LEDGER",
  TRUE ~ "INACTIVE_PENDING_ARCHIVE"
)
readr::write_csv(
  source_manifest,
  file.path(manifest_dir, "active_source_manifest.csv"),
  na = ""
)

output_paths <- relative_files("output")
output_paths <- output_paths[!grepl("^output/manifests/", output_paths)]
output_manifest <- hash_inventory(sort(output_paths))
output_manifest$artifact_class <- dplyr::case_when(
  grepl("/models?/|[.]rds$", output_manifest$path) ~ "model_or_r_object",
  grepl("/figures?/|[.](pdf|png|svg)$", output_manifest$path, ignore.case = TRUE) ~ "figure",
  grepl("audit|report", output_manifest$path, ignore.case = TRUE) ~ "audit_or_report",
  TRUE ~ "table_or_other_output"
)
output_manifest$migration_status <- "CURRENT_PATH_UNMOVED"
readr::write_csv(
  output_manifest,
  file.path(manifest_dir, "output_artifact_manifest.csv"),
  na = ""
)

message("Wrote ", nrow(source_manifest), " source rows and ",
        nrow(output_manifest), " output rows without moving artifacts.")
