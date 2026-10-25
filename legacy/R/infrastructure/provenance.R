# Module: R/infrastructure/provenance
# Responsibility: Record source, package, figure, and release provenance.
# Inputs: Source paths, Git state, registered figures, and retained originals.
# Outputs: Checksums, manifests, provenance and retained-figure checks.
# Side effects: Reads files and Git metadata; explicit copy helper writes retained copies.

# Provenance helpers for the stage-based targets workflow.

#' Sha256 File.
#'
#' @param path A project-relative or absolute filesystem path.
#' @return An R object described by the function name; its structure is validated before later stages use it.
sha256_file <- function(path) {
  if (!file.exists(path) || dir.exists(path)) {
    return(NA_character_)
  }
  output <- system2("shasum", c("-a", "256", path), stdout = TRUE, stderr = TRUE)
  if (!length(output)) {
    return(NA_character_)
  }
  strsplit(output[[1]], "  ", fixed = TRUE)[[1]][[1]]
}

#' Git Provenance Value.
#'
#' @param root A project-relative or absolute filesystem path.
#' @param arguments A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
git_provenance_value <- function(root, arguments) {
  output <- suppressWarnings(
    system2(
      "git", c("-C", root, arguments),
      stdout = TRUE, stderr = TRUE
    )
  )
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) return(NA_character_)
  paste(output, collapse = "\n")
}

#' Build Run Provenance.
#'
#' @param root A project-relative or absolute filesystem path.
#' @param release_gate A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
build_run_provenance <- function(root = getwd(), release_gate = NULL) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  if (!is.null(release_gate)) {
    assert_release_validation(release_gate)
  }

  status <- git_provenance_value(
    root, c("status", "--porcelain=v1", "--untracked-files=normal")
  )
  package_version_or_na <- function(package) {
    if (!requireNamespace(package, quietly = TRUE)) return(NA_character_)
    as.character(utils::packageVersion(package))
  }
  artifact <- function(path) sha256_file(file.path(root, path))

  values <- c(
    recorded_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    git_branch = git_provenance_value(root, c("branch", "--show-current")),
    git_commit = git_provenance_value(root, c("rev-parse", "HEAD")),
    git_worktree_dirty = as.character(isTRUE(nzchar(status))),
    r_version = R.version.string,
    r_platform = R.version$platform,
    operating_system = paste(Sys.info()[c("sysname", "release", "machine")], collapse = " "),
    targets_version = package_version_or_na("targets"),
    tarchetypes_version = package_version_or_na("tarchetypes"),
    quarto_r_package_version = package_version_or_na("quarto"),
    v_phylomaker2_version = package_version_or_na("V.PhyloMaker2"),
    renv_lock_sha256 = artifact("renv.lock"),
    targets_script_sha256 = artifact("_targets.R")
  )
  tibble::tibble(field = names(values), value = unname(values))
}

#' Relative Project Path.
#'
#' @param path A project-relative or absolute filesystem path.
#' @param root A project-relative or absolute filesystem path.
#' @return An R object described by the function name; its structure is validated before later stages use it.
relative_project_path <- function(path, root = getwd()) {
  normalized_root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  normalized_path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  prefix <- paste0(normalized_root, "/")
  ifelse(startsWith(normalized_path, prefix), substring(normalized_path, nchar(prefix) + 1L), normalized_path)
}

#' Figure Dimensions.
#'
#' @param path A project-relative or absolute filesystem path.
#' @return An R object described by the function name; its structure is validated before later stages use it.
figure_dimensions <- function(path) {
  extension <- tolower(tools::file_ext(path))
  width <- NA_real_
  height <- NA_real_
  unit <- NA_character_

  if (identical(extension, "pdf") && nzchar(Sys.which("pdfinfo"))) {
    info <- system2("pdfinfo", path, stdout = TRUE, stderr = FALSE)
    size_line <- grep("^Page size:", info, value = TRUE)
    if (length(size_line)) {
      values <- regmatches(size_line[[1]], gregexpr("[0-9.]+", size_line[[1]]))[[1]]
      if (length(values) >= 2L) {
        width <- as.numeric(values[[1]])
        height <- as.numeric(values[[2]])
        unit <- "points"
      }
    }
  }

  if (identical(extension, "png") && nzchar(Sys.which("sips"))) {
    info <- system2("sips", c("-g", "pixelWidth", "-g", "pixelHeight", path), stdout = TRUE, stderr = FALSE)
    width_line <- grep("pixelWidth:", info, value = TRUE)
    height_line <- grep("pixelHeight:", info, value = TRUE)
    if (length(width_line) && length(height_line)) {
      width <- as.numeric(sub(".*pixelWidth: *", "", width_line[[1]]))
      height <- as.numeric(sub(".*pixelHeight: *", "", height_line[[1]]))
      unit <- "pixels"
    }
  }

  data.frame(width = width, height = height, dimension_unit = unit)
}

#' Known Figure Producer.
#'
#' @param path A project-relative or absolute filesystem path.
#' @return An R object described by the function name; its structure is validated before later stages use it.
known_figure_producer <- function(path) {
  registry <- figure_registry()
  registered <- registry[registry$filename == path, , drop = FALSE]
  if (nrow(registered) == 1L) {
    return(data.frame(
      generating_script_if_known = registered$script[[1]],
      target_if_known = registered$target[[1]]
    ))
  }
  filename <- basename(path)
  script <- NA_character_
  target <- NA_character_

  target_map <- c(
    "compare_rel.abund_vs.sqrt.rel.abund.pdf" = "figure_sensitivity",
    "destination_climate_warming_gradient.pdf" = "figure_destination_climate_warming_gradient"
  )
  if (filename %in% names(target_map)) {
    target <- unname(target_map[[filename]])
    script <- dplyr::case_when(
      filename == "destination_climate_warming_gradient.pdf" ~
        "R/figures/climate.R",
      TRUE ~ "R/figures/consistency.R"
    )
  }
  additional_script_map <- c(
    "repeated_species_consistency.pdf" = "R/figures/consistency.R",
    "temporal_transition_matrices.pdf" = "R/figures/temporal.R"
  )
  if (filename %in% names(additional_script_map) && is.na(script)) {
    script <- unname(additional_script_map[[filename]])
  }

  if (grepl("^figures/manuscript/fig04_species_minmax_classification[.]pdf$", path)) {
    script <- "R/figures/species_context.R"
    target <- "figure04_species_minmax_classification_pdf"
  } else if (grepl("^figures/manuscript/fig04_species_minmax_classification[.]png$", path)) {
    script <- "R/figures/species_context.R"
    target <- "figure04_species_minmax_classification_png"
  } else if (grepl("^figures/manuscript/fig05_climate_only_fixed_effect_prediction_surfaces[.]pdf$", path)) {
    script <- "scripts/render_figure.R"
    target <- "figure05_climate_only_fixed_effect_pdf"
  } else if (grepl("^figures/manuscript/fig05_climate_only_fixed_effect_prediction_surfaces[.]png$", path)) {
    script <- "scripts/render_figure.R"
    target <- "figure05_climate_only_fixed_effect_png"
  } else if (grepl("^figures/manuscript/fig06_predictor_synthesis_support_and_heterogeneity[.]pdf$", path)) {
    script <- "R/figures/predictor_synthesis.R"
    target <- "figure06_predictor_synthesis_pdf"
  } else if (grepl("^figures/manuscript/fig06_predictor_synthesis_support_and_heterogeneity[.]png$", path)) {
    script <- "R/figures/predictor_synthesis.R"
    target <- "figure06_predictor_synthesis_png"
  } else if (grepl("^figures/manuscript/fig07_predictor_response_atlas[.]pdf$", path)) {
    script <- "scripts/render_figure.R"
    target <- "figure07_predictor_response_atlas_pdf"
  } else if (grepl("^figures/manuscript/fig07_predictor_response_atlas[.]png$", path)) {
    script <- "scripts/render_figure.R"
    target <- "figure07_predictor_response_atlas_png"
  } else if (grepl("^figures/manuscript/appendix_family_minmax_classification[.]pdf$", path)) {
    script <- "scripts/render_figure.R"
    target <- "appendix_family_minmax_pdf"
  } else if (grepl("^figures/manuscript/appendix_family_minmax_classification[.]png$", path)) {
    script <- "scripts/render_figure.R"
    target <- "appendix_family_minmax_png"
  } else if (grepl("^figures/manuscript/appendix_pev2_pev3_prediction_atlas[.]pdf$", path)) {
    script <- "scripts/render_figure.R"
    target <- "appendix_pev2_pev3_atlas_pdf"
  } else if (grepl("^figures/manuscript/appendix_pev2_pev3_prediction_atlas[.]png$", path)) {
    script <- "scripts/render_figure.R"
    target <- "appendix_pev2_pev3_atlas_png"
  } else if (grepl("^figures/manuscript/appendix_phylogeny_pcoa_summary[.]pdf$", path)) {
    script <- "scripts/render_figure.R"
    target <- "appendix_phylogeny_pcoa_summary_pdf"
  } else if (grepl("^figures/manuscript/appendix_phylogeny_pcoa_summary[.]png$", path)) {
    script <- "scripts/render_figure.R"
    target <- "appendix_phylogeny_pcoa_summary_png"
  } else if (grepl("^figures/manuscript/appendix_model_fixed_and_random_effects[.]pdf$", path)) {
    script <- "R/evaluation/model_effects.R"
    target <- "appendix_model_fixed_and_random_effects_pdf"
  } else if (grepl("^figures/manuscript/appendix_model_fixed_and_random_effects[.]png$", path)) {
    script <- "R/evaluation/model_effects.R"
    target <- "appendix_model_fixed_and_random_effects_png"
  } else if (grepl("^figures/manuscript/appendix_model_predictive_performance[.]pdf$", path)) {
    script <- "R/evaluation/model_performance.R"
    target <- "appendix_model_predictive_performance_pdf"
  } else if (grepl("^figures/manuscript/appendix_model_predictive_performance[.]png$", path)) {
    script <- "R/evaluation/model_performance.R"
    target <- "appendix_model_predictive_performance_png"
  } else if (grepl("^figures/manuscript/fig01_winners_losers[.]pdf$", path)) {
    script <- "R/figures/winners_losers.R"
    target <- "figure01_winners_losers_pdf"
  } else if (grepl("^figures/manuscript/fig01_winners_losers[.]png$", path)) {
    script <- "R/figures/winners_losers.R"
    target <- "figure01_winners_losers_png"
  } else if (grepl("^figures/manuscript/fig02_repeated_species_consistency[.]pdf$", path)) {
    script <- "R/figures/consistency.R"
    target <- "figure02_repeated_species_consistency_pdf"
  } else if (grepl("^figures/manuscript/fig02_repeated_species_consistency[.]png$", path)) {
    script <- "R/figures/consistency.R"
    target <- "figure02_repeated_species_consistency_png"
  } else if (grepl("^figures/manuscript/appendix_temporal_winning_losing_strength_correlations[.]pdf$", path)) {
    script <- "R/figures/temporal.R"
    target <- "appendix_temporal_winning_losing_strength_correlations_file"
  } else if (grepl("^figures/manuscript/appendix_temporal_winning_losing_strength_correlations[.]png$", path)) {
    script <- "R/figures/temporal.R"
    target <- "appendix_temporal_winning_losing_strength_correlations_png"
  } else if (grepl("^figures/manuscript/fig03_temporal_category_flows[.]pdf$", path)) {
    script <- "R/figures/temporal.R"
    target <- "figure03_temporal_category_flows_file"
  } else if (grepl("^figures/manuscript/fig03_temporal_category_flows[.]png$", path)) {
    script <- "R/figures/temporal.R"
    target <- "figure03_temporal_category_flows_png"
  } else if (grepl("^figures/manuscript/appendix_destination_climate_warming_gradient[.]pdf$", path)) {
    script <- "R/figures/climate.R"
    target <- "appendix_destination_climate_warming_gradient_file"
  } else if (grepl("^figures/manuscript/appendix_destination_climate_warming_gradient[.]png$", path)) {
    script <- "R/figures/climate.R"
    target <- "appendix_destination_climate_warming_gradient_png"
  } else if (grepl("^figures/manuscript/appendix_temporal_transition_matrices[.]pdf$", path)) {
    script <- "R/figures/temporal.R"
    target <- "appendix_temporal_transition_matrices_file"
  } else if (grepl("^figures/manuscript/appendix_temporal_transition_matrices[.]png$", path)) {
    script <- "R/figures/temporal.R"
    target <- "appendix_temporal_transition_matrices_png"
  }

  data.frame(generating_script_if_known = script, target_if_known = target)
}

#' Default Figure Status.
#'
#' @param path A project-relative or absolute filesystem path.
#' @return An R object described by the function name; its structure is validated before later stages use it.
default_figure_status <- function(path) {
  registry <- figure_registry()
  registered <- registry$status[registry$filename == path]
  if (length(registered) == 1L) return(registered[[1]])
  if (grepl(
    paste0(
      "^figures/manuscript/(fig01_winners_losers[.](pdf|png)|",
      "fig02_repeated_species_consistency[.](pdf|png)|",
      "fig03_temporal_category_flows[.](pdf|png)|",
      "fig04_species_minmax_classification[.](pdf|png)|",
      "fig05_climate_only_fixed_effect_prediction_surfaces[.](pdf|png)|",
      "fig06_predictor_synthesis_support_and_heterogeneity[.](pdf|png)|",
      "fig07_predictor_response_atlas[.](pdf|png)|",
      "appendix_family_minmax_classification[.](pdf|png)|",
      "appendix_pev2_pev3_prediction_atlas[.](pdf|png)|",
      "appendix_phylogeny_pcoa_summary[.](pdf|png)|",
      "appendix_model_predictive_performance[.](pdf|png)|",
      "appendix_temporal_winning_losing_strength_correlations[.](pdf|png)|",
      "appendix_destination_climate_warming_gradient[.](pdf|png)|",
      "appendix_temporal_transition_matrices[.](pdf|png))$"
    ),
    path
  )) {
    return("FINAL_MANUSCRIPT")
  }
  if (grepl("^figures/retained/", path)) {
    return("KEEP_FOR_MANUSCRIPT")
  }
  if (grepl("^figures/manuscript/", path)) {
    return("CURRENT_ANALYSIS_FIGURE")
  }
  if (grepl("^figures/candidates/", path)) {
    return("CURRENT_ANALYSIS_FIGURE")
  }
  "UNKNOWN_DO_NOT_TOUCH"
}

#' Figure Analysis Generation.
#'
#' @param path A project-relative or absolute filesystem path.
#' @return An R object described by the function name; its structure is validated before later stages use it.
figure_analysis_generation <- function(path) {
  registry <- figure_registry()
  if (path %in% registry$filename) {
    return("registered current-workflow figure")
  }
  if (grepl(
    paste0(
      "^figures/manuscript/(fig01_winners_losers[.](pdf|png)|",
      "fig02_repeated_species_consistency[.](pdf|png)|",
      "fig03_temporal_category_flows[.](pdf|png)|",
      "fig04_species_minmax_classification[.](pdf|png)|",
      "fig05_climate_only_fixed_effect_prediction_surfaces[.](pdf|png)|",
      "fig06_predictor_synthesis_support_and_heterogeneity[.](pdf|png)|",
      "fig07_predictor_response_atlas[.](pdf|png)|",
      "appendix_family_minmax_classification[.](pdf|png)|",
      "appendix_pev2_pev3_prediction_atlas[.](pdf|png)|",
      "appendix_phylogeny_pcoa_summary[.](pdf|png)|",
      "appendix_model_predictive_performance[.](pdf|png)|",
      "appendix_temporal_winning_losing_strength_correlations[.](pdf|png)|",
      "appendix_destination_climate_warming_gradient[.](pdf|png)|",
      "appendix_temporal_transition_matrices[.](pdf|png))$"
    ),
    path
  )) {
    return("explicitly approved current manuscript figure")
  }
  if (grepl("^figures/retained/", path)) {
    return("verified retained copy")
  }
  if (grepl("^figures/manuscript/", path)) {
    return("current current targets analysis")
  }
  if (grepl("^figures/candidates/", path)) {
    return("current current targets figure candidate")
  }
  "unknown"
}

#' Build Figure Provenance.
#'
#' @param root A project-relative or absolute filesystem path.
#' @param retained_config A validated workflow configuration or specification object.
#' @return An R object described by the function name; its structure is validated before later stages use it.
build_figure_provenance <- function(root = getwd(), retained_config = "config/retained_figures.csv") {
  candidates <- c(
    list.files(file.path(root, "figures"), pattern = "\\.(pdf|png|svg)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  )
  candidates <- sort(unique(candidates[file.exists(candidates)]))
  paths <- relative_project_path(candidates, root)

  retained <- if (file.exists(file.path(root, retained_config))) {
    read.csv(file.path(root, retained_config), stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    data.frame(figure_path = character(), expected_sha256 = character(), status = character(), reason = character())
  }

  rows <- lapply(seq_along(candidates), function(index) {
    path <- paths[[index]]
    dimensions <- figure_dimensions(candidates[[index]])
    producer <- known_figure_producer(path)
    override <- retained[retained$figure_path == path, , drop = FALSE]
    if (!nrow(override) && grepl("^figures/retained/", path)) {
      override <- retained[
        basename(retained$figure_path) == basename(path), , drop = FALSE
      ]
    }
    status <- if (nrow(override)) override$status[[1]] else default_figure_status(path)
    hash <- sha256_file(candidates[[index]])
    expected <- if (nrow(override)) override$expected_sha256[[1]] else NA_character_
    hash_ok <- if (is.na(expected) || !nzchar(expected)) NA else identical(hash, expected)
    proposed_destination <- if (identical(status, "KEEP_FOR_MANUSCRIPT") || identical(status, "KEEP_HISTORICAL")) {
      file.path("figures", "retained", basename(path))
    } else if (identical(status, "FINAL_MANUSCRIPT")) {
      path
    } else if (identical(status, "CURRENT_ANALYSIS_FIGURE")) {
      file.path("figures", "candidates", basename(path))
    } else {
      NA_character_
    }

    data.frame(
      figure_path = path,
      figure_name = basename(path),
      file_hash_sha256 = hash,
      expected_hash_sha256 = expected,
      expected_hash_matches = hash_ok,
      format = toupper(tools::file_ext(path)),
      width = dimensions$width,
      height = dimensions$height,
      dimension_unit = dimensions$dimension_unit,
      generating_script_if_known = producer$generating_script_if_known,
      target_if_known = producer$target_if_known,
      source_code_available = !is.na(producer$generating_script_if_known) &&
        file.exists(file.path(root, producer$generating_script_if_known)),
      scientific_purpose = if (nrow(override)) override$reason[[1]] else tools::file_path_sans_ext(gsub("[_-]+", " ", basename(path))),
      current_vs_old_analysis = figure_analysis_generation(path),
      manually_retained_possible = status %in% c("KEEP_FOR_MANUSCRIPT", "KEEP_HISTORICAL", "UNKNOWN_DO_NOT_TOUCH"),
      proposed_status = status,
      proposed_destination = proposed_destination,
      refinement_note = NA_character_,
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  result$duplicate_stem <- tolower(tools::file_path_sans_ext(result$figure_name))
  result$duplicate_format_relationships <- vapply(
    seq_len(nrow(result)),
    function(index) {
      peers <- result$figure_path[
        result$duplicate_stem == result$duplicate_stem[[index]] &
          result$figure_path != result$figure_path[[index]]
      ]
      if (length(peers)) paste(peers, collapse = "; ") else NA_character_
    },
    character(1)
  )
  result
}

#' Audit Retained Figure Copies.
#'
#' @param root A project-relative or absolute filesystem path.
#' @param retained_config A validated workflow configuration or specification object.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
audit_retained_figure_copies <- function(
    root = getwd(), retained_config = "config/retained_figures.csv") {
  retained <- read.csv(
    file.path(root, retained_config), stringsAsFactors = FALSE, check.names = FALSE
  )
  rows <- lapply(seq_len(nrow(retained)), function(index) {
    source_relative <- retained$figure_path[[index]]
    destination_relative <- file.path(
      "figures", "retained", basename(source_relative)
    )
    source <- file.path(root, source_relative)
    destination <- file.path(root, destination_relative)
    expected <- retained$expected_sha256[[index]]
    source_hash <- sha256_file(source)
    destination_hash <- sha256_file(destination)
    data.frame(
      source_path = source_relative,
      destination_path = destination_relative,
      source_sha256 = source_hash,
      destination_sha256 = destination_hash,
      expected_sha256 = expected,
      source_matches_expected = identical(source_hash, expected),
      copy_matches_source = identical(destination_hash, source_hash),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  if (any(!result$source_matches_expected | !result$copy_matches_source)) {
    stop("One or more retained figure copies failed checksum validation.",
         call. = FALSE)
  }
  result
}

#' Retained Figure Release Check.
#'
#' @param audit A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
retained_figure_release_check <- function(audit) {
  required <- c(
    "source_path", "destination_path", "source_matches_expected",
    "copy_matches_source"
  )
  if (!is.data.frame(audit) || !all(required %in% names(audit))) {
    stop("Retained-figure audit lacks its checksum validation columns.",
         call. = FALSE)
  }
  tibble::tibble(
    check = paste0("retained figure: ", audit$source_path),
    pass = as.logical(
      audit$source_matches_expected & audit$copy_matches_source
    ),
    detail = paste(audit$source_path, "->", audit$destination_path)
  )
}

#' Copy Retained Figures Verified.
#'
#' @param root A project-relative or absolute filesystem path.
#' @param retained_config A validated workflow configuration or specification object.
#' @return The path or paths written by the declared side effect.
copy_retained_figures_verified <- function(
    root = getwd(), retained_config = "config/retained_figures.csv") {
  retained <- read.csv(
    file.path(root, retained_config), stringsAsFactors = FALSE, check.names = FALSE
  )
  destinations <- vapply(seq_len(nrow(retained)), function(index) {
    source <- file.path(root, retained$figure_path[[index]])
    destination <- file.path(
      root, "figures", "retained", basename(retained$figure_path[[index]])
    )
    expected <- retained$expected_sha256[[index]]
    if (!identical(sha256_file(source), expected)) {
      stop("Retained source figure hash mismatch: ", source, call. = FALSE)
    }
    if (file.exists(destination) &&
        !identical(sha256_file(destination), expected)) {
      stop("Existing retained destination has unexpected content: ",
           destination, call. = FALSE)
    }
    if (!file.exists(destination)) {
      dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
      if (!file.copy(source, destination, overwrite = FALSE, copy.mode = TRUE)) {
        stop("Could not copy retained figure to ", destination, call. = FALSE)
      }
    }
    if (!identical(sha256_file(destination), expected)) {
      stop("Retained figure copy verification failed: ", destination,
           call. = FALSE)
    }
    normalizePath(destination, winslash = "/", mustWork = TRUE)
  }, character(1))
  unname(destinations)
}

#' Build Artifact Manifest.
#'
#' @param paths A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param root A project-relative or absolute filesystem path.
#' @param lineage A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
build_artifact_manifest <- function(paths, root = getwd(),
                                    lineage = "current targets") {
  producer <- names(paths)
  if (is.null(producer)) producer <- rep(NA_character_, length(paths))
  normalized <- normalizePath(paths, winslash = "/", mustWork = FALSE)
  if (anyDuplicated(normalized)) {
    stop("Artifact manifest contains duplicated file paths.", call. = FALSE)
  }
  exists <- file.exists(normalized) & !dir.exists(normalized)
  if (any(!exists)) {
    stop("Artifact manifest is missing: ",
         paste(paths[!exists], collapse = ", "), call. = FALSE)
  }
  information <- file.info(normalized)
  data.frame(
    artifact = unname(producer),
    path = relative_project_path(normalized, root),
    size_bytes = as.numeric(information$size),
    modified_utc = format(information$mtime, tz = "UTC", usetz = TRUE),
    sha256 = vapply(normalized, sha256_file, character(1)),
    lineage = lineage,
    status = "CURRENT_RELEASE_OUTPUT",
    stringsAsFactors = FALSE
  )
}

#' Build Scientific Source Inventory.
#'
#' @param paths A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param root A project-relative or absolute filesystem path.
#' @return An R object described by the function name; its structure is validated before later stages use it.
build_scientific_source_inventory <- function(paths, root = getwd()) {
  normalized <- normalizePath(paths, winslash = "/", mustWork = FALSE)
  if (anyDuplicated(normalized)) {
    stop("Scientific source inventory contains duplicated paths.", call. = FALSE)
  }
  exists <- file.exists(normalized) & !dir.exists(normalized)
  if (any(!exists)) {
    stop(
      "Scientific source inventory is missing: ",
      paste(paths[!exists], collapse = ", "), call. = FALSE
    )
  }
  information <- file.info(normalized)
  relative <- relative_project_path(normalized, root)
  data.frame(
    path = relative,
    source_role = dplyr::case_when(
      relative == "_targets.R" ~ "target_entry_point",
      startsWith(relative, "R/") & !startsWith(relative, "R/targets/") ~
        "active_function",
      startsWith(relative, "R/targets/") ~ "target_stage",
      startsWith(relative, "config/") ~ "configuration_or_decision_ledger",
      relative == "analysis_workflow.qmd" ~ "analysis_and_report_source",
      startsWith(relative, "scripts/") ~ "command_line_entry_or_maintenance",
      startsWith(relative, "tests/") ~ "test_or_regression_fixture",
      TRUE ~ "project_documentation_or_environment"
    ),
    size_bytes = as.numeric(information$size),
    sha256 = vapply(normalized, sha256_file, character(1)),
    stringsAsFactors = FALSE
  )
}

#' Phylogeny Package Provenance.
#'
#' @return An R object described by the function name; its structure is validated before later stages use it.
phylogeny_package_provenance <- function() {
  if (!requireNamespace("V.PhyloMaker2", quietly = TRUE)) {
    return(data.frame(package = "V.PhyloMaker2", installed = FALSE))
  }
  description <- utils::packageDescription("V.PhyloMaker2")
  fields <- c(
    "Package", "Version", "SourceURL", "RemoteType", "RemoteHost",
    "RemoteRepo", "RemoteUsername", "RemoteRef", "RemoteSha", "GithubRepo",
    "GithubUsername", "GithubRef", "GithubSHA1"
  )
  values <- vapply(fields, function(field) {
    if (identical(field, "SourceURL")) {
      return("https://github.com/jinyizju/V.PhyloMaker2")
    }
    value <- description[[field]]
    if (is.null(value)) NA_character_ else as.character(value)
  }, character(1))
  data.frame(package = "V.PhyloMaker2", installed = TRUE, field = fields, value = unname(values), stringsAsFactors = FALSE)
}

#' Validate Phylogeny Package Provenance.
#'
#' @param provenance A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param expected_version A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param expected_sha A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param expected_source A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
validate_phylogeny_package_provenance <- function(
    provenance,
    expected_version = "0.1.0",
    expected_sha = "7af3fb5152f691af2e4ec9d5e2e467d1b50505e9",
    expected_source = "https://github.com/jinyizju/V.PhyloMaker2") {
  values <- stats::setNames(provenance$value, provenance$field)
  checks <- data.frame(
    check = c(
      "V.PhyloMaker2 installed", "V.PhyloMaker2 version",
      "V.PhyloMaker2 Git source", "V.PhyloMaker2 Git commit"
    ),
    pass = c(
      nrow(provenance) > 0L && all(provenance$installed),
      identical(unname(values[["Version"]]), expected_version),
      identical(unname(values[["SourceURL"]]), expected_source) &&
        identical(unname(values[["RemoteUsername"]]), "jinyizju") &&
        identical(unname(values[["RemoteRepo"]]), "V.PhyloMaker2"),
      identical(unname(values[["RemoteSha"]]), expected_sha)
    ),
    detail = c(
      "Required package is available.",
      paste("expected", expected_version, "observed", values[["Version"]]),
      paste("expected", expected_source, "observed", values[["SourceURL"]]),
      paste("expected", expected_sha, "observed", values[["RemoteSha"]])
    ),
    stringsAsFactors = FALSE
  )
  if (any(!checks$pass)) {
    stop(
      "V.PhyloMaker2 provenance failed: ",
      paste(checks$check[!checks$pass], collapse = "; "), call. = FALSE
    )
  }
  checks
}
