# Step 19: resident phylogeny and response-independent PCoA -----------------
#
# This isolated step constructs phylogenetic information for the updated
# maximum-duration resident analysis set: species are recovered from the saved
# M_both object, then the declared downstream-only exclusions are applied. It does not read
# the AUC response, fit an AUC model, select axes using a response, or alter the
# project's TNRS-resolved ecological taxonomy.

source("config/analysis.R")
source("R/pipeline/utils.R")
source("R/pipeline/frequentist_hierarchical_auc.R")

required_packages <- c(
  "ape", "dplyr", "ggplot2", "ggrepel", "here", "readr", "tibble", "tidyr",
  "V.PhyloMaker2"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing packages: ", paste(missing_packages, collapse = ", "),
    ". V.PhyloMaker2 is available from github.com/jinyizju/V.PhyloMaker2.",
    call. = FALSE
  )
}

output_dir <- file.path(
  "output", "step_by_step_modeling", "resident_phylogeny_pcoa"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_files <- c(
  taxonomy_audit = file.path(output_dir, "tnrs_taxonomy_audit.csv"),
  crosswalk = file.path(output_dir, "tnrs_backbone_crosswalk.csv"),
  infraspecific_audit = file.path(
    output_dir, "infraspecific_placement_audit.csv"
  ),
  disagreements = file.path(output_dir, "taxonomy_disagreement_resolution.csv"),
  problematic = file.path(output_dir, "unresolved_problematic_taxa.csv"),
  tree_rds = file.path(output_dir, "resident_phylogeny.rds"),
  tree_newick = file.path(output_dir, "resident_phylogeny.newick"),
  distances = file.path(output_dir, "phylogenetic_distance_matrix.rds"),
  species_order = file.path(output_dir, "phylogeny_species_order.csv"),
  eigenvalues = file.path(output_dir, "pcoa_eigenvalue_summary.csv"),
  scores = file.path(output_dir, "species_pev_scores.csv"),
  extremes = file.path(output_dir, "pcoa_axis_extreme_species.csv"),
  tree_figure = file.path(output_dir, "resident_phylogenetic_tree.pdf"),
  scree_figure = file.path(output_dir, "pcoa_scree_broken_stick.pdf"),
  cumulative_figure = file.path(output_dir, "pcoa_cumulative_variation.pdf"),
  pev12_figure = file.path(output_dir, "pev1_vs_pev2.pdf"),
  pev13_figure = file.path(output_dir, "pev1_vs_pev3.pdf"),
  report = file.path(output_dir, "phylogeny_pcoa_audit_report.txt")
)

normalize_species_name <- function(x) {
  x <- gsub("_", " ", as.character(x), fixed = TRUE)
  trimws(gsub("[[:space:]]+", " ", x))
}

underscore_species_name <- function(x) {
  gsub(" ", "_", normalize_species_name(x), fixed = TRUE)
}

collapse_notes <- function(...) {
  values <- unlist(list(...), use.names = FALSE)
  values <- unique(values[!is.na(values) & nzchar(values)])
  if (length(values)) paste(values, collapse = "; ") else ""
}

write_csv_checked <- function(x, path) {
  readr::write_csv(x, path, na = "")
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Failed to write ", path, call. = FALSE)
  }
  invisible(path)
}

# -------------------------------------------------------------------------
# 1. Updated downstream species set and canonical TNRS taxonomy
# -------------------------------------------------------------------------

# Read the saved current M_both only to recover species identifiers, then apply
# the downstream scope filter. No response or fitted value is retained or used.
m_both_file <- file.path(
  "output", "step_by_step_modeling", "species_both",
  "auc_change_species_both_model.rds"
)
if (!file.exists(m_both_file)) {
  stop("The saved current M_both model is required: ", m_both_file,
       call. = FALSE)
}
m_both_object <- readRDS(m_both_file)
if (!inherits(m_both_object, "rma.mv") || is.null(m_both_object$data) ||
    !"species_id" %in% names(m_both_object$data)) {
  stop("The saved M_both object does not contain its analysis species_id.",
       call. = FALSE)
}
m_both_species <- m_both_object$data |>
  dplyr::transmute(species_id = as.character(.data$species_id)) |>
  filter_downstream_model_species("species_id") |>
  dplyr::distinct() |>
  dplyr::arrange(.data$species_id)
rm(m_both_object)

if (any(m_both_species$species_id %in% resident_downstream_model_exclusions())) {
  stop("A prespecified downstream-excluded taxon remains in the phylogeny set.",
       call. = FALSE)
}

project_taxonomy <- readr::read_csv(
  file.path("output", "tables", "analysis_dataset.csv"),
  show_col_types = FALSE
) |>
  dplyr::transmute(
    species_id = normalize_species_name(.data$accepted_name),
    tnrs_accepted_name = normalize_species_name(.data$accepted_name),
    tnrs_genus = trimws(as.character(.data$genus)),
    tnrs_family = trimws(as.character(.data$family))
  ) |>
  dplyr::distinct()

if (anyDuplicated(project_taxonomy$species_id)) {
  duplicate_rows <- project_taxonomy |>
    dplyr::count(.data$species_id, name = "n_taxonomy_rows") |>
    dplyr::filter(.data$n_taxonomy_rows > 1L)
  stop(
    "Canonical project taxonomy is non-unique for: ",
    paste(duplicate_rows$species_id, collapse = ", "), call. = FALSE
  )
}

taxonomy <- m_both_species |>
  dplyr::left_join(project_taxonomy, by = "species_id")
if (nrow(taxonomy) != nrow(m_both_species) ||
    any(is.na(taxonomy$tnrs_accepted_name))) {
  stop("The exact M_both species set did not map one-to-one to project TNRS taxonomy.",
       call. = FALSE)
}

binomial_pattern <- "^[A-Z][[:alpha:]-]+ [a-z][[:alpha:]-]+$"
infraspecific_pattern <- "\\b(subsp|ssp|var|f)\\.|\\bforma\\b"
authorship_pattern <- "[(),;]|\\b(ex|sensu|auct\\.)\\b|&|[[:space:]]{2,}"

taxonomy_audit <- taxonomy |>
  dplyr::mutate(
    duplicated_species_id = duplicated(.data$species_id) |
      duplicated(.data$species_id, fromLast = TRUE),
    duplicated_tnrs_accepted_name =
      duplicated(.data$tnrs_accepted_name) |
      duplicated(.data$tnrs_accepted_name, fromLast = TRUE),
    missing_genus = is.na(.data$tnrs_genus) | .data$tnrs_genus == "",
    missing_family = is.na(.data$tnrs_family) | .data$tnrs_family == "",
    infraspecific_name = grepl(
      infraspecific_pattern, .data$tnrs_accepted_name, perl = TRUE
    ),
    parent_binomial = dplyr::if_else(
      .data$infraspecific_name,
      vapply(strsplit(.data$tnrs_accepted_name, " ", fixed = TRUE),
             function(parts) paste(parts[seq_len(min(2L, length(parts)))],
                                   collapse = " "), character(1)),
      NA_character_
    ),
    authorship_or_formatting_artifact = grepl(
      authorship_pattern, .data$tnrs_accepted_name, perl = TRUE
    ) | grepl("_", .data$tnrs_accepted_name, fixed = TRUE),
    straightforward_binomial = grepl(
      binomial_pattern, .data$tnrs_accepted_name, perl = TRUE
    ),
    taxonomy_problem = .data$duplicated_species_id |
      .data$duplicated_tnrs_accepted_name | .data$missing_genus |
      .data$missing_family | .data$authorship_or_formatting_artifact |
      !.data$straightforward_binomial,
    audit_note = vapply(seq_len(dplyr::n()), function(i) {
      collapse_notes(
        if (duplicated_species_id[[i]]) "duplicated species_id" else "",
        if (duplicated_tnrs_accepted_name[[i]]) {
          "duplicated TNRS accepted name"
        } else "",
        if (missing_genus[[i]]) "missing genus" else "",
        if (missing_family[[i]]) "missing family" else "",
        if (infraspecific_name[[i]]) "infraspecific name" else "",
        if (authorship_or_formatting_artifact[[i]]) {
          "authorship or formatting artifact"
        } else "",
        if (!straightforward_binomial[[i]]) "not a straightforward binomial" else ""
      )
    }, character(1))
  )
write_csv_checked(taxonomy_audit, output_files[["taxonomy_audit"]])

if (any(taxonomy_audit$duplicated_species_id) ||
    any(taxonomy_audit$duplicated_tnrs_accepted_name) ||
    any(taxonomy_audit$missing_genus) || any(taxonomy_audit$missing_family)) {
  stop(
    "Canonical taxonomy has duplicate identifiers/names or missing genus/family; ",
    "the audit was saved and no phylogeny was built.", call. = FALSE
  )
}

# -------------------------------------------------------------------------
# 2. TNRS-to-LCVP crosswalk and deterministic scenario-3 placement
# -------------------------------------------------------------------------

backbone_tips <- V.PhyloMaker2::tips.info.LCVP |>
  dplyr::transmute(
    backbone_tip_label = as.character(.data$species),
    backbone_name = normalize_species_name(.data$species),
    backbone_genus = trimws(as.character(.data$genus)),
    backbone_family = trimws(as.character(.data$family))
  ) |>
  dplyr::distinct()
if (anyDuplicated(backbone_tips$backbone_name)) {
  stop("The installed LCVP backbone tip-name table is not unique.", call. = FALSE)
}

# The optional alternative-name candidates come only from the project's
# existing TNRS review. This is not a second name-resolution operation.
tnrs_review <- readr::read_csv(
  file.path("output", "tables", "tnrs_vegetation_review.csv"),
  show_col_types = FALSE
)
review_candidates <- tnrs_review |>
  dplyr::transmute(
    species_id = normalize_species_name(.data$accepted_name_reviewed),
    accepted_candidate = normalize_species_name(.data$accepted_name),
    matched_candidate = normalize_species_name(.data$matched_name),
    tnrs_score = suppressWarnings(as.numeric(.data$score))
  ) |>
  tidyr::pivot_longer(
    cols = c("accepted_candidate", "matched_candidate"),
    names_to = "candidate_source", values_to = "candidate_name"
  ) |>
  dplyr::filter(
    !is.na(.data$species_id), nzchar(.data$species_id),
    !is.na(.data$candidate_name), nzchar(.data$candidate_name)
  ) |>
  dplyr::inner_join(
    backbone_tips |>
      dplyr::select("backbone_name"),
    by = c("candidate_name" = "backbone_name")
  ) |>
  dplyr::mutate(
    source_priority = dplyr::if_else(
      .data$candidate_source == "accepted_candidate", 1L, 2L
    ),
    tnrs_score_sort = dplyr::coalesce(.data$tnrs_score, -Inf)
  ) |>
  dplyr::arrange(
    .data$species_id, .data$source_priority,
    dplyr::desc(.data$tnrs_score_sort), .data$candidate_name
  ) |>
  dplyr::group_by(.data$species_id) |>
  dplyr::slice(1L) |>
  dplyr::ungroup() |>
  dplyr::select(
    "species_id", alias_backbone_name = "candidate_name",
    "candidate_source"
  )

crosswalk_input <- taxonomy |>
  dplyr::left_join(
    taxonomy_audit |>
      dplyr::select("species_id", "infraspecific_name", "parent_binomial"),
    by = "species_id"
  ) |>
  dplyr::left_join(
    backbone_tips |>
      dplyr::select(
        tnrs_accepted_name = "backbone_name",
        exact_backbone_name = "backbone_name"
    ),
    by = "tnrs_accepted_name"
  ) |>
  dplyr::left_join(
    backbone_tips |>
      dplyr::select(
        parent_binomial = "backbone_name",
        parent_backbone_name = "backbone_name"
      ),
    by = "parent_binomial"
  ) |>
  dplyr::left_join(review_candidates, by = "species_id") |>
  dplyr::mutate(
    parent_in_backbone = .data$infraspecific_name &
      !is.na(.data$parent_backbone_name),
    pre_match_type = dplyr::case_when(
      .data$parent_in_backbone ~
        "infraspecific -> parent species placement",
      !is.na(.data$exact_backbone_name) ~ "exact species match",
      !is.na(.data$alias_backbone_name) ~ "synonym/backbone-name match",
      TRUE ~ "no species-level backbone match"
    ),
    chosen_backbone_name = dplyr::coalesce(
      .data$parent_backbone_name, .data$exact_backbone_name,
      .data$alias_backbone_name
    )
  )

# Never silently merge analysis species that happen to resolve to one backbone
# tip. For duplicated non-exact alternatives, retain each TNRS name and let S3
# make an auditable graft instead.
duplicate_backbone_names <- crosswalk_input |>
  dplyr::filter(!is.na(.data$chosen_backbone_name)) |>
  dplyr::count(.data$chosen_backbone_name, name = "n_analysis_species") |>
  dplyr::filter(.data$n_analysis_species > 1L) |>
  dplyr::pull(.data$chosen_backbone_name)

crosswalk_input <- crosswalk_input |>
  dplyr::mutate(
    duplicate_backbone_mapping = !is.na(.data$chosen_backbone_name) &
      .data$chosen_backbone_name %in% duplicate_backbone_names,
    alias_suppressed_to_avoid_merge =
      .data$duplicate_backbone_mapping &
      .data$pre_match_type == "synonym/backbone-name match",
    phylo_input_name = dplyr::if_else(
      .data$alias_suppressed_to_avoid_merge,
      .data$tnrs_accepted_name,
      dplyr::coalesce(.data$chosen_backbone_name, .data$tnrs_accepted_name)
    ),
    phylo_input_label = underscore_species_name(.data$phylo_input_name),
    shared_parent_placement = .data$parent_in_backbone &
      .data$duplicate_backbone_mapping,
    lcvp_synonym_found = !is.na(.data$alias_backbone_name) &
      .data$alias_backbone_name != .data$parent_binomial
  ) |>
  dplyr::group_by(.data$phylo_input_label) |>
  dplyr::mutate(
    shared_parent_placement = dplyr::n() > 1L &
      any(.data$parent_in_backbone)
  ) |>
  dplyr::ungroup()

invalid_duplicate_placements <- crosswalk_input |>
  dplyr::summarise(
    n_placement_taxa = dplyr::n(),
    any_parent_placement = any(.data$parent_in_backbone),
    .by = "phylo_input_label"
  ) |>
  dplyr::filter(.data$n_placement_taxa > 1L, !.data$any_parent_placement)
if (nrow(invalid_duplicate_placements)) {
  stop(
    "A duplicated phylogeny placement is not justified by an infraspecific ",
    "parent-species mapping. Inspect the taxonomy audit.", call. = FALSE
  )
}

sp_list <- crosswalk_input |>
  dplyr::transmute(
    species = .data$phylo_input_label,
    genus = .data$tnrs_genus,
    family = .data$tnrs_family
  ) |>
  dplyr::distinct(.data$species, .keep_all = TRUE) |>
  as.data.frame(stringsAsFactors = FALSE)

phylo_warnings <- character()
phylo_messages <- character()
phylo_result <- withCallingHandlers(
  V.PhyloMaker2::phylo.maker(
    sp.list = sp_list,
    tree = V.PhyloMaker2::GBOTB.extended.LCVP,
    nodes = V.PhyloMaker2::nodes.info.1.LCVP,
    output.sp.list = TRUE,
    output.tree = FALSE,
    scenarios = "S3"
  ),
  warning = function(condition) {
    phylo_warnings <<- c(phylo_warnings, conditionMessage(condition))
    invokeRestart("muffleWarning")
  },
  message = function(condition) {
    phylo_messages <<- c(phylo_messages, conditionMessage(condition))
    invokeRestart("muffleMessage")
  }
)

if (!all(c("scenario.3", "species.list") %in% names(phylo_result)) ||
    !inherits(phylo_result$scenario.3, "phylo")) {
  stop("V.PhyloMaker2 did not return the documented scenario.3 tree.",
       call. = FALSE)
}

placement <- tibble::as_tibble(phylo_result$species.list) |>
  dplyr::transmute(
    phylo_input_label = as.character(.data$species),
    returned_genus = as.character(.data$genus),
    returned_family = as.character(.data$family),
    backbone_status = as.character(.data$status)
  )

genus_lookup <- backbone_tips |>
  dplyr::filter(!is.na(.data$backbone_genus), nzchar(.data$backbone_genus)) |>
  dplyr::group_by(.data$backbone_genus) |>
  dplyr::summarise(
    genus_backbone_family = paste(
      sort(unique(.data$backbone_family)), collapse = " | "
    ),
    genus_in_backbone = TRUE,
    .groups = "drop"
  )
family_lookup <- backbone_tips |>
  dplyr::filter(!is.na(.data$backbone_family), nzchar(.data$backbone_family)) |>
  dplyr::distinct(backbone_family = .data$backbone_family) |>
  dplyr::mutate(family_in_backbone = TRUE)

crosswalk <- crosswalk_input |>
  dplyr::left_join(placement, by = "phylo_input_label") |>
  dplyr::left_join(
    backbone_tips,
    by = c("phylo_input_name" = "backbone_name")
  ) |>
  dplyr::left_join(
    genus_lookup, by = c("tnrs_genus" = "backbone_genus")
  ) |>
  dplyr::left_join(
    family_lookup, by = c("tnrs_family" = "backbone_family")
  ) |>
  dplyr::mutate(
    genus_in_backbone = dplyr::coalesce(.data$genus_in_backbone, FALSE),
    family_in_backbone = dplyr::coalesce(.data$family_in_backbone, FALSE),
    match_type = dplyr::case_when(
      .data$backbone_status == "prune" & .data$parent_in_backbone ~
        "infraspecific -> parent species placement",
      .data$backbone_status == "prune" &
        .data$pre_match_type == "exact species match" ~ "exact species match",
      .data$backbone_status == "prune" &
        .data$pre_match_type == "synonym/backbone-name match" ~
        "synonym/backbone-name match",
      .data$backbone_status == "bind" & .data$genus_in_backbone ~
        "grafted within genus",
      .data$backbone_status == "bind" & .data$family_in_backbone ~
        "grafted within family",
      TRUE ~ "unresolved"
    ),
    backbone_matched_name = dplyr::if_else(
      .data$match_type %in% c(
        "exact species match", "synonym/backbone-name match",
        "infraspecific -> parent species placement"
      ),
      .data$phylo_input_name, NA_character_
    ),
    backbone_genus_reported = dplyr::case_when(
      .data$match_type %in% c(
        "exact species match", "synonym/backbone-name match",
        "infraspecific -> parent species placement"
      ) ~ .data$backbone_genus,
      .data$match_type == "grafted within genus" ~ .data$tnrs_genus,
      TRUE ~ NA_character_
    ),
    backbone_family_reported = dplyr::case_when(
      .data$match_type %in% c(
        "exact species match", "synonym/backbone-name match",
        "infraspecific -> parent species placement"
      ) ~ .data$backbone_family,
      .data$match_type == "grafted within genus" ~
        .data$genus_backbone_family,
      .data$match_type == "grafted within family" ~ .data$tnrs_family,
      TRUE ~ NA_character_
    ),
    genus_disagreement = !is.na(.data$backbone_genus_reported) &
      nzchar(.data$backbone_genus_reported) &
      .data$tnrs_genus != .data$backbone_genus_reported,
    family_disagreement = !is.na(.data$backbone_family_reported) &
      nzchar(.data$backbone_family_reported) &
      !vapply(seq_len(dplyr::n()), function(i) {
        tnrs_family[[i]] %in% strsplit(
          backbone_family_reported[[i]], " \\| ", fixed = FALSE
        )[[1]]
      }, logical(1)),
    any_taxonomy_disagreement = .data$genus_disagreement |
      .data$family_disagreement,
    resolution_status = dplyr::case_when(
      !.data$any_taxonomy_disagreement ~ "not applicable",
      .data$match_type == "unresolved" ~ "outstanding",
      TRUE ~ "resolved"
    ),
    resolution_decision = dplyr::case_when(
      !.data$any_taxonomy_disagreement ~ "no disagreement",
      .data$match_type == "exact species match" ~
        "use exact LCVP species placement; retain TNRS ecology ID and taxonomy fields",
      .data$match_type == "synonym/backbone-name match" ~
        "use LCVP synonym/backbone nomenclature for placement; retain TNRS ecology ID",
      .data$match_type == "infraspecific -> parent species placement" ~
        "use LCVP parent species for placement; retain infraspecific ecology ID",
      .data$match_type %in% c("grafted within genus", "grafted within family") ~
        "retain TNRS genus/family for deterministic S3 placement",
      TRUE ~ "no defensible placement decision"
    ),
    resolution_reason = dplyr::case_when(
      !.data$any_taxonomy_disagreement ~ "TNRS and placement taxonomy agree",
      .data$match_type == "exact species match" ~
        "the accepted binomial is an exact backbone tip; disagreement is confined to TNRS versus LCVP classification fields",
      .data$match_type == "synonym/backbone-name match" ~
        "nomenclatural synonym or accepted-name difference documented by the existing TNRS review",
      .data$match_type == "infraspecific -> parent species placement" ~
        "the backbone contains the known parent binomial but not infraspecific divergence",
      .data$match_type == "grafted within genus" ~
        "species absent from backbone; TNRS genus exists in backbone",
      .data$match_type == "grafted within family" ~
        "species and TNRS genus absent from backbone; TNRS family exists",
      TRUE ~ "placement unavailable"
    ),
    notes = vapply(seq_len(dplyr::n()), function(i) {
      collapse_notes(
        if (match_type[[i]] == "grafted within genus") {
          paste0("S3 placement within TNRS genus ", tnrs_genus[[i]])
        } else "",
        if (match_type[[i]] == "grafted within family") {
          paste0("S3 placement within TNRS family ", tnrs_family[[i]])
        } else "",
        if (match_type[[i]] == "unresolved") "not represented by S3" else "",
        if (match_type[[i]] == "synonym/backbone-name match") {
          paste0(
            "alternative backbone name taken from existing TNRS review field ",
            candidate_source[[i]]
          )
        } else "",
        if (match_type[[i]] == "infraspecific -> parent species placement") {
          paste0("placed at LCVP parent species ", parent_binomial[[i]])
        } else "",
        if (alias_suppressed_to_avoid_merge[[i]]) {
          "alternative backbone name suppressed to prevent silent tip merge"
        } else "",
        if (duplicate_backbone_mapping[[i]]) {
          "multiple ecological taxa initially mapped to this backbone taxon"
        } else "",
        if (genus_disagreement[[i]]) "TNRS/backbone genus disagreement" else "",
        if (family_disagreement[[i]]) "TNRS/backbone family disagreement" else ""
      )
    }, character(1))
  ) |>
  dplyr::select(
    "species_id",
    `TNRS accepted name` = "tnrs_accepted_name",
    `TNRS genus` = "tnrs_genus",
    `TNRS family` = "tnrs_family",
    `backbone matched name` = "backbone_matched_name",
    `backbone genus` = "backbone_genus_reported",
    `backbone family` = "backbone_family_reported",
    "backbone_status", "match_type", "phylo_input_label",
    "infraspecific_name", "parent_binomial", "parent_in_backbone",
    "shared_parent_placement", "lcvp_synonym_found",
    "duplicate_backbone_mapping", "genus_disagreement",
    "family_disagreement", "any_taxonomy_disagreement",
    "resolution_status", "resolution_decision", "resolution_reason", "notes"
  )

write_csv_checked(crosswalk, output_files[["crosswalk"]])

infraspecific_placement_audit <- crosswalk |>
  dplyr::filter(.data$infraspecific_name) |>
  dplyr::mutate(
    final_phylogenetic_placement = dplyr::case_when(
      .data$match_type == "infraspecific -> parent species placement" ~
        .data$parent_binomial,
      .data$match_type == "grafted within genus" ~
        paste0("S3 graft within genus ", .data$`TNRS genus`),
      .data$match_type == "grafted within family" ~
        paste0("S3 graft within family ", .data$`TNRS family`),
      !is.na(.data$`backbone matched name`) ~ .data$`backbone matched name`,
      TRUE ~ "unresolved"
    ),
    reason = dplyr::case_when(
      .data$parent_in_backbone ~ paste0(
        "Parent binomial ", .data$parent_binomial,
        " is a backbone tip; no infraspecific divergence was invented."
      ),
      .data$match_type == "grafted within genus" ~ paste0(
        "Full name and parent binomial are absent as backbone tips; ",
        "deterministic S3 placement within TNRS genus ", .data$`TNRS genus`, "."
      ),
      .data$match_type == "grafted within family" ~ paste0(
        "Full name, parent binomial, and genus are absent as backbone tips; ",
        "deterministic S3 placement within TNRS family ",
        .data$`TNRS family`, "."
      ),
      TRUE ~ .data$notes
    )
  ) |>
  dplyr::transmute(
    `ecological species_id` = .data$species_id,
    `full TNRS name` = .data$`TNRS accepted name`,
    `parent binomial` = .data$parent_binomial,
    `Parent species tip present in GBOTB.extended.LCVP?` =
      dplyr::if_else(.data$parent_in_backbone, "Yes", "No"),
    `LCVP synonym found?` =
      dplyr::if_else(.data$lcvp_synonym_found, "Yes", "No"),
    `final phylogenetic placement` = .data$final_phylogenetic_placement,
    reason = .data$reason
  ) |>
  dplyr::arrange(.data$`ecological species_id`)
write_csv_checked(
  infraspecific_placement_audit, output_files[["infraspecific_audit"]]
)

disagreement_resolution <- crosswalk |>
  dplyr::filter(.data$any_taxonomy_disagreement) |>
  dplyr::select(
    "species_id", "TNRS accepted name", "TNRS genus", "TNRS family",
    "backbone matched name", "backbone genus", "backbone family",
    "match_type", "genus_disagreement", "family_disagreement",
    "resolution_status", "resolution_decision", "resolution_reason"
  )
write_csv_checked(disagreement_resolution, output_files[["disagreements"]])

problematic <- taxonomy_audit |>
  dplyr::select(
    "species_id", "tnrs_accepted_name", "tnrs_genus", "tnrs_family",
    dplyr::starts_with("duplicated_"), "missing_genus", "missing_family",
    "infraspecific_name", "parent_binomial", "authorship_or_formatting_artifact",
    "straightforward_binomial", "taxonomy_problem",
    taxonomy_audit_note = "audit_note"
  ) |>
  dplyr::left_join(
    crosswalk |>
      dplyr::select(
        "species_id", "backbone_status", "match_type",
        "parent_in_backbone", "shared_parent_placement",
        "duplicate_backbone_mapping", "genus_disagreement",
        "family_disagreement", "resolution_status", "resolution_decision",
        "resolution_reason", backbone_note = "notes"
      ),
    by = "species_id"
  ) |>
  dplyr::filter(
    .data$taxonomy_problem | .data$match_type == "unresolved" |
      .data$duplicate_backbone_mapping | .data$genus_disagreement |
      .data$family_disagreement
  )
write_csv_checked(problematic, output_files[["problematic"]])

unresolved_n <- sum(crosswalk$match_type == "unresolved")
unresolved_fraction <- unresolved_n / nrow(crosswalk)
if (unresolved_n > 10L || unresolved_fraction > 0.05) {
  stop(
    unresolved_n, " species (", scales::percent(unresolved_fraction),
    ") could not be placed. Crosswalk/problem tables were saved; stopping ",
    "because this is a substantial loss.", call. = FALSE
  )
}

# -------------------------------------------------------------------------
# 3. Final tree and patristic distance audit
# -------------------------------------------------------------------------

tree <- phylo_result$scenario.3
tip_map <- crosswalk |>
  dplyr::filter(.data$match_type != "unresolved") |>
  dplyr::select("phylo_input_label", "species_id")

# When several ecological identities share one known parent species, retain
# each ecological ID but do not invent divergence among them. Replace the one
# backbone tip with a zero-length polytomy so their patristic distances and
# downstream phylogenetic coordinates are identical.
tip_groups <- split(tip_map$species_id, tip_map$phylo_input_label)
duplicate_tip_groups <- tip_groups[lengths(tip_groups) > 1L]
for (placement_label in names(duplicate_tip_groups)) {
  where <- match(placement_label, tree$tip.label)
  if (is.na(where)) {
    stop("Shared parent placement tip is absent from the returned tree: ",
         placement_label, call. = FALSE)
  }
  daughter <- ape::stree(
    length(duplicate_tip_groups[[placement_label]]), type = "star",
    tip.label = duplicate_tip_groups[[placement_label]]
  )
  daughter$edge.length <- rep(0, nrow(daughter$edge))
  tree <- ape::bind.tree(tree, daughter, where = where, position = 0)
}

singleton_map <- tip_map |>
  dplyr::add_count(.data$phylo_input_label, name = "n_ecological_taxa") |>
  dplyr::filter(.data$n_ecological_taxa == 1L)
map_index <- match(tree$tip.label, singleton_map$phylo_input_label)
mapped_singletons <- !is.na(map_index)
tree$tip.label[mapped_singletons] <-
  singleton_map$species_id[map_index[mapped_singletons]]

expected_represented <- crosswalk$species_id[crosswalk$match_type != "unresolved"]
if (anyDuplicated(tree$tip.label) ||
    !setequal(tree$tip.label, expected_represented) ||
    length(tree$tip.label) != length(expected_represented)) {
  stop("Final tree tips do not exactly represent the auditable species set.",
       call. = FALSE)
}

saveRDS(tree, output_files[["tree_rds"]])
ape::write.tree(tree, file = output_files[["tree_newick"]])

# Audit tree validity before constructing the distance matrix.
edge_lengths <- as.numeric(tree$edge.length)
tip_depths <- ape::node.depth.edgelength(tree)[seq_len(ape::Ntip(tree))]
root_to_tip_discrepancy <- max(tip_depths) - min(tip_depths)
root_to_tip_relative_discrepancy <- root_to_tip_discrepancy /
  max(1, max(abs(tip_depths)))
tree_checks <- c(
  rooted = ape::is.rooted(tree),
  ultrametric = ape::is.ultrametric(tree, option = 1, tol = 1e-8),
  finite_edge_lengths = length(edge_lengths) > 0L && all(is.finite(edge_lengths)),
  nonnegative_edge_lengths = length(edge_lengths) > 0L &&
    all(edge_lengths >= 0),
  root_to_tip_discrepancy_acceptable = is.finite(root_to_tip_discrepancy) &&
    root_to_tip_relative_discrepancy <= 1e-8
)
if (!all(tree_checks)) {
  stop(
    "Final tree validation failed: ",
    paste(names(tree_checks)[!tree_checks], collapse = ", "), call. = FALSE
  )
}

species_order <- tibble::tibble(
  tree_tip_position = seq_along(tree$tip.label),
  species_id = tree$tip.label,
  newick_tip_label = underscore_species_name(tree$tip.label)
) |>
  dplyr::left_join(
    taxonomy |>
      dplyr::select(
        "species_id", "tnrs_accepted_name", "tnrs_genus", "tnrs_family"
      ),
    by = "species_id"
  )
write_csv_checked(species_order, output_files[["species_order"]])

phylogenetic_distance <- ape::cophenetic.phylo(tree)
distance_tolerance <- 1e-10 * max(1, max(abs(phylogenetic_distance)))
distance_checks <- c(
  square = nrow(phylogenetic_distance) == ncol(phylogenetic_distance),
  unique_taxa = !anyDuplicated(rownames(phylogenetic_distance)) &&
    !anyDuplicated(colnames(phylogenetic_distance)),
  row_column_names_identical = identical(
    rownames(phylogenetic_distance), colnames(phylogenetic_distance)
  ),
  row_order_matches_tree = identical(
    rownames(phylogenetic_distance), tree$tip.label
  ),
  symmetric = max(abs(
    phylogenetic_distance - t(phylogenetic_distance)
  )) <= distance_tolerance,
  negligible_diagonal = max(abs(diag(phylogenetic_distance))) <=
    distance_tolerance,
  complete = !anyNA(phylogenetic_distance) &&
    all(is.finite(phylogenetic_distance))
)
if (!all(distance_checks)) {
  stop(
    "Phylogenetic distance validation failed: ",
    paste(names(distance_checks)[!distance_checks], collapse = ", "),
    call. = FALSE
  )
}
saveRDS(phylogenetic_distance, output_files[["distances"]])

# -------------------------------------------------------------------------
# 4. Response-independent PCoA, negative-eigenvalue audit, and PEVs
# -------------------------------------------------------------------------

pcoa_uncorrected <- ape::pcoa(
  stats::as.dist(phylogenetic_distance), correction = "none"
)
uncorrected_eigenvalues <- as.numeric(
  pcoa_uncorrected$values$Eigenvalues
)
positive_sum <- sum(uncorrected_eigenvalues[
  uncorrected_eigenvalues > 0
])
negative_values <- uncorrected_eigenvalues[
  uncorrected_eigenvalues < 0
]
negative_ratio <- if (positive_sum > 0) {
  sum(abs(negative_values)) / positive_sum
} else {
  Inf
}

# A correction is considered genuinely required when negative inertia exceeds
# one part per million of positive inertia. Both the uncorrected and corrected
# summaries are retained whenever this transparent threshold is crossed.
meaningful_negative_threshold <- 1e-6
correction_required <- length(negative_values) > 0L &&
  is.finite(negative_ratio) && negative_ratio > meaningful_negative_threshold

if (correction_required) {
  pcoa_corrected <- ape::pcoa(
    stats::as.dist(phylogenetic_distance), correction = "lingoes"
  )
} else {
  pcoa_corrected <- NULL
}

broken_stick <- function(n_positive) {
  if (n_positive < 1L) return(numeric())
  vapply(seq_len(n_positive), function(j) {
    sum(1 / seq.int(j, n_positive)) / n_positive
  }, numeric(1))
}

make_eigen_summary <- function(pcoa_object, solution, corrected = FALSE) {
  values <- if (corrected) {
    as.numeric(pcoa_object$values$Corr_eig)
  } else {
    as.numeric(pcoa_object$values$Eigenvalues)
  }
  # ape::pcoa() has already omitted the analytically zero centering axis. Keep
  # every remaining strictly positive eigenvalue, including small terminal
  # axes, rather than imposing an analysis-specific numerical truncation.
  positive <- values > 0
  positive_values <- values[positive]
  positive_total <- sum(positive_values)
  proportion <- rep(NA_real_, length(values))
  cumulative <- rep(NA_real_, length(values))
  expectation <- rep(NA_real_, length(values))
  exceeds <- rep(NA, length(values))
  if (length(positive_values)) {
    proportion[positive] <- positive_values / positive_total
    cumulative[positive] <- cumsum(positive_values) / positive_total
    expectation[positive] <- broken_stick(length(positive_values))
    exceeds[positive] <- proportion[positive] > expectation[positive]
  }
  tibble::tibble(
    solution = solution,
    axis_number = seq_along(values),
    PEV = paste0("PEV", seq_along(values)),
    eigenvalue = values,
    positive_eigenvalue = positive,
    percent_positive_variation = 100 * proportion,
    cumulative_percent = 100 * cumulative,
    broken_stick_expectation_percent = 100 * expectation,
    exceeds_broken_stick = exceeds,
    correction = if (corrected) "Lingoes" else "none"
  )
}

eigen_summary <- make_eigen_summary(
  pcoa_uncorrected, "uncorrected", corrected = FALSE
)
if (correction_required) {
  eigen_summary <- dplyr::bind_rows(
    eigen_summary,
    make_eigen_summary(pcoa_corrected, "corrected", corrected = TRUE)
  )
}
write_csv_checked(eigen_summary, output_files[["eigenvalues"]])

selected_pcoa <- if (correction_required) pcoa_corrected else pcoa_uncorrected
selected_solution <- if (correction_required) "corrected" else "uncorrected"
selected_coordinates <- if (correction_required) {
  as.matrix(selected_pcoa$vectors.cor)
} else {
  as.matrix(selected_pcoa$vectors)
}
selected_summary <- eigen_summary |>
  dplyr::filter(
    .data$solution == .env$selected_solution,
    .data$positive_eigenvalue
  )

if (ncol(selected_coordinates) != nrow(selected_summary)) {
  stop(
    "PCoA coordinate count does not equal the number of positive axes.",
    call. = FALSE
  )
}
colnames(selected_coordinates) <- paste0("PEV", seq_len(ncol(selected_coordinates)))
if (is.null(rownames(selected_coordinates))) {
  rownames(selected_coordinates) <- rownames(phylogenetic_distance)
}
if (!setequal(rownames(selected_coordinates), tree$tip.label)) {
  stop("PCoA score identifiers do not match final tree tips.", call. = FALSE)
}
selected_coordinates <- selected_coordinates[tree$tip.label, , drop = FALSE]

# Enforce exact coordinate identity for ecological taxa represented by the
# same parent-species tip. Eigen decomposition can otherwise introduce tiny
# machine-level differences despite identical distance rows.
shared_coordinate_groups <- crosswalk |>
  dplyr::filter(.data$shared_parent_placement) |>
  dplyr::group_by(.data$`backbone matched name`) |>
  dplyr::summarise(species = list(.data$species_id), .groups = "drop")
for (species_group in shared_coordinate_groups$species) {
  reference_coordinates <- selected_coordinates[species_group[[1]], ]
  selected_coordinates[species_group, ] <- matrix(
    rep(reference_coordinates, each = length(species_group)),
    nrow = length(species_group), byrow = FALSE
  )
}

scores <- tibble::as_tibble(selected_coordinates, rownames = "species_id") |>
  dplyr::left_join(
    crosswalk |>
      dplyr::select(
        "species_id",
        tnrs_accepted_name = "TNRS accepted name",
        tnrs_genus = "TNRS genus",
        tnrs_family = "TNRS family",
        backbone_matched_name = "backbone matched name",
        "match_type", "infraspecific_name", "genus_disagreement",
        "family_disagreement", "resolution_status"
      ),
    by = "species_id"
  ) |>
  dplyr::select(
    "species_id", "tnrs_accepted_name", "tnrs_genus", "tnrs_family",
    "backbone_matched_name", "match_type", "infraspecific_name",
    "genus_disagreement", "family_disagreement", "resolution_status",
    dplyr::starts_with("PEV")
  )
write_csv_checked(scores, output_files[["scores"]])

extreme_axes <- intersect(c("PEV1", "PEV2", "PEV3"), names(scores))
axis_extremes <- scores |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(extreme_axes), names_to = "PEV", values_to = "score"
  ) |>
  dplyr::mutate(
    absolute_score = abs(.data$score),
    disagreement_flag = .data$genus_disagreement | .data$family_disagreement
  ) |>
  dplyr::group_by(.data$PEV) |>
  dplyr::slice_max(.data$absolute_score, n = 10L, with_ties = FALSE) |>
  dplyr::arrange(.data$PEV, dplyr::desc(.data$absolute_score)) |>
  dplyr::ungroup() |>
  dplyr::select(
    "PEV", "species_id", family = "tnrs_family", "score",
    "absolute_score", "match_type", "infraspecific_name",
    "disagreement_flag", "resolution_status"
  )
write_csv_checked(axis_extremes, output_files[["extremes"]])

# -------------------------------------------------------------------------
# 5. Five inspection-only figures
# -------------------------------------------------------------------------

grDevices::pdf(output_files[["tree_figure"]], width = 10, height = 10)
graphics::par(mar = c(0.5, 0.5, 2.5, 0.5))
ape::plot.phylo(
  tree, type = "fan", show.tip.label = FALSE, no.margin = TRUE,
  edge.color = "grey42", edge.width = 0.45
)
graphics::title(
  main = paste0(
    "Resident phylogenetic hypothesis (", length(tree$tip.label), " species)"
  ),
  line = -1
)
grDevices::dev.off()

plot_summary <- selected_summary |>
  dplyr::mutate(
    axis_label = factor(.data$PEV, levels = .data$PEV),
    axis_number = seq_len(dplyr::n())
  )

scree_plot <- ggplot2::ggplot(
  plot_summary,
  ggplot2::aes(x = .data$axis_number)
) +
  ggplot2::geom_line(
    ggplot2::aes(y = .data$percent_positive_variation),
    color = "#2B8CBE", linewidth = 0.7
  ) +
  ggplot2::geom_point(
    ggplot2::aes(y = .data$percent_positive_variation),
    color = "#2B8CBE", size = 0.8
  ) +
  ggplot2::geom_line(
    ggplot2::aes(y = .data$broken_stick_expectation_percent),
    color = "#D95F0E", linewidth = 0.7
  ) +
  ggplot2::geom_point(
    ggplot2::aes(y = .data$broken_stick_expectation_percent),
    color = "#D95F0E", size = 0.8
  ) +
  ggplot2::scale_x_continuous(
    trans = "log10",
    breaks = c(1, 2, 3, 5, 10, 20, 50, 100, 200, 500),
    minor_breaks = NULL
  ) +
  ggplot2::labs(
    x = "Phylogenetic eigenvector axis (log10 scale)",
    y = "Percent of positive phylogenetic variation",
    title = paste0("PCoA scree and broken-stick expectation (", selected_solution, ")"),
    subtitle = "Blue: PEV variation; orange: response-independent broken stick"
  ) +
  ggplot2::theme_minimal(base_size = 10)
ggplot2::ggsave(
  output_files[["scree_figure"]], scree_plot, width = 11, height = 6
)

cumulative_plot <- ggplot2::ggplot(
  plot_summary,
  ggplot2::aes(x = .data$axis_number, y = .data$cumulative_percent)
) +
  ggplot2::geom_hline(
    yintercept = c(50, 75, 90), color = "grey82", linewidth = 0.35
  ) +
  ggplot2::geom_line(color = "#2B8CBE", linewidth = 0.8) +
  ggplot2::geom_point(
    data = plot_summary |>
      dplyr::filter(.data$axis_number %in% c(1, 2, 3, 5, 10)),
    color = "#D95F0E", size = 2
  ) +
  ggplot2::scale_x_continuous(
    breaks = unique(c(1, 2, 3, 5, 10, 20, 50,
                      seq(100, nrow(plot_summary), by = 100))),
    minor_breaks = NULL
  ) +
  ggplot2::scale_y_continuous(limits = c(0, 100)) +
  ggplot2::labs(
    x = "Number of positive PEV axes retained",
    y = "Cumulative positive phylogenetic variation (%)",
    title = "Cumulative phylogenetic variation"
  ) +
  ggplot2::theme_minimal(base_size = 10)
ggplot2::ggsave(
  output_files[["cumulative_figure"]], cumulative_plot,
  width = 11, height = 6
)

match_colors <- c(
  "exact species match" = "#0072B2",
  "synonym/backbone-name match" = "#E69F00",
  "infraspecific -> parent species placement" = "#009E73",
  "grafted within genus" = "#CC79A7",
  "grafted within family" = "#D55E00",
  "unresolved" = "#666666"
)

make_pev_scatter <- function(y_axis, output_path) {
  if (!y_axis %in% names(scores)) {
    stop("Requested PEV scatter axis is unavailable: ", y_axis, call. = FALSE)
  }
  label_data <- scores |>
    dplyr::mutate(
      leverage = abs(.data$PEV1) / max(abs(.data$PEV1)) +
        abs(.data[[y_axis]]) / max(abs(.data[[y_axis]]))
    ) |>
    dplyr::slice_max(.data$leverage, n = 10L, with_ties = FALSE)
  plot <- ggplot2::ggplot(
    scores,
    ggplot2::aes(
      x = .data$PEV1, y = .data[[y_axis]], color = .data$match_type
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey86", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, color = "grey86", linewidth = 0.3) +
    ggplot2::geom_point(alpha = 0.78, size = 1.9) +
    ggrepel::geom_text_repel(
      data = label_data, ggplot2::aes(label = .data$species_id),
      size = 2.6, max.overlaps = Inf, box.padding = 0.25,
      min.segment.length = 0
    ) +
    ggplot2::scale_color_manual(values = match_colors, name = "Placement type") +
    ggplot2::labs(
      x = "PEV1", y = y_axis,
      title = paste("Resident phylogenetic coordinates:", "PEV1 vs", y_axis)
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1))
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      legend.text = ggplot2::element_text(size = 7.5),
      legend.title = ggplot2::element_text(size = 8),
      legend.key.height = grid::unit(0.28, "cm")
    )
  ggplot2::ggsave(output_path, plot, width = 11, height = 8)
}

make_pev_scatter("PEV2", output_files[["pev12_figure"]])
make_pev_scatter("PEV3", output_files[["pev13_figure"]])

# -------------------------------------------------------------------------
# 6. Concise, dynamic audit report
# -------------------------------------------------------------------------

match_counts <- crosswalk |>
  dplyr::count(.data$match_type, name = "n")
match_n <- function(type) {
  value <- match_counts$n[match_counts$match_type == type]
  if (length(value)) value[[1]] else 0L
}
match_line <- function(type, label) {
  n <- match_n(type)
  paste0(label, ": ", n, " (", sprintf("%.1f%%", 100 * n / nrow(crosswalk)), ")")
}

cumulative_at <- function(n) {
  if (!nrow(selected_summary)) return(NA_real_)
  selected_summary$cumulative_percent[[min(n, nrow(selected_summary))]]
}

broken_axes <- selected_summary$PEV[
  !is.na(selected_summary$exceeds_broken_stick) &
    selected_summary$exceeds_broken_stick
]
most_negative <- if (length(negative_values)) min(negative_values) else 0
disagreements_detected <- sum(crosswalk$any_taxonomy_disagreement)
disagreements_resolved <- sum(
  crosswalk$any_taxonomy_disagreement & crosswalk$resolution_status == "resolved"
)
disagreements_outstanding <- sum(
  crosswalk$any_taxonomy_disagreement & crosswalk$resolution_status == "outstanding"
)
extreme_problem_count <- axis_extremes |>
  dplyr::filter(
    grepl("grafted|unresolved", .data$match_type) |
      .data$disagreement_flag
  ) |>
  dplyr::summarise(n = dplyr::n_distinct(.data$species_id)) |>
  dplyr::pull(.data$n)
axis_extreme_lines <- unlist(lapply(extreme_axes, function(axis) {
  table <- axis_extremes |>
    dplyr::filter(.data$PEV == .env$axis) |>
    dplyr::mutate(
      summary = paste0(
        .data$species_id, " | ", .data$family, " | score=",
        format(.data$score, digits = 6), " | abs=",
        format(.data$absolute_score, digits = 6), " | ", .data$match_type,
        " | infra=", .data$infraspecific_name,
        " | disagreement=", .data$disagreement_flag
      )
    ) |>
    dplyr::pull(.data$summary)
  c(paste0(axis, " top 10 absolute scores:"), paste0("  ", table))
}), use.names = FALSE)
infraspecific_report_lines <- capture.output(print(
  infraspecific_placement_audit, n = Inf, width = Inf
))

report_lines <- c(
  "RESIDENT PHYLOGENY AND RESPONSE-INDEPENDENT PCoA AUDIT",
  "======================================================",
  "",
  "1. ANALYSIS SPECIES SET",
  "-----------------------",
  paste0("Updated downstream resident species set: ", nrow(taxonomy)),
  paste0("Excluded before phylogeny/downstream modeling: ",
         length(resident_downstream_model_exclusions()), " (",
         paste(resident_downstream_model_exclusions(), collapse = "; "), ")"),
  "The excluded taxa remain in upstream community and winner/loser calculations.",
  paste0("Species represented in final tree: ", length(tree$tip.label)),
  paste0("All phylogeny species exactly match updated ecological analysis: ",
         setequal(tree$tip.label, taxonomy$species_id)),
  "No AUC response was retained, modeled, inspected, or used to select axes.",
  "The existing TNRS accepted name, genus, and family remain canonical.",
  "",
  "2. TAXONOMY AND PLACEMENT",
  "-------------------------",
  paste0("V.PhyloMaker2 version: ",
         as.character(utils::packageVersion("V.PhyloMaker2"))),
  "Backbone: GBOTB.extended.LCVP with nodes.info.1.LCVP",
  "Placement: deterministic scenario 3 (S3)",
  paste0(
    "Verified call: V.PhyloMaker2::phylo.maker(sp.list, ",
    "tree = GBOTB.extended.LCVP, nodes = nodes.info.1.LCVP, ",
    "output.sp.list = TRUE, output.tree = FALSE, scenarios = \"S3\")"
  ),
  match_line("exact species match", "Exact species matches"),
  match_line("synonym/backbone-name match", "Synonym/backbone-name matches"),
  match_line("infraspecific -> parent species placement",
             "Infraspecific -> parent species placements"),
  match_line("grafted within genus", "Grafted within genus"),
  match_line("grafted within family", "Grafted within family"),
  match_line("unresolved", "Unresolved"),
  paste0("TNRS/backbone genus disagreements: ",
         sum(crosswalk$genus_disagreement)),
  paste0("TNRS/backbone family disagreements: ",
         sum(crosswalk$family_disagreement)),
  paste0("Genus/family disagreements detected: ", disagreements_detected),
  paste0("Resolved with recorded decision/reason: ", disagreements_resolved),
  paste0("Outstanding disagreements: ", disagreements_outstanding),
  "TNRS ecological identifiers were retained; placement decisions are recorded in taxonomy_disagreement_resolution.csv.",
  paste0("Ecological taxa initially mapping to duplicated backbone taxa: ",
         sum(crosswalk$duplicate_backbone_mapping)),
  "Duplicated backbone mappings were never silently merged; see the crosswalk.",
  "Newick convention encodes spaces as underscores; phylogeny_species_order.csv maps every Newick tip back to species_id.",
  paste0("V.PhyloMaker2 warnings captured: ", length(unique(phylo_warnings))),
  "",
  "Canonical-taxonomy audit",
  paste0("Duplicated species_id values: ",
         sum(taxonomy_audit$duplicated_species_id)),
  paste0("Duplicated TNRS accepted names: ",
         sum(taxonomy_audit$duplicated_tnrs_accepted_name)),
  paste0("Missing genera: ", sum(taxonomy_audit$missing_genus)),
  paste0("Missing families: ", sum(taxonomy_audit$missing_family)),
  paste0("Infraspecific names: ", sum(taxonomy_audit$infraspecific_name)),
  paste0("Authorship/formatting artifacts: ",
         sum(taxonomy_audit$authorship_or_formatting_artifact)),
  paste0("Names not straightforward binomials: ",
         sum(!taxonomy_audit$straightforward_binomial)),
  "",
  "Infraspecific placement audit",
  paste0(
    "Parent species tip present in GBOTB.extended.LCVP: Yes = ",
    sum(infraspecific_placement_audit[[
      "Parent species tip present in GBOTB.extended.LCVP?"
    ]] == "Yes"),
    "; No = ",
    sum(infraspecific_placement_audit[[
      "Parent species tip present in GBOTB.extended.LCVP?"
    ]] == "No"), "."
  ),
  paste0(
    "Absence from GBOTB.extended.LCVP indicates absence from the ",
    "phylogenetic backbone and does not imply that the taxon is unrecognized ",
    "by current LCVP/WCVP-based taxonomic resources."
  ),
  infraspecific_report_lines,
  "",
  "3. TREE DIAGNOSTICS",
  "-------------------",
  paste0("Rooted: ", tree_checks[["rooted"]]),
  paste0("Ultrametric: ", tree_checks[["ultrametric"]]),
  paste0("Negative branch lengths: ", sum(edge_lengths < 0)),
  paste0("Non-finite branch lengths: ", sum(!is.finite(edge_lengths))),
  paste0("Zero-length branches: ", sum(edge_lengths == 0)),
  paste0("Minimum branch length: ", format(min(edge_lengths), scientific = TRUE)),
  paste0("Maximum branch length: ", format(max(edge_lengths), scientific = TRUE)),
  paste0("Minimum root-to-tip distance: ", format(min(tip_depths), digits = 9)),
  paste0("Maximum root-to-tip distance: ", format(max(tip_depths), digits = 9)),
  paste0("Maximum root-to-tip discrepancy: ",
         format(root_to_tip_discrepancy, scientific = TRUE, digits = 7)),
  "",
  "Distance-matrix validation",
  paste0(names(distance_checks), ": ", unname(distance_checks)),
  "",
  "4. PCoA DIAGNOSTICS",
  "-------------------",
  paste0("Uncorrected negative eigenvalues: ", length(negative_values)),
  paste0("Most negative uncorrected eigenvalue: ",
         format(most_negative, scientific = TRUE, digits = 7)),
  paste0("Total negative magnitude / positive variation: ",
         format(negative_ratio, scientific = TRUE, digits = 7)),
  paste0("Meaningful-negative threshold: ",
         format(meaningful_negative_threshold, scientific = TRUE)),
  paste0("Meaningful negative eigenvalues present: ", correction_required),
  paste0("Correction applied: ",
         if (correction_required) "Lingoes" else "none"),
  paste0("PCoA solution used for inspection PEVs: ", selected_solution),
  paste0("Positive PCoA axes retained: ", nrow(selected_summary)),
  paste0("PEV1 positive variation (%): ",
         format(cumulative_at(1), digits = 6)),
  paste0("PEV1-2 cumulative positive variation (%): ",
         format(cumulative_at(2), digits = 6)),
  paste0("PEV1-3 cumulative positive variation (%): ",
         format(cumulative_at(3), digits = 6)),
  paste0("PEV1-5 cumulative positive variation (%): ",
         format(cumulative_at(5), digits = 6)),
  paste0("PEV1-10 cumulative positive variation (%): ",
         format(cumulative_at(10), digits = 6)),
  paste0("Axes exceeding descriptive broken-stick expectation: ",
         if (length(broken_axes)) paste(broken_axes, collapse = ", ") else "none"),
  "This is a response-independent description of eigenstructure and is not used as an axis-retention rule for ecological modeling.",
  "",
  "5. AXIS LEVERAGE AUDIT",
  "----------------------",
  axis_extreme_lines,
  paste0("Distinct grafted/disagreement taxa among the 30 axis-extreme rows: ",
         extreme_problem_count),
  "",
  "6. INTERPRETIVE BOUNDARY",
  "------------------------",
  paste0(
    "PEV orientation is arbitrary. These axes summarize phylogenetic ",
    "distances only and are not interpreted relative to Delta AUC. No choice ",
    "has been made about which, if any, PEVs should enter an ecological model."
  ),
  "Higher-order PEVs increasingly represent fine phylogenetic structure, some of which depends on deterministic grafting of taxa absent from the backbone.",
  "",
  "Saved outputs",
  paste0(seq_along(output_files), ". ", unname(output_files))
)
writeLines(report_lines, output_files[["report"]])

missing_outputs <- unname(output_files)[
  !file.exists(unname(output_files)) |
    file.info(unname(output_files))$size <= 0
]
if (length(missing_outputs)) {
  stop("Missing or empty requested outputs: ",
       paste(missing_outputs, collapse = ", "), call. = FALSE)
}

cat("Resident phylogeny/PCoA audit completed.\n")
cat("Requested species:", nrow(taxonomy), "\n")
cat("Represented species:", length(tree$tip.label), "\n")
cat("Positive PEV axes:", nrow(selected_summary), "\n")
cat("Correction:", if (correction_required) "Lingoes" else "none", "\n")
cat("Report:", output_files[["report"]], "\n")
