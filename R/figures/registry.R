# Module: R/figures/registry
# Responsibility: Dispatch registered export-based publication figures.
# Inputs: Figure key and existing target exports.
# Outputs: PDF/PNG manuscript figures.
# Side effects: Renderer calls write registered figure files.

# Semantic registry shared by the publication graph and figure-only CLI.

#' List standalone figure renderers.
#'
#' @return A data frame with stable renderer IDs, functions, and descriptions.
figure_renderer_registry <- function() {
  data.frame(
    id = c(
      "climate-response", "predictor-atlas", "predictor-synthesis",
      "family-context", "phylogeny-axes", "phylogeny-summary"
    ),
    renderer = c(
      "render_climate_response_figure",
      "render_predictor_response_atlas",
      "render_predictor_synthesis_figure",
      "render_family_context_appendix",
      "render_phylogeny_axis_atlas",
      "render_phylogeny_summary_appendix"
    ),
    description = c(
      "Resident climate fixed-effect response surfaces",
      "Trait, functional, and phylogenetic predictor-response atlas",
      "Predictor support and heterogeneity synthesis",
      "Family context-classification appendix",
      "Secondary phylogenetic-axis response atlas",
      "Phylogenetic-coordinate and eigenvalue summary"
    ),
    stringsAsFactors = FALSE
  )
}

#' Render one figure from existing exported artifacts.
#'
#' @param id Stable renderer ID from [figure_renderer_registry()].
#' @param arguments Renderer-specific positional arguments.
#' @return A character vector of generated artifact paths, invisibly.
render_figure_from_exports <- function(id, arguments = character()) {
  registry <- figure_renderer_registry()
  row <- registry[registry$id == id, , drop = FALSE]
  if (nrow(row) != 1L) {
    stop(
      "Unknown figure ID `", id, "`. Use `scripts/render_figure.R --list`.",
      call. = FALSE
    )
  }
  renderer <- get(row$renderer[[1]], mode = "function", inherits = TRUE)
  renderer(arguments)
}
