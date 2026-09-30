# Module: R/validation/dependency_audit
# Responsibility: Implement the named scientific workflow stage or shared infrastructure.
# Inputs: Function arguments and explicitly declared target dependencies.
# Outputs: Validated in-memory objects and declared workflow artifacts.
# Side effects: Limited to functions whose contracts explicitly write, copy, or render files.

# Executable audit of target ancestry at the current/history boundary.

#' Target Ancestors.
#'
#' @param target A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param edges A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
target_ancestors <- function(target, edges) {
  visited <- character()
  frontier <- target
  while (length(frontier)) {
    upstream <- unique(edges$from[edges$to %in% frontier])
    upstream <- setdiff(upstream, visited)
    if (!length(upstream)) break
    visited <- c(visited, upstream)
    frontier <- upstream
  }
  unique(visited)
}

#' Audit Legacy Fixture Dependencies.
#'
#' @param root A project-relative or absolute filesystem path.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
audit_legacy_fixture_dependencies <- function(root = NULL) {
  if (is.null(root)) {
    root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
    repeat {
      if (file.exists(file.path(root, "_targets.R"))) break
      parent <- dirname(root)
      if (identical(parent, root)) {
        stop("Could not locate _targets.R from the current directory.", call. = FALSE)
      }
      root <- parent
    }
  }
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  previous_directory <- setwd(root)
  on.exit(setwd(previous_directory), add = TRUE)
  network <- targets::tar_network(
    targets_only = TRUE,
    outdated = FALSE,
    callr_function = NULL,
    script = file.path(root, "_targets.R"),
    store = file.path(root, "_targets")
  )
  edges <- network$edges
  inventory <- pipeline_target_stage_table()
  historical_seed <- inventory$target[
    inventory$stage == "legacy" &
      grepl(
        "^(legacy_|frozen_release_reference)|equivalence|migration|reproduction|resident_consistency_checks",
        inventory$target
      )
  ]
  result <- dplyr::bind_rows(lapply(inventory$target, function(target) {
    direct <- sort(unique(edges$from[edges$to == target]))
    ancestors <- target_ancestors(target, edges)
    stage <- inventory$stage[inventory$target == target][[1]]
    observed_fixture_ancestry <- target %in% historical_seed ||
      any(ancestors %in% historical_seed)
    role <- if (stage == "legacy") "LEGACY_CERTIFICATION" else
      "CURRENT_WORKFLOW"
    tibble::tibble(
      target = target,
      stage = stage,
      role = role,
      direct_dependencies = paste(direct, collapse = "; "),
      fixture_ancestor = observed_fixture_ancestry,
      pass = stage == "legacy" || !observed_fixture_ancestry,
      detail = if (stage == "legacy") {
        "Historical ancestry is quarantined in the temporary certification factory."
      } else {
        "Current workflow targets must have no historical-reference ancestry."
      }
    )
  }))
  if (any(!result$pass)) {
    stop(
      "Current/history dependency boundary failed for: ",
      paste(result$target[!result$pass], collapse = ", "),
      call. = FALSE
    )
  }
  result
}
