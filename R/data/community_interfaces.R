# Module: R/data/community_interfaces
# Responsibility: Validate community and taxonomy tables before analysis joins.
# Inputs: Standardized community and reviewed taxonomy tables.
# Outputs: Fail-closed interface checks.
# Side effects: None.

# Focused validation functions for community and taxonomy targets.

#' Validate Community Interface.
#'
#' @param community A data frame whose required columns are checked before use.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
validate_community_interface <- function(community) {
  required <- c("source_id", "destination_id", "species_id")
  missing <- setdiff(required, names(community))
  if (length(missing)) {
    stop("Community data are missing: ", paste(missing, collapse = ", "))
  }
  tibble::tibble(rows = nrow(community), species = dplyr::n_distinct(community$species_id), pass = TRUE)
}

#' Validate Taxonomy Interface.
#'
#' @param taxonomy A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
validate_taxonomy_interface <- function(taxonomy) {
  required <- c("accepted_name", "family")
  missing <- setdiff(required, names(taxonomy))
  if (length(missing)) stop("Taxonomy data are missing: ", paste(missing, collapse = ", "))
  tibble::tibble(rows = nrow(taxonomy), duplicated_names = anyDuplicated(taxonomy$accepted_name), pass = !anyDuplicated(taxonomy$accepted_name))
}
