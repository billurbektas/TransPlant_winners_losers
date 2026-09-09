# Narrow validation interfaces for upstream community and taxonomy targets.

validate_canonical_community_interface <- function(community) {
  required <- c("source_id", "destination_id", "species_id")
  missing <- setdiff(required, names(community))
  if (length(missing)) {
    stop("Canonical community is missing: ", paste(missing, collapse = ", "))
  }
  tibble::tibble(rows = nrow(community), species = dplyr::n_distinct(community$species_id), pass = TRUE)
}

validate_canonical_taxonomy_interface <- function(taxonomy) {
  required <- c("accepted_name", "family")
  missing <- setdiff(required, names(taxonomy))
  if (length(missing)) stop("Canonical taxonomy is missing: ", paste(missing, collapse = ", "))
  tibble::tibble(rows = nrow(taxonomy), duplicated_names = anyDuplicated(taxonomy$accepted_name), pass = !anyDuplicated(taxonomy$accepted_name))
}
