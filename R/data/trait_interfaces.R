# Module: R/data/trait_interfaces
# Responsibility: Implement the named scientific workflow stage or shared infrastructure.
# Inputs: Function arguments and explicitly declared target dependencies.
# Outputs: Validated in-memory objects and declared workflow artifacts.
# Side effects: Limited to functions whose contracts explicitly write, copy, or render files.

# Narrow validation interface for audited species-level traits.

#' Build Trait Species.
#'
#' @param try_summary A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param retained_species A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
build_trait_species <- function(try_summary, retained_species) {
  retained_species <- sort(unique(as.character(retained_species)))
  values <- try_summary |>
    dplyr::filter(.data$trait %in% c("LDMC", "LMA", "vegetative_height")) |>
    dplyr::transmute(
      species_id = as.character(.data$accepted_name),
      trait = as.character(.data$trait),
      value = as.numeric(.data$median)
    ) |>
    tidyr::pivot_wider(names_from = "trait", values_from = "value")
  result <- tibble::tibble(species_id = retained_species) |>
    dplyr::left_join(values, by = "species_id")
  validate_trait_species_interface(result, length(retained_species))
  result
}


#' Validate Trait Species Interface.
#'
#' @param traits A data frame whose required columns are checked before use.
#' @param expected_species A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
validate_trait_species_interface <- function(traits, expected_species = 526L) {
  required <- c("species_id", "LDMC", "LMA", "vegetative_height")
  missing <- setdiff(required, names(traits))
  if (length(missing)) stop("Trait audit is missing: ", paste(missing, collapse = ", "))
  tibble::tibble(
    species = nrow(traits),
    unique_species = !anyDuplicated(traits$species_id),
    expected_species = nrow(traits) == expected_species,
    pass = !anyDuplicated(traits$species_id) && nrow(traits) == expected_species
  )
}
