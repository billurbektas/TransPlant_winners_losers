# Module: R/infrastructure/io
# Responsibility: Implement the named scientific workflow stage or shared infrastructure.
# Inputs: Function arguments and explicitly declared target dependencies.
# Outputs: Validated in-memory objects and declared workflow artifacts.
# Side effects: Limited to functions whose contracts explicitly write, copy, or render files.

#' Read a pipeline CSV without changing empty diagnostic text
#'
#' @param path Existing CSV path.
#' @details Missing diagnostic text becomes an empty string, matching the
#'   established CSV read-back used for resident-model inputs.
#' @return A tibble with the established CSV parsing behavior.
read_pipeline_csv <- function(path) {
  result <- readr::read_csv(path, show_col_types = FALSE)
  empty_text_fields <- intersect(
    names(result), c("warnings", "error", "CR1_warnings", "contrast_warnings")
  )
  for (field in empty_text_fields) {
    result[[field]][is.na(result[[field]])] <- ""
  }
  result
}

#' Normalize Name.
#'
#' @param x A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
normalize_name <- function(x) {
  x |>
    stringi::stri_trans_general("Latin-ASCII") |>
    stringr::str_squish() |>
    stringr::str_to_lower()
}

#' Repair Utf8.
#'
#' @param x A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
repair_utf8 <- function(x) {
  x <- as.character(x)
  invalid <- !stringi::stri_enc_isutf8(x)
  x[invalid] <- stringi::stri_encode(x[invalid], from = "windows-1252", to = "UTF-8")
  x
}

#' Assert Columns.
#'
#' @param x A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param required A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param object A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
assert_columns <- function(x, required, object = deparse(substitute(x))) {
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop(object, " is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(x)
}

#' Assert Unique.
#'
#' @param x A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param keys A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param object A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
assert_unique <- function(x, keys, object = deparse(substitute(x))) {
  assert_columns(x, keys, object)
  duplicated_keys <- x |>
    dplyr::count(dplyr::across(dplyr::all_of(keys)), name = "n") |>
    dplyr::filter(.data$n > 1L)
  if (nrow(duplicated_keys)) {
    stop(object, " is not unique by ", paste(keys, collapse = ", "), call. = FALSE)
  }
  invisible(x)
}

#' Write Rds Atomic.
#'
#' @param x A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param path A project-relative or absolute filesystem path.
#' @return The path or paths written by the declared side effect.
write_rds_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
  saveRDS(x, temporary)
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    stop("Could not atomically move output to ", path, call. = FALSE)
  }
  path
}

#' Write Csv Atomic.
#'
#' @param x A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param path A project-relative or absolute filesystem path.
#' @return The path or paths written by the declared side effect.
write_csv_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
  readr::write_csv(x, temporary, na = "")
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    stop("Could not atomically move output to ", path, call. = FALSE)
  }
  path
}

#' Read Config Csv.
#'
#' @param root A project-relative or absolute filesystem path.
#' @param filename A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param ... Additional arguments forwarded to the delegated operation.
#' @return An R object described by the function name; its structure is validated before later stages use it.
read_config_csv <- function(root, filename, ...) {
  readr::read_csv(file.path(root, "config", filename), show_col_types = FALSE, ...)
}

#' Safe Scale.
#'
#' @param x A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
safe_scale <- function(x) {
  if (all(is.na(x)) || stats::sd(x, na.rm = TRUE) == 0) return(rep(NA_real_, length(x)))
  as.numeric(scale(x))
}
