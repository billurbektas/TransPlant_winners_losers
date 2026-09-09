normalize_name <- function(x) {
  x |>
    stringi::stri_trans_general("Latin-ASCII") |>
    stringr::str_squish() |>
    stringr::str_to_lower()
}

repair_utf8 <- function(x) {
  x <- as.character(x)
  invalid <- !stringi::stri_enc_isutf8(x)
  x[invalid] <- stringi::stri_encode(x[invalid], from = "windows-1252", to = "UTF-8")
  x
}

assert_columns <- function(x, required, object = deparse(substitute(x))) {
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop(object, " is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(x)
}

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

read_config_csv <- function(root, filename, ...) {
  readr::read_csv(file.path(root, "config", filename), show_col_types = FALSE, ...)
}

safe_scale <- function(x) {
  if (all(is.na(x)) || stats::sd(x, na.rm = TRUE) == 0) return(rep(NA_real_, length(x)))
  as.numeric(scale(x))
}
