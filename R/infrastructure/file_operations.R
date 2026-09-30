# Module: R/infrastructure/file_operations
# Responsibility: Implement the named scientific workflow stage or shared infrastructure.
# Inputs: Function arguments and explicitly declared target dependencies.
# Outputs: Validated in-memory objects and declared workflow artifacts.
# Side effects: Limited to functions whose contracts explicitly write, copy, or render files.

# Explicit file-export helpers used only by file targets.

#' Write Text Atomic.
#'
#' @param text A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param path A project-relative or absolute filesystem path.
#' @return The path or paths written by the declared side effect.
write_text_atomic <- function(text, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(
    pattern = paste0(basename(path), "-"),
    tmpdir = dirname(path),
    fileext = ".tmp"
  )
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  writeLines(enc2utf8(as.character(text)), temporary, useBytes = TRUE)
  if (!file.rename(temporary, path)) {
    stop("Could not publish staged text file: ", path, call. = FALSE)
  }
  path
}
