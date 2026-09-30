# Module: R/infrastructure/file_operations
# Responsibility: Implement the named scientific workflow stage or shared infrastructure.
# Inputs: Function arguments and explicitly declared target dependencies.
# Outputs: Validated in-memory objects and declared workflow artifacts.
# Side effects: Limited to functions whose contracts explicitly write, copy, or render files.

# Explicit file-export helpers used only by file targets.

#' Require an Existing File Without Modifying It.
#'
#' @param path A project-relative or absolute file path.
#' @details This helper is for immutable inputs such as Phase 1 historical
#'   references. It has no side effects and fails closed when the file is
#'   missing or is a directory.
#' @return `path`, invisibly unchanged, after confirming it identifies an existing regular file.
require_existing_file <- function(path) {
  if (length(path) != 1L || is.na(path) || !nzchar(path) ||
      !file.exists(path) || dir.exists(path)) {
    stop("Required input file does not exist: ", path, call. = FALSE)
  }
  path
}

#' Copy File Atomic Verified.
#'
#' @param input A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param output A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return The path or paths written by the declared side effect.
copy_file_atomic_verified <- function(input, output) {
  if (length(input) > 1L) {
    candidates <- input[
      file.exists(input) & !dir.exists(input) &
        tolower(tools::file_ext(input)) %in% c("html", "pdf", "docx", "epub")
    ]
    if (length(candidates) != 1L) {
      stop("Expected exactly one rendered report in the input file vector.", call. = FALSE)
    }
    input <- candidates[[1]]
  }
  if (!file.exists(input) || dir.exists(input)) {
    stop("Input file does not exist: ", input, call. = FALSE)
  }
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(basename(output), "-"), tmpdir = dirname(output))
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  if (!file.copy(input, temporary, overwrite = TRUE, copy.mode = TRUE)) {
    stop("Could not stage file copy: ", input, call. = FALSE)
  }
  if (!identical(sha256_file(input), sha256_file(temporary))) {
    stop("Staged file checksum differs from source: ", input, call. = FALSE)
  }
  if (file.exists(output)) unlink(output)
  if (!file.rename(temporary, output)) {
    stop("Could not publish staged file: ", output, call. = FALSE)
  }
  output
}

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
