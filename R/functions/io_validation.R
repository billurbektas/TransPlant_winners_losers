# Explicit file-export helpers used only by file targets.

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
