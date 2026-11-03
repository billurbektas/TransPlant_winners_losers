# Caches keyed by their inputs (Bills, 3 Nov); used by steps 02, 05, 06, 06b, 07a, 07c,
# 07d and 07e. A cached fit is reused only when the key stored with it equals the key of
# today's inputs: the rows it is fitted to (only the columns it uses), its model
# specification and the script's cache version. When an upstream step changes those rows,
# the key changes and the fit is redone; fits whose inputs did not change are kept. The
# cache version is raised by hand when the fitting code itself changes.
cache_key <- function(...) digest::digest(list(...), algo = "xxhash64")

# the cached object if its key matches, NULL otherwise (missing, older or stale)
read_cache <- function(path, key) {
  if (!file.exists(path)) return(NULL)
  x <- readRDS(path)
  if (identical(attr(x, "cache_key"), key)) x else NULL
}

write_cache <- function(x, path, key) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  attr(x, "cache_key") <- key
  saveRDS(x, path)
  x
}
