if (!requireNamespace("renv", quietly = TRUE)) stop("Install renv first.", call. = FALSE)
if (file.exists("renv.lock")) renv::restore(prompt = FALSE)
if (!requireNamespace("targets", quietly = TRUE)) stop("targets is unavailable after renv restore.", call. = FALSE)

workers <- as.integer(Sys.getenv("TRANSPLANT_CHELSA_WORKERS", "2"))
message("Running targets pipeline; CHELSA worker limit: ", workers)
targets::tar_make()
