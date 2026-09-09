# Explicit dependency refresh. This is the only script that updates GitHub code
# or renv.lock; ordinary targets runs use the recorded lockfile.

if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv", repos = "https://cloud.r-project.org")

cran_packages <- c(
  "clubSandwich", "DescTools", "TNRS", "crew", "dplyr", "emmeans", "geometry",
  "ggalluvial", "ggplot2", "ggrepel", "here", "knitr", "lme4", "mclogit",
  "metafor", "patchwork", "purrr", "readr", "rmarkdown", "sandwich",
  "scales", "scatterpie", "sf", "stringi", "stringr", "targets", "terra",
  "testthat", "tibble", "tidyr"
)

renv::install(cran_packages, prompt = FALSE)
renv::install("billurbektas/tidyTRY", prompt = FALSE, rebuild = TRUE)

git_head <- function(repository) {
  line <- system2("git", c("ls-remote", repository, "HEAD"), stdout = TRUE, stderr = TRUE)
  if (!length(line) || !grepl("^[0-9a-f]{40}", line[[1]])) {
    stop("Could not resolve HEAD for ", repository, call. = FALSE)
  }
  strsplit(line[[1]], "[[:space:]]+")[[1]][[1]]
}

tidy_repository <- "https://github.com/billurbektas/tidyTRY"
chelsa_repository <- "https://github.com/billurbektas/get_chelsa"
tidy_sha <- git_head(tidy_repository)
chelsa_sha <- git_head(chelsa_repository)

dir.create("external", showWarnings = FALSE)
chelsa_url <- sprintf("https://raw.githubusercontent.com/billurbektas/get_chelsa/%s/R/get_chelsa.R", chelsa_sha)
download.file(chelsa_url, file.path("external", "get_chelsa.R"), mode = "wb", quiet = FALSE)

sources <- data.frame(
  component = c("tidyTRY", "get_chelsa"),
  repository = c(tidy_repository, chelsa_repository),
  git_sha = c(tidy_sha, chelsa_sha),
  recorded_date = as.character(Sys.Date())
)
write.csv(sources, file.path("config", "external_sources.csv"), row.names = FALSE, na = "")

renv::snapshot(
  packages = c("renv", cran_packages, "tidyTRY"),
  prompt = FALSE,
  force = TRUE
)

message("Dependency refresh complete. Review config/external_sources.csv and renv.lock before committing.")
