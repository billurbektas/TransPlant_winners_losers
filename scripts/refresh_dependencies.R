#!/usr/bin/env Rscript

# Explicit dependency refresh for the current targets workflow. Ordinary
# targets runs never update packages, GitHub sources, or renv.lock.

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}

cran_packages <- c(
  "clubSandwich", "DescTools", "TNRS", "crew", "dplyr", "emmeans",
  "geometry", "ggalluvial", "ggplot2", "ggrepel", "here", "knitr",
  "lme4", "mclogit", "metafor", "patchwork", "purrr", "quarto", "readr",
  "rmarkdown", "sandwich", "scales", "scatterpie", "sf", "stringi",
  "stringr", "tarchetypes", "targets", "terra", "testthat", "tibble",
  "tidyr"
)

vphylomaker_sha <- "7af3fb5152f691af2e4ec9d5e2e467d1b50505e9"
renv::install(cran_packages, prompt = FALSE)
renv::install("billurbektas/tidyTRY", prompt = FALSE, rebuild = TRUE)
renv::install(
  paste0("jinyizju/V.PhyloMaker2@", vphylomaker_sha),
  prompt = FALSE, rebuild = TRUE
)

git_head <- function(repository) {
  line <- system2(
    "git", c("ls-remote", repository, "HEAD"),
    stdout = TRUE, stderr = TRUE
  )
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
chelsa_url <- sprintf(
  "https://raw.githubusercontent.com/billurbektas/get_chelsa/%s/R/get_chelsa.R",
  chelsa_sha
)
download.file(
  chelsa_url, file.path("external", "get_chelsa.R"),
  mode = "wb", quiet = FALSE
)

sources <- data.frame(
  component = c("tidyTRY", "get_chelsa", "V.PhyloMaker2"),
  repository = c(
    tidy_repository, chelsa_repository,
    "https://github.com/jinyizju/V.PhyloMaker2"
  ),
  git_sha = c(tidy_sha, chelsa_sha, vphylomaker_sha),
  recorded_date = as.character(Sys.Date())
)
write.csv(
  sources, file.path("config", "external_sources.csv"),
  row.names = FALSE, na = ""
)

renv::snapshot(
  packages = c(
    "renv", cran_packages, "tidyTRY", "V.PhyloMaker2"
  ),
  prompt = FALSE, force = TRUE
)

message(
  "Dependency refresh complete. Review config/external_sources.csv and ",
  "renv.lock before committing."
)
