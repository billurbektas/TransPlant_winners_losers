# A data frame as a Markdown table, for the step reports (every analysis
# script writes one). Same logic as prep/R/report.R, which serves the prep
# pipeline.
md_table <- function(d) {
  d <- dplyr::mutate(d, dplyr::across(dplyr::everything(), ~ ifelse(is.na(.x), "", as.character(.x))))
  c(paste("|", paste(names(d), collapse = " | "), "|"),
    paste("|", paste(rep("---", ncol(d)), collapse = " | "), "|"),
    apply(d, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}
