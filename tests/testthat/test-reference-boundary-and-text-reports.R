test_that("frozen human-readable reports are exact verified copies", {
  copies <- tibble::tribble(
    ~copy, ~sha256,
    "reports/reference_text/step22_downstream_526_report.txt",
    "a84a4d1aeda3c10ec2943d7253f419d22f092c0b9664f0148ebbf72903f3724c",
    "reports/reference_text/step21_trait_data_audit_report.txt",
    "4f8cf56f554cbb42844f4f74664fba5db8f20c0438f7af06e74084a0b9af9853",
    "reports/reference_text/step19_phylogeny_pcoa_audit_report.txt",
    "2da6484d4da1e3ee7524b65a41ccf360f03b4efc557357206dc25506fc6382e5"
  )
  purrr::pwalk(copies, function(copy, sha256) {
    copy <- project_path(copy)
    expect_true(file.exists(copy))
    expect_identical(sha256_file(copy), sha256)
  })
})

test_that("current scientific targets have no frozen-fixture ancestry", {
  audit <- audit_legacy_fixture_dependencies()
  scientific <- dplyr::filter(
    audit, .data$role == "CURRENT_WORKFLOW"
  )
  validation <- dplyr::filter(audit, .data$role == "LEGACY_CERTIFICATION")
  expect_true(all(scientific$pass))
  expect_false(any(scientific$fixture_ancestor))
  expect_true(all(validation$pass))
  expect_true(any(validation$fixture_ancestor))

  current_terminal_targets <- c(
    "release_validation", "release_artifact_manifest",
    "report_data_foundation", "report_response_estimation",
    "report_resident_climate_model", "report_family_context_model",
    "report_predictor_models", "report_model_synthesis",
    "report_release_validation"
  )
  terminals <- dplyr::filter(
    scientific, .data$target %in% current_terminal_targets
  )
  expect_setequal(terminals$target, current_terminal_targets)
  expect_false(any(terminals$fixture_ancestor))

  phylogeny <- dplyr::filter(
    scientific, .data$target == "phylogeny_model_sequence"
  )
  expect_match(
    phylogeny$direct_dependencies,
    "phylogeny_predictor_specification",
    fixed = TRUE
  )
  expect_false(grepl("reference", phylogeny$direct_dependencies, fixed = TRUE))
})

test_that("the current text report is assembled from active objects", {
  source_text <- paste(
    readLines(
      project_path("R", "targets", "definitions", "publication.R"),
      warn = FALSE
    ),
    collapse = "\n"
  )
  start <- regexpr("current_model_report,", source_text, fixed = TRUE)[[1]]
  end <- regexpr("current_model_results_file,", source_text, fixed = TRUE)[[1]]
  expect_gt(start, 0L)
  expect_gt(end, start)
  block <- substring(source_text, start, end - 1L)
  expect_match(block, "resident_climate_model_results", fixed = TRUE)
  expect_match(block, "heterogeneity_synthesis", fixed = TRUE)
  expect_match(block, "current_scientific_release_validation", fixed = TRUE)
  expect_false(grepl("legacy_transition_reference", block, fixed = TRUE))
  expect_false(grepl("step22_reference", block, fixed = TRUE))
  expect_false(grepl("archive/", block, fixed = TRUE))
})
