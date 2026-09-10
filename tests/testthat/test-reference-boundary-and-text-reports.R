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

test_that("canonical scientific targets have no frozen-fixture ancestry", {
  audit <- audit_canonical_fixture_dependencies()
  scientific <- dplyr::filter(
    audit, .data$role == "CANONICAL_SCIENTIFIC_COMPUTATION"
  )
  validation <- dplyr::filter(audit, .data$role == "REFERENCE_VALIDATION")
  expect_true(all(scientific$pass))
  expect_false(any(scientific$fixture_ancestor))
  expect_true(all(validation$pass))
  expect_true(all(validation$fixture_ancestor))

  phylogeny <- dplyr::filter(
    scientific, .data$target == "phylogeny_model_sequence_526"
  )
  expect_match(
    phylogeny$direct_upstream_targets,
    "phylogeny_predictor_specification_526",
    fixed = TRUE
  )
  expect_false(grepl("reference", phylogeny$direct_upstream_targets, fixed = TRUE))
})

test_that("the current text report is assembled from active canonical objects", {
  source_text <- paste(
    readLines(project_path("R", "targets", "publication_targets.R"), warn = FALSE),
    collapse = "\n"
  )
  start <- regexpr("current_downstream_model_report,", source_text, fixed = TRUE)[[1]]
  end <- regexpr("current_downstream_model_results_file,", source_text, fixed = TRUE)[[1]]
  expect_gt(start, 0L)
  expect_gt(end, start)
  block <- substring(source_text, start, end - 1L)
  expect_match(block, "mboth_526_refit_bundle", fixed = TRUE)
  expect_match(block, "heterogeneity_synthesis", fixed = TRUE)
  expect_match(block, "active_downstream_release_validation", fixed = TRUE)
  expect_false(grepl("transition_reference", block, fixed = TRUE))
  expect_false(grepl("step22_reference", block, fixed = TRUE))
  expect_false(grepl("archive/", block, fixed = TRUE))
})
