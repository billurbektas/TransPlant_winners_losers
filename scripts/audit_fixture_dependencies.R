#!/usr/bin/env Rscript

root <- normalizePath(getwd())
source(file.path(root, "R", "functions", "dependency_audit.R"))
source(file.path(root, "R", "functions", "io_validation.R"))

audit <- audit_canonical_fixture_dependencies()
output <- file.path(root, "output", "audits", "canonical_fixture_dependency_audit.csv")
write_csv_atomic(audit, output)
print(audit, n = Inf, width = Inf)
cat("\nWrote:", output, "\n")
