# ==============================================================================
# Resident continuous-AUC frequentist hierarchical analysis
#
# This runner is deliberately separate from _targets.R and from the existing
# winner/loser/colonizer analyses. It consumes their preserved continuous AUC
# contrasts but never uses significance categories as the response.
# ==============================================================================

source("config/analysis.R")
source("R/pipeline/utils.R")
source("R/pipeline/frequentist_hierarchical_auc.R")
source("R/pipeline/frequentist_hierarchical_auc_reporting.R")

primary_horizon <- Sys.getenv("TRANSPLANT_AUC_PRIMARY_HORIZON", "max")
comparison_horizon <- Sys.getenv("TRANSPLANT_AUC_COMPARISON_HORIZON", "3")
species_display_min <- as.integer(
  Sys.getenv("TRANSPLANT_AUC_SPECIES_DISPLAY_MIN", "3")
)
run_leave_one_out <- !tolower(
  Sys.getenv("TRANSPLANT_AUC_SKIP_LOEO", "false")
) %in% c("1", "true", "yes")
render_visual_report <- !tolower(
  Sys.getenv("TRANSPLANT_AUC_SKIP_VISUAL_REPORT", "false")
) %in% c("1", "true", "yes")

cat("Running resident continuous-AUC hierarchical analysis...\n")
result <- run_frequentist_hierarchical_auc(
  primary_horizon = primary_horizon,
  comparison_horizon = comparison_horizon,
  species_display_min = species_display_min,
  output_dir = file.path("output", "frequentist_hierarchical_auc"),
  run_leave_one_out = run_leave_one_out,
  render_visual_report = render_visual_report
)
cat("Completed. Report: ", result$report, "\n", sep = "")
if (isTRUE(render_visual_report)) {
  cat("Visual report: ", result$visual_report, "\n", sep = "")
}
