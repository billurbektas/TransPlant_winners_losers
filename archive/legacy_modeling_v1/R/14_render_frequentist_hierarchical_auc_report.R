# Render the resident continuous-AUC reports from saved outputs without refits.

source("config/analysis.R")
source("R/pipeline/utils.R")
source("R/pipeline/frequentist_hierarchical_auc.R")
source("R/pipeline/frequentist_hierarchical_auc_reporting.R")

output_dir <- file.path("output", "frequentist_hierarchical_auc")
audit_report <- render_frequentist_hierarchical_auc_data_audit(output_dir)
visual_report <- render_frequentist_hierarchical_auc_visual_report(output_dir)

cat("Data audit: ", audit_report, "\n", sep = "")
cat("Visual report: ", visual_report, "\n", sep = "")
