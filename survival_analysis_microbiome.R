# ============================================================
# Survival analysis of baseline alpha diversity
# Kaplan-Meier, log-rank test and univariable Cox regression
# ============================================================


# Description:
#   Evaluate the association between baseline alpha diversity
#   indices and progression-free survival (PFS).
#
# Analyses:
#   1. Kaplan-Meier curves using median-based groups
#   2. Log-rank test to compare survival curves
#   3. Univariable Cox regression using continuous indices
#   4. Proportional hazards assumption using Schoenfeld residuals
#
# Required input:
#   ps_heces.rds or a previously generated df_div object
#
# Required variables:
#   Tiempo.x, CodigoPaciente.x, PFS_time, Evento_PFS,
#   Shannon, Chao1, InvSimpson
#
# IMPORTANT:
#   Confirm that PFS_time is measured in days and that
#   Evento_PFS is coded as 1 = event and 0 = censored.
# ============================================================


# 1. Load packages --------------------------------------------

required_packages <- c(
  "dplyr",
  "survival",
  "survminer",
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running the script: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(dplyr)
library(survival)
library(survminer)
library(ggplot2)


# 2. Define input and output paths ----------------------------

# Option A: Load a previously generated alpha-diversity table.
# This file must contain one row per sample and the required
# clinical and alpha-diversity variables.

input_file <- "df_div.csv"

output_dir <- "results/survival_analysis"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# 3. Load data ------------------------------------------------

if (!file.exists(input_file)) {
  stop(
    "Input file not found: ", input_file,
    "\nSet input_file to the path of your df_div.csv file."
  )
}

df_div <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# 4. Check required variables ---------------------------------

required_variables <- c(
  "Tiempo.x",
  "CodigoPaciente.x",
  "PFS_time",
  "Evento_PFS",
  "Shannon",
  "Chao1",
  "InvSimpson"
)

missing_variables <- setdiff(
  required_variables,
  colnames(df_div)
)

if (length(missing_variables) > 0) {
  stop(
    "Missing required variables: ",
    paste(missing_variables, collapse = ", ")
  )
}


# 5. Select baseline samples ----------------------------------

df_baseline <- df_div %>%
  filter(Tiempo.x == "PreC1")

if (nrow(df_baseline) == 0) {
  stop("No baseline samples found for Tiempo.x == 'PreC1'.")
}

# Ensure that each patient contributes only one baseline row.
if (anyDuplicated(df_baseline$CodigoPaciente.x)) {
  stop(
    "Multiple baseline rows detected for at least one patient. ",
    "Resolve duplicate baseline samples before survival analysis."
  )
}


# 6. Prepare survival variables -------------------------------

# Convert variables to numeric when necessary.
numeric_variables <- c(
  "PFS_time",
  "Evento_PFS",
  "Shannon",
  "Chao1",
  "InvSimpson"
)

for (variable in numeric_variables) {
  df_baseline[[variable]] <- suppressWarnings(
    as.numeric(as.character(df_baseline[[variable]]))
  )
}

# Check event coding.
observed_events <- unique(
  na.omit(df_baseline$Evento_PFS)
)

if (!all(observed_events %in% c(0, 1))) {
  stop(
    "Evento_PFS must be coded as 0 = censored and 1 = event."
  )
}

if (length(observed_events) < 2) {
  warning(
    "Only one event category is observed. ",
    "Survival comparisons may not be informative."
  )
}

# Check survival times.
if (any(df_baseline$PFS_time <= 0, na.rm = TRUE)) {
  stop(
    "PFS_time contains zero or negative values. ",
    "Check the time origin and survival-time calculation."
  )
}

if (all(is.na(df_baseline$PFS_time))) {
  stop("PFS_time contains no valid numeric observations.")
}


# 7. Define alpha-diversity indices ---------------------------

diversity_indices <- c(
  "Shannon",
  "Chao1",
  "InvSimpson"
)


# 8. Create complete results tables ---------------------------

cox_results <- data.frame()
logrank_results <- data.frame()
ph_results <- data.frame()


# 9. Run survival analyses ------------------------------------

for (index in diversity_indices) {
  
  message("Analysing index: ", index)
  
  # ----------------------------------------------------------
  # 9.1 Prepare complete-case data for this index
  # ----------------------------------------------------------
  
  df_model <- df_baseline %>%
    filter(
      !is.na(PFS_time),
      !is.na(Evento_PFS),
      !is.na(.data[[index]])
    )
  
  if (nrow(df_model) < 5) {
    warning(
      "Too few complete observations for ", index,
      ". Skipping this index."
    )
    next
  }
  
  if (length(unique(df_model[[index]])) < 2) {
    warning(
      "Insufficient variation in ", index,
      ". Skipping this index."
    )
    next
  }
  
  # ----------------------------------------------------------
  # 9.2 Kaplan-Meier and log-rank test
  # ----------------------------------------------------------
  
  # Exploratory grouping based on the median.
  median_value <- median(
    df_model[[index]],
    na.rm = TRUE
  )
  
  group_variable <- paste0(index, "_group")
  
  df_model[[group_variable]] <- factor(
    ifelse(
      df_model[[index]] > median_value,
      "High",
      "Low"
    ),
    levels = c("Low", "High")
  )
  
  # Check that both groups are represented.
  if (length(unique(df_model[[group_variable]])) < 2) {
    warning(
      "Only one group is represented for ", index,
      ". Kaplan-Meier analysis skipped."
    )
  } else {
    
    # Kaplan-Meier curves.
    km_formula <- as.formula(
      paste0(
        "Surv(PFS_time, Evento_PFS) ~ ",
        group_variable
      )
    )
    
    km_fit <- survfit(
      km_formula,
      data = df_model
    )
    
    # Explicit log-rank test.
    logrank_test <- survdiff(
      km_formula,
      data = df_model,
      rho = 0
    )
    
    logrank_p <- pchisq(
      logrank_test$chisq,
      df = length(logrank_test$n) - 1,
      lower.tail = FALSE
    )
    
    logrank_results <- rbind(
      logrank_results,
      data.frame(
        Index = index,
        N = nrow(df_model),
        N_low = sum(df_model[[group_variable]] == "Low"),
        N_high = sum(df_model[[group_variable]] == "High"),
        Median_cutoff = median_value,
        Chi_square = unname(logrank_test$chisq),
        Degrees_freedom = length(logrank_test$n) - 1,
        P_value = logrank_p
      )
    )
    
    # Plot Kaplan-Meier curves.
    km_plot <- ggsurvplot(
      km_fit,
      data = df_model,
      pval = TRUE,
      conf.int = TRUE,
      risk.table = TRUE,
      risk.table.height = 0.25,
      legend.title = index,
      legend.labs = c("Low", "High"),
      xlab = "Time",
      ylab = "Progression-free survival probability",
      title = paste(
        "Baseline", index,
        "and progression-free survival"
      ),
      ggtheme = theme_bw()
    )
    
    # Save the plot, including the risk table.
    pdf_file <- file.path(
      output_dir,
      paste0("KM_", index, ".pdf")
    )
    
    pdf(
      pdf_file,
      width = 8,
      height = 8
    )
    
    print(km_plot)
    dev.off()
    
    # Save the plot as a PNG as well.
    png_file <- file.path(
      output_dir,
      paste0("KM_", index, ".png")
    )
    
    png(
      png_file,
      width = 2400,
      height = 2400,
      res = 300
    )
    
    print(km_plot)
    dev.off()
  }
  
  
  # ----------------------------------------------------------
  # 9.3 Univariable Cox regression
  # ----------------------------------------------------------
  
  # Alpha-diversity index entered as a continuous predictor.
  cox_formula <- as.formula(
    paste0(
      "Surv(PFS_time, Evento_PFS) ~ ",
      index
    )
  )
  
  cox_model <- coxph(
    cox_formula,
    data = df_model,
    ties = "efron"
  )
  
  cox_summary <- summary(cox_model)
  
  cox_results <- rbind(
    cox_results,
    data.frame(
      Index = index,
      N = nrow(df_model),
      Events = sum(df_model$Evento_PFS == 1),
      HR = unname(cox_summary$coefficients[1, "exp(coef)"]),
      CI_lower = unname(cox_summary$conf.int[1, "lower .95"]),
      CI_upper = unname(cox_summary$conf.int[1, "upper .95"]),
      P_value = unname(cox_summary$coefficients[1, "Pr(>|z|)"])
    )
  )
  
  
  # ----------------------------------------------------------
  # 9.4 Check proportional hazards assumption
  # ----------------------------------------------------------
  
  ph_test <- cox.zph(cox_model)
  
  ph_table <- as.data.frame(ph_test$table)
  ph_table$Term <- rownames(ph_table)
  
  ph_table$Index <- index
  
  ph_results <- rbind(
    ph_results,
    ph_table
  )
  
  # Save diagnostic plot.
  ph_plot_file <- file.path(
    output_dir,
    paste0("PH_assumption_", index, ".pdf")
  )
  
  pdf(
    ph_plot_file,
    width = 7,
    height = 5
  )
  
  plot(ph_test)
  dev.off()
}


# 10. Multiple-testing correction -----------------------------

# Apply Benjamini-Hochberg correction separately to the
# three Cox tests and the three log-rank tests.

if (nrow(cox_results) > 0) {
  cox_results$P_adjusted_BH <- p.adjust(
    cox_results$P_value,
    method = "BH"
  )
}

if (nrow(logrank_results) > 0) {
  logrank_results$P_adjusted_BH <- p.adjust(
    logrank_results$P_value,
    method = "BH"
  )
}


# 11. Export results -------------------------------------------

write.csv(
  cox_results,
  file.path(output_dir, "cox_regression_results.csv"),
  row.names = FALSE
)

write.csv(
  logrank_results,
  file.path(output_dir, "logrank_results.csv"),
  row.names = FALSE
)

write.csv(
  ph_results,
  file.path(output_dir, "cox_ph_assumption_results.csv"),
  row.names = FALSE
)


# 12. Print results --------------------------------------------

cat("\n==============================\n")
cat("COX REGRESSION RESULTS\n")
cat("==============================\n")
print(cox_results)

cat("\n==============================\n")
cat("LOG-RANK TEST RESULTS\n")
cat("==============================\n")
print(logrank_results)

cat("\n==============================\n")
cat("PROPORTIONAL HAZARDS TEST\n")
cat("==============================\n")
print(ph_results)

cat(
  "\nAnalysis completed. Results saved in:",
  output_dir,
  "\n"
)