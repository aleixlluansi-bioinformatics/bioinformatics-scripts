# ============================================================
# Logistic regression analyses
# Univariate GLM, Firth logistic regression and multivariable GLM
# ============================================================

library(broom)
library(logistf)
library(writexl)

# Dataset
data <- metadatanovaV23_2.2

# Output directory
dir.create("results", showWarnings = FALSE)

# ============================================================
# UNIVARIATE LOGISTIC REGRESSION
# Outcome: Binary_PWG23_2
# ============================================================

predictors <- c(
  "SEX", "ATB_0", "SMK", "SMK_O", "ALC", "HSF",
  "HCH_P", "SES", "WGT_0", "FDM_0",
  "Binary_Breastfeeding", "NSP_1D", "FMF", "DCA",
  "PVC60dias_newprova", "FVF", "AOM", "URI", "GII",
  "ATB", "CST", "URI_N3", "URI_N2", "AGE_G", "AGE_M",
  "HCH_N", "FVF_N", "HCH_NS", "AOM_N", "GII_N",
  "ATB_N", "CST_N", "CST_N2", "CST_N3"
)

outcome <- "Binary_PWG23_2"

if (!outcome %in% names(data)) {
  stop("The outcome variable is missing from the dataset.")
}

missing_predictors <- setdiff(predictors, names(data))

if (length(missing_predictors) > 0) {
  warning(
    "The following predictors are missing: ",
    paste(missing_predictors, collapse = ", ")
  )
}

predictors <- intersect(predictors, names(data))

univariate_results <- list()

for (predictor in predictors) {
  
  formula_model <- reformulate(
    predictor,
    response = outcome
  )
  
  logistic_model <- tryCatch(
    glm(
      formula_model,
      family = binomial(link = "logit"),
      data = data
    ),
    error = function(e) {
      message("Model failed for ", predictor, ": ", e$message)
      NULL
    }
  )
  
  if (is.null(logistic_model)) next
  
  model_results <- broom::tidy(
    logistic_model,
    conf.int = TRUE,
    exponentiate = TRUE
  )
  
  # Exclude the intercept
  model_results <- model_results[
    model_results$term != "(Intercept)",
    ,
    drop = FALSE
  ]
  
  if (nrow(model_results) == 0) next
  
  model_results$Predictor <- predictor
  
  univariate_results[[predictor]] <- model_results
}

univariate_results_df <- do.call(
  rbind,
  univariate_results
)

# Adjust p-values across all tested coefficients
univariate_results_df$p_adjusted_Bonferroni <- p.adjust(
  univariate_results_df$p.value,
  method = "bonferroni"
)

write_xlsx(
  univariate_results_df,
  "results/univariate_logistic_regression.xlsx"
)

print(univariate_results_df)


# ============================================================
# FIRTH PENALISED LOGISTIC REGRESSION
# Outcome: Binary_PWG23_2
# ============================================================

firth_results <- list()

for (predictor in predictors) {
  
  formula_model <- reformulate(
    predictor,
    response = outcome
  )
  
  firth_model <- tryCatch(
    logistf(
      formula_model,
      data = data
    ),
    error = function(e) {
      message("Firth model failed for ", predictor, ": ", e$message)
      NULL
    }
  )
  
  if (is.null(firth_model)) next
  
  coefficients <- firth_model$coef
  
  coef_data <- data.frame(
    term = names(coefficients),
    estimate = exp(unname(coefficients)),
    conf.low = exp(unname(firth_model$ci.lower)),
    conf.high = exp(unname(firth_model$ci.upper)),
    p.value = unname(firth_model$prob),
    stringsAsFactors = FALSE
  )
  
  # Exclude the intercept
  coef_data <- coef_data[
    coef_data$term != "(Intercept)",
    ,
    drop = FALSE
  ]
  
  if (nrow(coef_data) == 0) next
  
  coef_data$Predictor <- predictor
  
  firth_results[[predictor]] <- coef_data
}

if (length(firth_results) > 0) {
  
  firth_results_df <- do.call(rbind, firth_results)
  
  # Adjust p-values across all tested coefficients
  firth_results_df$p_adjusted_Bonferroni <- p.adjust(
    firth_results_df$p.value,
    method = "bonferroni"
  )
  
  write_xlsx(
    firth_results_df,
    "results/univariate_firth_logistic_regression.xlsx"
  )
  
  print(firth_results_df)
}



# ============================================================
# MULTIVARIABLE LOGISTIC REGRESSION
# Outcome: Binary2_RPWG4
# ============================================================

multivariable_formula <- Binary2_RPWG4 ~
  Breastfeeding + AGE_G + AGE_M + WGT_0 + GII_N

multivariable_model <- glm(
  multivariable_formula,
  family = binomial(link = "logit"),
  data = data
)

# Model summary
summary(multivariable_model)

# Odds ratios and 95% confidence intervals
multivariable_results <- broom::tidy(
  multivariable_model,
  conf.int = TRUE,
  exponentiate = TRUE
)

# Exclude intercept
multivariable_results <- multivariable_results[
  multivariable_results$term != "(Intercept)",
  ,
  drop = FALSE
]

write_xlsx(
  multivariable_results,
  "results/multivariable_logistic_regression.xlsx"
)

print(multivariable_results)