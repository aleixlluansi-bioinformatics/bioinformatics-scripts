# ============================================================
# Linear mixed-effects model for repeated measurements
# Outcome: Chao1 diversity
# Exposure: Previous antibiotic use
# Random effect: Subject ID
# ============================================================

# Required package
library(lmerTest)

# ------------------------------------------------------------
# 1. Check data
# ------------------------------------------------------------

required_variables <- c(
  "Chao1",
  "ATB_previous_N",
  "PIC"
)

if (!all(required_variables %in% names(metadatanovaV23_2.2))) {
  stop("One or more required variables are missing from the dataset.")
}

# Keep complete observations for the model
model_data <- metadatanovaV23_2.2[
  complete.cases(
    metadatanovaV23_2.2[, required_variables]
  ),
]

# ------------------------------------------------------------
# 2. Fit the mixed-effects model
# ------------------------------------------------------------

model_chao1 <- lmer(
  Chao1 ~ ATB_previous_N + (1 | PIC),
  data = model_data,
  REML = TRUE
)

# ------------------------------------------------------------
# 3. Display results
# ------------------------------------------------------------

summary(model_chao1)

# Fixed-effect estimates and p-values
model_coefficients <- coef(summary(model_chao1))

print(model_coefficients)

# Extract coefficients associated with previous antibiotic use
antibiotic_results <- model_coefficients[
  grepl("^ATB_previous_N", rownames(model_coefficients)),
  ,
  drop = FALSE
]

print(antibiotic_results)