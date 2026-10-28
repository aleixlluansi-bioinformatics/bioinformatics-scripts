# ============================================================
# PAIRED COMPARISON OF ALPHA DIVERSITY
# Outcome: Chao1
# Comparison: Baseline (Time = 0) vs EOS (Time = 1)
# Groups: Bread
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# ------------------------------------------------------------
# 1. Prepare the data
# ------------------------------------------------------------

# Replace ID_PARTICIPANTE with the actual participant ID column.
id_variable <- "ID_PARTICIPANTE"

required_variables <- c(
  id_variable,
  "Bread",
  "Time",
  "Chao1"
)

if (!all(required_variables %in% names(sampledata.phylum))) {
  stop("One or more required variables are missing.")
}

chao1_data <- sampledata.phylum %>%
  select(all_of(required_variables)) %>%
  filter(
    Time %in% c("0", "1"),
    !is.na(Chao1),
    !is.na(Bread),
    !is.na(.data[[id_variable]])
  ) %>%
  mutate(
    Time = as.character(Time),
    Bread = as.factor(Bread)
  )

# ------------------------------------------------------------
# 2. Keep participants with both measurements
# ------------------------------------------------------------

paired_data <- chao1_data %>%
  group_by(across(all_of(c(id_variable, "Bread")))) %>%
  filter(
    n_distinct(Time) == 2,
    all(c("0", "1") %in% Time)
  ) %>%
  ungroup()

# Check that each participant has one observation per time point.
duplicate_measurements <- paired_data %>%
  count(across(all_of(c(id_variable, "Bread", "Time")))) %>%
  filter(n > 1)

if (nrow(duplicate_measurements) > 0) {
  stop(
    "Some participants have multiple measurements per time point. ",
    "Resolve these duplicates before proceeding."
  )
}

# ------------------------------------------------------------
# 3. Paired tests by group
# ------------------------------------------------------------

run_paired_test <- function(data, outcome, group_name) {
  
  wide_data <- data %>%
    select(
      all_of(id_variable),
      Time,
      all_of(outcome)
    ) %>%
    pivot_wider(
      names_from = Time,
      values_from = all_of(outcome),
      names_prefix = "Time_"
    ) %>%
    filter(
      !is.na(Time_0),
      !is.na(Time_1)
    )
  
  if (nrow(wide_data) < 2) {
    warning("Insufficient paired observations in group: ", group_name)
    return(NULL)
  }
  
  # Difference: EOS - Baseline
  differences <- wide_data$Time_1 - wide_data$Time_0
  
  # Test normality of paired differences
  shapiro_p <- if (
    length(differences) >= 3 &&
    length(differences) <= 5000 &&
    length(unique(differences)) > 1
  ) {
    shapiro.test(differences)$p.value
  } else {
    NA_real_
  }
  
  # Select paired test
  if (!is.na(shapiro_p) && shapiro_p >= 0.05) {
    
    test <- t.test(
      wide_data$Time_1,
      wide_data$Time_0,
      paired = TRUE
    )
    
    test_name <- "Paired t-test"
    
  } else {
    
    test <- wilcox.test(
      wide_data$Time_1,
      wide_data$Time_0,
      paired = TRUE,
      exact = FALSE
    )
    
    test_name <- "Wilcoxon signed-rank test"
  }
  
  data.frame(
    Outcome = outcome,
    Group = group_name,
    N_pairs = nrow(wide_data),
    Mean_baseline = mean(wide_data$Time_0),
    Mean_EOS = mean(wide_data$Time_1),
    Mean_difference = mean(differences),
    Shapiro_p = shapiro_p,
    Test = test_name,
    P_value = test$p.value
  )
}

chao1_results <- lapply(
  levels(droplevels(paired_data$Bread)),
  function(group_name) {
    
    group_data <- paired_data %>%
      filter(Bread == group_name)
    
    run_paired_test(
      data = group_data,
      outcome = "Chao1",
      group_name = group_name
    )
  }
)

chao1_results <- bind_rows(chao1_results)

# Correct for testing both groups
chao1_results$P_adjusted_BH <- p.adjust(
  chao1_results$P_value,
  method = "BH"
)

print(chao1_results)

# Export results
dir.create("results", showWarnings = FALSE)

write.csv(
  chao1_results,
  "results/paired_chao1_results.csv",
  row.names = FALSE
)




# ============================================================
# PAIRED BOXPLOT: CHAO1 BY GROUP AND TIME
# ============================================================

plot_data <- paired_data %>%
  mutate(
    Time = factor(
      Time,
      levels = c("0", "1"),
      labels = c("Baseline W0", "EOS W8")
    )
  )

chao1_plot <- ggpaired(
  plot_data,
  x = "Time",
  y = "Chao1",
  id = id_variable,
  color = "Time",
  line.color = "black",
  line.size = 0.4,
  palette = c("#F8766D", "#619CFF"),
  facet.by = "Bread"
) +
  labs(
    title = "Chao1 index",
    x = "Time",
    y = "Chao1 index"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(
      size = 11,
      color = "black",
      face = "bold"
    ),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

print(chao1_plot)

ggsave(
  "results/paired_chao1_plot.png",
  plot = chao1_plot,
  width = 8,
  height = 5,
  dpi = 300
)