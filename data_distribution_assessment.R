# =============================================================================
# Data Distribution Assessment
# =============================================================================
# Purpose:
#   Explore numeric variables, summarize their distributions, create diagnostic
#   plots, and run Shapiro-Wilk tests where appropriate.
#
# Requirements:
#   Base R only.
#
# Usage:
#   1. Run the script in R / RStudio.
#   2. Select a CSV file when prompted.
#   3. Select the numeric variables to assess.
#   4. Results are written to a new output folder beside the input file.
#
# Important interpretation notes:
#   - Shapiro-Wilk tests whether data are compatible with a normal distribution;
#     it does not prove normality when p >= 0.05.
#   - In large samples, small deviations can produce significant p-values.
#   - In linear models, normality assumptions usually concern model residuals,
#     not necessarily the raw outcome or predictors.
#   - Repeated measures and clustered data are not independent; assess them
#     using a design-appropriate model and diagnostic approach.
#   - Outlier flags are exploratory and should not be used to remove values
#     automatically.
# =============================================================================

# ---------- Select input file ------------------------------------------------

input_file <- file.choose()
if (!file.exists(input_file)) stop("Input file does not exist.")

data <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(data) == 0L) stop("The selected file has no rows.")
if (ncol(data) == 0L) stop("The selected file has no columns.")

cat("\nLoaded file:", basename(input_file), "\n")
cat("Rows:", nrow(data), "| Columns:", ncol(data), "\n")

# ---------- Select numeric variables ----------------------------------------

numeric_cols <- names(data)[vapply(data, is.numeric, logical(1))]
if (length(numeric_cols) == 0L) {
  stop("No numeric columns were detected. Check the CSV column types.")
}

cat("\nNumeric columns detected:\n")
for (i in seq_along(numeric_cols)) {
  cat(sprintf("  %d. %s\n", i, numeric_cols[i]))
}
cat("Enter variable numbers separated by commas, or press Enter for all numeric columns.\n")

selection <- trimws(readline("Variables to assess: "))
if (selection == "") {
  selected_cols <- numeric_cols
} else {
  ids <- suppressWarnings(as.integer(trimws(strsplit(selection, ",")[[1]])))
  if (anyNA(ids) || any(!ids %in% seq_along(numeric_cols))) {
    stop("Invalid selection. Run the script again and enter valid variable numbers.")
  }
  selected_cols <- unique(numeric_cols[ids])
}

# ---------- Output folder ----------------------------------------------------

input_dir <- dirname(normalizePath(input_file))
input_stem <- tools::file_path_sans_ext(basename(input_file))
output_dir <- file.path(input_dir, paste0(input_stem, "_distribution_report"))
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- Descriptive statistics ------------------------------------------

safe_skewness <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3L) return(NA_real_)
  s <- sd(x)
  if (!is.finite(s) || s == 0) return(NA_real_)
  mean(((x - mean(x)) / s)^3)
}

safe_shapiro <- function(x) {
  x <- x[is.finite(x)]
  # R's shapiro.test requires between 3 and 5000 non-missing observations.
  if (length(x) < 3L) return(list(statistic = NA_real_, p_value = NA_real_,
                                  note = "Fewer than 3 non-missing values"))
  if (length(x) > 5000L) return(list(statistic = NA_real_, p_value = NA_real_,
                                     note = "More than 5000 values; test not run"))
  if (length(unique(x)) < 2L) return(list(statistic = NA_real_, p_value = NA_real_,
                                          note = "No variation"))
  result <- tryCatch(shapiro.test(x), error = function(e) NULL)
  if (is.null(result)) {
    return(list(statistic = NA_real_, p_value = NA_real_,
                note = "Test could not be calculated"))
  }
  list(statistic = unname(result$statistic),
       p_value = result$p.value,
       note = "")
}

summaries <- lapply(selected_cols, function(variable) {
  x_all <- data[[variable]]
  x <- x_all[is.finite(x_all)]
  n_total <- length(x_all)
  n_valid <- length(x)
  n_missing <- sum(is.na(x_all) | !is.finite(x_all))
  sh <- safe_shapiro(x)

  q <- if (n_valid > 0) {
    quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1), names = FALSE)
  } else {
    rep(NA_real_, 5)
  }

  data.frame(
    variable = variable,
    n_total = n_total,
    n_valid = n_valid,
    n_missing_or_nonfinite = n_missing,
    mean = if (n_valid > 0) mean(x) else NA_real_,
    sd = if (n_valid > 1) sd(x) else NA_real_,
    median = if (n_valid > 0) median(x) else NA_real_,
    q1 = q[2],
    q3 = q[4],
    min = q[1],
    max = q[5],
    skewness = safe_skewness(x),
    shapiro_w = sh$statistic,
    shapiro_p = sh$p_value,
    shapiro_note = sh$note,
    stringsAsFactors = FALSE
  )
})

summary_table <- do.call(rbind, summaries)
write.csv(summary_table, file.path(output_dir, "descriptive_statistics.csv"),
          row.names = FALSE, na = "")

# ---------- Plot each variable -----------------------------------------------

for (variable in selected_cols) {
  x_all <- data[[variable]]
  x <- x_all[is.finite(x_all)]
  if (length(x) == 0L) next

  safe_name <- gsub("[^A-Za-z0-9_-]", "_", variable)
  pdf(file.path(output_dir, paste0(safe_name, "_distribution.pdf")),
      width = 8, height = 6)
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  par(mfrow = c(1, 2))
  hist(x, breaks = "FD", main = paste("Histogram:", variable),
       xlab = variable, col = "grey85", border = "white")
  if (length(unique(x)) > 1L) {
    qqnorm(x, main = paste("Normal Q-Q plot:", variable))
    qqline(x, col = "red", lwd = 2)
  } else {
    plot.new()
    title(main = paste("Normal Q-Q plot:", variable))
    text(0.5, 0.5, "No variation: Q-Q plot not informative")
  }
  dev.off()
}

# ---------- Exploratory outlier flags using the 1.5 x IQR rule ---------------

outlier_rows <- list()
for (variable in selected_cols) {
  x <- data[[variable]]
  valid <- is.finite(x)
  if (sum(valid) < 4L) next

  q <- quantile(x[valid], probs = c(0.25, 0.75), na.rm = TRUE)
  iqr_value <- q[2] - q[1]
  if (iqr_value == 0) next

  low <- q[1] - 1.5 * iqr_value
  high <- q[2] + 1.5 * iqr_value
  idx <- which(valid & (x < low | x > high))

  if (length(idx) > 0L) {
    outlier_rows[[length(outlier_rows) + 1L]] <- data.frame(
      row_number = idx,
      variable = variable,
      value = x[idx],
      lower_fence = low,
      upper_fence = high,
      stringsAsFactors = FALSE
    )
  }
}

if (length(outlier_rows) > 0L) {
  outliers <- do.call(rbind, outlier_rows)
} else {
  outliers <- data.frame(
    row_number = integer(0),
    variable = character(0),
    value = numeric(0),
    lower_fence = numeric(0),
    upper_fence = numeric(0)
  )
}
write.csv(outliers, file.path(output_dir, "exploratory_outlier_flags.csv"),
          row.names = FALSE, na = "")

# ---------- Plain-language interpretation guide -----------------------------

guide <- c(
  "DATA DISTRIBUTION ASSESSMENT — INTERPRETATION GUIDE",
  "",
  "1. Review histograms and Q-Q plots alongside descriptive statistics.",
  "2. Shapiro-Wilk p < 0.05 indicates evidence against a normal distribution;",
  "   p >= 0.05 does not establish that the data are normal.",
  "3. Consider sample size, the scientific context, influential observations,",
  "   and the robustness of the intended statistical method.",
  "4. For linear models, assess residuals and model assumptions rather than",
  "   requiring every raw variable to be normally distributed.",
  "5. Outlier flags use the 1.5 x IQR rule and are exploratory, not automatic",
  "   grounds for deleting observations.",
  "6. For repeated-measures, clustered, or longitudinal data, use diagnostics",
  "   appropriate to the model and dependence structure.",
  "",
  paste("Input file:", basename(input_file)),
  paste("Variables assessed:", paste(selected_cols, collapse = ", "))
)
writeLines(guide, file.path(output_dir, "interpretation_guide.txt"))

cat("\nAssessment complete.\n")
cat("Output folder:\n", output_dir, "\n\n")
cat("Files created:\n")
cat("- descriptive_statistics.csv\n")
cat("- one PDF per variable (histogram and Q-Q plot)\n")
cat("- exploratory_outlier_flags.csv\n")
cat("- interpretation_guide.txt\n")
