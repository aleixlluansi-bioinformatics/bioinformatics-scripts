# =============================================================================
# Sample Size Calculator
# =============================================================================
# Purpose:
#   Calculate approximate sample sizes for common research designs.
#
# Supported designs:
#   1. Estimate one proportion
#   2. Estimate one mean
#   3. Compare two independent proportions
#   4. Compare two independent means
#   5. Compare paired means (paired t-test approximation)
#   6. Detect a Pearson correlation
#
# Important:
#   - These calculations depend on the assumptions and effect sizes entered.
#   - They are planning estimates, not substitutes for a protocol-specific
#     statistical analysis plan.
#   - Clustered, repeated-measures, survival, non-inferiority, equivalence,
#     diagnostic-accuracy, and complex regression designs may require dedicated
#     methods.
#   - For two-group comparisons, the default allocation ratio is 1:1.
#
# Requirements:
#   Base R only.
#
# Usage:
#   Run this file in R or RStudio and follow the prompts.
# =============================================================================

# ---------- Utility functions -------------------------------------------------

ask_numeric <- function(prompt, default = NULL, lower = -Inf, upper = Inf,
                        allow_blank = FALSE) {
  repeat {
    suffix <- if (!is.null(default)) paste0(" [default: ", default, "]") else ""
    raw <- readline(paste0(prompt, suffix, ": "))
    if (identical(raw, "") && !is.null(default)) return(as.numeric(default))
    if (identical(raw, "") && allow_blank) return(NA_real_)

    value <- suppressWarnings(as.numeric(raw))
    if (length(value) == 1L && is.finite(value) &&
        value >= lower && value <= upper) {
      return(value)
    }
    cat("Please enter a valid number between ", lower, " and ", upper, ".\n",
        sep = "")
  }
}

ask_choice <- function(prompt, choices) {
  cat("\n", prompt, "\n", sep = "")
  for (i in seq_along(choices)) {
    cat(sprintf("  %d. %s\n", i, choices[i]))
  }
  repeat {
    x <- suppressWarnings(as.integer(readline("Select an option: ")))
    if (!is.na(x) && x %in% seq_along(choices)) return(x)
    cat("Please enter one of the listed option numbers.\n")
  }
}

ask_alpha <- function() {
  ask_numeric("Two-sided significance level alpha", default = 0.05,
              lower = 0.000001, upper = 0.20)
}

ask_power <- function() {
  ask_numeric("Desired statistical power (e.g. 0.80 or 0.90)",
              default = 0.80, lower = 0.50, upper = 0.999)
}

add_dropout <- function(n, dropout_pct) {
  if (dropout_pct < 0 || dropout_pct >= 100) {
    stop("Expected loss percentage must be >= 0 and < 100.")
  }
  ceiling(n / (1 - dropout_pct / 100))
}

print_result <- function(title, n_per_group = NULL, total_n = NULL,
                         details = character()) {
  cat("\n", strrep("=", 68), "\n", title, "\n", strrep("=", 68), "\n",
      sep = "")
  if (!is.null(n_per_group)) {
    cat("Required evaluable participants per group: ", n_per_group, "\n",
        sep = "")
  }
  if (!is.null(total_n)) {
    cat("Required total evaluable participants: ", total_n, "\n", sep = "")
  }
  if (length(details)) cat(paste0(details, collapse = "\n"), "\n")
  cat("\n")
}

# ---------- Select design -----------------------------------------------------

designs <- c(
  "Estimate one population proportion",
  "Estimate one population mean",
  "Compare two independent proportions",
  "Compare two independent means",
  "Compare paired means (before-after / matched pairs)",
  "Detect a Pearson correlation"
)

choice <- ask_choice("Choose the study design:", designs)

# ---------- 1. One proportion -------------------------------------------------

if (choice == 1) {
  p <- ask_numeric("Expected population proportion (e.g. 0.30)",
                   lower = 0.000001, upper = 0.999999)
  margin <- ask_numeric("Desired absolute margin of error (e.g. 0.05)",
                        lower = 0.000001, upper = 0.50)
  conf <- ask_numeric("Confidence level (e.g. 0.95)", default = 0.95,
                      lower = 0.80, upper = 0.9999)
  dropout <- ask_numeric("Expected loss / non-response (%)", default = 10,
                         lower = 0, upper = 99.9)

  z <- qnorm(1 - (1 - conf) / 2)
  n <- ceiling(z^2 * p * (1 - p) / margin^2)
  n_adjusted <- add_dropout(n, dropout)

  print_result(
    "Sample size to estimate one proportion",
    total_n = n_adjusted,
    details = c(
      paste0("Evaluable sample size before losses: ", n),
      paste0("Confidence level: ", conf),
      paste0("Absolute margin of error: ", margin),
      paste0("Expected proportion: ", p),
      paste0("Loss allowance: ", dropout, "%"),
      "Method: normal-approximation formula; use exact methods for small samples."
    )
  )
}

# ---------- 2. One mean -------------------------------------------------------

if (choice == 2) {
  sigma <- ask_numeric("Expected population standard deviation (SD)",
                       lower = 0.000001)
  margin <- ask_numeric("Desired margin of error in original units",
                        lower = 0.000001)
  conf <- ask_numeric("Confidence level (e.g. 0.95)", default = 0.95,
                      lower = 0.80, upper = 0.9999)
  dropout <- ask_numeric("Expected loss / non-response (%)", default = 10,
                         lower = 0, upper = 99.9)

  z <- qnorm(1 - (1 - conf) / 2)
  n <- ceiling((z * sigma / margin)^2)
  n_adjusted <- add_dropout(n, dropout)

  print_result(
    "Sample size to estimate one mean",
    total_n = n_adjusted,
    details = c(
      paste0("Evaluable sample size before losses: ", n),
      paste0("Expected SD: ", sigma),
      paste0("Desired margin of error: ", margin),
      paste0("Confidence level: ", conf),
      paste0("Loss allowance: ", dropout, "%"),
      "Method: normal approximation using an externally estimated SD."
    )
  )
}

# ---------- 3. Two independent proportions ----------------------------------

if (choice == 3) {
  p1 <- ask_numeric("Expected proportion in group 1 (e.g. 0.30)",
                    lower = 0.000001, upper = 0.999999)
  p2 <- ask_numeric("Expected proportion in group 2 (e.g. 0.50)",
                    lower = 0.000001, upper = 0.999999)
  alpha <- ask_alpha()
  power <- ask_power()
  ratio <- ask_numeric("Allocation ratio n2 / n1", default = 1,
                       lower = 0.01, upper = 100)
  dropout <- ask_numeric("Expected loss / non-evaluable (%)", default = 10,
                         lower = 0, upper = 99.9)

  # Two-sided normal approximation for independent proportions.
  p_bar <- (p1 + ratio * p2) / (1 + ratio)
  z_alpha <- qnorm(1 - alpha / 2)
  z_power <- qnorm(power)
  numerator <- (z_alpha * sqrt(p_bar * (1 - p_bar) * (1 + 1 / ratio)) +
                  z_power * sqrt(p1 * (1 - p1) + p2 * (1 - p2) / ratio))^2
  n1 <- ceiling(numerator / (p1 - p2)^2)
  n2 <- ceiling(ratio * n1)
  n1_adj <- add_dropout(n1, dropout)
  n2_adj <- add_dropout(n2, dropout)

  print_result(
    "Comparison of two independent proportions",
    n_per_group = NULL,
    total_n = n1_adj + n2_adj,
    details = c(
      paste0("Group 1 evaluable n: approximately ", n1,
             " (recruit ", n1_adj, ")"),
      paste0("Group 2 evaluable n: approximately ", n2,
             " (recruit ", n2_adj, ")"),
      paste0("Expected proportions: ", p1, " vs ", p2),
      paste0("Absolute difference: ", abs(p1 - p2)),
      paste0("Alpha: ", alpha, "; power: ", power),
      paste0("Allocation ratio n2/n1: ", ratio),
      paste0("Loss allowance: ", dropout, "%"),
      "Method: two-sided normal approximation; confirm with a validated exact or score-based method."
    )
  )
}

# ---------- 4. Two independent means ----------------------------------------

if (choice == 4) {
  mean1 <- ask_numeric("Expected mean in group 1")
  mean2 <- ask_numeric("Expected mean in group 2")
  sd <- ask_numeric("Assumed common within-group SD", lower = 0.000001)
  alpha <- ask_alpha()
  power <- ask_power()
  ratio <- ask_numeric("Allocation ratio n2 / n1", default = 1,
                       lower = 0.01, upper = 100)
  dropout <- ask_numeric("Expected loss / non-evaluable (%)", default = 10,
                         lower = 0, upper = 99.9)

  delta <- abs(mean1 - mean2)
  if (delta == 0) stop("The expected difference must be greater than zero.")
  # Approximate two-sample t-test calculation using the noncentral t distribution.
  n1 <- 2
  repeat {
    n2 <- ceiling(ratio * n1)
    df <- n1 + n2 - 2
    ncp <- delta / (sd * sqrt(1 / n1 + 1 / n2))
    crit <- qt(1 - alpha / 2, df = df)
    achieved <- pt(-crit, df = df, ncp = ncp) +
      1 - pt(crit, df = df, ncp = ncp)
    if (achieved >= power || n1 > 1e6) break
    n1 <- n1 + 1
  }
  n1_adj <- add_dropout(n1, dropout)
  n2_adj <- add_dropout(n2, dropout)

  print_result(
    "Comparison of two independent means",
    total_n = n1_adj + n2_adj,
    details = c(
      paste0("Group 1 evaluable n: ", n1, " (recruit ", n1_adj, ")"),
      paste0("Group 2 evaluable n: ", n2, " (recruit ", n2_adj, ")"),
      paste0("Expected difference: ", delta),
      paste0("Assumed common SD: ", sd),
      paste0("Standardized effect (Cohen's d): ", round(delta / sd, 3)),
      paste0("Alpha: ", alpha, "; target power: ", power,
             "; approximate achieved power: ", round(achieved, 3)),
      paste0("Allocation ratio n2/n1: ", ratio),
      paste0("Loss allowance: ", dropout, "%"),
      "Assumptions: independent groups, common SD, two-sided test, approximately normal outcome."
    )
  )
}

# ---------- 5. Paired means ---------------------------------------------------

if (choice == 5) {
  delta <- ask_numeric("Expected mean within-person difference (absolute value)",
                       lower = 0.000001)
  sd_diff <- ask_numeric("SD of within-person differences",
                         lower = 0.000001)
  alpha <- ask_alpha()
  power <- ask_power()
  dropout <- ask_numeric("Expected loss / non-evaluable pairs (%)",
                         default = 10, lower = 0, upper = 99.9)

  # Approximate paired t-test calculation via noncentral t distribution.
  n <- 3
  repeat {
    df <- n - 1
    ncp <- (delta / sd_diff) * sqrt(n)
    crit <- qt(1 - alpha / 2, df = df)
    achieved <- pt(-crit, df = df, ncp = ncp) +
      1 - pt(crit, df = df, ncp = ncp)
    if (achieved >= power || n > 1e6) break
    n <- n + 1
  }
  n_adj <- add_dropout(n, dropout)

  print_result(
    "Paired means / within-person change",
    total_n = n_adj,
    details = c(
      paste0("Required evaluable pairs: ", n),
      paste0("Pairs to recruit after loss allowance: ", n_adj),
      paste0("Expected mean difference: ", delta),
      paste0("SD of paired differences: ", sd_diff),
      paste0("Standardized paired effect (dz): ", round(delta / sd_diff, 3)),
      paste0("Alpha: ", alpha, "; target power: ", power,
             "; approximate achieved power: ", round(achieved, 3)),
      paste0("Loss allowance: ", dropout, "%"),
      "Use the SD of within-person differences, not the SD at either time point."
    )
  )
}

# ---------- 6. Pearson correlation ------------------------------------------

if (choice == 6) {
  r <- ask_numeric("Smallest correlation to detect (e.g. 0.30)",
                   lower = -0.999, upper = 0.999)
  if (r == 0) stop("A correlation of zero cannot be distinguished from zero using this approach.")
  alpha <- ask_alpha()
  power <- ask_power()
  dropout <- ask_numeric("Expected loss / non-evaluable (%)", default = 10,
                         lower = 0, upper = 99.9)

  # Fisher z approximation for testing H0: rho = 0.
  z_r <- atanh(r)
  n <- ceiling(3 + ((qnorm(1 - alpha / 2) + qnorm(power)) / z_r)^2)
  n_adj <- add_dropout(n, dropout)

  print_result(
    "Sample size to detect a Pearson correlation",
    total_n = n_adj,
    details = c(
      paste0("Evaluable sample size before losses: ", n),
      paste0("Correlation to detect: ", r),
      paste0("Alpha: ", alpha, "; target power: ", power),
      paste0("Loss allowance: ", dropout, "%"),
      "Method: Fisher z approximation for testing a correlation against zero.",
      "Assumptions: independent observations and a suitable bivariate-normal approximation."
    )
  )
}

cat("Reminder: document the assumptions and justify the effect size using prior evidence or pilot data.\n")
