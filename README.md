# R Statistical Tools

A small collection of base-R scripts for common research tasks.

## Scripts

### `sample_size_calculator.R`
Interactive sample-size estimates for:
- Estimation of one proportion
- Estimation of one mean
- Comparison of two independent proportions
- Comparison of two independent means
- Paired means
- Pearson correlation

The calculations are planning approximations. Select assumptions from prior evidence or pilot data and confirm protocol-specific calculations with a statistician. This script does not cover every study design.

### `data_distribution_assessment.R`
Imports a CSV file and creates:
- Descriptive statistics
- Histograms and normal Q-Q plots
- Shapiro-Wilk results when the sample size is within the test's supported range
- Exploratory outlier flags using the 1.5 × IQR rule
- A plain-language interpretation guide

## Requirements

- R (base installation; no additional packages required)

## How to use

1. Download the `.R` file(s) or clone this repository.
2. Open a script in RStudio or another R environment.
3. Run the script and follow the prompts.
4. For the distribution assessment, select a CSV file containing numeric columns.

## Important notes

- These scripts are educational and research-support tools, not substitutes for a prespecified statistical analysis plan.
- Sample-size results depend on the chosen effect size, variability, alpha, power, allocation ratio, and loss assumptions.
- Normality tests should not be interpreted in isolation. Inspect plots and consider the assumptions of the intended model.
- Outlier flags are exploratory and do not imply that observations should be removed.
- Do not use identifiable or confidential participant data in public repositories or example files.

## Author

Add your name and contact/profile information if desired.

## License

Choose and add a license before inviting reuse. The appropriate license may depend on institutional or project requirements.
