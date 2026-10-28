# Shiny Applications

A collection of original interactive applications developed in R using the Shiny framework. These applications aim to make statistical analyses and data exploration more accessible through user-friendly graphical interfaces, without requiring users to write R code.

## Repository Structure

```text
shiny-apps/
│
├── README.md
│
└── data-distribution-assessment/
    ├── app.R
    └── README.md
```

*Note: The repository structure may evolve as new applications are developed and added.*

## Available Applications

### 1. Data Distribution Assessment

**Directory:** `data-distribution-assessment/`

An interactive application for exploring numerical variables in a dataset and assessing their distributions.

**Features:**
- Upload datasets in CSV format.
- Preview uploaded data.
- Select numerical variables for analysis.
- Calculate descriptive statistics.
- Visualize distributions using histograms and boxplots.
- Assess normality using the Shapiro–Wilk test.
- Generate Q–Q plots for visual assessment of normality.
- Explore missing data.

**Input:** A CSV file containing the variables to analyze.

**Output:** Descriptive statistics, tables, and interactive visualizations displayed in the application.

Refer to the README in the application directory for detailed instructions and input requirements.

## Getting Started

### Requirements

- [R](https://www.r-project.org/)
- [RStudio](https://posit.co/download/rstudio-desktop/) (recommended)
- The R packages required by each application, such as `shiny` and `ggplot2`.

### Run an Application Locally

1. Download or clone this repository.
2. Open R or RStudio.
3. Navigate to the directory containing the application.
4. Install the required packages if they are not already installed.
5. Run the application.

For example:

```r
install.packages(c("shiny", "ggplot2"))
shiny::runApp("data-distribution-assessment")
```

The application will open in your default web browser.

## Online Access

Applications can also be deployed to a compatible hosting service, allowing users to access them through a web browser without installing R or running code locally.

Online availability depends on whether an application has been deployed. The presence of its source code in this repository does not necessarily mean that a public web version is available.

## Development and Reproducibility

Applications are developed in R and documented to support transparency, reproducibility, and reuse.

Each application should include:
- A description of its purpose and functionality.
- Input data requirements.
- Package dependencies.
- Instructions for running the application.
- Relevant methodological considerations and limitations.

## Disclaimer

These applications are provided for research and educational purposes. Results should be interpreted in the context of the dataset, statistical assumptions, and study design. Users are responsible for validating the results before drawing scientific conclusions.

## License

Please refer to the repository's license information for terms governing the use, modification, and redistribution of these applications.


