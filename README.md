# Bioinformatics and Statistical Analysis Scripts

A collection of reusable scripts for bioinformatics, statistical analysis, and microbiome research. This repository aims to provide practical and reproducible tools for processing biological data, performing statistical analyses, and visualizing research results.

## Repository Structure

The repository is organized into four main categories according to the purpose of each script.

```text
bioinformatics-scripts/
│
├── README.md
│
├── statistical-analysis/
│   ├── README.md
│   └── data-distribution-assessment/
│       ├── data_distribution_assessment.R
│       ├── input_file.csv
│       └── README.md
│
├── bioinformatics-analysis/
│   └── README.md
│
├── microbiome-analysis/
│   ├── README.md
│   ├── alpha-diversity/
│   ├── beta-diversity/
│   │   ├── pcoa_permanova.R
│   │   └── README.md
│   ├── differential-abundance/
│   ├── clinical-associations/
│   │   ├── survival_analysis_microbiome.R
│   │   └── README.md
│   └── longitudinal-analysis/
│
└── shiny-apps/
    └── data-distribution-assessment/
        ├── app.R
        └── README.md
```

*Note: The directory tree describes the intended organization of the repository. Some directories and scripts may be added as the repository grows.*

## 1. Statistical Analysis

**Path:** `statistical-analysis/`

Scripts for general statistical analysis and exploratory data analysis, applicable to a wide range of research datasets.

### Data Distribution Assessment

**Path:** `statistical-analysis/data-distribution-assessment/`

Tools for exploring data distributions and assessing assumptions before performing statistical analyses.

The planned functionality includes:

- Descriptive statistics.
- Assessment of missing data.
- Distribution visualization using histograms and boxplots.
- Normality assessment using the Shapiro–Wilk test and Q–Q plots.

The example file `input_file.csv` illustrates the expected input format. Check the corresponding README for the exact input requirements.

## 2. Bioinformatics Analysis

**Path:** `bioinformatics-analysis/`

Scripts for common bioinformatics workflows, particularly those involving sequencing data and computational processing of biological datasets.

Potential areas of development include:

- Sequencing quality control.
- FASTQ file processing.
- Sequencing data management.
- Preparation and submission of sequencing data to public repositories, such as the European Nucleotide Archive (ENA).

The availability and requirements of each workflow are described in its corresponding directory.

## 3. Microbiome Analysis

**Path:** `microbiome-analysis/`

Scripts for the statistical analysis and visualization of microbiome data, with a focus on microbial community composition, diversity, clinical associations, and longitudinal study designs.

### Alpha Diversity

**Path:** `microbiome-analysis/alpha-diversity/`

Analysis of within-sample microbial diversity, including metrics such as Shannon diversity, Chao1 richness, and the inverse Simpson index.

### Beta Diversity

**Path:** `microbiome-analysis/beta-diversity/`

Analysis of differences in microbial community composition between samples.

Planned and available workflows may include:

- Principal Coordinates Analysis (PCoA).
- PERMANOVA for assessing associations between community composition and explanatory variables.
- Bray–Curtis dissimilarity.
- Aitchison distance based on centered log-ratio (CLR) transformation.
- Visualization of microbial community structure.

### Differential Abundance

**Path:** `microbiome-analysis/differential-abundance/`

Scripts for identifying microbial taxa associated with clinical variables, experimental conditions, or study groups.

### Clinical Associations

**Path:** `microbiome-analysis/clinical-associations/`

Scripts for investigating associations between microbiome features and clinical outcomes.

Examples include:

- Survival analysis using Kaplan–Meier curves.
- Log-rank tests.
- Cox proportional hazards regression.
- Assessment of associations between microbial diversity and progression-free survival or other time-to-event outcomes.

### Longitudinal Analysis

**Path:** `microbiome-analysis/longitudinal-analysis/`

Tools for studying changes in microbial diversity and community composition over time, accounting for repeated measurements when appropriate.

Possible approaches include mixed-effects models and within-subject comparisons.

## 4. Shiny Applications

**Path:** `shiny-apps/`

Interactive applications built with R and Shiny to make selected analyses accessible through a graphical interface.

### Data Distribution Assessment

**Path:** `shiny-apps/data-distribution-assessment/`

An interactive application for uploading a dataset and exploring its numerical variables through descriptive statistics, distribution plots, missing-data summaries, and normality assessment.

The application can be run locally using R and Shiny. Web access requires deployment to a hosting service.

## Getting Started

### Requirements

Requirements depend on the specific script or application. Most scripts are written in R and may require packages such as:

- `ggplot2`
- `dplyr`
- `vegan`
- `phyloseq`
- `microbiome`
- `survival`
- `survminer`
- `shiny`

Not every script requires all these packages. Consult the README in the relevant directory for the exact dependencies.

### Running an R Script

1. Clone or download the repository.
2. Install the required R packages.
3. Prepare the input data in the format specified by the script.
4. Open the `.R` file in RStudio or another R environment.
5. Update input and output paths if necessary.
6. Run the script.

### Running a Shiny Application

Open R in the application directory and run:

```r
install.packages("shiny")  # If not already installed
shiny::runApp()
```

Install any additional packages required by the application before running it.

## Reproducibility and Data Management

To support reproducible research:

- Each analysis should document its input format, dependencies, execution steps, and expected outputs.
- Example datasets should contain synthetic or appropriately anonymized data.
- File paths and analysis parameters should be clearly identified.
- Statistical methods and relevant assumptions should be documented.
- Study-specific scripts should be adapted and validated before being applied to new datasets.

**Do not upload identifiable patient data, confidential research data, credentials, or other sensitive information to a public repository.**

## Contributing

Suggestions, bug reports, and improvements are welcome. When contributing a new script, please:

1. Place it in the most appropriate category.
2. Use a descriptive filename.
3. Include a README explaining its purpose and requirements.
4. Provide a small example dataset when appropriate and safe to share.
5. Document important parameters, assumptions, and limitations.

## Disclaimer

These scripts are provided for research and educational purposes. They are not guaranteed to be suitable for every dataset or study design. Users are responsible for checking data quality, verifying assumptions, and validating results before drawing scientific conclusions.

## License

A license has not yet been specified. Please check the repository's license information before reusing, modifying, or redistributing its contents.

## Author
Aleix Lluansí

Aleix Lluansi


