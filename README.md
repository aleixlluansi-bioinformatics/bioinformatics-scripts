# Bioinformatics and Statistical Analysis Scripts

A collection of reusable scripts developed for bioinformatics, statistical analysis, and microbiome research. This repository brings together practical workflows for processing sequencing data, analysing microbial communities, performing statistical modelling, and visualising research results.

The scripts are intended to support reproducible research and provide adaptable starting points for biomedical and biological data analysis.

## Repository Structure

The repository is organised into four main areas: bioinformatics analysis, microbiome analysis, statistical analysis, and Shiny applications.

```text
bioinformatics-scripts/
│
├── README.md
│
├── bioinformatics-analysis/
│   ├── amplicon-sequencing/
│   │   ├── 03.phyloseq_run1_4_H.R
│   │   └── dada2_amplicon_pipeline.R
│   │
│   └── dpipeline_shotgun_BSC/
│       ├── 1. fastqc.slurm
│       ├── 2.multiqc.slurm
│       ├── 3.trimgalore.slurm
│       ├── 4.fastqc_trimmed.sh
│       ├── 5.multiqc_trimmed.sh
│       ├── 6.bowtie2.sh
│       ├── 7.read_stats.sh
│       └── 8.filter_samples.sh
│
├── microbiome-analysis/
│   ├── Aldex2_protocol_2.R
│   ├── Coda4Microbiome_pipeline.R
│   ├── Interactive ordination.R
│   ├── PCA Aitchison distances with permanova.R
│   ├── PCOA unweighted unifrac with permanova.R
│   ├── PICRUST2 pipeline.R
│   ├── ancombc.R
│   ├── SparCC_correlations.R
│   ├── unsupervisedPCA.R
│   └── venndiagram.R
│
├── statistical-analysis/
│   ├── ML_prediction/
│   │   ├── CatboostPipeline.R
│   │   ├── RandomforestPipeline.R
│   │   ├── SeleccionPredictores.R
│   │   └── XGboostPipeline.R
│   │
│   ├── Anova_Tukey_Grafics.R
│   ├── Descriptive_tables.R
│   ├── Kaplan_meier.R
│   ├── Linear mixed-effects model.R
│   ├── Univariate_Multivariate_Analyses.R
│   ├── data_distribution_assessment.R
│   ├── sample_size_calculator.R
│   └── survival_analysis_microbiome.R
│
└── shiny-apps/
    └── data-distribution-assessment/
        ├── app.R
        └── README.md
```

*Note: This structure reflects the current organisation described in this README. Update the tree whenever files are added, renamed, or moved. The presence of a script does not necessarily mean that it has been tested or validated for general use.*

## 1. Bioinformatics Analysis

**Directory:** `bioinformatics-analysis/`

Scripts and workflows for processing biological sequencing data, including quality control, read filtering, sequence processing, and preparation of sequencing datasets.

### Amplicon Sequencing

**Directory:** `bioinformatics-analysis/amplicon-sequencing/`

Workflows for analysing amplicon sequencing data, including:

- Read-quality assessment and filtering.
- Sequence denoising and amplicon sequence variant (ASV) inference using DADA2.
- Paired-end read merging and chimera removal.
- Taxonomic assignment.
- Construction and processing of phyloseq objects.

The exact parameters, reference databases, and processing steps should be adapted to the sequencing platform, target region, primers, and experimental design.

### Shotgun Metagenomics

**Directory:** `bioinformatics-analysis/dpipeline_shotgun_BSC/`

Shell scripts and SLURM job scripts for shotgun metagenomic sequencing workflows.

The available scripts cover steps such as:

- Raw-read quality control using FastQC.
- Quality-report aggregation using MultiQC.
- Read trimming using Trim Galore.
- Quality assessment after trimming.
- Host-read filtering using Bowtie2.
- Read-count summaries and sample filtering.

These workflows may require access to a high-performance computing (HPC) environment. Check the scripts for their software dependencies, reference databases, input-file conventions, and resource requirements before execution.

## 2. Microbiome Analysis

**Directory:** `microbiome-analysis/`

Scripts for statistical analysis and visualisation of microbiome datasets, focusing on microbial community structure, differential abundance, associations, and compositional data analysis.

### Community Structure and Ordination

Available scripts include workflows involving:

- Principal component analysis (PCA).
- Principal coordinates analysis (PCoA).
- Aitchison distances and centred log-ratio (CLR) transformations.
- Unweighted UniFrac distances.
- Interactive ordination visualisations.
- Permutational multivariate analysis of variance (PERMANOVA).

Distance metrics and statistical methods should be selected according to the data characteristics and research question. Where samples are repeated within individuals, the permutation scheme or statistical model should account for the study design when appropriate.

### Differential Abundance

Scripts include approaches based on:

- ANCOM-BC.
- ALDEx2.
- Other methods implemented in the available analysis workflows.

Differential-abundance methods make different assumptions and may produce different results. Review their documentation and assess whether each method is appropriate for the dataset, sampling design, and research question.

### Compositional Data Analysis and Associations

Available workflows include:

- CoDA4Microbiome.
- SparCC-based microbial association analysis.
- Compositional transformations and multivariate analyses.

Interpret associations carefully, particularly when analysing relative-abundance data. Associations between microbial taxa do not necessarily indicate direct biological interactions or causal relationships.

### Functional Profiling

**Script:** `PICRUST2 pipeline.R`

A workflow for predicting potential microbial functional profiles from marker-gene sequencing data using PICRUSt2.

Predicted functional profiles are computational inferences, not direct measurements of microbial genes, gene expression, or metabolic activity.

## 3. Statistical Analysis

**Directory:** `statistical-analysis/`

Reusable scripts for descriptive analysis, hypothesis testing, regression, survival analysis, mixed-effects modelling, sample-size calculations, and machine learning.

### Descriptive and Exploratory Analysis

Relevant scripts:

- `Descriptive_tables.R`
- `data_distribution_assessment.R`

These scripts support exploratory data analysis, descriptive summaries, and assessment of data distributions.

Depending on the implementation, analyses may include missing-data summaries, descriptive statistics, distribution plots, and normality assessments.

### Hypothesis Testing

**Script:** `Anova_Tukey_Grafics.R`

A workflow for analysis of variance (ANOVA), post-hoc comparisons using Tukey's method, and graphical representation of results.

The suitability of ANOVA depends on the study design and the assumptions of the fitted model.

### Regression Analysis

**Script:** `Univariate_Multivariate_Analyses.R`

Scripts for univariable and multivariable analyses, including regression-based approaches.

The specific model family should be selected according to the outcome variable and study design. Multivariable models can be used to adjust for measured covariates, although adjustment does not by itself eliminate confounding.

### Mixed-Effects Models

**Script:** `Linear mixed-effects model.R`

Tools for fitting linear mixed-effects models when observations are clustered or repeated within individuals, centres, or other grouping units.

These models can help account for dependence between observations and quantify variation at different levels of the study design.

### Survival Analysis

Relevant scripts:

- `Kaplan_meier.R`
- `survival_analysis_microbiome.R`

Workflows for time-to-event analysis, potentially including:

- Kaplan–Meier survival curves.
- Log-rank tests.
- Cox proportional hazards regression.
- Associations between microbiome features and clinical outcomes.

The event definition, time origin, censoring rules, and model assumptions should be checked for each dataset. Cox models require assessment of the proportional-hazards assumption.

### Sample-Size Calculation

**Script:** `sample_size_calculator.R`

A utility for sample-size calculations. The appropriate calculation depends on the primary outcome, target effect size, statistical power, significance level, study design, and expected data loss.

Check the script's documentation to identify the specific calculations it supports.

### Machine Learning and Prediction

**Directory:** `statistical-analysis/ML_prediction/`

Scripts for predictive modelling and predictor selection, including:

- CatBoost.
- Random Forest.
- XGBoost.
- Predictor-selection workflows.

Machine-learning performance should be assessed using an appropriate validation strategy. When multiple observations come from the same individual, data splitting should prevent information leakage between training and test sets. Feature selection and hyperparameter tuning should also be performed within the training or resampling procedure.

## 4. Shiny Applications

**Directory:** `shiny-apps/`

Interactive applications developed with R and Shiny to make selected analyses accessible through a graphical interface.

### Data Distribution Assessment

**Directory:** `shiny-apps/data-distribution-assessment/`

An application for exploring uploaded datasets and assessing numerical variables.

Depending on the implemented functionality, it may provide:

- Descriptive statistics.
- Missing-data summaries.
- Histograms and boxplots.
- Normality assessment.
- Quantile–quantile (Q–Q) plots.

Consult the application-specific README for its dependencies, input requirements, and usage instructions.

## Getting Started

### Requirements

Requirements vary by script. Most R workflows require R and the relevant CRAN or Bioconductor packages. Bioinformatics workflows may additionally require command-line tools, reference databases, and access to an HPC environment.

Examples of R packages used across different workflows include:

- `ggplot2`
- `dplyr`
- `vegan`
- `phyloseq`
- `microbiome`
- `survival`
- `survminer`
- `shiny`

These packages are examples, not a universal dependency list. Install only the packages required by the workflow you intend to run.

### Running an R Script

1. Clone or download the repository.
2. Read the relevant script and its accompanying documentation.
3. Install the required dependencies.
4. Prepare the input data in the expected format.
5. Update input and output paths and any study-specific parameters.
6. Run the script in R or RStudio.
7. Check warnings, intermediate outputs, and diagnostic plots before interpreting the results.

Where possible, use project-relative paths and record software and package versions to improve reproducibility.

### Running a Shiny Application

From the application's directory, open R and run:

```r
install.packages("shiny")  # Only if Shiny is not already installed
shiny::runApp()
```

Install any additional packages required by the application. The application runs locally unless it is deployed to a hosting service.

### Running Shell or SLURM Workflows

Before running a shell script or submitting a SLURM job:

1. Check the required software and environment modules.
2. Review input and output paths.
3. Confirm reference-database locations.
4. Verify resource requests, including CPUs, memory, and runtime.
5. Test the workflow on a small subset of data.

Do not assume that scripts written for one HPC environment will run unchanged on another.

## Reproducibility and Data Management

The following practices are recommended:

- Document each script's purpose, input format, dependencies, parameters, and expected outputs.
- Keep example datasets synthetic or appropriately de-identified.
- Avoid hard-coded personal or institution-specific paths where possible.
- Record relevant software and package versions.
- Document statistical assumptions and limitations.
- Validate workflows before applying them to new datasets.
- Use version control to track changes to scripts and documentation.

**Do not commit identifiable patient data, confidential research data, credentials, passwords, or other sensitive information to a public repository.**

## Contributing

Suggestions, bug reports, and improvements are welcome.

When adding or modifying a script:

1. Place it in the most appropriate directory.
2. Use a clear, descriptive filename.
3. Document its purpose, inputs, outputs, dependencies, and parameters.
4. Include a small example dataset when appropriate and safe to share.
5. Explain important assumptions and limitations.
6. Test the script on suitable example data before describing it as validated.

## Disclaimer

These scripts are provided for research and educational purposes. They are not guaranteed to be suitable for every dataset, experimental design, or research question.

Users are responsible for checking data quality, verifying assumptions, validating computational workflows, and interpreting results appropriately before drawing scientific conclusions.

## License

No license is specified in this README. Unless a license is added to the repository, reuse and redistribution may be subject to copyright restrictions. Add an appropriate `LICENSE` file if you intend to permit others to reuse or modify the scripts.

## Author

**Aleix Lluansí**

