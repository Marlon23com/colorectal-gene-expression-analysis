# Colorectal Cancer Differential Gene Expression Analysis

This project focuses on differential expression analysis using RNA-seq count data for cancer versus normal tissue samples. The workflow uses DESeq2 and R to analyze gene expression changes and identify statistically significant differentially expressed genes (DEGs).

## Project overview

The analysis compares:
- cancer samples
- normal tissue samples

The workflow includes data preprocessing, normalization, DESeq2 modeling, visualization, and DEG identification.

## Included files
- `DESeq2_baseR_analysis.R` — base R implementation
- `DESeq2_tidyverse_analysis.R` — tidyverse-based workflow
- `DESeq2_analysis_explained.Rmd` — annotated explanatory analysis file
- `baseRworkingscript.R` — working script with additional exploration
- `sample_descriptive_stats_raw.csv` — descriptive stats for raw data
- `sample_descriptive_stats_normalized.csv` — normalized sample metrics
- `all_DEGs.csv` — filtered differentially expressed genes

## Tools used
- R
- DESeq2
- Bioconductor
- ggplot2 / tidyverse

## How to run
1. Open RStudio.
2. Install required packages.
3. Run the analysis script for the chosen workflow.
4. Review output tables and plots.

## Objective
This project demonstrates gene expression analysis, RNA-seq preprocessing, differential expression testing, and biological interpretation of significant genes in a cancer-related dataset.
