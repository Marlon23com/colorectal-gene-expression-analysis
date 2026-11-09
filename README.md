# Colorectal Cancer Differential Gene Expression Analysis

This project analyzes RNA-seq expression data from colorectal cancer and normal tissue samples to identify genes that are significantly altered between conditions. The workflow uses DESeq2 in R to perform differential expression analysis and generate biologically relevant gene lists and summary outputs.

## Project overview

The analysis follows a standard bioinformatics workflow, including:
- sample preparation and condition labeling
- normalization and exploratory data quality checks
- DESeq2 differential expression modeling
- identification of significant genes and effect sizes
- output tables and summary statistics for downstream interpretation

This project demonstrates core computational biology skills in RNA-seq data analysis, statistical testing, and interpretation of transcriptomic changes associated with cancer biology.

## Included files
- `DESeq2_baseR_analysis.R` — base R differential expression workflow
- `DESeq2_tidyverse_analysis.R` — tidyverse-based DESeq2 analysis
- `DESeq2_analysis_explained.Rmd` — annotated explanation of the analysis pipeline
- `baseRworkingscript.R` — exploratory working script
- `all_DEGs.csv` — differentially expressed genes output
- `all_DEGs_tidy.csv` — tidy-formatted DEG results
- `sample_descriptive_stats_raw.csv` — raw sample summary statistics
- `sample_descriptive_stats_normalized.csv` — normalized sample summary statistics
