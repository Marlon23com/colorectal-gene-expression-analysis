# Colorectal Cancer Differential Gene Expression Analysis

This project analyzes RNA-seq count data comparing colorectal cancer and normal tissue samples. The workflow uses DESeq2 in R to detect differentially expressed genes and summarize the key biological findings from the dataset.

## Project overview

The analysis compares cancer and normal samples using a standard differential expression pipeline, including:
- data preparation and sample labeling
- normalization and quality checks
- DESeq2 modeling for differential expression
- identification of significant genes
- output tables and exploratory visualizations

This project demonstrates practical bioinformatics skills in RNA-seq analysis and interpretation of gene expression changes associated with cancer.

## Included files
- `DESeq2_baseR_analysis.R` — base R differential expression workflow
- `DESeq2_tidyverse_analysis.R` — tidyverse-based DESeq2 analysis
- `DESeq2_analysis_explained.Rmd` — annotated explanation of the analysis
- `baseRworkingscript.R` — exploratory working script
- `all_DEGs.csv` — differentially expressed genes output
- `all_DEGs_tidy.csv` — tidy-formatted DEG output
- `sample_descriptive_stats_raw.csv` — raw sample summary statistics
- `sample_descriptive_stats_normalized.csv` — normalized sample summary statistics
