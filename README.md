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

## Core principles used

This project follows the standard RNA-seq differential expression workflow used in bioinformatics:

- count matrix normalization
- filtering lowly expressed genes
- metadata design for condition labels
- DESeq2 modeling for gene expression comparison
- multiple testing correction
- visualization of distribution and differential expression

The main statistical principle in this analysis is that RNA-seq counts are not directly comparable across samples because library sizes and technical variation differ. DESeq2 models count data with a negative binomial distribution and estimates variance before comparing cancer and normal groups. This is appropriate for sequencing count data because the counts are discrete, non-normal, and often vary in sequencing depth across samples.

The analytical thought process behind the workflow was to first ensure the data were suitable for comparison. I filtered low-count genes to remove noise and reduce false positives, then built a design matrix with cancer and normal labels so the model could test the effect of disease status on gene expression. After fitting the DESeq2 model, the analysis focused on identifying genes with consistent differences in expression while accounting for biological and technical variability.

The multiple testing correction is essential because thousands of genes are tested simultaneously. If we simply used ordinary p-values without adjustment, a large number of false positives would appear. The script therefore compares results using two correction methods: Benjamini-Hochberg (BH) and Holm.

### Benjamini-Hochberg (BH) correction
The BH procedure controls the false discovery rate (FDR), meaning it limits the proportion of false positives among the genes called significant. In practical terms, it ranks p-values from smallest to largest, adjusts each value based on its rank, and identifies a set of genes that are significant while controlling overall false discoveries. This is commonly preferred in RNA-seq analyses because it balances sensitivity and false-positive control when thousands of genes are tested.

### Holm correction
The Holm procedure is a family-wise error rate (FWER) method. Unlike BH, which controls the expected proportion of false discoveries, Holm is more conservative and aims to control the probability of making at least one false positive across the whole gene list. It adjusts p-values in a stepwise way, giving stronger protection against any false positive but often reducing power compared with BH.

In this project, using both corrections demonstrates a more rigorous statistical approach: BH is useful when the priority is identifying likely biologically relevant genes while controlling FDR, whereas Holm provides a stricter safeguard against false positives across the full set of tests.

## Intended result

The intended result is to identify differentially expressed genes between cancer tissue and normal tissue. In practical terms, the output should produce a list of genes with:

- log2 fold changes
- p-values
- adjusted p-values
- biological significance

The project also aims to generate diagnostic and exploratory outputs, including:
- PCA
- MA plots
- volcano plots
- heatmaps
- sample distribution summaries

These outputs help assess whether the dataset behaves as expected and whether the identified genes are biologically meaningful rather than merely statistical artifacts.

## Included files
- `DESeq2_baseR_analysis.R` — base R differential expression workflow
- `DESeq2_tidyverse_analysis.R` — tidyverse-based DESeq2 analysis
- `DESeq2_analysis_explained.Rmd` — annotated explanation of the analysis pipeline
- `baseRworkingscript.R` — exploratory working script
- `all_DEGs.csv` — differentially expressed genes output
- `all_DEGs_tidy.csv` — tidy-formatted DEG results
- `sample_descriptive_stats_raw.csv` — raw sample summary statistics
- `sample_descriptive_stats_normalized.csv` — normalized sample summary statistics
