# Colorectal Cancer Differential Gene Expression Analysis

This project analyzes RNA-seq expression data from colorectal cancer and normal tissue samples to identify genes significantly altered between conditions. The workflow uses DESeq2 in R to perform differential expression analysis and generate biologically relevant gene lists and diagnostic summaries for interpretation.

## Project overview

The analysis follows a standard bioinformatics workflow, including:
- sample preparation and condition labeling
- normalization and exploratory data quality checks
- DESeq2 differential expression modeling
- identification of significant genes and effect sizes
- output tables and summary statistics for downstream interpretation

This project demonstrates core computational biology skills in RNA-seq data analysis, statistical testing, and interpretation of transcriptomic changes associated with cancer biology.

## Thought process and analytical rationale

The central challenge in RNA-seq analysis is that raw counts cannot be interpreted directly as expression levels because sequencing depth, library size, and technical variability differ between samples. A gene with a higher count in one sample may not necessarily be biologically more active if that sample simply received more sequencing reads. For this reason, the first major step in the workflow was normalization, which adjusts for technical differences so that count values are more comparable across samples.

I then filtered lowly expressed genes to remove weakly detected features that add noise and reduce statistical power. This helps ensure that the model is based on genes with reliable measurable expression. Once the expression matrix was cleaned, metadata was constructed to define the biological comparison: cancer versus normal samples.

DESeq2 was chosen because it is specifically designed for count-based transcriptomic data. Unlike standard linear models applied to continuous variables, DESeq2 models counts with a negative binomial distribution, which better reflects the variance structure of RNA-seq data. This allows the analysis to estimate biological variability and determine which genes show significant expression changes between groups.

The comparison between cancer and normal tissue is biologically meaningful because it allows us to ask which genes are differentially expressed in disease relative to healthy tissue. In a cancer context, this can reveal pathways, biomarkers, or regulatory programs that are activated or suppressed during tumor progression.

## Multiple testing correction: BH and Holm

A key statistical issue in RNA-seq analysis is the problem of multiple testing. In a gene-expression experiment, we are testing thousands of genes simultaneously. If we used an unadjusted significance threshold such as p < 0.05 for each gene, we would expect many false positives simply by chance.

This is why multiple-testing correction is essential. In this analysis, two correction procedures were considered: Benjamini-Hochberg (BH) and Holm.

### Benjamini-Hochberg (BH)
The BH procedure controls the False Discovery Rate (FDR), which is the expected proportion of false positives among the rejected hypotheses. It is widely used in transcriptomics because the goal is often to identify a manageable set of genes while allowing some controlled false discoveries.

In practice, the BH method ranks the p-values from smallest to largest, adjusts them according to their rank, and identifies genes that remain significant after this FDR control. This is particularly useful in RNA-seq because we expect large numbers of genes to be tested, and we want a balance between sensitivity and specificity.

### Holm correction
The Holm correction is a family-wise error rate (FWER) method, which is generally more conservative than BH. Instead of controlling the proportion of false positives in the set of discoveries, it controls the probability of making at least one false rejection across the entire family of hypotheses.

This means Holm is stricter and usually reduces the number of significant findings compared with BH. It is useful when the aim is to minimize the chance of even one false positive, but it can reduce sensitivity when many tests are being performed.

In this project, comparing BH and Holm was valuable because it showed how the interpretation of the results changes depending on the statistical criterion chosen. BH is often preferred in expression studies when the focus is on discovering candidate genes, while Holm is more conservative and suitable when strict control over false positives is required.

## Intended result

The intended result is to identify differentially expressed genes between cancer tissue and normal tissue. In practical terms, the output should produce a list of genes with:
- log2 fold changes
- p-values
- adjusted p-values
- biological significance

The project also aims to generate diagnostic and exploratory outputs such as:
- PCA
- MA plots
- volcano plots
- heatmaps
- sample distribution summaries

These analyses help determine whether the data are behaving as expected and whether the identified genes are biologically meaningful. PCA and clustering visualizations assess sample separation by condition, while MA and volcano plots help summarize effect size and significance. Together, these outputs provide both statistical evidence and biological context for the expression differences observed.

## Included files
- `DESeq2_baseR_analysis.R` — base R differential expression workflow
- `DESeq2_tidyverse_analysis.R` — tidyverse-based DESeq2 analysis
- `DESeq2_analysis_explained.Rmd` — annotated explanation of the analysis pipeline
- `baseRworkingscript.R` — exploratory working script
- `all_DEGs.csv` — differentially expressed genes output
- `all_DEGs_tidy.csv` — tidy-formatted DEG results
- `sample_descriptive_stats_raw.csv` — raw sample summary statistics
- `sample_descriptive_stats_normalized.csv` — normalized sample summary statistics
