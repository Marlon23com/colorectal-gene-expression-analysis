# DESeq2 Differential Expression Analysis (tidyverse Version)
# Author: [Your Name]
# Date: [Today's Date]

# 1. Load Required Packages ---------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("reactome.db", quietly = TRUE)) BiocManager::install("reactome.db")
if (!requireNamespace("gplots", quietly = TRUE)) install.packages("gplots")
if (!requireNamespace("genefilter", quietly = TRUE)) BiocManager::install("genefilter")
if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")

library(DESeq2)
library(org.Hs.eg.db)
library(reactome.db)
library(gplots)
library(genefilter)
library(RColorBrewer)
library(tidyverse)

# 2. Data Import --------------------------------------------------------------
setwd("C:/Users/hp/Desktop/2nd Kota/Gene Structure & Function/Project/Dataset")
design <- read_tsv('E-GEOD-50760-experiment-design.tsv') %>% slice(1:37)
raw.counts <- read_tsv('E-GEOD-50760-raw-counts.tsv')

# 3. Data Preprocessing -------------------------------------------------------
raw.counts <- raw.counts %>% select(`Gene Name`, all_of(design$Run))
# Remove rows with missing 'Gene Name' before making rownames
raw.counts <- raw.counts %>% filter(!is.na(`Gene Name`) & `Gene Name` != "")
raw.counts <- raw.counts %>% mutate(`Gene Name` = make.unique(as.character(`Gene Name`)))
count_data <- raw.counts %>% column_to_rownames("Gene Name") %>% as.matrix()

# Filter low-expressed genes
keep_genes <- rowSums(count_data > 10) >= 5
filtered_counts <- count_data[keep_genes, ]

design <- design %>% mutate(type = c(rep('cancer', 18), rep('normal', 19)))
metadata <- design %>% select(Run, type) %>% column_to_rownames("Run")
metadata$type <- as.factor(metadata$type)

# 4. Exploratory Data Analysis (EDA) ------------------------------------------
boxplot(log2(filtered_counts + 1), las = 2, col = "lightblue", main = "Boxplot: Raw Counts (log2)", ylab = "log2(count + 1)")
top_genes <- head(order(rowMeans(filtered_counts), decreasing = TRUE), 20)
matplot(t(filtered_counts[top_genes, ]), type = "l", lty = 1, col = rainbow(20), main = "Top 20 Highly Expressed Genes (Raw)", ylab = "Counts", xlab = "Sample")
legend("topright", legend = rownames(filtered_counts)[top_genes], col = rainbow(20), lty = 1, cex = 0.5)

descr_stats <- function(mat) {
  data.frame(
    Min = apply(mat, 2, min),
    Q1 = apply(mat, 2, quantile, 0.25),
    Median = apply(mat, 2, median),
    Mean = apply(mat, 2, mean),
    Q3 = apply(mat, 2, quantile, 0.75),
    Max = apply(mat, 2, max),
    SD = apply(mat, 2, sd)
  )
}
raw_stats <- descr_stats(filtered_counts)
write.csv(raw_stats, "sample_descriptive_stats_raw.csv")
barplot(raw_stats$Mean, main = "Sample Means (Raw)", col = "skyblue", las = 2, ylab = "Mean Counts", names.arg = colnames(filtered_counts), cex.names = 0.7)
barplot(raw_stats$SD, main = "Sample SDs (Raw)", col = "orange", las = 2, ylab = "SD Counts", names.arg = colnames(filtered_counts), cex.names = 0.7)

# 5. DESeq2 Analysis & Normalization ------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = filtered_counts, colData = metadata, design = ~ type)
dds <- dds[rowSums(counts(dds)) > 0, ]
dds <- DESeq(dds)

# Get results with different p-adjustment methods
desq2_res_bh <- results(dds, contrast = c("type", "cancer", "normal"), pAdjustMethod = "BH")
desq2_res_holm <- results(dds, contrast = c("type", "cancer", "normal"), pAdjustMethod = "holm")

desq2_res_bh_df <- as.data.frame(desq2_res_bh)
desq2_res_bh_df$gene <- rownames(count_data)[match(rownames(desq2_res_bh_df), rownames(dds))]

print(head(desq2_res_bh_df, 10))

norm_counts <- counts(dds, normalized = TRUE)
norm_stats <- descr_stats(norm_counts)
write.csv(norm_stats, "sample_descriptive_stats_normalized.csv")
barplot(norm_stats$Mean, main = "Sample Means (Normalized)", col = "lightgreen", las = 2, ylab = "Mean Norm. Counts", names.arg = colnames(norm_counts), cex.names = 0.7)
barplot(norm_stats$SD, main = "Sample SDs (Normalized)", col = "purple", las = 2, ylab = "SD Norm. Counts", names.arg = colnames(norm_counts), cex.names = 0.7)

# Combined boxplot of raw and normalized counts (side by side for comparison)
par(mar = c(8, 4, 4, 2) + 0.1)  # Increase bottom margin for sample names
boxplot(log2(filtered_counts + 1),
        las = 2, col = "lightblue", main = "Boxplot: Raw Counts (log2)", ylab = "log2(count + 1)",
        names = colnames(filtered_counts), cex.axis = 0.7)
boxplot(log2(norm_counts + 1),
        las = 2, col = "lightgreen", main = "Boxplot: Normalized Counts (log2)", ylab = "log2(count + 1)",
        names = colnames(norm_counts), cex.axis = 0.7)
par(mar = c(5, 4, 4, 2) + 0.1)  # Reset margins

boxplot(list(Raw = as.vector(log2(filtered_counts + 1)), Normalized = as.vector(log2(norm_counts + 1))),
        col = c("lightblue", "lightgreen"),
        main = "Boxplot: Raw vs Normalized Counts (log2)",
        ylab = "log2(count + 1)")

# Violin plots for top 20 highly expressed genes (raw and normalized)
if (!requireNamespace("vioplot", quietly = TRUE)) install.packages("vioplot")
library(vioplot)

par(mar = c(8, 4, 4, 2) + 0.1)
vioplot::vioplot(
  as.list(as.data.frame(filtered_counts[top_genes, ])),
  names = colnames(filtered_counts),
  col = "lightblue",
  main = "Violin Plot: Top 20 Highly Expressed Genes (Raw)",
  ylab = "Counts",
  las = 2,
  cex.axis = 0.7
)

vioplot::vioplot(
  as.list(as.data.frame(norm_counts[top_genes, ])),
  names = colnames(norm_counts),
  col = "lightgreen",
  main = "Violin Plot: Top 20 Highly Expressed Genes (Normalized)",
  ylab = "Norm. Counts",
  las = 2,
  cex.axis = 0.7
)
par(mar = c(5, 4, 4, 2) + 0.1)

# Volcano Plot (EnhancedVolcano)
if (requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  library(EnhancedVolcano)
  EnhancedVolcano(desq2_res_bh_df,
                  lab = desq2_res_bh_df$gene,
                  x = 'log2FoldChange',
                  y = 'pvalue',
                  pCutoff = 1e-5,
                  FCcutoff = 1.333,
                  title = 'Volcano Plot (EnhancedVolcano)',
                  subtitle = 'cancer vs normal',
                  caption = 'BH adjusted p-values')
} else {
  message('EnhancedVolcano package not installed. Skipping EnhancedVolcano plot.')
}

# 6. QC & Diagnostic Plots ----------------------------------------------------
# MA Plot using ggplot2 (tidyverse style)
res_tbl <- as_tibble(desq2_res_bh_df, rownames = "ensembl")
ggplot(res_tbl, aes(x = baseMean, y = log2FoldChange, color = padj < 0.1)) +
  geom_point(alpha = 0.5) +
  scale_x_log10() +
  labs(title = "MA Plot", x = "Mean Expression", y = "Log2 Fold Change") +
  theme_minimal()

plotDispEsts(dds, ylim = c(1e-6, 1e1), main = "Dispersion Plot")
ggplot(res_tbl, aes(x = pvalue)) +
  geom_histogram(bins = 20, fill = "grey", color = "black") +
  labs(title = "P-value Histogram", x = "P-value", y = "Count") +
  theme_minimal()

rld <- rlog(dds)
meanSdPlot <- function(mat, main = "Mean-SD Plot") {
  means <- rowMeans(mat)
  sds <- apply(mat, 1, sd)
  plot(means, sds, pch = 20, col = rgb(0,0,1,0.3), main = main, xlab = "Mean", ylab = "SD")
  abline(lm(sds ~ means), col = "red")
}
meanSdPlot(assay(rld), main = "Mean-SD Plot (rlog)")

# 7. Sample & Gene Clustering -------------------------------------------------
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
colours = colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
heatmap.2(sampleDistMatrix, trace = "none", col = colours, main = "Sample Distance Heatmap")

# PCA plot (tidyverse style, ggplot2)
pcaData <- plotPCA(rld, intgroup = "type", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = type, label = name)) +
  geom_point(size = 3) +
  geom_text(vjust = 1.5, size = 2.5) +
  labs(title = "PCA Plot", x = paste0("PC1 (", percentVar[1], "% variance)"),
       y = paste0("PC2 (", percentVar[2], "% variance)")) +
  theme_minimal()

topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 35)
heatmap.2(assay(rld)[topVarGenes, ], scale = "row", trace = "none", dendrogram = "column", col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255), main = "Top Variable Genes")

# 8. Reactome Pathway Enrichment (tidyverse) ----------------------------------
res2 <- res_tbl %>% filter(!is.na(padj), !is.na(entrezid), entrezid %in% keys(reactome.db, "ENTREZID"))
reactomeTable <- AnnotationDbi::select(reactome.db, keys = as.character(res2$entrezid), keytype = "ENTREZID", columns = c("ENTREZID", "REACTOMEID"))
incm <- do.call(rbind, with(reactomeTable, tapply(ENTREZID, factor(REACTOMEID), function(x) res2$entrezid %in% x)))
colnames(incm) <- res2$entrezid
within <- function(x, lower, upper) (x >= lower & x <= upper)
incm <- incm[within(rowSums(incm), lower = 20, upper = 80), ]

testCategory <- function(reactomeID) {
  isMember <- incm[reactomeID, ]
  tibble(
    reactomeID = reactomeID,
    numGenes = sum(isMember),
    avgLFC = mean(res2$log2FoldChange[isMember]),
    sdLFC = sd(res2$log2FoldChange[isMember]),
    pvalue = t.test(res2$log2FoldChange[isMember])$p.value
  )
}
reactomeResult <- map_dfr(rownames(incm), testCategory) %>%
  mutate(padjust = p.adjust(pvalue, "BH"))

# 9. Annotation ---------------------------------------------------------------
convertIDs <- function(ids, fromKey, toKey, db, ifMultiple = c("putNA", "useFirst")) {
  stopifnot(inherits(db, "AnnotationDb"))
  ifMultiple <- match.arg(ifMultiple)
  suppressWarnings(selRes <- AnnotationDbi::select(db, keys = ids, keytype = fromKey, columns = c(fromKey, toKey)))
  if (ifMultiple == "putNA") {
    duplicatedIds <- selRes[duplicated(selRes[, 1]), 1]
    selRes <- selRes[!selRes[, 1] %in% duplicatedIds, ]
  }
  return(selRes[match(ids, selRes[, 1]), 2])
}
res_tbl <- as_tibble(as.data.frame(res), rownames = "ensembl") %>%
  mutate(hgnc_symbol = convertIDs(ensembl, "ENSEMBL", "SYMBOL", org.Hs.eg.db),
         entrezid = convertIDs(ensembl, "ENSEMBL", "ENTREZID", org.Hs.eg.db))

# 10. Export Results ----------------------------------------------------------
write_csv(res_tbl, "all_DEGs_tidy.csv")
write_csv(reactomeResult, "reactome_enrichment_tidy.csv")

# 11. Gene Expression Profiling with geneProfileR ----------------------------------
# This section demonstrates how to use geneProfileR for gene expression visualization.
if (!requireNamespace("geneProfileR", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_github("bcm-uga/geneProfileR")
}
library(geneProfileR)

# Example: Plot expression profiles for selected genes (replace with your genes of interest)
selected_genes <- rownames(norm_counts)[1:5]  # Or use gene symbols if available
# Ensure metadata rownames match norm_counts column names
stopifnot(all(colnames(norm_counts) %in% rownames(metadata)))
geneProfilePlot(
  expressionMatrix = norm_counts[selected_genes, ],
  sampleAnnot = metadata,
  group = "type",
  main = "Gene Expression Profiles (geneProfileR)"
)

# End of Script
