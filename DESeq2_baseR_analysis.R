# DESeq2 Differential Expression Analysis (Base R Version)
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

library(DESeq2)
library(org.Hs.eg.db)
library(reactome.db)
library(gplots)
library(genefilter)
library(RColorBrewer)

# 2. Data Import --------------------------------------------------------------
setwd("C:/Users/hp/Desktop/2nd Kota/Gene Structure & Function/Project/Dataset")
design <- read.delim('E-GEOD-50760-experiment-design.tsv', sep = '\t')
design <- design[1:37, ]
raw.counts <- read.delim('E-GEOD-50760-raw-counts.tsv')

# 3. Data Preprocessing -------------------------------------------------------
# Subset columns to match design
raw.counts <- raw.counts[, c("Gene.Name", design$Run)]
rownames(raw.counts) <- make.unique(as.character(raw.counts$Gene.Name))
count_data <- as.matrix(raw.counts[, -1])
rownames(count_data) <- rownames(raw.counts)

# Filter low-expressed genes
keep_genes <- rowSums(count_data > 10) >= 5
filtered_counts <- count_data[keep_genes, ]

# Prepare metadata
metadata <- data.frame(id = design$Run, type = c(rep('cancer', 18), rep('normal', 19)))
rownames(metadata) <- metadata$id
metadata$type <- as.factor(metadata$type)

# 4. Exploratory Data Analysis (EDA) ------------------------------------------
# Descriptive statistics for each sample (raw)
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
res_bh <- results(dds, contrast = c("type", "cancer", "normal"), pAdjustMethod = "BH")
res_holm <- results(dds, contrast = c("type", "cancer", "normal"), pAdjustMethod = "holm")

# Convert to data frame and add gene names
res_bh_df <- as.data.frame(res_bh)
res_bh_df$gene <- rownames(count_data)[match(rownames(res_bh_df), rownames(dds))]

# Print head of results for quick inspection
print(head(res_bh_df, 10))

# Normalized counts
norm_counts <- counts(dds, normalized = TRUE)

# Descriptive statistics for each sample (normalized)
norm_stats <- descr_stats(norm_counts)
write.csv(norm_stats, "sample_descriptive_stats_normalized.csv")
barplot(norm_stats$Mean, main = "Sample Means (Normalized)", col = "lightgreen", las = 2, ylab = "Mean Norm. Counts", names.arg = colnames(norm_counts), cex.names = 0.7)
barplot(norm_stats$SD, main = "Sample SDs (Normalized)", col = "purple", las = 2, ylab = "SD Norm. Counts", names.arg = colnames(norm_counts), cex.names = 0.7)

# Combined boxplot of raw and normalized counts (side by side for comparison)
par(mar = c(8, 4, 4, 2) + 0.1) # Increase bottom margin for sample names
boxplot(log2(filtered_counts + 1),
        las = 2, col = "lightblue", main = "Boxplot: Raw Counts (log2)", ylab = "log2(count + 1)",
        names = colnames(filtered_counts), cex.axis = 0.7)

boxplot(log2(norm_counts + 1),
        las = 2, col = "lightgreen", main = "Boxplot: Normalized Counts (log2)", ylab = "log2(count + 1)",
        names = colnames(norm_counts), cex.axis = 0.7)
par(mar = c(5, 4, 4, 2) + 0.1) # Reset margins

# Violin plots for top 20 highly expressed genes (raw and normalized)
if (!requireNamespace("vioplot", quietly = TRUE)) install.packages("vioplot")
library(vioplot)

# Raw counts violin plot
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

# Normalized counts violin plot
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

# 6. QC & Diagnostic Plots ----------------------------------------------------
plotMA(res_bh, ylim = c(-3, 3), main = "MA Plot")
plotDispEsts(dds, ylim = c(1e-6, 1e1), main = "Dispersion Plot")
hist(res_bh$pvalue, breaks = 20, col = "grey", main = "P-value Histogram")

# Volcano Plot (EnhancedVolcano)
if (requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  library(EnhancedVolcano)
  EnhancedVolcano(res_bh_df,
                  lab = res_bh_df$gene,
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

# Volcano Plot (Base R)
volcano_data <- data.frame(log2FC = res_bh$log2FoldChange, negLogP = -log10(res_bh$pvalue),
                           padj = res_bh$padj, gene = rownames(res_bh))
with(volcano_data, plot(log2FC, negLogP, pch = 20, col = ifelse(padj < 0.05 & abs(log2FC) > 1, "red", "grey"),
                        main = "Volcano Plot", xlab = "log2 Fold Change", ylab = "-log10(p-value)"))
abline(v = c(-1, 1), col = "blue", lty = 2)
abline(h = -log10(0.05), col = "green", lty = 2)
legend("topright", legend = c("Significant", "Not significant"), col = c("red", "grey"), pch = 20)

# rlog variance plot
rld <- rlog(dds)
meanSdPlot <- function(mat, main = "Mean-SD Plot") {
  means <- rowMeans(mat)
  sds <- apply(mat, 1, sd)
  plot(means, sds, pch = 20, col = rgb(0,0,1,0.3), main = main, xlab = "Mean", ylab = "SD")
  abline(lm(sds ~ means), col = "red")
}
par(mfrow = c(1, 2))
meanSdPlot(log2(filtered_counts + 1), main = "Mean-SD Plot (Raw)")
meanSdPlot(assay(rld), main = "Mean-SD Plot (rlog)")
par(mfrow = c(1, 1))

# 7. Sample & Gene Clustering -------------------------------------------------
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
colours = colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
heatmap.2(sampleDistMatrix, trace = "none", col = colours, main = "Sample Distance Heatmap")

plotPCA(rld, intgroup = "type")
title("PCA Plot")

topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 35)
heatmap.2(assay(rld)[topVarGenes, ], scale = "row", trace = "none", dendrogram = "both", 
          Colv = TRUE, Rowv = TRUE, col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255), 
          main = "Top Variable Genes Heatmap (Gene & Sample Clustering)",
          labCol = colnames(assay(rld)), labRow = rownames(assay(rld))[topVarGenes],
          margins = c(8, 8))

# 8. Reactome Pathway Enrichment (Base R) -------------------------------------
res2 <- res_bh[res_bh$entrezid %in% keys(reactome.db, "ENTREZID") & !is.na(res_bh$padj), ]
reactomeTable <- AnnotationDbi::select(reactome.db, keys = as.character(res2$entrezid), keytype = "ENTREZID", columns = c("ENTREZID", "REACTOMEID"))
incm <- do.call(rbind, with(reactomeTable, tapply(ENTREZID, factor(REACTOMEID), function(x) res2$entrezid %in% x)))
colnames(incm) <- res2$entrezid
within <- function(x, lower, upper) (x >= lower & x <= upper)
incm <- incm[within(rowSums(incm), lower = 20, upper = 80), ]

testCategory <- function(reactomeID) {
  isMember <- incm[reactomeID, ]
  data.frame(
    reactomeID = reactomeID,
    numGenes = sum(isMember),
    avgLFC = mean(res2$log2FoldChange[isMember]),
    sdLFC = sd(res2$log2FoldChange[isMember]),
    pvalue = t.test(res2$log2FoldChange[isMember])$p.value,
    stringsAsFactors = FALSE
  )
}
reactomeResult <- do.call(rbind, lapply(rownames(incm), testCategory))
reactomeResult$padjust <- p.adjust(reactomeResult$pvalue, "BH")

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
res_bh$hgnc_symbol <- convertIDs(rownames(res_bh), "ENSEMBL", "SYMBOL", org.Hs.eg.db)
res_bh$entrezid <- convertIDs(rownames(res_bh), "ENSEMBL", "ENTREZID", org.Hs.eg.db)

# 10. Export Results ----------------------------------------------------------
write.csv(as.data.frame(res_bh), "all_DEGs.csv")
write.csv(as.data.frame(reactomeResult), "reactome_enrichment.csv")

# End of Script
