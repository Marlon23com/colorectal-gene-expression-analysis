# Load required Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("DESeq2", "EnhancedVolcano", "RColorBrewer", "pheatmap", 
                       "clusterProfiler", "org.Hs.eg.db", "edgeR", 
                       "SummarizedExperiment", "limma"))

BiocManager::install('org.Hs.eg.db')
BiocManager::install('GO.db', force = TRUE)

# Load libraries
library(DESeq2)
library(EnhancedVolcano)
library(RColorBrewer)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(edgeR)
library(SummarizedExperiment)
library(limma)

# Set working directory
setwd("C:/Users/hp/Desktop/2nd Kota/Gene Structure & Function/Project/Dataset")

# 1. Data Loading & Preprocessing
#-----------------------------------------------------
# Load experimental design
design <- read.delim('E-GEOD-50760-experiment-design.tsv', sep = '\t')
design <- design[1:37, ]
design$label <- c(rep('cancer', 18), rep('normal', 19))  # Add sample labels

# Load raw counts and subset to match design
raw.counts <- read.delim('E-GEOD-50760-raw-counts.tsv')

raw.counts <- raw.counts[, c("Gene.Name", design$Run)]

# View data
View(raw.counts)

# Prepare metadata
metadata <- data.frame(id = design$Run, type = design$label)
rownames(metadata) <- metadata$id
View(metadata)

# Check if data aligns
cat("Count data samples:", ncol(raw.counts) - 1, "\n")  # Subtract 1 for Gene.Name column
cat("Metadata samples:  ", nrow(metadata), "\n")

# 2. Quality Control & Filtering
#-----------------------------------------------------
# Remove Gene.Name column and convert to matrix
# Store gene names before removing the column
gene_names <- raw.counts$Gene.Name

# Create count matrix (now you'll have matching gene names)
count_data <- as.matrix(raw.counts[, -which(colnames(raw.counts) == "Gene.Name")])
rownames(count_data) <- gene_names  # Important: attach gene names as row names

# Filter low-expressed genes (at least 5 samples with >10 counts)
keep_genes <- rowSums(count_data > 10) >= 5  
filtered_counts <- count_data[keep_genes, ]
View(filtered_counts)

# Ensure column names in count_data match row names in metadata
stopifnot(ncol(filtered_counts) == nrow(metadata))
stopifnot(all(colnames(filtered_counts) %in% rownames(metadata)))

# 3. DESeq2 Analysis
#-----------------------------------------------------

# 2. Convert your metadata columns to factors (important!)
metadata$type <- as.factor(metadata$type)


# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = filtered_counts,
  colData = metadata,
  design = ~ type
)

# Normalize and run DESeq2
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)

# Get results
res_bh <- results(dds, contrast = c("type", "cancer", "normal"), pAdjustMethod = "BH")
res_holm <- results(dds, contrast = c("type", "cancer", "normal"), pAdjustMethod = "holm")
view(res_bh)

# Get results with gene name for plotting

# Convert to data frame and add gene names
res_bh_df <- as.data.frame(res_bh)
res_bh_df$gene <- rownames(count_data)[match(rownames(res_bh_df), rownames(dds))]

# view results
head(res_bh_df, 2)  # Should show log2FC, pvalue columns
rownames(res_bh_df)  # Should return gene IDs

# Check if rownames exist
if(is.null(rownames(res_bh_df))) {
  # If no rownames, create a dummy ID column
  res_df$ID <- paste0("Gene_", 1:nrow(res_bh_df))
  lab_col <- "ID"
} else {
  lab_col <- rownames(res_bh_df)
}

view(lab_col)


##### ====================================================
# Apply significance thresholds
sig_genes <- res_bh[
  which(res_bh$padj <= 0.01 & abs(res_bh$log2FoldChange) >= 1 & res_bh$baseMean >= 20),
]
sig_genes <- sig_genes[order(sig_genes$padj), ]

# 4. Visualization
#-----------------------------------------------------
# PCA plot
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "type")

EnhancedVolcano(
  res_bh_df,
  lab = res_bh_df$gene,  # Use explicit gene column
  x = 'log2FoldChange',
  y = 'pvalue',
  pCutoff = 1e-5,
  FCcutoff = 1.333
)


# Heatmap
top_genes <- rownames(sig_genes)[1:min(50, nrow(sig_genes))]
pheatmap(assay(vsd)[top_genes, ], annotation_col = metadata)

# MA plot from DESeq2 results
plotMA(res_bh, main="MA Plot: Cancer vs Normal", ylim=c(-5,5))
abline(h=c(-1,1), col="red", lty=2)  # Add log2FC cutoff lines

library(ggplot2)

# Create MA plot data
ma_data <- data.frame(
  mean = res_bh$baseMean,
  lfc = res_bh$log2FoldChange,
  significant = res_bh$padj < 0.01 & abs(res_bh$log2FoldChange) > 1
)

ggplot(ma_data, aes(x=log10(mean), y=lfc)) +
  geom_point(aes(color=significant), alpha=0.6, size=1) +
  scale_color_manual(values=c("grey", "red")) +
  geom_hline(yintercept=c(-1,1), linetype="dashed", color="blue") +
  labs(x="Mean of normalized counts (log10)", 
       y="Log2 fold change",
       title="MA Plot") +
  theme_minimal()

# 5. Functional Enrichment
#-----------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)

# Prepare significant genes (using adjusted p-value < 0.01 and |log2FC| > 1)
sig_genes <- subset(res_bh_df, padj < 0.01 & abs(log2FoldChange) > 1)
gene_list <- sig_genes$gene  # Your gene symbols

# Convert gene symbols to ENTREZ IDs
entrez_ids <- bitr(gene_list, 
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
view(entrez_ids)
# ========veriFYING GENES===========
# Check how many genes were successfully converted to ENTREZ IDs
print(paste("Input genes:", length(gene_list)))
print(paste("Mapped genes:", nrow(entrez_ids)))

# Verify the first few conversions
head(entrez_ids, 3)

print(paste("Number of significant genes:", nrow(sig_genes)))
# =========looks good!==========


# ==============verify database contents================
# Check gene list
print(head(gene_list, 10))

# Verify database contents
columns(org.Hs.eg.db)  # Available ID types
keytypes(org.Hs.eg.db) # Mappable types

# Test with known human genes
test_genes <- c("TP53", "BRCA1", "ACTB")
bitr(test_genes, "SYMBOL", "ENTREZID", org.Hs.eg.db)
# ==========looks GOOD====================

# GO Enrichment Analysis
ego <- enrichGO(gene ,
                OrgDb = org.Hs.eg.db,
                ont = "BP",  # Biological Process
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2)

# View top enriched terms
head(ego, 3)
view(ego)

# Plot GO enrichment results
dotplot(ego, showCategory=15, title="GO Biological Process Enrichment")

# KEGG Pathway Analysis (if needed)
kk <- enrichKEGG(gene = entrez_ids$ENTREZID,
                 organism = 'hsa',
                 pvalueCutoff = 0.05)
barplot(kk, showCategory=15, title="KEGG Pathway Enrichment")

# 6. Export Results
#-----------------------------------------------------
write.csv(as.data.frame(sig_genes), "significant_DEGs.csv")
write.csv(as.data.frame(res_bh), "all_DEGs_bh.csv")
write.csv(as.data.frame(res_holm), "all_DEGs_holm.csv")