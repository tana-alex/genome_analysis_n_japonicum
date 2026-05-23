#!/usr/bin/env Rscript
# 02_visualization.R

#==============================================================================
# PACKAGE LOADING
#==============================================================================

library("DESeq2")
library("ggplot2")
library("pheatmap")
library("ggrepel")

#==============================================================================
# LOAD DATA FROM PREVIOUS ANALYSIS
#==============================================================================

results_dir <- "results"

# Create plots directory
plots_dir <- paste(results_dir, "plots", sep="/")
if (!dir.exists(plots_dir)){
  dir.create(plots_dir)
}

# Load DESeq2 objects
load(paste(results_dir, "deseq2_objects.RData", sep="/"))

#==============================================================================
# PCA PLOTS
#==============================================================================

# Variance stabilizing transformation for plotting
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

# Basic PCA plot
pca_plot <- plotPCA(vsd, intgroup = "condition")
print(pca_plot)

# Enhanced PCA with labels
pca_enhanced <- pca_plot + 
  geom_text(aes(label = name), vjust = -1) +
  theme_bw() +
  ggtitle("PCA - Control vs Heat Treated") +
  scale_color_manual(values = c("Control" = "blue", "Heat_treated" = "red"))

print(pca_enhanced)

# Save PCA plot
ggsave(paste(plots_dir, "pca_plot.pdf", sep="/"), 
       plot = pca_enhanced, width = 8, height = 6)

#==============================================================================
# VOLCANO PLOT
#==============================================================================

# Convert results to dataframe
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# Add significance category
res_df$category <- "Not significant"
res_df$category[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Upregulated"
res_df$category[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Downregulated"

# Create volcano plot
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, 
                                   y = -log10(padj),
                                   color = category,
                                   label = gene)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_text(data = subset(res_df, padj < 0.05),
            vjust = -1, size = 3) +
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "Not significant" = "grey")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  theme_bw() +
  labs(title = "Heat Treated vs Control",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value")

print(volcano_plot)

# Save volcano plot
ggsave(paste(plots_dir, "volcano_plot.pdf", sep="/"), 
       plot = volcano_plot, width = 10, height = 8)

#==============================================================================
# HEATMAPS
#==============================================================================

# Heatmap of top 20 most significant genes
top_genes <- head(order(res$padj), 20)

heatmap_top20 <- pheatmap(assay(vsd)[top_genes, ],
                          annotation_col = sample_info,
                          scale = "row",
                          show_rownames = TRUE,
                          cluster_cols = FALSE,
                          main = "Top 20 Most Significant Genes")

# Save heatmap
pdf(paste(plots_dir, "heatmap_top20.pdf", sep="/"))
print(heatmap_top20)
dev.off()

#==============================================================================
# MA PLOT
#==============================================================================

# MA plot
pdf(paste(plots_dir, "ma_plot.pdf", sep="/"))
plotMA(res, ylim = c(-5, 5))
dev.off()

cat("All plots saved in:", plots_dir, "\n")
