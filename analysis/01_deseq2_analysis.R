#!/usr/bin/env Rscript
# 01_deseq2_analysis.R


#==============================================================================
# PACKAGE LOADING
#==============================================================================

library("DESeq2")
library("ggplot2")
library("dplyr")

#==============================================================================
# SETUP DIRECTORIES
#==============================================================================

in_dir <- "data"
results_dir <- "results"

if (!dir.exists(results_dir)){
  dir.create(results_dir)
}

#==============================================================================
# DATA LOADING AND PREPROCESSING  
#==============================================================================

# Load count matrix
data <- read.table(paste(in_dir, "gene_counts.txt", sep="/"), 
                   header = TRUE, check.names = FALSE)

head(data)
dim(data)

genes <- data$Geneid

# Extract count columns
count_matrix <- data[, c("Control_1.bam", "Control_2.bam", "Control_3.bam",
                         "Heat_treated_1.bam", "Heat_treated_2.bam", "Heat_treated_3.bam")]

# Set gene IDs as row names
rownames(count_matrix) <- genes

# Clean column names
colnames(count_matrix) <- c("Control_1", "Control_2", "Control_3",
                            "Heat_treated_1", "Heat_treated_2", "Heat_treated_3")

head(count_matrix)

# Create sample information table
sample_info <- data.frame(
  condition = factor(c("Control", "Control", "Control",
                       "Heat_treated", "Heat_treated", "Heat_treated")),
  row.names = colnames(count_matrix)
)

#==============================================================================
# DESEQ2 ANALYSIS
#==============================================================================

# Set up dataset
dds <- DESeqDataSetFromMatrix(countData = round(count_matrix), 
                              colData = sample_info, 
                              design = ~ condition)

# Run DESeq
dds <- DESeq(dds)
normed <- assay(rlog(dds, blind=FALSE))
resultsNames(dds)

# Save normalized expression
write.table(normed, paste(results_dir, "normalized_expression.txt", sep='/'),
            quote=FALSE, sep="\t", row.names = TRUE)

# Filter lowly expressed genes
means <- apply(normed, 1, mean)
expThresh <- 3.4  # log2 of 10 reads
expressed <- means >= expThresh
normed_expr <- normed[expressed,]

# Get results
res <- results(dds, contrast = c("condition", "Heat_treated", "Control"))

# Summary
summary(res)

#==============================================================================
# DIFFERENTIAL EXPRESSION FILTERING
#==============================================================================

# Significant genes
sig_genes <- subset(res, 
                    padj < 0.05 & 
                    abs(log2FoldChange) > 1)

cat("Total significant genes:", nrow(sig_genes), "\n")

# Upregulated genes
up_genes <- subset(res,
                   padj < 0.05 &
                   log2FoldChange > 1)

# Downregulated genes
down_genes <- subset(res,
                     padj < 0.05 &
                     log2FoldChange < -1)

cat("Upregulated genes:", nrow(up_genes), "\n")
cat("Downregulated genes:", nrow(down_genes), "\n")

#==============================================================================
# SAVE RESULTS
#==============================================================================

# Save DESeq2 results
write.csv(as.data.frame(res), paste(results_dir, "deseq2_results.csv", sep="/"), 
          row.names = TRUE)

# Save significant genes
write.csv(as.data.frame(up_genes), paste(results_dir, "upregulated_genes.csv", sep="/"), 
          row.names = TRUE)

write.csv(as.data.frame(down_genes), paste(results_dir, "downregulated_genes.csv", sep="/"), 
          row.names = TRUE)

# Save R objects for next scripts
save(dds, res, up_genes, down_genes, sample_info, 
     file = paste(results_dir, "deseq2_objects.RData", sep="/"))

cat("DESeq2 analysis complete! Results saved in:", results_dir, "\n")
