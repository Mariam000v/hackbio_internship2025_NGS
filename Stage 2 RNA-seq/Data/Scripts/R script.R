# ============================================================
# Differential Expression Analysis (DESeq2) + Gene Annotation
# Multi-Comparison Script for BD-II Microglia RNA-seq
# ============================================================

# --- 1. Setup Environment ---
setwd("C:/meme plans/Bioinformatics internships/hackbio project2/R2")

# --- 2. Load Required Libraries ---
library(readxl)
library(DESeq2)
library(pheatmap)
library(biomaRt)

# --- 3. Import Data ---
gene_counts <- read.delim("gene_counts.txt", header = TRUE)
metadata <- read_excel("metadata.xlsx", col_names = TRUE)

# --- 4. Prepare Metadata ---
metadata$Group <- factor(metadata$Group, levels = c("SHC","FHC","SBD","FBD"))
metadata$State <- factor(metadata$State, levels = c("Healthy control", "Bipolar disorder"))

# --- 5. Prepare Counts Data ---
counts <- gene_counts[, 7:14]  # adjust columns if needed
rownames(counts) <- gene_counts$Geneid

# --- 6. Function to Run DESeq2, Filter DEGs, Volcano & Annotation ---
run_DE_analysis <- function(counts, metadata, group1, group2, outprefix) {
  
  # Subset metadata for the two groups
  subset_meta <- metadata[metadata$Group %in% c(group1, group2), ]
  
  # Subset counts using SampleID (safer than rownames)
  subset_counts <- counts[, subset_meta$SampleID]
  
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(
    countData = subset_counts,
    colData = subset_meta,
    design = ~ Group
  )
  
  # Run DESeq2
  dds <- DESeq(dds)
  res <- results(dds)
  
  # Normalize counts
  norm_counts <- counts(dds, normalized = TRUE)
  
  # Filter significant DEGs
  res <- na.omit(res)
  sig_res <- res[res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]
  upregulated <- sig_res[sig_res$log2FoldChange > 0, ]
  downregulated <- sig_res[sig_res$log2FoldChange < 0, ]
  
  # Volcano Plot
  sig_genes <- res$padj < 0.05 & abs(res$log2FoldChange) >= 1
  colors <- ifelse(sig_genes,
                   ifelse(res$log2FoldChange > 1, "red", "blue"),
                   "grey")
  plot(res$log2FoldChange, -log10(res$padj),
       col = colors, pch = 19, cex = 0.5,
       xlab = "Log2 Fold Change",
       ylab = "-log10(Adjusted P-value)",
       main = paste("Volcano Plot:", group1, "vs", group2))
  abline(v = c(-1, 1), h = -log10(0.05), lty = 2, col = "black")
  legend("bottomright",
         legend = c("Upregulated", "Downregulated", "Not significant"),
         col = c("red", "blue", "grey"), pch = 15, cex = 0.7)
  
  # Heatmap (Scaled DEGs)
  if (nrow(sig_res) > 1) {
    degs <- norm_counts[rownames(sig_res), ]
    degs_scaled <- t(scale(t(degs)))
    annotation_col <- data.frame(Group = subset_meta$Group)
    rownames(annotation_col) <- subset_meta$SampleID
    
    heat_colors <- colorRampPalette(c("blue", "white", "red"))(100)
    
    pheatmap(degs_scaled,
             cluster_rows = TRUE, cluster_cols = TRUE,
             show_rownames = FALSE, show_colnames = TRUE,
             annotation_col = annotation_col,
             color = heat_colors,
             main = paste("Heatmap of DEGs:", group1, "vs", group2))
  }
  
  # Gene Annotation via biomaRt
  up_gene_ids <- gsub("\\..*", "", rownames(upregulated))
  down_gene_ids <- gsub("\\..*", "", rownames(downregulated))
  
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "uswest")
  
  annot_data <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "description"),
    filters = "ensembl_gene_id",
    values = unique(c(up_gene_ids, down_gene_ids)),
    mart = ensembl
  )
  
  # Merge annotations
  up_annot <- merge(
    data.frame(ensembl_gene_id = up_gene_ids, upregulated),
    annot_data,
    by = "ensembl_gene_id",
    all.x = TRUE
  )
  down_annot <- merge(
    data.frame(ensembl_gene_id = down_gene_ids, downregulated),
    annot_data,
    by = "ensembl_gene_id",
    all.x = TRUE
  )
  
  # Save annotated DEGs
  write.csv(up_annot, paste0(outprefix, "_upregulated_genes.csv"), row.names = FALSE)
  write.csv(down_annot, paste0(outprefix, "_downregulated_genes.csv"), row.names = FALSE)
  
  cat("\n✅ DEGs for", group1, "vs", group2, "saved successfully!\n")
}

# --- 7. Run Analyses for All Comparisons ---
run_DE_analysis(counts, metadata, "SBD", "FBD", "SBD_vs_FBD")
run_DE_analysis(counts, metadata, "SBD", "SHC", "SBD_vs_SHC")
run_DE_analysis(counts, metadata, "FBD", "FHC", "FBD_vs_FHC")
