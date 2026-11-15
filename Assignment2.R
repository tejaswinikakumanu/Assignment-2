# Assignment 2 - SARS-CoV-2 miRNA Expression Analysis

# Software:
# R 4.5.1
# Bioconductor 3.22
# DESeq2
# clusterProfiler
# org.Hs.eg.db
# GO.db

# Input expression matrix
counts <- read.csv(
  "results/miRNA_expression_counts.csv",
  row.names = 1,
  check.names = FALSE
)

# Sample metadata
metadata <- read.csv(
  "results/metadata.csv",
  row.names = 1,
  check.names = FALSE
)

# Differential-expression results
deg_mock_sars <- read.csv(
  "results/DEG_Mock_vs_SARSCoV2.csv"
)

deg_sars_time <- read.csv(
  "results/DEG_SARSCoV2_72H_vs_24H.csv"
)

# GO enrichment results
go_mock_sars <- read.csv(
  "results/GO_BP_Mock_vs_SARSCoV2.csv"
)

go_sars_time <- read.csv(
  "results/GO_BP_SARSCoV2_72H_vs_24H.csv"
)

# Summary
cat("Expression matrix:", nrow(counts), "genes x",
    ncol(counts), "samples\n")

cat("Mock vs SARS-CoV-2 significant DEGs:",
    sum(!is.na(deg_mock_sars$padj) &
        deg_mock_sars$padj < 0.05), "\n")

cat("SARS-CoV-2 72H vs 24H significant DEGs:",
    sum(!is.na(deg_sars_time$padj) &
        deg_sars_time$padj < 0.05), "\n")

cat("GO analysis files loaded successfully.\n")
