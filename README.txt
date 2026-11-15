
ASSIGNMENT 2 - SARS-CoV-2 miRNA EXPRESSION ANALYSIS

============================================================
1. DATASET
============================================================

BioProject:
PRJNA901149

Samples:
- Mock 24H: 3 samples
- Mock 72H: 3 samples
- SARS-CoV-2 24H: 3 samples
- SARS-CoV-2 72H: 3 samples

Total samples: 12


============================================================
2. REFERENCE GENOME
============================================================

Human reference genome:
GRCh38


============================================================
3. ANALYSIS WORKFLOW
============================================================

The analysis workflow was:

FASTQ data
-> quality control
-> read trimming
-> STAR alignment to GRCh38
-> featureCounts quantification
-> miRNA expression matrix
-> DESeq2 differential expression
-> GO Biological Process enrichment


============================================================
4. R / BIOCONDUCTOR VERSION
============================================================

Command:

R --version

Output:

R version 4.5.1 (2025-06-13) -- "Great Square Root"


Command:

R -q -e 'BiocManager::version()'

Output:

[1] "3.22"


============================================================
5. METADATA
============================================================

Command:

wc -l results/metadata.csv

Output:

13 results/metadata.csv

The file contains one header plus 12 samples.


============================================================
6. DESEQ2 INSTALLATION
============================================================

Command:

BiocManager::install("DESeq2", ask=FALSE, update=FALSE)

Successful installation:

* DONE (DESeq2)


DESeq2 loading check:

Command:

R -q -e 'library(DESeq2); cat("DESeq2 OK\n")'

Output:

DESeq2 OK


============================================================
7. FINAL miRNA EXPRESSION MATRIX
============================================================

The final expression matrix was constructed from the featureCounts
gene-count file and the miRNA gene-ID list.

Command:

/N/soft/rhel8/r/gnu/4.5.1/bin/R -q -e 'fc<-read.delim("counts/gene_counts.txt",comment.char="#",check.names=FALSE); ids<-readLines("results/miRNA_gene_ids.txt"); m<-fc[fc[["Geneid"]]%in%ids,c(1,7:18)]; colnames(m)<-c("GeneID","SRR22269872","SRR22269873","SRR22269874","SRR22269875","SRR22269876","SRR22269877","SRR22269878","SRR22269879","SRR22269880","SRR22269881","SRR22269882","SRR22269883"); write.csv(m,"results/miRNA_expression_counts.csv",row.names=FALSE); cat("DONE:",nrow(m),"genes x",ncol(m)-1,"samples\n")'

Output:

DONE: 1879 genes x 12 samples


Missing-value check:

Command:

/N/soft/rhel8/r/gnu/4.5.1/bin/R -q -e 'x<-read.csv("results/miRNA_expression_counts.csv",check.names=FALSE); cat("NA values:",sum(is.na(x)),"\n")'

Output:

NA values: 0


Final matrix:
- 1,879 miRNA genes
- 12 samples
- 0 NA values


============================================================
8. DESEQ2 COMPARISON 1
============================================================

Comparison:

Mock vs SARS-CoV-2

Design:

~ Condition

Contrast:

SARS-CoV-2 vs Mock

Command:

LD_PRELOAD=/N/soft/rhel8/gcc/14.2.0/lib64/libstdc++.so.6 R_LIBS_USER="$HOME/R/4.5/library" /N/soft/rhel8/r/gnu/4.5.1/bin/R -q -e 'library(DESeq2); cts<-read.csv("results/miRNA_expression_counts.csv",row.names=1,check.names=FALSE); meta<-read.csv("results/metadata.csv",row.names=1,check.names=FALSE); meta<-meta[colnames(cts),,drop=FALSE]; cts<-as.matrix(cts); storage.mode(cts)<-"integer"; cts<-cts[rowSums(cts)>0,]; meta$Condition<-factor(meta$Condition,levels=c("Mock","SARS-CoV-2")); dds<-DESeqDataSetFromMatrix(countData=cts,colData=meta,design=~Condition); dds<-DESeq(dds); res<-as.data.frame(results(dds,contrast=c("Condition","SARS-CoV-2","Mock"))); res$GeneID<-rownames(res); res<-res[order(res$padj),]; write.csv(res,"results/DEG_Mock_vs_SARSCoV2.csv",row.names=FALSE); sig<-res[!is.na(res$padj)&res$padj<0.05,]; write.csv(sig,"results/DEG_Mock_vs_SARSCoV2_significant.csv",row.names=FALSE); cat("DONE:",nrow(sig),"significant DEGs\n")'

Output:

DONE; 1 significant DEGs


Result:
1 significant DEG


============================================================
9. DESEQ2 COMPARISON 2
============================================================

Comparison:

SARS-CoV-2 72H vs SARS-CoV-2 24H

Only the six SARS-CoV-2 samples were used.

Design:

~ Time

Contrast:

72H vs 24H

Before the analysis, the infected-sample count totals were checked.

Output:

SRR22269872  511251
SRR22269873  229300
SRR22269874  923227
SRR22269878  363015
SRR22269879  723137
SRR22269880  490715


Successful DESeq2 command:

LD_PRELOAD=/N/soft/rhel8/gcc/14.2.0/lib64/libstdc++.so.6 R_LIBS_USER="$HOME/R/4.5/library" /N/soft/rhel8/r/gnu/4.5.1/bin/R -q -e 'library(DESeq2); cts<-read.csv("results/miRNA_expression_counts.csv",row.names=1,check.names=FALSE); meta<-read.csv("results/metadata.csv",row.names=1,check.names=FALSE); meta<-meta[colnames(cts),,drop=FALSE]; keep<-meta$Condition=="SARS-CoV-2"; cts<-as.matrix(cts[,keep,drop=FALSE]); storage.mode(cts)<-"integer"; meta<-meta[keep,,drop=FALSE]; meta$Time<-factor(as.character(meta$Time),levels=c("24","72")); stopifnot(all(colnames(cts)==rownames(meta))); dds<-DESeqDataSetFromMatrix(countData=cts,colData=meta,design=~Time); dds<-DESeq(dds); res<-as.data.frame(results(dds,contrast=c("Time","72","24"))); res$GeneID<-rownames(res); res<-res[order(res$padj),]; write.csv(res,"results/DEG_SARSCoV2_72H_vs_24H.csv",row.names=FALSE); sig<-res[!is.na(res$padj)&res$padj<0.05,]; write.csv(sig,"results/DEG_SARSCoV2_72H_vs_24H_significant.csv",row.names=FALSE); cat("DONE:",nrow(sig),"significant DEGs\n")'

Output:

DONE: 25 significant DEGs


Result:
25 significant DEGs


============================================================
10. GO PACKAGE SETUP
============================================================

clusterProfiler check:

Command:

LD_PRELOAD=/N/soft/rhel8/gcc/14.2.0/lib64/libstdc++.so.6 R_LIBS_USER="$HOME/R/4.5/library" /N/soft/rhel8/r/gnu/4.5.1/bin/R -q -e 'cat("clusterProfiler:",requireNamespace("clusterProfiler",quietly=TRUE),"\n")'

Output:

clusterProfiler: TRUE


org.Hs.eg.db installation:

Command:

LD_PRELOAD=/N/soft/rhel8/gcc/14.2.0/lib64/libstdc++.so.6 R_LIBS_USER="$HOME/R/4.5/library" /N/soft/rhel8/r/gnu/4.5.1/bin/R -q -e 'BiocManager::install("org.Hs.eg.db",lib=Sys.getenv("R_LIBS_USER"),ask=FALSE,update=FALSE,force=TRUE)'

Output:

* DONE (org.Hs.eg.db)


GO.db installation:

Command:

LD_PRELOAD=/N/soft/rhel8/gcc/14.2.0/lib64/libstdc++.so.6 R_LIBS_USER="$HOME/R/4.5/library" /N/soft/rhel8/r/gnu/4.5.1/bin/R -q -e 'BiocManager::install("GO.db",lib=Sys.getenv("R_LIBS_USER"),ask=FALSE,update=FALSE,force=TRUE)'

Output:

* DONE (GO.db)


Ensembl annotation check:

Command:

/N/soft/rhel8/r/gnu/4.5.1/bin/R -q -e 'library(org.Hs.eg.db); print("ENSEMBL" %in% keytypes(org.Hs.eg.db))'

Output:

[1] TRUE


============================================================
11. GO ENRICHMENT
============================================================

GO Biological Process enrichment was performed for both significant
DEG sets.

Method:
- org.Hs.eg.db
- GO.db
- Ensembl identifiers
- Hypergeometric test
- Benjamini-Hochberg correction
- Adjusted p-value < 0.05


Final GO output:

GO COMPLETE: 0 significant terms; comparison 2 = 0 significant terms


GO result files:

GO_BP_Mock_vs_SARSCoV2.csv
GO_BP_Mock_vs_SARSCoV2_significant.csv
GO_BP_SARSCoV2_72H_vs_24H.csv
GO_BP_SARSCoV2_72H_vs_24H_significant.csv


The GO result file was checked and contained GO terms, including:

GO:0035278
miRNA-mediated gene silencing by inhibition of translation

GO:0150128
negative regulation of insulin production

These terms did not pass the adjusted p-value threshold.


============================================================
12. FINAL RESULTS
============================================================

Expression matrix:
1,879 genes x 12 samples

Missing values:
0

Mock vs SARS-CoV-2:
1 significant DEG

SARS-CoV-2 72H vs 24H:
25 significant DEGs

Significant GO terms:
0 for Mock vs SARS-CoV-2
0 for SARS-CoV-2 72H vs 24H


============================================================
13. SOFTWARE
============================================================

- FastQC
- STAR
- featureCounts
- R 4.5.1
- Bioconductor 3.22
- DESeq2
- clusterProfiler
- org.Hs.eg.db
- GO.db


============================================================
14. FINAL OUTPUT FILES
============================================================

- Expression_Levels_All_Genes_Per_Sample.csv
- DEGs_Mock_vs_SARS.csv
- DEGs_SARS_24H_vs_72H.csv
- GO_Enrichment_Mock_vs_SARS.csv
- GO_Enrichment_SARS_24H_vs_72H.csv
- Summary_Mock_vs_SARS.csv
- Summary_SARS_24H_vs_72H.csv
- Assignment2.R
- README.txt
- README.md

