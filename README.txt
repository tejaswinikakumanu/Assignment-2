Assignment 2 - SARS-CoV-2 miRNA Expression Analysis

Dataset:
12 samples total:
- Mock: 24H (3 samples) and 72H (3 samples)
- SARS-CoV-2 infected: 24H (3 samples) and 72H (3 samples)

Reference:
- Human genome build: GRCh38
- GTF annotation used with the GRCh38 reference

Preprocessing and alignment:
1. FastQC was used for sequence quality control.
2. Adapter/quality trimming was performed to generate 12 trimmed FASTQ files.
3. A STAR genome index was generated using the GRCh38 reference and GTF annotation.
4. All 12 trimmed samples were aligned with STAR.
5. Alignment produced 12 BAM files.
6. featureCounts was used to generate the gene count table.
7. A miRNA gene-ID list was used to extract miRNA counts.
8. The final miRNA expression matrix contained 1,879 miRNA genes and 12 samples.
9. The final expression matrix contained 0 missing/NA values.

Statistical analysis:
Software environment:
- R version 4.5.1
- Bioconductor version 3.22
- DESeq2
- clusterProfiler
- org.Hs.eg.db
- GO.db

Sample metadata:
- Mock 24H: SRR22269881, SRR22269882, SRR22269883
- Mock 72H: SRR22269875, SRR22269876, SRR22269877
- SARS-CoV-2 24H: SRR22269878, SRR22269879, SRR22269880
- SARS-CoV-2 72H: SRR22269872, SRR22269873, SRR22269874

Differential expression analysis:
Comparison 1:
- Mock vs SARS-CoV-2
- DESeq2 design: ~ Condition
- Contrast: SARS-CoV-2 vs Mock
- Significant threshold: adjusted p-value (padj) < 0.05
- Significant DEGs: 1

Comparison 2:
- SARS-CoV-2 72H vs SARS-CoV-2 24H
- Analysis used the six infected samples
- DESeq2 design: ~ Time
- Contrast: 72H vs 24H
- Significant threshold: adjusted p-value (padj) < 0.05
- Significant DEGs: 25

GO enrichment:
- GO biological-process enrichment was performed for both significant DEG sets.
- Human annotation was obtained using org.Hs.eg.db and GO.db.
- GO enrichment was evaluated using a hypergeometric test.
- Multiple-testing adjustment used the Benjamini-Hochberg method.
- Significant GO terms were defined using adjusted p-value < 0.05.
- No GO biological-process terms passed the adjusted p-value cutoff for either comparison.

Output files:
- Expression_Levels_All_Genes_Per_Sample.csv
- DEGs_Mock_vs_SARS.csv
- DEGs_SARS_24H_vs_72H.csv
- GO_Enrichment_Mock_vs_SARS.csv
- GO_Enrichment_SARS_24H_vs_72H.csv
- Summary_Mock_vs_SARS.csv
- Summary_SARS_24H_vs_72H.csv

Notes:
The analysis used STAR alignment and Ensembl/GRCh38-associated miRNA gene identifiers. The GO analysis was performed using GO annotations associated with the analyzed gene identifiers. Results are exploratory and depend on the reference annotation and analysis choices described above.
