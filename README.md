# Assignment 2 - SARS-CoV-2 miRNA Expression Analysis

## Overview

This project analyzes miRNA expression in human respiratory samples from Mock Control and SARS-CoV-2 infected conditions at 24 hours and 72 hours.

### Samples

- Mock 24H: 3 samples
- Mock 72H: 3 samples
- SARS-CoV-2 24H: 3 samples
- SARS-CoV-2 72H: 3 samples
- Total: 12 samples

## Reference and preprocessing

- Human reference genome: GRCh38
- Quality control: FastQC
- Adapter/quality trimming: performed on all samples
- Alignment: STAR
- Quantification: featureCounts
- miRNA gene identifiers were used to construct the expression matrix.

The final expression matrix contains:

- 1,879 miRNA genes
- 12 samples
- 0 missing/NA values

## Differential expression analysis

### Mock vs SARS-CoV-2

DESeq2 design:

```text
~ Condition
```

Contrast:

```text
SARS-CoV-2 vs Mock
```

Significance threshold:

```text
adjusted p-value (padj) < 0.05
```

Result:

- 1 significant DEG

### SARS-CoV-2 72H vs 24H

Only the six SARS-CoV-2 samples were used.

DESeq2 design:

```text
~ Time
```

Contrast:

```text
72H vs 24H
```

Significance threshold:

```text
adjusted p-value (padj) < 0.05
```

Result:

- 25 significant DEGs

## GO Biological Process enrichment

GO Biological Process enrichment was performed for both significant DEG sets.

The analysis used:

- org.Hs.eg.db
- GO.db
- Ensembl gene identifiers
- Hypergeometric enrichment testing
- Benjamini-Hochberg multiple-testing correction
- adjusted p-value < 0.05 as the significance threshold

### Mock vs SARS-CoV-2

No GO Biological Process terms passed the adjusted p-value threshold.

### SARS-CoV-2 72H vs 24H

GO Biological Process terms were identified, but none passed the adjusted p-value threshold.

## Software

- R 4.5.1
- Bioconductor 3.22
- DESeq2
- clusterProfiler
- org.Hs.eg.db
- GO.db

## Output files

- `Assignment2.R`
- `README.txt`
- `Expression_Levels_All_Genes_Per_Sample.csv`
- `DEGs_Mock_vs_SARS.csv`
- `DEGs_SARS_24H_vs_72H.csv`
- `GO_Enrichment_Mock_vs_SARS.csv`
- `GO_Enrichment_SARS_24H_vs_72H.csv`
- `Summary_Mock_vs_SARS.csv`
- `Summary_SARS_24H_vs_72H.csv`

## Notes

The analysis used STAR alignment and GRCh38-associated Ensembl miRNA gene identifiers. Differential-expression and GO results reflect the reference annotation, filtering, statistical design, and analysis choices used in this workflow.
