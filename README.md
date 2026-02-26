# DESeq2 sex-biased expression analysis in R

This repository contains an R script that runs a complete DESeq2
workflow to identify differentially expressed genes between female
ovary and male testis samples, converts normalized counts to TPM, and
generates several diagnostic plots (MA plot, PCA, per-gene counts,
heatmap of top genes).

The code is designed as a template you can adapt to other two-condition
RNA-seq comparisons.

## Files

- `deseq2_sex_diff_expression.R`  
  Main script implementing:
  - DESeq2 differential expression (Male vs Female).  
  - TPM conversion from normalized counts.  
  - Extraction of the top 150 up- and down-regulated genes.  
  - MA plot, PCA, gene-level counts plot, and heatmap of top genes.

- `counts.reads.txt` (not included)  
  Example input file expected by the script; you should provide your own
  counts table in this format (see below).

## Input format

`counts.reads.txt` should be a tab-delimited text file with:

- First column: gene IDs (used as row names).  
- One column per RNA-seq sample (e.g., 6 columns for 3 female ovary and
  3 male testis samples). Column names in the script are:

  - `Female.ovary_1Aligned.sortedByCoord.out.bam`  
  - `Female.ovary_2Aligned.sortedByCoord.out.bam`  
  - `Female.ovary_3Aligned.sortedByCoord.out.bam`  
  - `Male.testis_1Aligned.sortedByCoord.out.bam`  
  - `Male.testis_2Aligned.sortedByCoord.out.bam`  
  - `Male.testis_3Aligned.sortedByCoord.out.bam`

  You can change these names in the script to match your own file.

- An additional column called `Length` with gene lengths in base pairs,
  used to compute TPM.

Example (header + a couple of rows):

```text
gene_id    Female.ovary_1Aligned.sortedByCoord.out.bam  ...  Male.testis_3Aligned.sortedByCoord.out.bam  Length
GENE1      120                                         ...  305                                         1450
GENE2      0                                           ...  25                                          2100
...
