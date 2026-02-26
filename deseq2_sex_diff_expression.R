# deseq2_sex_diff_expression.R
# DESeq2 workflow for differential expression (Male vs Female)
# with TPM conversion and top-gene visualization.

#-----------------------------
# 1. Setup
#-----------------------------

# Install DESeq2 once (run interactively, not inside scripts for reproducibility):
# if (!require("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# BiocManager::install("DESeq2")

library(DESeq2)
library(RColorBrewer)
library(pheatmap)

#-----------------------------
# 2. Read raw counts
#-----------------------------

# Expected input: counts.reads.txt
# - tab-delimited
# - first column: gene IDs
# - one column per sample
# - an additional column named "Length" with gene lengths (bp) for TPM

count_file <- "counts.reads.txt"

countData <- read.table(
  file      = count_file,
  header    = TRUE,
  row.names = 1,
  sep       = "\t",
  check.names = FALSE
)

# Inspect dimensions
cat("Raw countData dimensions:", dim(countData), "\n")

# Select the 6 RNA-seq samples used here
# Adjust these column names to match your file
sample_cols <- c(
  "Female.ovary_1Aligned.sortedByCoord.out.bam",
  "Female.ovary_2Aligned.sortedByCoord.out.bam",
  "Female.ovary_3Aligned.sortedByCoord.out.bam",
  "Male.testis_1Aligned.sortedByCoord.out.bam",
  "Male.testis_2Aligned.sortedByCoord.out.bam",
  "Male.testis_3Aligned.sortedByCoord.out.bam"
)

new_countData <- countData[, sample_cols]
cat("Selected countData dimensions:", dim(new_countData), "\n")

#-----------------------------
# 3. Sample metadata
#-----------------------------

condition <- factor(
  c("Female", "Female", "Female", "Male", "Male", "Male"),
  levels = c("Female", "Male")
)

colData <- data.frame(
  row.names = colnames(new_countData),
  condition = condition
)

print(colData)

#-----------------------------
# 4. Run DESeq2
#-----------------------------

dds <- DESeqDataSetFromMatrix(
  countData = new_countData,
  colData   = colData,
  design    = ~ condition
)

# Prefilter low counts
dds <- dds[rowSums(counts(dds)) > 1, ]

dds <- DESeq(dds)

# Male vs Female (Male / Female)
res <- results(dds, contrast = c("condition", "Male", "Female"), alpha = 0.05)
res <- res[complete.cases(res), ]

cat("Result dimensions:", dim(res), "\n")
print(head(res))

# Export DE results
write.table(
  res,
  file      = "DESeq2_Male_vs_Female_results.txt",
  sep       = "\t",
  quote     = FALSE,
  col.names = NA
)

# Summary
summary(res)

#-----------------------------
# 5. Counts -> TPM
#-----------------------------

# Assuming `countData` has a column named "Length" in bp
if (!"Length" %in% colnames(countData)) {
  stop("counts.reads.txt must contain a 'Length' column for TPM calculation.")
}

norm_counts <- counts(dds, normalized = TRUE)

counts_to_tpm <- function(countMat, geneLengths) {
  rpk <- countMat / (geneLengths / 1000)           # reads per kilobase
  scalingFactors <- colSums(rpk, na.rm = TRUE) / 1e6
  t(t(rpk) / scalingFactors)
}

TPM <- counts_to_tpm(norm_counts, countData$Length)
write.table(
  TPM,
  file      = "TPM_normalized_counts.txt",
  sep       = "\t",
  quote     = FALSE,
  col.names = NA
)

#-----------------------------
# 6. Top genes by adjusted p-value
#-----------------------------

n <- 150
resOrdered <- res[order(res$padj), ]

top_up   <- resOrdered[resOrdered$log2FoldChange >  1, ][1:n, ]
top_down <- resOrdered[resOrdered$log2FoldChange < -1, ][1:n, ]

topResults <- rbind(top_up, top_down)
topResults <- topResults[complete.cases(topResults), ]

write.table(
  topResults,
  file      = "topResults_150.txt",
  sep       = "\t",
  quote     = FALSE,
  col.names = NA
)

cat("Top results preview:\n")
print(topResults[c(1:5, (nrow(topResults)-4):nrow(topResults)),
                 c("baseMean", "log2FoldChange", "padj")])

#-----------------------------
# 7. Plots
#-----------------------------

# MA plot
plotMA(res,
       main = "DESeq2: Male testis vs Female ovary",
       ylim = c(-15, 15))

# PCA on rlog-transformed counts
rld <- rlogTransformation(dds, blind = TRUE)
plotPCA(rld, intgroup = "condition")

# Counts for the most significant gene
plotCounts(dds,
           gene    = which.min(res$padj),
           intgroup = "condition",
           pch      = 19)

# Heatmap of top genes
hmcol <- brewer.pal(11, "RdBu")
norm_counts_all <- counts(dds, normalized = TRUE)

mat_top <- norm_counts_all[rownames(topResults), , drop = FALSE]

pheatmap(
  log(mat_top + 0.5),
  color         = colorRampPalette(rev(hmcol))(255),
  cluster_rows  = TRUE,
  cluster_cols  = TRUE,
  scale         = "none",
  show_rownames = FALSE,
  main          = "Top differentially expressed genes"
)
