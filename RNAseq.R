# Install packages
path <- "/Users/salemo/RNAseq_workshop/airway2"
devtools::install_local(path)

packages <- c(
  "airway2",
  "EnhancedVolcano", "ExploreModelMatrix",
  "apeglm", "pheatmap", "iSEE", "iSEEu"
)
for (package in packages) {
  if (!require(package, quietly = TRUE)) {
    BiocManager::install(package)
  }
}


# Import packages
suppressPackageStartupMessages({
  library(airway2)
  library(tximeta)
  library(DESeq2)
  library(org.Hs.eg.db)
  library(SummarizedExperiment)
  library(ExploreModelMatrix)
  library(apeglm)
  library(pheatmap)
  library(iSEE)
  library(iSEEu)
})

# Explore airway2
dir <- system.file("extdata", package = "airway2")
list.files(dir)
list.files(file.path(dir, "quants"))

# Read the metadata
coldata <- read.delim(file.path(dir, "SraRunTable.txt"))
colnames(coldata)
coldata <- coldata[, c("Run", "cell_line", "treatment")]
colnames(coldata)[colnames(coldata) == "Run"] <- "names"
coldata$files <- file.path(dir, "quants", coldata$names, "quant.sf")
head(coldata)
all(file.exists(coldata$files))


# Import Salmon quantifications
## Import quantifications on the transcript level
st <- tximeta(coldata = coldata, type = "salmon", dropInfReps = TRUE)

column_sums_counts <- colSums(assay(st, "counts"))
columns_sums_abundance <- colSums(assay(st, "abundance"))
print(column_sums_counts)
print(columns_sums_abundance)

# Explore rows
rowRanges(st)
rowData(st)

# Summarizing on gene level
sg <- tximeta::summarizeToGene(st)
sg <- tximeta::addIds(sg, "SYMBOL", gene = TRUE)
rowData(sg)
colData(sg)

# Relevel the treatments and set the control as reference
colData(sg)$treatment <- factor(colData(sg)$treatment)
colData(sg)$treatment <- relevel(colData(sg)$treatment, ref = "Untreated")


vd <- VisualizeDesign(
  sampleData = colData(sg),
  designFormula = ~ cell_line + treatment
)
vd$plotlist

# Generate a DESeqDataSet object from a SummarizedExperiment object
dds <- DESeqDataSet(sg, design = ~ cell_line + treatment)

# Filtering lowly expressed genes remove genes that have zero or single reads
# across the samples
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep, ]
dim(dds)

# Mean-variance stabilization using vst
vsd <- DESeq2::vst(dds, blind = FALSE)
# blind = FALSE, which means that differences between cell lines and treatment
# (the variables in the design) will not contribute to the expected
# variance-mean trend of the experiment.

# Plot the mean-variance relationship
vsn::meanSdPlot(assay(vsd), ranks = FALSE)

# Principal component analysis
DESeq2::plotPCA(vsd, intgroup = "treatment")

# Differential expression analysis
dds <- DESeq2::DESeq(dds)

# Plot the estimated dispersions
DESeq2::plotDispEsts(dds)

# Extract the results
res <- DESeq2::results(dds)
head(res)
summary(res)

## Remove the genes that were filtered out in the independent filtering
hist(res$pvalue[!is.na(res$padj)])

# We also add a couple of extra columns that will be useful for the interactive
# visualization later
rowData(dds)$log10Dispersion <- log10(rowData(dds)$dispersion)

restmp <- DataFrame(res)
restmp$log10BaseMean <- log10(restmp$baseMean)
restmp$mlog10PValue <- -log10(restmp$pvalue)
colnames(restmp) <- paste0("DESeq2_dex_vs_untrt_", colnames(restmp))
rowData(dds) <- cbind(rowData(dds), restmp)

# How do we filter the results based on LFC and or the FDR ?  FDR
res_fdr <- results(dds, alpha = 0.05, contrast = c(
  "treatment", "Dexamethasone",
  "Untreated"
))
table(res_fdr$padj < 0.05)

# LFC
res_lfc <- results(dds, lfcThreshold = 1, contrast = c(
  "treatment", "Dexamethasone",
  "Untreated"
))
summary(res_lfc)


# Plotting
DESeq2::plotMA(res, ylim = c(-5, 5))

DESeq2::resultsNames(dds)
resape <- DESeq2::lfcShrink(dds,
  coef = "treatment_Dexamethasone_vs_Untreated",
  type = "apeglm"
)

DESeq2::plotMA(resape, ylim = c(-5, 5))


# Heatmap
stopifnot(rownames(vsd) == rownames(res))
mat <- assay(vsd)
rownames(mat) <- ifelse(!is.na(rowData(vsd)$SYMBOL),
  rowData(vsd)$SYMBOL, rownames(vsd)
)
mat <- mat[head(order(res$padj), 30), ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(vsd)[, c("treatment"), drop = FALSE])
pheatmap(mat, annotation_col = df)

# Volcano plot
EnhancedVolcano(res,
  lab = rowData(vsd)$SYMBOL,
  x = "log2FoldChange",
  y = "pvalue"
)

# Interactive visualization
dds <- iSEEu::registerAveAbPatterns(dds, "log10BaseMean")
dds <- iSEEu::registerLogFCPatterns(dds, "log2FoldChange")
dds <- iSEEu::registerPValuePatterns(dds, "pvalue")
app <- iSEE(dds, initial = list(
  MAPlot(),
  VolcanoPlot(),
  RowDataTable(),
  FeatureAssayPlot()
))
shiny::runApp(app)

# Export list of DE genes
stopifnot(all(rownames(res) == rownames(dds)))
res$symbol <- rowData(dds)$SYMBOL

res_ordered <- res[order(res$padj), ]
head(res_ordered)

res_ordered_df <- as.data.frame(res_ordered)[seq_len(100), ]
write.table(cbind(id = rownames(res_ordered_df), res_ordered_df),
  file = "results.txt", quote = FALSE, sep = "\t",
  row.names = FALSE
)

# What's next ?
# - Pathway analysis with (clusterProfiler)
# - Gene set enrichment analysis (fgsea)
