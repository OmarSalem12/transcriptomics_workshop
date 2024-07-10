# Install missing packages
packages <- c(
  "Seurat",
  "dplyr",
  "patchwork"
)
for (package in packages) {
  if (!require(package, quietly = TRUE)) {
    install.packages(package)
  }
}

# Import libraries
library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
# After downloading the dataset from:
# https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
# Unzip and untar the file and set the path to the unzipped folder
pbmc_data <- Read10X(
  data.dir = "/Users/salemo/RNAseq_workshop/filtered_gene_bc_matrices/hg19"
)

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(
  counts = pbmc_data, project = "pbmc3k", min.cells = 3, min.features = 200
)

# Let's look at the counts
pbmc.data[1:50, 1:30]

# Calculate QC metrics and add them to metadata
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3
)

## I don't recommend filltering out large number of features
pbmc <- subset(pbmc,
  subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5
)

# Normalize the data
pbmc <- NormalizeData(pbmc,
  normalization.method = "LogNormalize", scale.factor = 10000
)

# FInd HVGs
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scale the data
# Default is HVGs, but you can specify features
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes, vars.to.regress = "percent.mt")

# Dimensionality reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Check the variance captured by each PC
ElbowPlot(pbmc)

# Cluster
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# UMAP (non-linear dimensionality reduction)
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Visualize
DimPlot(pbmc, reduction = "umap", label = TRUE)


# Cell annotations
# Find cluster markers
pbmc_markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# Here we can extract the top 10 markers for each cluster and plot them to
# help annotate the clusters
# Alternatively (or usually in addition), we can plot some canonical markers

cell_markers <- list(
  "Naive CD4+ T" = c("IL7R", "CCR7"),
  "CD14+ Mono" = c("CD14", "LYZ"),
  "Memory CD4+" = c("S100A4"),
  "B" = c("MS4A1"),
  "CD8+ T" = c("CD8A"),
  "FCGR3A+ Mono" = c("FCGR3A", "MS4A7"),
  "NK" = c("GNLY", "NKG7"),
  "DC" = c("FCER1A", "CST3"),
  "Platelet" = c("PPBP")
)

pdf("Canonical_markers.pdf", width = 15)
Seurat::DotPlot(pbmc, features = cell_markers)
dev.off()

# We can also use some functions to Visualize individuak markers of interest

# Violin plot
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
# Feature plot
FeaturePlot(pbmc, features = c(
  "MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
  "CD8A"
))

# From the above steps we should have an idea on the identity of each cell type
new_cluster_ids <- c(
  "Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
  "NK", "DC", "Platelet"
)
names(new_cluster_ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new_cluster_ids)
DimPlot(pbmc, reduction = "umap", label = TRUE) + NoLegend()

# Now that we have annotated the clusters, what's next ?
## Answer the biological questions!! Can you think of any ?
