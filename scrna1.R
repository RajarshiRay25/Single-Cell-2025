# library import

library(dplyr)
library(Seurat)
library(patchwork)


# Load counts dataset

pbmc.data <- Read10X(data.dir = "D:/Single Cell Tutorial/PBMC1/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19")

# Create Seurat object

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc", min.cells = 3, min.features = 200)
pbmc

# Standard preprocessing and QC 

## Adding mitochondrial count percentage within samples within metadata

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

## View meta data

View(pbmc@meta.data)

## Create violin plots to observe the QC data

VlnPlot(pbmc,features = c("nCount_RNA","nFeature_RNA","percent.mt"),ncol = 3)

# We filter cells that have unique feature counts over 2,500 or less than 200
# We filter cells that have >5% mitochondrial counts

## Develop linear plots to assess the correlation analysis

plot_1 <- FeatureScatter(pbmc,feature1 = "nCount_RNA",feature2 = "percent.mt")
plot_2 <- FeatureScatter(pbmc,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")

plot_1 + plot_2


## Perform filtering

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc


# Normalise the data

pbmc <- NormalizeData(pbmc)
View(pbmc)
# Determine the Most variable genes which are potentially significant expressed high or low in the cells

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst",nfeatures = 2000)

# Find the top 10 variable genes

variable.genes.10 <- head(VariableFeatures(pbmc),10)

# Plot the Variable genes 

plot.1 <- VariableFeaturePlot(pbmc)
LabelPoints(plot = plot.1,points = variable.genes.10,repel = TRUE)

# Scale the data to range the expression level between 0 and 1

all.genes <- rownames(pbmc)
all.genes

pbmc <- ScaleData(pbmc,features = all.genes)

# Linear Dimensional Reduction - PCA

pbmc <- RunPCA(pbmc,features = VariableFeatures(object = pbmc))


VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)


# Elbow method to obtain the dimensionality

ElbowPlot(pbmc)
