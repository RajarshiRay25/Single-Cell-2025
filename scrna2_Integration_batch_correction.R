# This script will ensure the smooth analysis of scRNA datasets across several experimental setups and analytical mediums.

## Libraries

library(Seurat)
library(SeuratData)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(gridExtra)

## Access data files - GSE180665

dirs <- list.dirs(path = "data/",recursive = F,full.names = F)
dirs

for (x in dirs) {
  name <- gsub('_filtered_feature_bc_matrix','', x)
  
  cts <- ReadMtx(mtx = paste0('data/',x,'/matrix.mtx.gz'),
                 features = paste0('data/',x,'/features.tsv.gz'),
                 cells = paste0('data/',x,'/barcodes.tsv.gz'))
  
  # create seurat objects
  assign(name, CreateSeuratObject(counts = cts))
}

name

## merge datasets

merged_seurat <- merge(HB17_background, y = c(HB17_PDX, HB30_PDX),
                       add.cell.ids = ls()[15:17],
                       project = 'HB')


merged_seurat
View(merged_seurat)


# Perform QC

View(merged_seurat@meta.data)

## Include sample name column within meta data from the rownames of metadata 

merged_seurat$sample <- rownames(merged_seurat@meta.data)

## Split sample column into separate sample information

merged_seurat@meta.data <- separate(merged_seurat@meta.data,col="sample",into=c("Patient","Type","Barcode"),sep="_")


## Calculate the mitochondrial percentage  for QC

merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat,pattern = '^MT-')

View(merged_seurat@meta.data)

## Perform filtering of QC parameters as per convention to retain only significant counts

merged_seurat_filter <- subset(merged_seurat, subset = nCount_RNA > 800 & nFeature_RNA > 500 & mitoPercent < 10)


View(merged_seurat_filter@meta.data)

## Perform the standard scRNA analysis workflow to prepare the data
# Normalise the counts
merged_seurat_filter <- NormalizeData(object = merged_seurat_filter) 

# Find the variable genes/features
merged_seurat_filter <- FindVariableFeatures(object = merged_seurat_filter)

# Scale the data
merged_seurat_filter <- ScaleData(object = merged_seurat_filter)

View(merged_seurat_filter@assays$RNA$counts.1)

# Perform dimensionality reduction test PCA 

merged_seurat_filter <- RunPCA(object = merged_seurat_filter)

ElbowPlot(merged_seurat_filter)

DimHeatmap(merged_seurat_filter, dims = 1:5, cells = 500, balanced = TRUE)

## Cluster the cells

merged_seurat_filter <- FindNeighbors(merged_seurat_filter, dims = 1:5)
merged_seurat_filter <- FindClusters(merged_seurat_filter)

merged_seurat_filter <- RunUMAP(object = merged_seurat_filter, dims = 1:5)

## Run the visual plots to analyse 

plot1 <- DimPlot(merged_seurat_filter,reduction = "umap")
DimPlot(merged_seurat_filter,reduction = "umap")

# plot -- We see here that same cell type has different clusters which shouldnt happen unless technical variation is there
p1 <- DimPlot(merged_seurat_filter, reduction = 'umap', group.by = 'Patient')
p2 <- DimPlot(merged_seurat_filter, reduction = 'umap', group.by = 'Type',
              cols = c('red','green','blue'))

grid.arrange(p1, p2, ncol = 2, nrow = 2)

# Perform integration and correct batch effects

obj.list <- SplitObject(merged_seurat_filter,split.by = "Patient")
View(obj.list)


for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}

# Select the features for integration

features <- SelectIntegrationFeatures(object.list = obj.list)

# Determine the integration anchors where the clustering similiarity is possible


anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)

# integrate data
seurat.integrated <- IntegrateData(anchorset = anchors)

# Scale data, run PCA and UMAP and visualize integrated data
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:5)


p3 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Patient')
p4 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Type',
              cols = c('red','green','blue'))


grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)

# Saving the data

saveRDS(seurat.integrated, file = "seurat_integrated.rds")


# Load the history

seurat.integrated <- readRDS("seurat_integrated.rds")
