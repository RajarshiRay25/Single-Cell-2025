 # Import libraries

library(Seurat)
library(tidyverse)

# Data loading

sc.data <- readRDS('seurat_integrated.rds')
View(sc.data)
View(sc.data@meta.data)
# Visualise Data clusters

clusters <- DimPlot(sc.data,reduction = 'umap',group.by = 'seurat_clusters',label = T)
clusters

background <- DimPlot(sc.data,reduction = 'umap',group.by = 'Type')
background

clusters | background
