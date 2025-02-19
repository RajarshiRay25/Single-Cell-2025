# Define libraries

library(Seurat)
library(SeuratData)
library(tidyverse)

# Define dataset

ifnb <- readRDS('ifnb_harmony.rds')
View(ifnb@meta.data)

# visualize data
clusters <- DimPlot(ifnb, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE)
condition <- DimPlot(ifnb, reduction = 'umap', group.by = 'stim')

condition|clusters

# Find all the DEG biomarkers in the dataset across each cell types

FindAllMarkers(ifnb,
               logfc.threshold = 0.25,
               min.pct = 0.1,
               only.pos = TRUE,
               test.use = 'DESeq2',
               slot = 'counts')

# Find conserved markers across condition types across cell type

marker.feature.cluster.3 <- FindConservedMarkers(ifnb,
                                                 ident.1 = 3,
                                                 grouping.var = 'stim')

View(marker.feature.cluster.3)

# Visualise the top marker plot

FeaturePlot(ifnb,features = c('FCGR3A'),min.cutoff = 'q10')

# Now we move on to annotate the cell types based on meta data

View(ifnb@meta.data)

## To see the name of cell barcodes

Idents(ifnb)

## Set the names to seurat annotations

Idents(ifnb) <- ifnb@meta.data$seurat_annotations
View(Idents(ifnb))


DimPlot(ifnb,reduction = 'umap', label = T)

# Find marker genes between conditions

ifnb$celltype.cond <- paste0(ifnb$seurat_annotations,'_',ifnb$stim)

ifnb$celltype.cond

## Now incorporate the newly condition based annotation in the meta data

Idents(ifnb) <- ifnb$celltype.cond

View(ifnb@meta.data)

DimPlot(ifnb, reduction = 'umap', label = TRUE)

cd16.response <- FindMarkers(ifnb, ident.1 = 'CD16 Mono_STIM', ident.2 = 'CD16 Mono_CTRL')

head(cd16.response)

# plotting conserved features vs DE features between conditions
head(marker.feature.cluster.3)


FeaturePlot(ifnb, features = c('FCGR3A', 'AIF1', 'IFIT1'), split.by = 'stim', min.cutoff = 'q10')









