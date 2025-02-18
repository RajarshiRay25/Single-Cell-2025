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

merged_seurat <- merge(HB17_background, y = c(HB17_PDX, HB17_tumor, HB30_PDX, HB30_tumor, HB53_background,
                                              HB53_tumor),
                       add.cell.ids = ls()[15:21],
                       project = 'HB')


merged_seurat
View(merged_seurat)


# Perform QC

View(merged_seurat@meta.data)

## Include sample name column within meta data

merged_seurat$sample <- rownames(merged_seurat@meta.data)

## Split sample column into separate sample information

merged_seurat@meta.data <- separate(merged_seurat@meta.data,col="sample",into=c("Patient","Type","Barcode"),sep="_")


# Calculate the mitochondrial percentage  for QC

merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat,pattern = '^MT-')

View(merged_seurat@meta.data)



