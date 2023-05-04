#!/usr/bin/env Rscript

# The goal of this script is to merge the celltype and sample annotations into the seuratobject.

library(Seurat)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(paste0("Script needs 4 arguments. Current input is:", args))
}

seurat_rds <- args[1]
celltypes_csv <- args[2]
sample_metadata_xlsx <- args[3]
seurat_annotated_rds <- args[4]

seuratObject <- readRDS(seurat_rds)
celltypes <- read.csv(celltypes_csv)
sample_metadata <- readxl::read_excel(sample_metadata_xlsx) %>%
  dplyr::mutate(Sample_ID = paste0("S", Sample_ID))

seuratObject_metadata <- seuratObject@meta.data %>%
  dplyr::mutate(CellID = rownames(.)) %>%
  dplyr::left_join(celltypes, by = "CellID") %>%
  dplyr::left_join(sample_metadata, by = c("sampleID" = "Sample_ID")) %>%
  data.frame(., row.names = .$CellID)

seuratObject <- seuratObject[,rownames(seuratObject_metadata)]
seuratObject@meta.data <- seuratObject_metadata

seuratObject <- DietSeurat(seuratObject, 
                        counts = TRUE, 
                        data = TRUE, 
                        scale.data = FALSE)

seuratObject <- SCTransform(seuratObject, verbose = FALSE, conserve.memory = TRUE)

seuratObject <- RunPCA(object = seuratObject, npcs = 100, seed.use = 79804123)
ElbowPlot(seuratObject, ndims = 100)

seuratObject <- FindNeighbors(seuratObject, reduction = "pca", dims = 1:31)
seuratObject <- FindClusters(seuratObject, resolution = 0.5, verbose = FALSE)
seuratObject <- RunUMAP(seuratObject, dims = 1:31, seed.use = 512352)

saveRDS(seuratObject, seurat_annotated_rds, compress = "gzip")

sessionInfo()