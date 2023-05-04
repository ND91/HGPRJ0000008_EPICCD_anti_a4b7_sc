#!/usr/bin/env Rscript

# The goal of this script is to cluster the data in an unsupervised fashion.

library(Seurat)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

seurat_rds <- args[1]
seurat_clustered_rds <- args[2]
seurat_umap_csv <- args[3]

seuratObject <- readRDS(seurat_rds)

seuratObject <- DietSeurat(seuratObject, 
                           counts = TRUE, 
                           data = TRUE, 
                           scale.data = FALSE)

seuratObject <- SCTransform(seuratObject, verbose = FALSE, conserve.memory = TRUE)
seuratObject <- RunPCA(object = seuratObject, npcs = 100, seed.use = 879667)
seuratObject <- FindNeighbors(seuratObject, reduction = "pca", dims = 1:43)
seuratObject <- FindClusters(seuratObject, resolution = 0.5, verbose = FALSE)
seuratObject <- RunUMAP(seuratObject, dims = 1:43, seed.use = 4321)

saveRDS(seuratObject, seurat_clustered_rds, compress = "gzip")

umap_df <- data.frame(CB = colnames(seuratObject),
                      Embeddings(seuratObject[["umap"]]),
                      seuratObject@meta.data)

write.csv(umap_df, seurat_umap_csv)

sessionInfo()