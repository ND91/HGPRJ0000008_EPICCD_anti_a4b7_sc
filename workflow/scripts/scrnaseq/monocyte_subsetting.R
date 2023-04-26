#!/usr/bin/env Rscript

# The goal of this script is to subset the monocytes, cluster the data in an unsupervised fashion.

library(Seurat)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

seurat_rds <- args[1]
seurat_monocyte_rds <- args[2]
seurat_umap_csv <- args[3]

seuratObject <- readRDS(seurat_rds)

seuratObject_monocyte <- seuratObject[,which(seuratObject@meta.data$manual_l2 == "Monocyte")]

seuratObject_monocyte <- DietSeurat(seuratObject_monocyte, 
                                    counts = TRUE, 
                                    data = TRUE, 
                                    scale.data = FALSE)

seuratObject_monocyte <- SCTransform(seuratObject_monocyte, verbose = FALSE, conserve.memory = TRUE)
seuratObject_monocyte <- RunPCA(object = seuratObject_monocyte, npcs = 100, seed.use = 67894123)
seuratObject_monocyte <- FindNeighbors(seuratObject_monocyte, reduction = "pca", dims = 1:15)
seuratObject_monocyte <- FindClusters(seuratObject_monocyte, resolution = 0.5, verbose = FALSE)
seuratObject_monocyte <- RunUMAP(seuratObject_monocyte, dims = 1:15, seed.use = 512352)

saveRDS(seuratObject_monocyte, seurat_monocyte_rds, compress = "gzip")

umap_df <- data.frame(CB = colnames(seuratObject_monocyte),
                      Embeddings(seuratObject_monocyte[["umap"]]),
                      seuratObject_monocyte@meta.data)

write.csv(umap_df, seurat_umap_csv)

sessionInfo()