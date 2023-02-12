#!/usr/bin/env bash

# The goal of this script is to perform pseudobulk differential expression analysis on the different celltypes.

library(Seurat)
library(muscat)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

seurat_rds <- args[1]
degs_rds <- args[2]
pb_rds <- args[3]

seuratObject <- readRDS(seurat_rds)

seuratObject@meta.data$Response <- ifelse(seuratObject@meta.data$Response == "Non-responder", "Nonresponder", "Responder")

seuratObject <- seuratObject[,-which(seuratObject@meta.data$manual_l3 %in% names(which(table(seuratObject@meta.data$manual_l3)<=10)))]

sce <- as.SingleCellExperiment(seuratObject)

sce <- prepSCE(sce, 
               kid = "manual_l3",
               gid = "Response",
               sid = "SampleID",
               drop = F)

pb <- aggregateData(sce,
                    assay = "counts", 
                    fun = "sum",
                    by = c("cluster_id", "sample_id"))

degs <- pbDS(pb, method = "DESeq2", filter = "none")
degs_results <- degs$table$Responder
degs_results <- lapply(degs_results, function(celltype){celltype[order(celltype$p_val),]})

saveRDS(degs, degs_rds, compress = "gzip")
saveRDS(pb, pb_rds, compress = "gzip")

sessionInfo()