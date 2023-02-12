#!/usr/bin/env Rscript

# The goal of this script is to merge the demultiplexing, celltype, and sample annotations into the seuratobject.

library(Seurat)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(paste0("Script needs 5 arguments. Current input is:", args))
}

seurat_rds <- args[1]
demultiplexed_csv <- args[2]
celltypes_csv <- args[3]
sample_metadata_csv <- args[4]
seurat_annotated_rds <- args[5]

seuratObject <- readRDS(seurat_rds)
demultiplexed <- read.csv(demultiplexed_csv)
celltypes <- read.csv(celltypes_csv)
sample_metadata <- readxl::read_excel(sample_metadata_csv)

seuratObject_metadata <- seuratObject@meta.data %>%
  dplyr::left_join(demultiplexed, by = "CellID") %>%
  dplyr::filter(!SampleID %in% c("Negative", "Multiplets")) %>%
  dplyr::left_join(celltypes, by = "CellID") %>%
  dplyr::left_join(sample_metadata, by = c("SampleID" = "Sample_ID")) %>%
  data.frame(., row.names = .$CellID)

seuratObject <- seuratObject[,rownames(seuratObject_metadata)]
seuratObject@meta.data <- seuratObject_metadata

saveRDS(seuratObject, seurat_annotated_rds, compress = "gzip")

sessionInfo()