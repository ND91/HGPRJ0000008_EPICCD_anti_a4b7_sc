#!/usr/bin/env Rscript
# The goal of this script is to import and normalize single cell experiment sample

suppressPackageStartupMessages(library(Seurat))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(paste0("Script needs 4 arguments. Current input is:", args))
}

count_path <- args[1]
seurat_rds_path <- args[2]
runid <- args[3]
lowerbound_numis <- args[4]

cnts <- Read10X(count_path)
expr_cnts <- cnts$`Gene Expression`
abc_cnts <- cnts$`Antibody Capture`

seuratObject <- CreateSeuratObject(counts = expr_cnts, min.cells = 3)
seuratObject[["HTO"]] <- CreateAssayObject(counts = abc_cnts)

# Select cells that have equal to or more than 600 UMIs
seuratObject <- seuratObject[,seuratObject@meta.data$nFeature_RNA>=600]

seuratObject <- NormalizeData(seuratObject)

# Append the run ID into the cellbarcode to prevent any collisions later on.
runid_rfriendly <- paste0("S", gsub("-", "_", runid))
seuratObject@meta.data$RunID <- runid
seuratObject@meta.data$Study <- "Current study"
seuratObject <- RenameCells(object = seuratObject, add.cell.id = runid_rfriendly)
seuratObject@meta.data$CellID <- rownames(seuratObject@meta.data)

saveRDS(seuratObject, seurat_rds_path, compress = "gzip")

sessionInfo()
