#!/usr/bin/env bash

# The goal of this script is to perform differential abundance analysis on the different celltypes relative to the total population of PBMCs.

library(Seurat)
library(speckle)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

seurat_rds <- args[1]
dacs_involvedvuninvolved_l1rl0_csv <- args[2]
dacs_involvedvuninvolved_immune_l3rl1_csv <- args[3]

seuratObject <- readRDS(seurat_rds)

# l1rl0

dacs_l1rl0 <- propeller(clusters = seuratObject@meta.data$manual_l0, 
                        sample = seuratObject@meta.data$sampleID, 
                        group = seuratObject@meta.data$Phenotype, 
                        robust = T)

write.csv(dacs_l1rl0, dacs_involvedvuninvolved_l1rl0_csv)

# Immune: l3rl1

seuratObject_immune <- seuratObject[,seuratObject@meta.data$manual_l0 == "Immune"]

dacs_immune_l3rl1 <- propeller(clusters = seuratObject_immune@meta.data$manual_l3, 
                               sample = seuratObject_immune@meta.data$sampleID, 
                               group = seuratObject_immune@meta.data$Phenotype, 
                               robust = T)

write.csv(dacs_immune_l3rl1, dacs_involvedvuninvolved_immune_l3rl1_csv)

sessionInfo()