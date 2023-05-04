#!/usr/bin/env Rscript
# The goal of this script is to perform UMAP dimension reduction on the masscytometry data.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(paste0("Script needs 4 arguments. Current input is:", args))
}

library(SingleCellExperiment)
library(scater)
library(dplyr)
library(uwot)

sce_annotated_rds <- args[1]
nsubsample <- args[2]
sce_dimred_ss_rds <- args[3]
umap_ss_csv <- args[4]

sce_annotated <- readRDS(sce_annotated_rds)

sce_annotated_ol <- sce_annotated[,which(colData(sce_annotated)$PBMC_scrnaseq == "Yes")]

set.seed(796318)
cells_ss <- sample(ncol(sce_annotated_ol), size = nsubsample, replace = F)

sce_dimred_ss <- sce_annotated_ol[,cells_ss]
umap_df <- DataFrame(scater::calculateUMAP(sce_dimred_ss, exprs_values = "ncounts") %>%
                       data.frame() %>%
                       dplyr::rename(UMAP_1 = 1,
                                     UMAP_2 = 2))
reducedDims(sce_dimred_ss) <- list(UMAP = umap_df)

saveRDS(sce_dimred_ss, sce_dimred_ss_rds)

umap_ss_df <- data.frame(CB = colnames(sce_dimred_ss),
                         reducedDims(sce_dimred_ss)[[1]],
                         colData(sce_dimred_ss))

write.csv(umap_ss_df, umap_ss_csv)

sessionInfo()