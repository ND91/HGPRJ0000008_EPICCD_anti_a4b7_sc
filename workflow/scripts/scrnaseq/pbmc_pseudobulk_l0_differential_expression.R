#!/usr/bin/env Rscript

# The goal of this script is to create a pseudobulk DDS object for the entire PBMC dataset.

library(Seurat)
library(dplyr)
library(DESeq2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(paste0("Script needs 4 arguments. Current input is:", args))
}

seurat_rds <- args[1]
dds_rds <- args[2]
rld_rds <- args[3]
degs_rvnr_csv <- args[4]

seuratObject <- readRDS(seurat_rds)

pb_sample_metadata <- data.frame(SampleID = seuratObject@meta.data$SampleID,
                                 Sex = seuratObject@meta.data$Sex,
                                 Age = seuratObject@meta.data$Age,
                                 Response = seuratObject@meta.data$Response) %>%
  unique() %>%
  dplyr::mutate(Response_recoded = ifelse(Response == "Responder", "R", "NR"))
rownames(pb_sample_metadata) <- pb_sample_metadata$SampleID

pb_counts <- GetAssayData(seuratObject, slot = "counts") %*% model.matrix(~0+SampleID, data = seuratObject@meta.data)
colnames(pb_counts) <- gsub("SampleID", "", colnames(pb_counts))

pb_counts <- pb_counts[which(Matrix::rowSums(pb_counts) != 0),]

pb_dds <- DESeqDataSetFromMatrix(countData = pb_counts,
                                 colData = pb_sample_metadata[colnames(pb_counts),],
                                 design = ~Response_recoded+Sex+Age)

pb_dds <- DESeq(pb_dds)
pb_rld <- rlog(pb_dds)

degs_rvnr <- results(pb_dds, contrast = c("Response_recoded", "R", "NR")) %>%
  data.frame() %>%
  tibble::rownames_to_column("GeneID") %>%
  dplyr::arrange(pvalue)

saveRDS(pb_dds, dds_rds, compress = "gzip")
saveRDS(pb_rld, rld_rds, compress = "gzip")
write.csv(degs_rvnr, degs_rvnr_csv)

sessionInfo()