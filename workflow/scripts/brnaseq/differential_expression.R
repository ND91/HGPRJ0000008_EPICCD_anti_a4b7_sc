#!/usr/bin/env Rscript

# The goal of this script is to perform differential expression analysis

library(dplyr)
library(SummarizedExperiment)
library(GenomicRanges)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(DESeq2)
library(readxl)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(paste0("Script needs 4 arguments. Current input is:", args))
}

se_rds <- args[1]
dds_rds <- args[2]
rld_rds <- args[3]
degs_csv <- args[4]

se <- readRDS(se_rds)
colData(se)$Response_coded <- factor(ifelse(colData(se)$Response == "Responder", "R", "NR"), levels = c("NR", "R"))

se <- se[,colData(se)$PBMC_scrnaseq == "Yes"]

dds <- DESeqDataSet(se, design = ~Response_coded+Sex+Age)
rld <- rlog(dds)

dds <- DESeq(dds)

degs_rvnr <- results(dds, contrast = c("Response_coded", "R", "NR"), independentFiltering = T) %>%
  data.frame() %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::arrange(pvalue) %>%
  dplyr::mutate(encode = rownames(.),
                ensembl = gsub("\\.[0-9]+$", "", encode),
                symbol = mapIds(org.Hs.eg.db,
                                keys = ensembl,
                                keytype="ENSEMBL",
                                column="SYMBOL", 
                                multiVals = "first"),
                entrez = mapIds(org.Hs.eg.db,
                                keys = ensembl,
                                keytype="ENSEMBL",
                                column="ENTREZID", 
                                multiVals = "first"))
  
saveRDS(dds, dds_rds, compress = "gzip")
saveRDS(rld, rld_rds, compress = "gzip")
write.csv(degs_rvnr, degs_csv)

sessionInfo()