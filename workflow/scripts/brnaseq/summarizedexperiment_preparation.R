#!/usr/bin/env Rscript

# The goal of this script is to cleanup the counts coming from featureCounts, annotate it with the necessary feature and sample metadata and subsequently remove feature genes that present no measurable expression in any sample.

library(dplyr)
library(SummarizedExperiment)
library(GenomicRanges)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(readxl)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

raw_counts_txt <- args[1]
sample_metadata_xlsx <- args[2]
cleaned_counts_se_rds <- args[3]

raw_counts <- read.csv(raw_counts_txt, sep = "\t", skip = 1)

# Counts

counts <- raw_counts %>%
  dplyr::select(Geneid, matches("bam")) %>%
  dplyr::rename_all(function(i){
    stringr::str_replace_all(i, '.+\\.bam\\.(.+)\\..+\\.bam', '\\1')
  }) %>%
  data.frame(., row.names = .$Geneid) %>%
  dplyr::select(-Geneid) %>%
  as.matrix()

# Feature metadata

feature_metadata <- raw_counts %>%
  dplyr::select(Geneid, Chr, Start, End, Strand, Length)

feature_metadata_list <- mapply(function(chrom, startpos, endpos, strand){
  cbind(chr = unlist(chrom),
        start = unlist(startpos),
        end = unlist(endpos),
        strand = unlist(strand))
}, chrom = strsplit(feature_metadata$Chr, ";"),
startpos = strsplit(feature_metadata$Start, ";"),
endpos = strsplit(feature_metadata$End, ";"),
strand = strsplit(feature_metadata$Strand, ";")
)

names(feature_metadata_list) <- feature_metadata$Geneid

feature_metadata_df <- do.call(rbind, lapply(names(feature_metadata_list), function(geneid){
  data.frame(feature_metadata_list[[geneid]], 
             Geneid = geneid,
             Genelength = feature_metadata[which(feature_metadata$Geneid == geneid), "Length"])
})) %>%
  dplyr::mutate(ensembl = gsub("(ENSG[0-9]+)\\.[0-9]+", "\\1", Geneid),
                symbol = mapIds(x = org.Hs.eg.db,
                                keys=ensembl,
                                column="SYMBOL",
                                keytype="ENSEMBL",
                                multiVals="first"),
                entrez = mapIds(x = org.Hs.eg.db,
                                keys=ensembl,
                                column="ENTREZID",
                                keytype="ENSEMBL",
                                multiVals="first"))

feature_metadata_grangeslist <- makeGRangesListFromDataFrame(feature_metadata_df,
                                                             split.field = "Geneid", 
                                                             seqnames.field="chr",
                                                             start.field="start",
                                                             end.field="end",
                                                             strand.field="strand",
                                                             keep.extra.columns=T)

feature_metadata_grangeslist <- feature_metadata_grangeslist[rownames(counts)]

# Sample metadata

sample_metadata <- readxl::read_excel(sample_metadata_xlsx) %>%
  dplyr::filter(Sample_ID %in% colnames(counts)) %>%
  data.frame(., row.names = .$Sample_ID)

sample_metadata <- sample_metadata[colnames(counts),]

# Create summarizedExperiment

se_obj <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=counts),
                                                     rowRanges=feature_metadata_grangeslist, 
                                                     colData=sample_metadata)

se_obj <- se_obj[which(rowSums(assay(se_obj)) != 0),]

saveRDS(se_obj, cleaned_counts_se_rds, compress = "gzip")

sessionInfo()