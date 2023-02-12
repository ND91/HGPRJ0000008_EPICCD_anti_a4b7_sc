#!/usr/bin/env Rscript
# The goal of this script is to import the masscytometry data from Jan Verhoeff, convert it to a SCE object, harmonize the nomenclature with the scRNAseq experiment, and add sample and donor metadata.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(paste0("Script needs 4 arguments. Current input is:", args))
}

library(Matrix)
library(dplyr)
library(SingleCellExperiment)

masscytometry_rds <- args[1]
files_pbmc_masscytometry_xlsx <- args[2]
sample_metadata_xlsx <- args[3]
feature_metadata_xlsx <- args[4]
sce_annotated_rds <- args[5]

masscytometry <- readRDS(masscytometry_rds)
files_pbmc_masscytometry <- readxl::read_excel(files_pbmc_masscytometry_xlsx)
sample_metadata <- readxl::read_excel(sample_metadata_xlsx)
feature_metadata <- readxl::read_excel(feature_metadata_xlsx) %>%
  data.frame(., row.names = .$`Target protein`)

# Merge all the annotations into the event metadata
sample_metadata_full <- sample_metadata %>%
  dplyr::left_join(files_pbmc_masscytometry, by = "Sample_ID") %>%
  dplyr::mutate(Run_ID_short = gsub("^(c[0-9]+)_.+$", "\\1", Run_ID))

event_metadata_annotated <- masscytometry %>%
  dplyr::select(SOM, Phenotype, file, Lineage) %>%
  dplyr::mutate(
    manual_l4 = case_when(
      Phenotype %in% c("CD27- B Cells") ~ "B naive",
      Phenotype %in% c("CD27+ Memory B Cells") ~ "B memory",
      Phenotype %in% c("cDCs") ~ "CDC2",
      Phenotype %in% c("Classical Monocytes") ~ "Classical monocyte",
      Phenotype %in% c("Intermediate Monocytes") ~ "Intermediate monocyte",
      Phenotype %in% c("Non-classical Monocytes") ~ "Non-classical monocyte",
      Phenotype %in% c("DN TCells") ~ "DNT",
      Phenotype %in% c("Naive CD4 TCells") ~ "CD4 naive",
      Phenotype %in% c("Naive CD8 TCells") ~ "CD8 naive",
      Phenotype %in% c("NK CD16- CD27+") ~ "NK CD16-CD27+",
      Phenotype %in% c("pDCs") ~ "PDC",
      Phenotype %in% c("Doublets", "HLA-DR+ CD3+ CD14dim - possible doublets") ~ "Multiplet",
      Phenotype %in% c("Lin negative - CD45RO+ CD25+") ~ "Lin-CD45RO+CD25+",
      Phenotype %in% c("Lin negative - CD69++ CCR4+ a4b7+") ~ "Lin-CD69++CCR4+a4b7+",
      TRUE ~ Phenotype),
    manual_l3 = case_when(
      manual_l4 %in% c("Lin-CD45RO+CD25+", "Lin-CD69++CCR4+a4b7+") ~ "Lin-",
      manual_l4 %in% c("NK CD16-CD27+", "NK CD16-") ~ "NK CD16-",
      TRUE ~ manual_l4),
    manual_l2 = case_when(
      manual_l3 %in% c("Classical monocyte", "Intermediate monocyte", "Non-classical monocyte") ~ "Monocyte",
      manual_l3 %in% c("CD4 CTL", "CD4 naive", "CD4 proliferating", "CD4 TCM", "CD4 TEM", "CD4 Treg") ~ "CD4T",
      manual_l3 %in% c("CD8 CTL", "CD8 naive", "CD8 proliferating", "CD8 TCM", "CD8 TEM") ~ "CD8T",
      manual_l3 %in% c("CDC2") ~ "CDC",
      manual_l3 %in% c("DNT") ~ "other T",
      manual_l3 %in% c("NK CD16-", "NK CD16+") ~ "NK",
      TRUE ~ manual_l3),
    manual_l1 = case_when(
      manual_l2 %in% c("CDC", "Platelet") ~ "Myeloid",
      manual_l2 %in% c("B memory", "B naive") ~ "B",
      manual_l2 %in% c("Monocyte", "CDC") ~ "Myeloid",
      manual_l2 %in% c("CD4T", "CD8T", "other T") ~ "T",
      manual_l2 %in% c("Multiplet") ~ "Multiplet",
      TRUE ~ manual_l2)) %>%
  dplyr::left_join(sample_metadata_full, by = c("file" = "Run_ID_short"))

masscytometry_expr <- masscytometry %>%
  dplyr::select(-c(SOM, Phenotype, file, Lineage)) %>%
  dplyr::rename(ITGA4 = CD49d,
                ITGAL = CD11a,
                CD8A = CD8a,
                IL2RA = CD25,
                TNR4 = CD134,
                TNR6 = CD95,
                HAVR2 = TIM.3,
                PDCD1 = PD.1,
                KLRB1 = CD161,
                TNR9 = `4.1BB`,
                HLADR = HLA.DR,
                CES1 = CES.1,
                IL7RA = CD127)

masscytometry_expr <- Matrix::Matrix(t(masscytometry_expr), sparse = T)

sce_annotated <- SingleCellExperiment(list(ncounts = masscytometry_expr),
                                      colData = DataFrame(event_metadata_annotated),
                                      rowData = DataFrame(feature_metadata))

colnames(sce_annotated) <- paste0("E", formatC(1:ncol(masscytometry_expr), width = nchar(ncol(masscytometry_expr)), format = "d", flag = "0"))

saveRDS(sce_annotated, sce_annotated_rds, compress = "gzip")

sessionInfo()
