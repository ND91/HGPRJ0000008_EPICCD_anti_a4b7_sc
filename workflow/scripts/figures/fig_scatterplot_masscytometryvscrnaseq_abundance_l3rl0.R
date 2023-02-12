#!/usr/bin/env Rscript

# The goal of this script is to create a scatterplot of l3 relative to all PBMCs for all individuals with the scRNAseq data on the x-axis and the mass cytometry data on the y-axis colored by l1.

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrastr)
library(ggrepel)
library(scales)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

seurat_rds <- args[1]
masscytometry_rds <- args[2]
scatterplot_pdf <- args[3]

seuratObject <- readRDS(seurat_rds)
masscytometry <- readRDS(masscytometry_rds)

scrnaseq_l3rl0_df <- seuratObject@meta.data %>%
  dplyr::mutate(manual_l3 = as.factor(manual_l3)) %>%
  dplyr::group_by(SampleID, Response, manual_l3, .drop = F) %>%
  dplyr::summarize(Nl3sample = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(SampleID, Response) %>%
  dplyr::mutate(Ncells = sum(Nl3sample),
                Ncellprop = Nl3sample/Ncells,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop)) %>%
  dplyr::rename(Sample_ID = SampleID,
                Nl3sample_scrnaseq = Nl3sample,
                Ncells_scrnaseq = Ncells,
                Ncellprop_scrnaseq = Ncellprop)

masscytometry_l3rl0_df <- masscytometry %>%
        dplyr::mutate(manual_l3 = as.factor(manual_l3)) %>%
        dplyr::group_by(Sample_ID, Response, manual_l3, .drop = F) %>%
        dplyr::summarize(Nl3sample = n()) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(Sample_ID, Response) %>%
        dplyr::mutate(Ncells = sum(Nl3sample),
                      Ncellprop = Nl3sample/Ncells,
                      Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop)) %>%
  dplyr::rename(Nl3sample_masscytometry = Nl3sample,
                Ncells_masscytometry = Ncells,
                Ncellprop_masscytometry = Ncellprop)

scrnaseq_masscytometry_l3rl0_df <- scrnaseq_l3rl0_df %>%
  dplyr::full_join(masscytometry_l3rl0_df, by = c("Sample_ID", "manual_l3", "Response"))

cor_stats <- summary(lm(Ncellprop_scrnaseq ~ Ncellprop_masscytometry, data = scrnaseq_masscytometry_l3rl0_df))

plotobj <- scrnaseq_masscytometry_l3rl0_df %>%
  ggplot(aes(x = Ncellprop_scrnaseq*100, y = Ncellprop_masscytometry*100)) +
  geom_smooth(method=lm, level=0.99) +
  geom_point() +
  labs(title = "Correlation scRNAseq and CyTOF",
       subtitle = "% PBMCs",
       x = "scRNAseq",
       y = "CyTOF") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom")

pdf(width = 5, height = 5, file = scatterplot_pdf)
print(plotobj)
dev.off()

