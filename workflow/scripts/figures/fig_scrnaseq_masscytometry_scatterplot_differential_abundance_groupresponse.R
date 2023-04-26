#!/usr/bin/env Rscript

# The goal of this script is to create a scatterplot of the effect sizes when comparing R with NR for scRNAseq and CyTOF.

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrastr)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

scrnaseq_dacs_csv <- args[1]
cytof_dacs_csv <- args[2]
scatterplot_pdf <- args[3]

scrnaseq_dacs <- read.csv(scrnaseq_dacs_csv) %>%
  dplyr::select(BaselineProp.clusters, PropRatio, Tstatistic, P.Value, FDR) %>%
  dplyr::rename(Celltype = BaselineProp.clusters,
                propratio_scrnaseq = PropRatio,
                t_scrnaseq = Tstatistic,
                pvalue_scrnaseq = P.Value,
                fdr_scrnaseq = FDR) %>%
  dplyr::mutate(log2_propratio_scrnaseq = log2(propratio_scrnaseq))
cytof_dacs <- read.csv(cytof_dacs_csv) %>%
  dplyr::select(BaselineProp.clusters, PropRatio, Tstatistic, P.Value, FDR) %>%
  dplyr::rename(Celltype = BaselineProp.clusters,
                propratio_cytof = PropRatio,
                t_cytof = Tstatistic,
                pvalue_cytof = P.Value,
                fdr_cytof = FDR) %>%
  dplyr::mutate(log2_propratio_cytof = log2(propratio_cytof))

dacs_overlapping <- scrnaseq_dacs %>%
  dplyr::inner_join(cytof_dacs, by = "Celltype") %>%
  dplyr::mutate(concordance = case_when(
                    log2_propratio_scrnaseq < 0 & log2_propratio_cytof < 0 ~ "concordant",
                    log2_propratio_scrnaseq > 0 & log2_propratio_cytof > 0 ~ "concordant",
                    log2_propratio_scrnaseq < 0 & log2_propratio_cytof > 0 ~ "discordant",
                    log2_propratio_scrnaseq > 0 & log2_propratio_cytof < 0 ~ "discordant",
                  ),
                concordance = factor(concordance, levels = c("discordant", "concordant")),
                label = ifelse(concordance == "concordant", Celltype, NA))

plotobj <- dacs_overlapping %>%
  ggplot(aes(x = log2_propratio_scrnaseq, y = log2_propratio_cytof, alpha = concordance)) +
  geom_vline(xintercept = 0, col = "#d3d3d3") +
  geom_hline(yintercept = 0, col = "#d3d3d3") +
  geom_point(show.legend = F) +
  geom_label_repel(aes(label = Celltype), show.legend = F) +
  labs(title = "Differential abundance R vs NR",
       subtitle = bquote(log[2]~"(R/NR)"),
       x = "scRNAseq",
       y = "CyTOF") +
  scale_alpha_manual(values = c("discordant" = 0.25,
                                "concordant" = 1)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom")

pdf(width = 7, height = 7, file = scatterplot_pdf)
print(plotobj)
dev.off()

sessionInfo()