#!/usr/bin/env Rscript

# The goal of this script is to create a boxplot of l1 relative to all PBMCs grouped by response and facetted by lineage.

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrastr)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

seurat_rds_path <- args[1]
boxplot_pdf <- args[3]

plotobj <- seuratObject@meta.data %>%
  dplyr::mutate(manual_l1 = as.factor(manual_l1)) %>%
  dplyr::group_by(SampleID, Response, manual_l1, .drop = F) %>%
  dplyr::summarize(Nl1sample = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(SampleID, Response) %>%
  dplyr::mutate(Ncells = sum(Nl1sample),
                Ncellprop = Nl1sample/Ncells,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop)) %>%
  ggplot(aes(x = forcats::fct_reorder(manual_l1, -Ncellprop), y = Ncellprop, col = Response)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
  labs(subtitle = "Proportion relative to all PBMCs",
       y = "Proportion") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(width = 5, height = 5, file = boxplot_pdf)
print(plotobj)
dev.off()

png(width = 4, height = 4, units = "in", res = 120, file = boxplot_png)
print(plotobj)
dev.off()
