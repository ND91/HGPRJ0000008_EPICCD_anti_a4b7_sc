#!/usr/bin/env Rscript

# The goal of this script is to create a boxplot of ITGA4 and A4B7 by manual_l3.

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrastr)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

sce_rds <- args[1]
comparison <- args[2]
boxplot_pdf <- args[3]

sce <- readRDS(sce_rds)

expression_df <- assay(sce)[c("ITGA4", "A4B7"),] %>%
  data.frame(., ProteinID = rownames(.)) %>%
  tidyr::pivot_longer(-ProteinID, names_to = "cellID", values_to = "expression") %>%
  dplyr::left_join(colData(sce) %>%
                     data.frame(., cellID = rownames(.)),
                   by = "cellID")

boxplot_ggplot2 <- expression_df %>%
  ggplot(aes(x = forcats::fct_reorder(manual_l3, -expression), y = expression)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_point_rast(position = position_dodge(width=0.75), alpha = 0.5) +
  labs(y = "Expression") +
  facet_wrap(~ProteinID, nrow = 2, scales = "free_y") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(width = 20, height = 10, file = boxplot_pdf)
print(boxplot_ggplot2)
dev.off()

sessionInfo()