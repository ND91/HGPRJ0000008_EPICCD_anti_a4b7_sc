#!/usr/bin/env Rscript

# The goal of this script is to create a boxplot of l3 relative to all PBMCs grouped by response and facetted by lineage.

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrastr)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

proportion_xlsx <- args[1]
boxplot_pdf <- args[2]

proportion_dcs <- readxl::read_excel("facs_pdcs.xlsx")

round(wilcox.test(pDC~Response,data=proportion_dcs)$p.value, 3)
wilcox.test(cDC~Response,data=proportion_dcs)

plotobj <- proportion_dcs %>%
  tidyr::pivot_longer(-c(DonorID, Response), names_to = "Celltype", values_to = "Proportion") %>%
  dplyr::mutate(label = ifelse(Celltype == "pDC", 
                               paste0("PDC\np-value = ", round(wilcox.test(pDC~Response,data=proportion_dcs)$p.value, 3)),
                               paste0("CDC2\np-value = ", round(wilcox.test(cDC~Response,data=proportion_dcs)$p.value, 3)))) %>%
  ggplot(aes(x = Response, y = Proportion, col = Response)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
  labs(subtitle = "Proportion relative to input PBMCs",
       y = "Proportion") +
  facet_wrap(~label, ncol = 2, scales = "free_y") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(width = 6, height = 3, file = boxplot_pdf)
print(plotobj)
dev.off()

png(width = 8, height = 4, units = "in", res = 120, file = boxplot_png)
print(plotobj)
dev.off()
