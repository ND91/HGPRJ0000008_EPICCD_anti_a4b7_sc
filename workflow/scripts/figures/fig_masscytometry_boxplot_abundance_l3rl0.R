#!/usr/bin/env Rscript

# The goal of this script is to create a boxplot of l3 relative to all PBMCs grouped by response and facetted by lineage.

library(dplyr)
library(ggplot2)
library(ggrastr)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

masscytometry_rds <- args[1]
boxplot_pdf <- args[2]

masscytometry <- readRDS(masscytometry_rds)

boxplot_ggobj <- masscytometry %>%
  dplyr::mutate(manual_l3 = as.factor(manual_l3)) %>%
  dplyr::group_by(Sample_ID, Response, manual_l3, .drop = F) %>%
  dplyr::summarize(Nl3sample = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Sample_ID, Response) %>%
  dplyr::mutate(Ncells = sum(Nl3sample),
                Ncellprop = Nl3sample/Ncells,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop)) %>%
  ggplot(aes(x = forcats::fct_reorder(manual_l3, -Ncellprop), y = Ncellprop, col = Response)) +
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

pdf(width = 10, height = 5, file = boxplot_pdf)
print(boxplot_ggobj)
dev.off()