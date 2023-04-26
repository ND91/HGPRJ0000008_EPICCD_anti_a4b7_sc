#!/usr/bin/env Rscript

# The goal of this script is to create a arrowplot of all DEGs.

library(dplyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

degs_rds <- args[1]
manual_l3_order_xlsx <- args[2]
arrowplot_degs_pdf <- args[3]

degs <- readRDS(degs_rds)

manual_l3_order <- readxl::read_excel(manual_l3_order_xlsx, col_names = F) %>% pull(...1)

global_sig_degs <- unique(unlist(lapply(degs$table$Responder, function(celltype){celltype[which(celltype$p_adj.glb<0.05),]$gene})))

global_sig_degs_df <- do.call(rbind, lapply(degs$table$Responder, function(celltype){celltype %>%
    data.frame() %>%
    dplyr::filter(gene %in% global_sig_degs) %>%
    dplyr::mutate(Direction = factor(ifelse(stat < 0, "Down", "Up"), levels = c("Down", "Up")),
                  Significant = ifelse(p_adj.glb < 0.05, "Significant", "NS"),
                  Significant = ifelse(is.na(Significant), "NS", Significant))
})) %>%
  dplyr::mutate(cluster_id = factor(cluster_id, levels = manual_l3_order))

plotobj <- ggplot(global_sig_degs_df, aes(x = cluster_id, 
                               y = gene, 
                               size = -log10(p_val),
                               #colour = Direction,
                               shape = Direction,
                               fill = logFC,
                               alpha = Significant)) +
  geom_point() +
  labs(y = "Gene",
       x = "Celltype",
       title = "R relative to NR") +
  scale_alpha_discrete(range = c(0.25, 1)) +
  scale_shape_manual(values = c(Up = 24, Down = 25)) +
  scale_fill_gradient2(high = "coral3", low = "deepskyblue3", mid = "white") +
  scale_color_manual(values = c(Up = "coral3", Down = "deepskyblue3")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.pos = "right")

pdf(width = 7.5, height = 8, file = arrowplot_degs_pdf)
print(plotobj)
dev.off()

sessionInfo()