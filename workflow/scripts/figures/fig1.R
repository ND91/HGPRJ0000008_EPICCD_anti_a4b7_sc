#!/usr/bin/env Rscript

# The goal of this script is to create figure 1 for the manuscript.
# Legend: Figure 1 Integrin a4b7 expression in PBMCs obtained from patients with CD. Uniform manifold approximation and projection (UMAP) visualization of the PBMCs as obtained through (A) single-cell RNA-sequencing (scRNAseq) and (B) mass cytometry by time of flight (CyTOF) from CD patients on VDZ that respond (R; N = 4) and that do not respond (NR; N = 4) colored by the cell identity. Expression of the markers used to annotate the PBMCs at the level of (C) gene expression through a dotplot where size and color intensity represent the percentage cells with measurable expression and the median expression, respectively, and a (D) protein expression through a heatmap with the color representing the median expression. (E) Scatterplot representing the percentage celltypes relative to all PBMCs for scRNAseq on the x-axis and CyTOF on the y-axis colored by lineage. (F) Visualization of ITGA4 and ITGB7 gene expression as visualized through UMAP (left) and a jitterplot superposed on a boxplot with celltype on the x-axis and normalized gene expression on the y-axis. (G) Visualization of integrin α4 and integrin α4β7 protein expression as visualized through UMAP (left) and a jitterplot superposed on a boxplot with celltype on the x-axis and normalized gene expression on the y-axis.

library(dplyr)
library(ggplot2)
library(ggrastr)
library(ggrepel)
library(ggpubr)
library(Seurat)
library(SingleCellExperiment)
library(viridis)
library(ComplexHeatmap)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

seurat_rds <- args[1]
sce_rds <- args[2]
manual_l3_order_xlsx <- args[3]
fig2_pdf <- args[4]

scrnaseq_seuratObject <- readRDS(seurat_rds)
masscytometry_sce <- readRDS(sce_rds)

manual_l3_order <- readxl::read_excel(manual_l3_order_xlsx, col_names = T)  %>%
  dplyr::mutate(number = as.numeric(factor(Celltype, levels = Celltype)),
                celltype_number = paste0(number, ". ", Celltype),
                celltype_number = factor(celltype_number, levels = celltype_number))

manual_l3_colors <- manual_l3_order$Color
names(manual_l3_colors) <- manual_l3_order$celltype_number

scrnaseq_umap_df <- data.frame(CB = colnames(scrnaseq_seuratObject),
                               Embeddings(scrnaseq_seuratObject[["umap"]]),
                               scrnaseq_seuratObject@meta.data) %>%
  dplyr::mutate(manual_l3 = factor(manual_l3, levels = manual_l3_order$Celltype),
                manual_l3_number = as.numeric(manual_l3),
                manual_l3_w_number = factor(paste0(as.numeric(manual_l3), ". ", manual_l3), levels = levels(manual_l3_order$celltype_number)))

masscytometry_umap_df <- data.frame(CB = colnames(masscytometry_sce),
                                    reducedDims(masscytometry_sce)[[1]],
                                    colData(masscytometry_sce))  %>%
  dplyr::mutate(manual_l3 = factor(manual_l3, levels = manual_l3_order$Celltype),
                manual_l3_number = as.numeric(manual_l3),
                manual_l3_w_number = factor(paste0(as.numeric(manual_l3), ". ", manual_l3), levels = levels(manual_l3_order$celltype_number)))

# Fig2A

fig2A <- ggplot(scrnaseq_umap_df, aes(x = UMAP_1, y = UMAP_2, col = manual_l3_w_number)) +
  geom_point_rast(show.legend = T, size = 0.5, col = "black") +
  geom_point_rast(show.legend = T, size = 0.25) +
  labs(y = "",
       x = "",
       title = "scRNAseq",
       subtitle = paste0(nrow(scrnaseq_umap_df), " cells")) +
  geom_label_repel(data = scrnaseq_umap_df %>%
                     dplyr::group_by(manual_l3_number, manual_l3_w_number) %>%
                     summarize(x = median(x = UMAP_1),
                               y = median(x = UMAP_2)),
                   mapping = aes(label = manual_l3_number, x = x, y = y, fill = manual_l3_w_number),
                   alpha = 0.75, 
                   show.legend = F,
                   col = "black") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(values = manual_l3_colors, drop = F) +
  scale_fill_manual(values = manual_l3_colors, drop = F) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(face = "bold"))

# Fig2B

fig2B <- ggplot(masscytometry_umap_df, aes(x = UMAP_1, y = UMAP_2, col = manual_l3_w_number)) +
  geom_point_rast(show.legend = T, size = 0.5, col = "black") +
  geom_point_rast(show.legend = T, size = 0.25) +
  labs(y = "",
       x = "",
       title = "CyTOF",
       subtitle = paste0("subsampled to ", nrow(masscytometry_umap_df), " cells")) +
  geom_label_repel(data = masscytometry_umap_df %>%
                     dplyr::group_by(manual_l3_number, manual_l3_w_number) %>%
                     summarize(x = median(x = UMAP_1),
                               y = median(x = UMAP_2)),
                   mapping = aes(label = manual_l3_number, x = x, y = y, fill = manual_l3_w_number),
                   alpha = 0.75, 
                   show.legend = F,
                   col = "black") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(values = manual_l3_colors, drop = F) +
  scale_fill_manual(values = manual_l3_colors, drop = F) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(face = "bold"))

# Fig2C

scrnaseq_marker_genes <- c("PTPRC", "MKI67", "CD3D", "CD4", "CD8A", "SELL", "CCR5", "S100A4", "FOXP3", "TIGIT", "CTLA4", "TRDC", "TRGC1", "GNLY", "KLRB1", "NCAM1", "NKG7", "GZMA", "GZMB", "GZMK", "GZMH", "IL7R", "CD2", "BANK1", "MS4A1", "CD27", "IGHD", "IGHM", "JCHAIN", "HLA-DRA", "CST3", "CD14", "FCGR3A", "CD1C", "FCER1A", "SIGLEC6", "AXL", "IRF4", "IRF8", "CLEC4C", "KLF4", "PPBP", "CD34")

scrnaseq_marker_expr <- GetAssayData(scrnaseq_seuratObject, assay = "RNA")[which(rownames(GetAssayData(scrnaseq_seuratObject, assay = "RNA")) %in% scrnaseq_marker_genes), ]

scrnaseq_marker_expr_median <- data.frame(GeneID = rownames(scrnaseq_marker_expr), scrnaseq_marker_expr) %>%
  tidyr::pivot_longer(-GeneID, names_to = "CB", values_to = "nUMIs") %>%
  dplyr::mutate(CB = stringr::str_replace(CB, "\\.", "-")) %>%
  dplyr::inner_join(data.frame(CB = rownames(scrnaseq_seuratObject@meta.data), 
                               Celltype = scrnaseq_seuratObject@meta.data$manual_l3), 
                    by = "CB") %>%
  dplyr::mutate(Celltype = factor(Celltype), 
                GeneID = factor(GeneID)) %>%
  dplyr::group_by(Celltype, GeneID, .drop = F) %>%
  dplyr::summarize(Median = median(log1p(nUMIs)),
                   Percentage = mean(nUMIs>0)*100) %>%
  dplyr::mutate(Median = ifelse(is.na(Median), 0, Median),
                Percentage = ifelse(is.na(Percentage), 0, Percentage),
                Celltype = factor(Celltype, manual_l3_order$Celltype),
                GeneID = factor(GeneID, scrnaseq_marker_genes))

fig2C <- scrnaseq_marker_expr_median %>% 
  ggplot(aes(x = GeneID, y = Celltype, col = Median)) +
  geom_tile(alpha = 0, col = "grey") +
  geom_point(aes(size = Percentage)) +
  scale_color_viridis() +
  theme_bw() +
  theme(legend.pos = "bottom", 
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing.x=unit(0, "lines"))

# Fig2D

masscytometry_marker_proteins <- c("CD45", "ITGA4", "ITGAL", "CD3", "CD4", "HAVR2", "CD44", "IL2RA", "CD8A", "CCR7", "CD7", "IL7RA", "PDCD1", "CTLA4", "CD28", "CD69", "CD45RA", "CD45RO", "CXCR3", "CCR4", "CCR9", "CD5", "CCR5", "TNR4", "TNR6", "TNR9", "CD9", "CD2", "CD57", "KLRB1", "CD16", "CCR10", "CD27", "HLADR", "CES1", "CD14")

masscytometry_marker_expr <- assay(masscytometry_sce)[which(rownames(masscytometry_sce) %in% masscytometry_marker_proteins),]

masscytometry_marker_expr_median <- data.frame(protein = rownames(masscytometry_marker_expr), as.matrix(masscytometry_marker_expr)) %>%
  tidyr::pivot_longer(-protein, names_to = "EventID", values_to = "expression") %>%
  dplyr::inner_join(data.frame(EventID = colnames(masscytometry_sce), 
                               Celltype = colData(masscytometry_sce)$manual_l3), 
                    by = "EventID") %>%
  dplyr::mutate(Celltype = factor(Celltype), 
                protein = factor(protein)) %>%
  dplyr::group_by(Celltype, protein, .drop = F) %>%
  dplyr::summarize(Median = mean(expression),
                   Percentage = mean(expression>0)*100) %>%
  dplyr::mutate(Median = ifelse(is.na(Median), 0, Median),
                Percentage = ifelse(is.na(Percentage), 0, Percentage),
                Celltype = factor(Celltype, manual_l3_order$Celltype),
                protein = factor(protein, masscytometry_marker_proteins))

masscytometry_marker_expr_median_wide <- masscytometry_marker_expr_median %>%
  tidyr::pivot_wider(-Percentage, names_from = Celltype, values_from = Median)
masscytometry_marker_expr_median_wide_df <- data.frame(masscytometry_marker_expr_median_wide[,-1], row.names = masscytometry_marker_expr_median_wide$protein)
colnames(masscytometry_marker_expr_median_wide_df) <- colnames(masscytometry_marker_expr_median_wide)[-1]

masscytometry_overlapping_celltypes <- manual_l3_order$Celltype[manual_l3_order$Celltype %in% colnames(masscytometry_marker_expr_median_wide_df)]

masscytometry_marker_expr_median_wide_df <- t(masscytometry_marker_expr_median_wide_df[,masscytometry_overlapping_celltypes])
masscytometry_marker_expr_median_wide_df <- masscytometry_marker_expr_median_wide_df[nrow(masscytometry_marker_expr_median_wide_df):1,]

fig2D <- ComplexHeatmap::pheatmap(masscytometry_marker_expr_median_wide_df, 
                                  scale = "column",
                                  color = viridis(1000),
                                  cluster_rows = F,
                                  row_names_side = "left",
                                  name = "arcsinh(5)")

# Fig2E

scrnaseq_l3rl0_df <- scrnaseq_seuratObject@meta.data %>%
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

masscytometry_l3rl0_df <- colData(masscytometry_sce) %>%
  data.frame() %>%
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

scatterplot_ggobj <- scrnaseq_masscytometry_l3rl0_df %>%
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

# Fig2F

scrnaseq_itga4_itgb7_expr_df <- data.frame(CB = colnames(scrnaseq_seuratObject),
                      Embeddings(scrnaseq_seuratObject[["umap"]]),
                      expr = GetAssayData(scrnaseq_seuratObject)[c("ITGA4"),],
                      gene = "ITGA4",
                      scrnaseq_seuratObject@meta.data) %>%
  dplyr::bind_rows(data.frame(CB = colnames(scrnaseq_seuratObject),
                              Embeddings(scrnaseq_seuratObject[["umap"]]),
                              expr = GetAssayData(scrnaseq_seuratObject)["ITGB7",],
                              gene = "ITGB7",
                              scrnaseq_seuratObject@meta.data)) %>%
  dplyr::left_join(manual_l3_order, by = c("manual_l3" = "Celltype")) %>%
  dplyr::mutate(manual_l3 = factor(manual_l3, levels = manual_l3_order$Celltype),
                expr_rank = rank(expr, ties.method="first"),
                label = paste0("scRNAseq: ", gene))

fig2F_left <- ggplot(scrnaseq_itga4_itgb7_expr_df, aes(x = UMAP_1, y = UMAP_2, order = expr_rank)) +
  geom_point_rast(show.legend = F, size = 0.5, col = "#d3d3d3") +
  geom_point_rast(data = scrnaseq_itga4_itgb7_expr_df %>%
                    dplyr::filter(expr>0), aes(col = expr), show.legend = T, size = 0.5) +
  labs(y = "",
       x = "") +
  guides(colour = guide_legend(override.aes = list(size = 3),
                               title = "nUMIs")) +
  facet_wrap(~label, nrow = 2, ncol = 1) +
  scale_color_viridis() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "left",
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(face = "bold"))

fig2F_right <- ggplot(scrnaseq_itga4_itgb7_expr_df, aes(x = manual_l3, y = expr)) +
  geom_jitter_rast(col = "#d3d3d3") +
  geom_boxplot(alpha = 0.9, outlier.shape = NA) +
  facet_wrap(~label, nrow = 2, ncol = 1, scales = "free_y") +
  labs(y = "nUMIs") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "none",
        axis.title.x = element_blank(),
        strip.text.x = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

fig2F <- ggarrange(fig2F_left, fig2F_right, nrow = 1, ncol = 2, widths = c(0.75, 1), align = "hv")

# Fig2G

masscytometry_itga4_itgb7_expr_df <- data.frame(expr = assay(masscytometry_sce)[c("ITGA4"),], protein = "ITGA4", cellID = colnames(masscytometry_sce)) %>%
  dplyr::left_join(data.frame(colData(masscytometry_sce), reducedDims(masscytometry_sce)[[1]], cellID = colnames(masscytometry_sce)),
                   by = "cellID") %>%
  dplyr::left_join(manual_l3_order, by = c("manual_l3" = "Celltype")) %>%
  dplyr::mutate(manual_l3 = factor(manual_l3, levels = manual_l3_order$Celltype),
                expr_rank = rank(expr, ties.method="first"),
                label = paste0("Mass cytometry: integrin a4"),
                protein = factor(protein))

fig2G_left <- ggplot(masscytometry_itga4_itgb7_expr_df, aes(x = UMAP_1, y = UMAP_2, order = expr_rank)) +
  geom_point_rast(show.legend = F, size = 0.5, col = "#d3d3d3") +
  geom_point_rast(data = masscytometry_itga4_itgb7_expr_df %>%
                    dplyr::filter(expr>0), aes(col = expr), show.legend = T, size = 0.5) +
  labs(y = "",
       x = "") +
  guides(colour = guide_legend(override.aes = list(size = 3),
                               title = "expr")) +
  facet_wrap(~label, nrow = 1, ncol = 1) +
  scale_color_viridis() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "left",
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(face = "bold"))

fig2G_right <- ggplot(masscytometry_itga4_itgb7_expr_df, aes(x = manual_l3, y = expr)) +
  geom_jitter_rast(col = "#d3d3d3") +
  geom_boxplot(alpha = 0.9, outlier.shape = NA) +
  facet_wrap(~label, nrow = 1, ncol = 1, scales = "free_y") +
  labs(y = "nUMIs") +
  scale_x_discrete(drop=FALSE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "none",
        axis.title.x = element_blank(),
        strip.text.x = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

fig2G <- ggarrange(fig2G_left, fig2G_right, nrow = 1, ncol = 2, widths = c(0.75, 1), align = "hv")

# Compiled

fig2AB <- ggarrange(fig2A, fig2B, nrow = 2, ncol = 1, labels = c("A", "B"))
fig2CD <- ggarrange(fig2C, fig2D, nrow = 2, ncol = 1, labels = c("C", "D"))
fig2ABCDE <- ggarrange(fig2AB, fig2CD, fig2E, nrow = 1, ncol = 3, widths = c(0.5, 1, 1), labels = c("", "", "E"))
fig2FG <- ggarrange(fig2F, fig2G, nrow = 2, ncol = 1, labels = c("F", "G"), heights = c(1, 0.5), align = "hv")

fig2 <- ggarrange(fig2ABCDE, fig2FG, nrow = 2, ncol = 1)

pdf(width = 11.7, height = 8.3, file = umap_pdf)
print(plotobj)
dev.off()

sessionInfo()