#!/usr/bin/env Rscript

# The goal of this script is to create figure 5 for the manuscript.

library(dplyr)
library(ggplot2)
library(ggrastr)
library(ggrepel)
library(ggpubr)
library(ggforce)
library(Seurat)
library(SingleCellExperiment)
library(viridis)
library(scales)
library(DESeq2)
library(stringr)
library(msigdbr)
library(cowplot)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9) {
  stop(paste0("Script needs 9 arguments. Current input is:", args))
}

seurat_rds <- args[1]
monocyte_seurat_rds <- args[2]
scrnaseq_degs_l3_rds <- args[3]
fgsea_list_rds <- args[4]
brnaseq_degs_csv <- args[5]
brnaseq_rld_rds <- args[6]
ccrl_mono_xlsx <- args[7]
manual_l3_order_xlsx <- args[8]
fig5_pdf <- args[9]

seuratObject <- readRDS(seurat_rds)
monocyte_seuratObject <- readRDS(monocyte_seurat_rds)
scrnaseq_degs_l3 <- readRDS(scrnaseq_degs_l3_rds)
fgsea_list <- readRDS(fgsea_list_rds)
brnaseq_degs <- read.csv(brnaseq_degs_csv)
brnaseq_rld <- readRDS(brnaseq_rld_rds)
ccrl_mono <- readxl::read_excel(ccrl_mono_xlsx) %>%
  dplyr::filter(!is.na(Gene) & !is.na(Partner))
manual_l3_order <- readxl::read_excel(manual_l3_order_xlsx, col_names = T)  %>%
  dplyr::mutate(number = as.numeric(factor(Celltype, levels = Celltype)),
                celltype_number = paste0(number, ". ", Celltype),
                celltype_number = factor(celltype_number, levels = celltype_number))

manual_l3_colors <- manual_l3_order$Color
names(manual_l3_colors) <- manual_l3_order$celltype_number

response_colors <- c(`Responder` = "#009E73",
                     `Non-responder` = "#D55E00")

# fig5A

monocyte_umap_df <- data.frame(CB = colnames(monocyte_seuratObject),
                               Embeddings(monocyte_seuratObject[["umap"]]),
                               monocyte_seuratObject@meta.data) %>%
  dplyr::mutate(manual_l3 = factor(manual_l3, levels = manual_l3_order$Celltype),
                manual_l3_number = as.numeric(manual_l3),
                manual_l3_w_number = factor(paste0(as.numeric(manual_l3), ". ", manual_l3), levels = levels(manual_l3_order$celltype_number)))


fig5A_left <- ggplot(monocyte_umap_df, aes(x = UMAP_1, y = UMAP_2, col = Response)) +
  geom_point_rast(show.legend = T, size = 0.75, col = "black") +
  geom_point_rast(show.legend = T, size = 0.5) +
  #geom_density2d() +
  labs(y = "",
       x = "",
       title = "Monocytes") +
  geom_label_repel(data = monocyte_umap_df %>%
                     dplyr::group_by(manual_l3, manual_l3_number, manual_l3_w_number) %>%
                     summarize(x = median(x = UMAP_1),
                               y = median(x = UMAP_2)),
                   mapping = aes(label = manual_l3, x = x, y = y),
                   alpha = 0.75, 
                   show.legend = F,
                   col = "black", 
                   max.overlaps = 100) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(values = response_colors, drop = F) +
  scale_fill_manual(values = response_colors, drop = F) +
  theme_bw() +
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

fig5A_right <- ggplot(monocyte_umap_df, aes(x = UMAP_1, y = UMAP_2, col = Response)) +
  geom_point_rast(show.legend = T, size = 0.5, alpha = 0.25) +
  geom_density2d() +
  labs(y = "",
       x = "",
       title = "") +
  geom_label_repel(data = monocyte_umap_df %>%
                     dplyr::group_by(manual_l3, manual_l3_number, manual_l3_w_number) %>%
                     summarize(x = median(x = UMAP_1),
                               y = median(x = UMAP_2)),
                   mapping = aes(label = manual_l3, x = x, y = y),
                   alpha = 0.75, 
                   show.legend = F,
                   col = "black", 
                   max.overlaps = 100) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(values = response_colors, drop = F) +
  scale_fill_manual(values = response_colors, drop = F) +
  theme_bw() +
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

# fig5B

fig5B <- scrnaseq_degs_l3$`Classical monocyte`$degs %>%
  data.frame() %>%
  dplyr::mutate(Significance = ifelse(padj<0.05, "Significant", "NS"),
                label = ifelse(gene %in% c("CFD", "MSR1", "CCL3", "CCL4", "CXCL2", "VSTM1", "RIPK2", "GPR183"), gene, NA)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_vline(xintercept = 0, col = "#d3d3d3") +
  geom_hline(yintercept = 0, col = "#d3d3d3") +
  geom_point_rast(aes(col = Significance)) +
  geom_label_repel(aes(label = label)) +
  labs(title = "Classical monocyte differential expression",
       subtitle = "Responders relative to non-responders",
       y = bquote('-'~log[10]~'(p-value)'),
       x = bquote('-'~log[2]~'(fold change)')) +
  scale_color_manual(values = c(`Significant` = "#000000", `NS` = "#d3d3d3"), drop = FALSE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom")

# fig5C 

classical_monocytes_brnaseq_scrnaseq <- scrnaseq_degs_l3$`Classical monocyte`$degs %>%
  data.frame() %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::rename_with(.cols = baseMean:padj, function(x){paste0("scrnaseq_", x)}) %>%
  dplyr::inner_join(brnaseq_degs %>%
                      data.frame() %>%
                      dplyr::rename_with(.cols = baseMean:padj, function(x){paste0("brnaseq_", x)}), by = c("gene" = "symbol")) %>%
  dplyr::select(gene, encode, scrnaseq_stat, scrnaseq_pvalue, brnaseq_stat, brnaseq_pvalue) %>%
  dplyr::filter(sign(scrnaseq_stat) == sign(brnaseq_stat))

genes_of_interest <- data.frame(symbol = c("CFD", "MSR1"),
                                encode = c("ENSG00000197766.8", "ENSG00000038945.15"))

fig5C <- data.frame(monocyte_seuratObject@meta.data[,c("CellID", "SampleID", "Response")], t(GetAssayData(monocyte_seuratObject)[genes_of_interest$symbol,])) %>%
  tidyr::pivot_longer(-c(CellID, SampleID, Response), names_to = "Gene", values_to = "Expression") %>%
  dplyr::mutate(SampleID = factor(SampleID, levels = c("F213", "F285", "F291", "F334", "F223", "F200", "F216", "F140"))) %>%
  ggplot(aes(x = SampleID, y = Expression, col = Response)) +
  geom_jitter_rast() +
  geom_boxplot(outlier.shape = NA) +
  scale_color_manual(values = response_colors, drop = F) +
  facet_wrap(~Gene, ncol = 1, scales = "free_y") +
  labs(title = "",
       subtitle = "",
       y = "nUMIs") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.text.x = element_text(face = "bold"))

# fig5D

scrnaseq_classical_monocyte_rld <- rlog(scrnaseq_degs_l3$`Classical monocyte`$dds)
scrnaseq_classical_monocyte_rld <- scrnaseq_classical_monocyte_rld[,c("F213", "F285", "F291", "F334", "F223", "F200", "F216", "F140")]

fig5D_left <- data.frame(colData(scrnaseq_classical_monocyte_rld)[,c("SampleID", "Response")], t(assay(scrnaseq_classical_monocyte_rld)[rownames(scrnaseq_classical_monocyte_rld) %in% genes_of_interest$symbol,])) %>%
  tidyr::pivot_longer(-c(SampleID, Response), names_to = "symbol", values_to = "Expression") %>%
  dplyr::left_join(scrnaseq_degs_l3$`Classical monocyte`$degs %>%
                     data.frame(),
                   by = c("symbol" = "gene")) %>%
  dplyr::mutate(label = paste0(symbol, "\np-value = ", formatC(pvalue, digits = 2, format = "e"))) %>%
  ggplot(aes(x = Response, y = Expression, col = Response)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter_rast() +
  scale_color_manual(values = response_colors, drop = F) +
  facet_wrap(~label, scales = "free_y", ncol = 1) +
  labs(title = "Classical monocytes",
       subtitle = "scRNAseq (pseudobulk)",
       y = bquote(log[2]~'(expression)')) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_text(face = "bold"))

fig5D_right <- data.frame(colData(brnaseq_rld)[,c("Sample_ID", "Response")], t(assay(brnaseq_rld)[rownames(brnaseq_rld) %in% genes_of_interest$encode,])) %>%
  tidyr::pivot_longer(-c(Sample_ID, Response), names_to = "encode", values_to = "Expression") %>%
  dplyr::left_join(genes_of_interest, by = "encode") %>%
  dplyr::select(Sample_ID, Response, symbol, Expression) %>%
  dplyr::rename(SampleID = Sample_ID) %>% 
  dplyr::left_join(brnaseq_degs %>%
                     data.frame(),
                   by = "symbol") %>%
  dplyr::mutate(label = paste0(symbol, "\np-value = ", formatC(pvalue, digits = 2, format = "e"))) %>%
  ggplot(aes(x = Response, y = Expression, col = Response)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter_rast() +
  scale_color_manual(values = response_colors, drop = F) +
  facet_wrap(~label, scales = "free_y", ncol = 1) +
  labs(subtitle = "Bulk RNAseq",
       y = bquote(log[2]~'(expression)')) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_text(face = "bold"))

# fig5E

scrnaseq_degs_l3_classical_monocytes_rld <- rlog(scrnaseq_degs_l3$`Classical monocyte`$dds)

fig5E <- data.frame(SampleID = colData(scrnaseq_degs_l3_classical_monocytes_rld)$SampleID,
                    Response = colData(scrnaseq_degs_l3_classical_monocytes_rld)$Response,
                    t(assay(scrnaseq_degs_l3_classical_monocytes_rld)[c("ITGA4", "ITGB7"),])) %>%
  tidyr::pivot_longer(-c(SampleID, Response), names_to = "Gene", values_to = "Expression") %>%
  dplyr::left_join(scrnaseq_degs_l3$`Classical monocyte`$degs %>% 
                     data.frame() %>%
                     dplyr::filter(gene %in% c("ITGA4", "ITGB7")), 
                   by = c("Gene" = "gene")) %>%
  dplyr::mutate(label = paste0(Gene, "\np-value = ", round(pvalue, 2))) %>%
  ggplot(aes(x = Response, y = Expression, col = Response)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter_rast() +
  scale_color_manual(values = response_colors) +
  facet_wrap(~label, scales = "free_y", nrow = 2) +
  labs(y = bquote(log[2]~'(expression)')) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_text(face = "bold"))

# fig5F

fig5F <- fgsea_list$`Classical monocyte` %>%
  data.frame() %>%
  dplyr::select(-leadingEdge) %>%
  dplyr::mutate(Significance = ifelse(padj < 0.05, "Significant", "NS"),
                label = ifelse(pathway %in% c("KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
                                              "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
                                              "KEGG_CHEMOKINE_SIGNALING_PATHWAY",
                                              "KEGG_JAK_STAT_SIGNALING_PATHWAY",
                                              "KEGG_TGF_BETA_SIGNALING_PATHWAY",
                                              "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION"),
                               pathway,
                               NA),
                label = str_to_sentence(gsub("KEGG ", "", gsub("_", " ", label))),
                label = gsub("Jak stat", "JAK STAT", gsub("_", " ", label)),
                label = gsub("Tgf beta", "TGFB", gsub("_", " ", label))) %>%
  ggplot(aes(x = NES, y = -log10(pval))) +
  geom_vline(xintercept = 0, col = "#d3d3d3") +
  geom_hline(yintercept = 0, col = "#d3d3d3") +
  geom_point_rast(aes(col = Significance)) +
  geom_label_repel(aes(label = label),
                   box.padding = 0.5, 
                   max.overlaps = Inf,
                   nudge_x = -2) +
  geom_point_rast(data = . %>%
                    dplyr::filter(!is.na(label)),
                  col = "red") +
  xlim(-3, 3) +
  #geom_label_repel(aes(label = label)) +
  labs(title = "Differentally enriched KEGG pathways",
       subtitle = "Responders relative to non-responders",
       y = bquote('-'~log[10]~'(p-value)'),
       x = "NES") +
  scale_color_manual(values = c(`Significant` = "#000000", `NS` = "#d3d3d3"), drop = FALSE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom")

# fig5G

degs_classical_monocytes_nomsig <- scrnaseq_degs_l3$`Classical monocyte`$degs %>%
  data.frame() %>%
  dplyr::filter(pvalue<0.05)

cytokine_cytokine_receptor_interaction <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>%
  dplyr::filter(gs_name == "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
                gene_symbol %in% degs_classical_monocytes_nomsig$gene)

cc_nomsig_classical_monocytes <- rownames(scrnaseq_classical_monocyte_rld)[which(rownames(scrnaseq_classical_monocyte_rld) %in% cytokine_cytokine_receptor_interaction$gene_symbol)]

fig5G <- grid::grid.grabExpr(pheatmap::pheatmap(assay(scrnaseq_classical_monocyte_rld)[cc_nomsig_classical_monocytes, ],
                                                scale = "row",
                                                main = "Cytokine cytokine receptor interaction genes", 
                                                annotation_col = data.frame(Response = colData(scrnaseq_classical_monocyte_rld)$Response, 
                                                                            row.names = colData(scrnaseq_classical_monocyte_rld)$SampleID),
                                                cluster_cols = F,
                                                annotation_colors = list(Response = response_colors)))

# fig5H

fig5H <- do.call(rbind, lapply(names(scrnaseq_degs_l3), function(celltype){
  scrnaseq_degs_l3[[celltype]]$degs %>%
    data.frame() %>%
    dplyr::filter(gene %in% ccrl_mono$Partner) %>%
    dplyr::mutate(Celltype = celltype,
                  Group = factor(ifelse(stat > 0, "Responder", "Non-responder"), levels = c("Non-responder", "Responder")),
                  Significance = ifelse(pvalue<0.05, "Significant", "NS"))
})) %>%
  dplyr::mutate(gene = factor(gene, levels = rev(unique(ccrl_mono$Partner))),
                Celltype = factor(Celltype, levels = manual_l3_order$Celltype)) %>%
  ggplot(aes(x = Celltype, y = gene)) +
  geom_point(aes(size = -log10(pvalue), col = Group, fill = stat, alpha = Significance)) +
  scale_alpha_discrete(range = c(0.25, 1)) +
  scale_fill_gradient2(high = "#009E73",
                       low = "#D55E00",
                       mid = "white") +
  scale_color_manual(values = response_colors) +
  labs(title = "Cytokine-ligand/receptor differential expression") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y =  element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"))

# Compiled

fig5A <- ggarrange(plotlist = list(fig5A_left, fig5A_right), nrow = 1, ncol = 2, labels = c("A", ""), align = "hv", legend = "bottom", common.legend = T)
fig5DE <- ggarrange(fig5D_left, fig5D_right, fig5E, nrow = 1, ncol = 3, align = "hv", common.legend = T, legend = "bottom", labels = c("D", "", "E"))
fig5BCDE <- ggarrange(plotlist = list(fig5B, fig5C, fig5DE), nrow = 1, ncol = 3, labels = c("B", "C", "", ""), widths = c(1, 0.5, 1), align = "hv")
fig5FGH <- plot_grid(plotlist = list(fig5F, fig5G, fig5H), nrow = 1, ncol = 3, labels = c("F", "G", "H"), rel_widths = c(1.25, 2, 2), align = "hv")

fig5 <- plot_grid(plotlist = list(fig5A, fig5BCDE, fig5FGH), nrow = 3, align = "hv")

pdf(width = 17, height = 17, file = fig5_pdf)
print(fig5)
dev.off()

sessionInfo()