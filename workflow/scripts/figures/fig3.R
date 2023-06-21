#!/usr/bin/env Rscript

# The goal of this script is to create figure 3 for the manuscript.

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

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9) {
  stop(paste0("Script needs 9 arguments. Current input is:", args))
}

seurat_rds <- args[1]
scrnaseq_dacs_csv <- args[2]
sce_rds <- args[3]
masscytometry_dacs_csv <- args[4]
scrnaseq_degs_l3_rds <- args[5]
scrnaseq_fgsea_list_rds <- args[6]
manual_l1_order_xlsx <- args[7]
manual_l3_order_xlsx <- args[8]
fig3_pdf <- args[9]

scrnaseq_seuratObject <- readRDS(seurat_rds)
scrnaseq_dacs <- read.csv(scrnaseq_dacs_csv)
sce_masscytometry <- readRDS(sce_rds)
masscytometry_dacs <- read.csv(masscytometry_dacs_csv)
scrnaseq_degs_l3 <- readRDS(scrnaseq_degs_l3_rds)
scrnaseq_fgsea_list <- readRDS(scrnaseq_fgsea_list_rds)

manual_l1_order <- readxl::read_excel(manual_l1_order_xlsx, col_names = T)  %>%
  dplyr::mutate(number = as.numeric(factor(Celltype, levels = Celltype)),
                celltype_number = paste0(number, ". ", Celltype),
                celltype_number = factor(celltype_number, levels = celltype_number))

manual_l1_colors <- manual_l1_order$Color
names(manual_l1_colors) <- manual_l1_order$Celltype

manual_l3_order <- readxl::read_excel(manual_l3_order_xlsx, col_names = T)  %>%
  dplyr::mutate(number = as.numeric(factor(Celltype, levels = Celltype)),
                celltype_number = paste0(number, ". ", Celltype),
                celltype_number = factor(celltype_number, levels = celltype_number))

manual_l3_number_colors <- manual_l3_order$Color
names(manual_l3_number_colors) <- manual_l3_order$celltype_number

response_colors <- c(`Responder` = "#009E73",
                     `Non-responder` = "#D55E00")

hierarchy <- unique(rbind(scrnaseq_seuratObject@meta.data[,c("manual_l1", "manual_l2", "manual_l3")], 
                          data.frame(colData(sce_masscytometry))[,c("manual_l1", "manual_l2", "manual_l3")]))

# Fig3A

fig3A <- scrnaseq_seuratObject@meta.data %>%
  dplyr::mutate(manual_l3 = as.factor(manual_l3)) %>%
  dplyr::group_by(SampleID, Response, manual_l3, .drop = F) %>%
  dplyr::summarize(Nl3sample = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(SampleID, Response) %>%
  dplyr::mutate(Ncells = sum(Nl3sample),
                Ncellprop = Nl3sample/Ncells,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop)) %>%
  dplyr::mutate(manual_l3 = factor(manual_l3, levels = manual_l3_order$Celltype)) %>%
  ggplot(aes(x = manual_l3, y = Ncellprop*100, col = Response)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
  labs(title = "scRNAseq",
       y = "%PBMCs") +
  scale_color_manual(values = response_colors) +
  scale_x_discrete(drop=FALSE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Fig3B

fig3B <- colData(sce_masscytometry) %>%
  data.frame() %>%
  dplyr::mutate(manual_l3 = as.factor(manual_l3)) %>%
  dplyr::group_by(Sample_ID, Response, manual_l3, .drop = F) %>%
  dplyr::summarize(Nl3sample = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Sample_ID, Response) %>%
  dplyr::mutate(Ncells = sum(Nl3sample),
                Ncellprop = Nl3sample/Ncells,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop)) %>%
  dplyr::mutate(manual_l3 = factor(manual_l3, levels = manual_l3_order$Celltype)) %>%
  ggplot(aes(x = manual_l3, y = Ncellprop*100, col = Response)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
  labs(title = "CyTOF",
       y = "%PBMCs") +
  scale_color_manual(values = response_colors) +
  scale_x_discrete(drop=FALSE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Fig3C

scrnaseq_dacs_df <- scrnaseq_dacs %>%
  dplyr::select(BaselineProp.clusters, PropRatio, Tstatistic, P.Value, FDR) %>%
  dplyr::rename(Celltype = BaselineProp.clusters,
                propratio_scrnaseq = PropRatio,
                t_scrnaseq = Tstatistic,
                pvalue_scrnaseq = P.Value,
                fdr_scrnaseq = FDR) %>%
  dplyr::mutate(log2_propratio_scrnaseq = -log2(propratio_scrnaseq))
masscytometry_dacs_df <- masscytometry_dacs %>%
  dplyr::select(BaselineProp.clusters, PropRatio, Tstatistic, P.Value, FDR) %>%
  dplyr::rename(Celltype = BaselineProp.clusters,
                propratio_cytof = PropRatio,
                t_cytof = Tstatistic,
                pvalue_cytof = P.Value,
                fdr_cytof = FDR) %>%
  dplyr::mutate(log2_propratio_cytof = -log2(propratio_cytof))

dacs_overlapping <- scrnaseq_dacs_df %>%
  dplyr::full_join(masscytometry_dacs_df, by = "Celltype") %>%
  dplyr::left_join(hierarchy, by = c("Celltype" = "manual_l3")) %>%
  dplyr::mutate(propratio_cytof = ifelse(is.na(propratio_cytof), 0, propratio_cytof),
                t_cytof = ifelse(is.na(t_cytof), 0, t_cytof),
                log2_propratio_cytof = ifelse(is.na(log2_propratio_cytof) | log2_propratio_cytof == Inf, 0, log2_propratio_cytof), 
                pvalue_cytof = ifelse(is.na(pvalue_cytof), 1, pvalue_cytof),
                fdr_cytof = ifelse(is.na(fdr_cytof), 1, fdr_cytof), 
                propratio_scrnaseq = ifelse(is.na(propratio_scrnaseq), 0, propratio_scrnaseq),
                t_scrnaseq = ifelse(is.na(t_scrnaseq), 0, t_scrnaseq),
                log2_propratio_scrnaseq = ifelse(is.na(log2_propratio_scrnaseq) | log2_propratio_scrnaseq == Inf, 0, log2_propratio_scrnaseq),
                pvalue_scrnaseq = ifelse(is.na(pvalue_scrnaseq), 1, pvalue_scrnaseq),
                fdr_scrnaseq = ifelse(is.na(fdr_scrnaseq), 1, fdr_scrnaseq),
                mean_log2_propratio_cytof_scrnaseq = (log2_propratio_cytof+log2_propratio_scrnaseq)/2,
                Significance = ifelse(pvalue_scrnaseq < 0.05 | pvalue_cytof < 0.05, "Significant", "NS"))
  

fig3C <- dacs_overlapping %>%
  dplyr::mutate(Lineage = factor(manual_l1, levels = manual_l1_order$Celltype)) %>%
  ggplot(aes(x = log2_propratio_scrnaseq, y = log2_propratio_cytof, col = Lineage), size = 5) +
  geom_vline(xintercept = 0, col = "#d3d3d3") +
  geom_hline(yintercept = 0, col = "#d3d3d3") +
  geom_point(show.legend = T) +
  geom_label_repel(aes(label = Celltype), show.legend = F) +
  labs(title = "Differential abundance R vs NR",
       subtitle = bquote(log[2]~"(R/NR)"),
       x = "scRNAseq",
       y = "CyTOF") +
  # scale_alpha_manual(values = c("NS" = 0.5,
  #                               "Significant" = 1)) +
  #scale_color_manual(values = manual_l1_colors) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom")

# Fig3D
tbnkmyeloid_overlapping_split <- dacs_overlapping %>% 
  dplyr::filter(manual_l1 %in% c("T", "NK", "B", "Myeloid")) %>%
  split(., .$manual_l1)

lapply(tbnkmyeloid_overlapping_split, function(manual_l1){
  t.test(x = manual_l1$mean_log2_propratio_cytof_scrnaseq, mu = 0, alternative = "two.sided")
})

fig3D <- dacs_overlapping %>%
  dplyr::filter(manual_l1 %in% c("T", "NK", "B", "Myeloid")) %>%
  dplyr::mutate(manual_l1 = factor(manual_l1, levels = c("T", "NK", "B", "Myeloid"))) %>%
  ggplot(aes(x = manual_l1 , y = mean_log2_propratio_cytof_scrnaseq)) +
  geom_hline(yintercept = 0, col = "#d3d3d3") +
  geom_boxplot(outlier.shape = NA, alpha = 0.25) +
  geom_jitter_rast() +
  labs(title = 'Mean ratio R/NR',
       y = bquote('mean'~log[2]~'(R/NR)')) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),
        axis.title.x = element_blank())

# Fig3E

fig3E_left <- colData(sce_masscytometry) %>%
  data.frame() %>%
  dplyr::mutate(manual_l3 = as.factor(manual_l3)) %>%
  dplyr::group_by(Sample_ID, Response, manual_l3, .drop = F) %>%
  dplyr::summarize(Nl3sample = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Sample_ID, Response) %>%
  dplyr::mutate(Ncells = sum(Nl3sample),
                Ncellprop = Nl3sample/Ncells,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop)) %>%
  dplyr::filter(manual_l3 == "CD8 TCM") %>%
  ggplot(aes(x = Response, y = Ncellprop*100, col = Response)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, show.legend = F) +
  geom_jitter(alpha = 0.5, show.legend = F) +
  labs(title = "CD8 TCM (CyTOF)",
       subtitle = paste0("p-value = ", round(masscytometry_dacs$P.Value[masscytometry_dacs$BaselineProp.clusters == "CD8 TCM"], 3)),
       y = "%PBMCs") +
  scale_color_manual(values = response_colors) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),)

fig3E_right <- scrnaseq_seuratObject@meta.data %>%
  dplyr::mutate(manual_l3 = as.factor(manual_l3)) %>%
  dplyr::group_by(SampleID, Response, manual_l3, .drop = F) %>%
  dplyr::summarize(Nl3sample = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(SampleID, Response) %>%
  dplyr::mutate(Ncells = sum(Nl3sample),
                Ncellprop = Nl3sample/Ncells,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop)) %>%
  dplyr::filter(manual_l3 == "CD8 TCM") %>%
  ggplot(aes(x = Response, y = Ncellprop*100, col = Response)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, show.legend = F) +
  geom_jitter(alpha = 0.5, show.legend = F) +
  labs(title = "CD8 TCM (scRNAseq)",
       subtitle = paste0("p-value = ", round(scrnaseq_dacs$P.Value[scrnaseq_dacs$BaselineProp.clusters == "CD8 TCM"], 3)),
       y = "%PBMCs") +
  scale_color_manual(values = response_colors) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Fig3F

fig3F <- scrnaseq_seuratObject@meta.data %>%
  dplyr::mutate(manual_l3 = as.factor(manual_l3)) %>%
  dplyr::group_by(SampleID, Response, manual_l3, .drop = F) %>%
  dplyr::summarize(Nl3sample = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(SampleID, Response) %>%
  dplyr::mutate(Ncells = sum(Nl3sample),
                Ncellprop = Nl3sample/Ncells,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop)) %>%
  dplyr::filter(manual_l3 == "MAIT") %>%
  ggplot(aes(x = Response, y = Ncellprop*100, col = Response)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, show.legend = F) +
  geom_jitter(alpha = 0.5, show.legend = F) +
  labs(title = "MAIT",
       subtitle = paste0("p-value = ", round(scrnaseq_dacs$P.Value[scrnaseq_dacs$BaselineProp.clusters == "MAIT"], 3)),
       y = "%PBMCs") +
  scale_color_manual(values = response_colors) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Fig3G

fig3G <- scrnaseq_seuratObject@meta.data %>%
  dplyr::mutate(manual_l3 = as.factor(manual_l3)) %>%
  dplyr::group_by(SampleID, Response, manual_l3, .drop = F) %>%
  dplyr::summarize(Nl3sample = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(SampleID, Response) %>%
  dplyr::mutate(Ncells = sum(Nl3sample),
                Ncellprop = Nl3sample/Ncells,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop)) %>%
  dplyr::filter(manual_l3 == "PDC") %>%
  ggplot(aes(x = Response, y = Ncellprop*100, col = Response)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, show.legend = F) +
  geom_jitter(alpha = 0.5, show.legend = F) +
  labs(title = "PDC",
       subtitle = paste0("p-value = ", round(scrnaseq_dacs$P.Value[scrnaseq_dacs$BaselineProp.clusters == "PDC"], 3)),
       y = "%PBMCs") +
  scale_color_manual(values = response_colors) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Fig3H

tcells <- c("CD4 naive", "CD4 TCM", "CD4 TEM", "CD4 CTL", "CD4 Treg", "CD4 proliferating", "CD8 naive", "CD8 TCM", "CD8 TEM", "MAIT", "NKT", "GDT", "DNT")

tcell_degs <- do.call(rbind, lapply(tcells, function(tcell){
  scrnaseq_degs_l3[[tcell]]$degs %>%
    data.frame() %>%
    dplyr::slice_head(n = 2) %>%
    dplyr::mutate(Celltype = tcell)
}))

tcell_goi <- c(unique(tcell_degs$gene), "ITGA4", "ITGB7")

tcell_de_goi <- do.call(rbind, lapply(tcells, function(tcell){
  scrnaseq_degs_l3[[tcell]]$degs %>%
    data.frame() %>%
    dplyr::filter(gene %in% tcell_goi) %>%
    dplyr::mutate(Celltype = tcell)
}))

fig3H <- tcell_de_goi %>%
  dplyr::mutate(Significance = ifelse(pvalue<0.05, "Significant", "NS"),
                Group = ifelse(stat<0, "Non-responder", "Responder"),
                gene = factor(gene, levels = tcell_goi),
                Celltype = factor(Celltype, levels = rev(tcells))) %>%
  ggplot(aes(x = gene, y = Celltype)) +
  geom_point(aes(size = -log10(pvalue), col = Group, fill = stat, alpha = Significance)) +
  scale_alpha_discrete(range = c(0.25, 1)) +
  #scale_shape_manual(values = c(Responder = 24, 
  #                              `Non-responder` = 25)) +
  scale_fill_gradient2(high = "#009E73",
                       low = "#D55E00",
                       mid = "white") +
  scale_color_manual(values = response_colors) +
  labs(title = "T cell DEGs") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y =  element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"))

# Fig3I

tcell_depws <- do.call(rbind, lapply(tcells, function(tcell){
  scrnaseq_fgsea_list[[tcell]] %>%
    data.frame() %>%
    dplyr::filter(padj<0.05)
}))

tcell_depwoi <- names(which(table(tcell_depws$pathway)>2))

tcell_de_pwoi <- do.call(rbind, lapply(tcells, function(tcell){
  scrnaseq_fgsea_list[[tcell]] %>%
    data.frame() %>%
    dplyr::select(-8) %>%
    dplyr::filter(pathway %in% tcell_depwoi) %>%
    dplyr::mutate(Celltype = tcell)
}))

fig3I <- tcell_de_pwoi %>%
    dplyr::mutate(Group = factor(ifelse(NES > 0, "Responder", "Non-responder"), levels = c("Non-responder", "Responder")),
                  Significance = ifelse(pval<0.05, "Significant", "NS")) %>%
  dplyr::mutate(pathway = factor(pathway, levels = rev(unique(tcell_depwoi))),
                Celltype = factor(Celltype, levels = manual_l3_order$Celltype)) %>%
  ggplot(aes(x = pathway, y = Celltype)) +
  geom_point(aes(size = -log10(pval), col = Group, fill = NES, alpha = Significance)) +
  scale_alpha_discrete(range = c(0.25, 1)) +
  #scale_shape_manual(values = c(Responder = 24, 
  #                              `Non-responder` = 25)) +
  scale_fill_gradient2(high = "#009E73",
                       low = "#D55E00",
                       mid = "white") +
  scale_color_manual(values = response_colors) +
  labs(title = "T cell KEGG GSEAs") +
  theme_bw() +
  guides(size = guide_legend(override.aes = list(shape=17))) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y =  element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"))

# Compiled

fig3AB <- ggarrange(fig3A, fig3B, nrow = 2, ncol = 1, labels = c("A", "B"), common.legend = T, legend = "bottom", align = "hv", widths = c(3, 0.5))
fig3CD <- ggarrange(fig3C, fig3D, nrow = 1, ncol = 2, labels = c("C", "D"), align = "hv", widths = c(1, 0.5))

fig3ABCD <- ggarrange(fig3AB, fig3CD, nrow = 1, ncol = 2, align = "hv", widths = c(1, 0.75))
fig3EFG <- ggarrange(fig3E_left, fig3E_right, fig3F, fig3G, nrow = 2, ncol = 2, align = "hv", labels = c("E", "", "F", "G"))
fig3EFGH <- ggarrange(fig3EFG, fig3H, nrow = 1, ncol = 2, align = "hv", labels = c("", "H"), widths = c(0.5, 1))

fig3 <- ggarrange(fig3ABCD, fig3EFGH, nrow = 2, ncol = 1, align = "hv", heights = c(1, 0.75))

pdf(width = 16.5, height = 15, file = fig3_pdf)
print(fig3)
dev.off()

sessionInfo()