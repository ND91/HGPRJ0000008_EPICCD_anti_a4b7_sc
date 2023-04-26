rule fig_response_characteristics:
  input:
    sample_metadata_xlsx=expand("{sample_metadata}", sample_metadata=config['sample_metadata']),
  output:
    boxplot_response_characteristics_pdf="{basedir}/output/figures/boxplot_response_characteristics.pdf",
  threads: 
    1
  conda:
    "../envs/r-tidyverse.yaml"
  log:
    "{basedir}/output/figures/response_characteristics.log",
  benchmark:
    "{basedir}/output/figures/response_characteristics_benchmark.txt",
  resources:
    mem_mb=47000,
  shell:
    """
    Rscript workflow/scripts/figures/fig_response_characteristics.R "{input.sample_metadata_xlsx}" "{output.boxplot_response_characteristics_pdf}" &> "{log}"
    """

rule fig_scrnaseq_umap_pbmc:
  input:
    umap_pbmc_csv="{basedir}/output/scrnaseq/clustered/umap.csv",
    manual_l3_order_xlsx=config['manual_l3_order'],
  output:
    umap_pbmc_pdf="{basedir}/output/figures/scrnaseq_umap_pbmc.pdf",
  threads: 
    1
  conda:
    "../envs/r-tidyverse.yaml"
  log:
    "{basedir}/output/figures/scrnaseq_umap_pbmc.log",
  benchmark:
    "{basedir}/output/figures/scrnaseq_umap_pbmc_benchmark.txt",
  resources:
    mem_mb=47000,
  shell:
    """
    Rscript workflow/scripts/figures/fig_scrnaseq_umap_pbmc.R "{input.umap_pbmc_csv}" "{input.manual_l3_order_xlsx}" "{output.umap_pbmc_pdf}" &> "{log}"
    """

rule fig_scrnaseq_dotplot_markers_l3:
  input:
    seurat_rds="{basedir}/output/scrnaseq/clustered/clustered_SeuratObject.Rds",
    manual_l3_order_xlsx=config['manual_l3_order'],
  output:
    dotplot_markers_pdf="{basedir}/output/figures/scrnaseq_dotplot_markers_l3.pdf",
  threads:
    1
  conda:
    "../envs/r-tidyverse.yaml"
  log:
    "{basedir}/output/figures/scrnaseq_dotplot_markers_l3.log",
  benchmark:
    "{basedir}/output/figures/scrnaseq_dotplot_markers_l3_benchmark.txt",
  resources:
    mem_mb=47000,
  shell:
    """
    Rscript workflow/scripts/figures/fig_scrnaseq_dotplot_expression_markers_l3.R {input.seurat_rds} {input.manual_l3_order_xlsx} {output.dotplot_markers_pdf} &> {log}
    """

rule fig_masscytometry_umap_pbmc:
  input:
    umap_pbmc_csv="{basedir}/output/masscytometry/dimred/umap_ss.csv",
    manual_l3_order_xlsx=config['manual_l3_order'],
  output:
    umap_pbmc_pdf="{basedir}/output/figures/masscytometry_umap_pbmc.pdf",
  threads: 
    1
  conda:
    "../envs/r-tidyverse.yaml"
  log:
    "{basedir}/output/figures/masscytometry_umap_pbmc.log",
  benchmark:
    "{basedir}/output/figures/masscytometry_umap_pbmc_benchmark.txt",
  resources:
    mem_mb=47000,
  shell:
    """
    Rscript workflow/scripts/figures/fig_masscytometry_umap_pbmc.R "{input.umap_pbmc_csv}" "{input.manual_l3_order_xlsx}" "{output.umap_pbmc_pdf}" &> "{log}"
    """

rule fig_masscytometry_heatmap_markers_l3:
  input:
    masscytometry_annotated_rds="{basedir}/output/masscytometry/cell_metadata/sce_annotated.Rds",
    manual_l3_order_xlsx=config['manual_l3_order'],
  output:
    heatmap_markers_pdf="{basedir}/output/figures/masscytometry_heatmap_markers_l3.pdf",
  threads:
    1
  conda:
    "../envs/r-complexheatmap.yaml",
  log:
    "{basedir}/output/figures/masscytometry_heatmap_markers_l3.log",
  benchmark:
    "{basedir}/output/figures/masscytometry_heatmap_markers_l3_benchmark.txt",
  resources:
    mem_mb=47000,
  shell:
    """
    Rscript workflow/scripts/figures/fig_masscytometry_heatmap_expression_markers_l3.R "{input.masscytometry_annotated_rds}" "{input.manual_l3_order_xlsx}" "{output.heatmap_markers_pdf}" &> "{log}"
    """

rule fig_scrnaseq_pbmc_ITGA4_ITGB7:
  input:
    seurat_rds="{basedir}/output/scrnaseq/clustered/clustered_SeuratObject.Rds",
    manual_l3_order_xlsx=config['manual_l3_order'],
  output:
    umap_pbmc_ITGA4_ITGB7_pdf="{basedir}/output/figures/scrnaseq_umap_pbmc_ITGA4_ITGB7.pdf",
    boxplot_pbmc_ITGA4_ITGB7_pdf="{basedir}/output/figures/scrnaseq_boxplot_pbmc_ITGA4_ITGB7.pdf",
  threads: 
    1
  conda:
    "../envs/r-seurat.yaml"
  log:
    "{basedir}/output/figures/scrnaseq_pbmc_ITGA4_ITGB7.log",
  benchmark:
    "{basedir}/output/figures/scrnaseq_pbmc_ITGA4_ITGB7_benchmark.txt",
  resources:
    mem_mb=47000,
  shell:
    """
    Rscript workflow/scripts/figures/fig_scrnaseq_pbmc_ITGA4_ITGB7_expression.R "{input.seurat_rds}" "{input.manual_l3_order_xlsx}" "{output.umap_pbmc_ITGA4_ITGB7_pdf}" "{output.boxplot_pbmc_ITGA4_ITGB7_pdf}" &> "{log}"
    """

rule fig_masscytometry_pbmc_ITGA4_A4B7:
  input:
    sce_ss_rds="{basedir}/output/masscytometry/dimred/sce_dimred_ss.Rds",
    sce_full_rds="{basedir}/output/masscytometry/cell_metadata/sce_annotated.Rds",
    manual_l3_order_xlsx=config['manual_l3_order'],
  output:
    umap_pbmc_ITGA4_A4B7_pdf="{basedir}/output/figures/masscytometry_umap_pbmc_ITGA4_A4B7.pdf",
    boxplot_pbmc_ITGA4_A4B7_pdf="{basedir}/output/figures/masscytometry_boxplot_pbmc_ITGA4_A4B7.pdf",
  threads: 
    1
  conda:
    "../envs/r-singlecellexperiment.yaml"
  log:
    "{basedir}/output/figures/masscytometry_pbmc_ITGA4_A4B7.log",
  benchmark:
    "{basedir}/output/figures/masscytometry_pbmc_ITGA4_A4B7_benchmark.txt",
  resources:
    mem_mb=47000,
  shell:
    """
    Rscript workflow/scripts/figures/fig_masscytometry_pbmc_ITGA4_A4B7_expression.R "{input.sce_ss_rds}" "{input.sce_full_rds}" "{input.manual_l3_order_xlsx}" "{output.umap_pbmc_ITGA4_A4B7_pdf}" "{output.boxplot_pbmc_ITGA4_A4B7_pdf}" &> "{log}"
    """
    
rule fig_scrnaseq_masscytometry_scatterplot_abundance_l3rl0_colmanuall1:
  input:
    seurat_rds="{basedir}/output/scrnaseq/clustered/clustered_SeuratObject.Rds",
    masscytometry_annotated_rds="{basedir}/output/masscytometry/cell_metadata/sce_annotated.Rds",
  output:
    scatterplot_pdf="{basedir}/output/figures/scrnaseq_masscytometry_scatterplot_abundance_l3rl0_colmanuall1.pdf",
  threads:
    1
  conda:
    "../envs/r-seurat.yaml"
  log:
    "{basedir}/output/figures/scrnaseq_masscytometry_scatterplot_abundance_l3rl0_colmanuall1.log",
  benchmark:
    "{basedir}/output/figures/scrnaseq_masscytometry_scatterplot_abundance_l3rl0_colmanuall1.txt",
  resources:
    mem_mb=47000,
  shell:
    """
    Rscript workflow/scripts/figures/fig_scatterplot_masscytometryvscrnaseq_abundance_l3rl0.R "{input.seurat_rds}" "{input.masscytometry_annotated_rds}" "{output.scatterplot_pdf}" &> "{log}"
    """
 
rule fig_scrnaseq_boxplot_abundance_l3rl0_groupresponse:
  input:
    seurat_rds="{basedir}/output/scrnaseq/clustered/clustered_SeuratObject.Rds",
  output:
    boxplot_pdf="{basedir}/output/figures/scrnaseq_boxplot_abundance_l3rl0_groupresponse.pdf",
  threads:
    1
  conda:
    "../envs/r-tidyverse.yaml"
  log:
    "{basedir}/output/figures/scrnaseq_boxplot_abundance_l3rl0_groupresponse.log",
  benchmark:
    "{basedir}/output/figures/scrnaseq_boxplot_abundance_l3rl0_groupresponse_benchmark.txt",
  resources:
    mem_mb=47000,
  shell:
    """
    Rscript workflow/scripts/figures/fig_scrnaseq_boxplot_abundance_l3rl0_groupresponse.R "{input.seurat_rds}" "{output.boxplot_pdf}" &> "{log}"
    """

rule fig_masscytometry_boxplot_abundance_l3rl0_groupresponse:
  input:
    masscytometry_annotated_rds="{basedir}/output/masscytometry/cell_metadata/masscytometry_annotated.Rds",
  output:
    boxplot_pdf="{basedir}/output/figures/masscytometry_boxplot_abundance_l3rl0_groupresponse.pdf",
  threads:
    1
  conda:
    "../envs/r-tidyverse.yaml",
  log:
    "{basedir}/output/figures/masscytometry_boxplot_abundance_l3rl0_groupresponse.log",
  benchmark:
    "{basedir}/output/figures/masscytometry_boxplot_abundance_l3rl0_groupresponse_benchmark.txt",
  resources:
    mem_mb=47000,
  shell:
    """
    Rscript workflow/scripts/figures/fig_masscytometry_boxplot_abundance_l3rl0_groupresponse.R "{input.masscytometry_annotated_rds}" "{output.boxplot_pdf}" &> "{log}"
    """
    
rule fig_scrnaseq_masscytometry_scatterplot_differential_abundance_groupresponse:
  input:
    dacs_scrnaseq_csv="{basedir}/output/scrnaseq/da/dacs.csv",
    dacs_masscytometry_csv="{basedir}/output/masscytometry/da/dacs.csv",
  output:
    scatterplot_pdf="{basedir}/output/figures/scrnaseq_masscytometry_scatterplot_differential_abundance_groupresponse.pdf",
  threads:
    1
  conda:
    "../envs/r-seurat.yaml",
  log:
    "{basedir}/output/figures/scrnaseq_masscytometry_scatterplot_differential_abundance_groupresponse.log",
  benchmark:
    "{basedir}/output/figures/scrnaseq_masscytometry_scatterplot_differential_abundance_groupresponse_benchmark.txt",
  resources:
    mem_mb=47000,
  shell:
    """
    Rscript workflow/scripts/figures/fig_scrnaseq_masscytometry_scatterplot_differential_abundance_groupresponse.R "{input.dacs_scrnaseq_csv}" "{input.dacs_masscytometry_csv}" "{output.scatterplot_pdf}" &> "{log}"
    """

rule fig_scrnaseq_monocyte_umap:
  input:
    monocyte_umap_csv="{basedir}/output/scrnaseq/monocyte/umap.csv",
    manual_l3_order_xlsx=config['manual_l3_order'],
  output:
    monocyte_umap_pdf="{basedir}/output/figures/scrnaseq_monocyte_umap.pdf",
  threads:
    1
  conda:
    "../envs/r-tidyverse.yaml"
  log:
    "{basedir}/output/figures/scrnaseq_monocyte_umap.log",
  benchmark:
    "{basedir}/output/figures/scrnaseq_monocyte_umap_benchmark.txt",
  resources:
    mem_mb=47000,
  shell:
    """
    Rscript workflow/scripts/figures/fig_scrnaseq_arrowplot_degs.R "{input.degs_rds}" "{input.manual_l3_order_xlsx}" "{output.arrowplot_degs_pdf}" &> "{log}"
    """

rule fig_scrnaseq_arrowplot_degs:
  input:
    degs_rds="{basedir}/output/scrnaseq/de/degs.Rds",
    manual_l3_order_xlsx=config['manual_l3_order'],
  output:
    arrowplot_degs_pdf="{basedir}/output/figures/scrnaseq_arrowplot_degs.pdf",
  threads:
    1
  conda:
    "../envs/r-tidyverse.yaml"
  log:
    "{basedir}/output/figures/scrnaseq_arrowplot_degs.log",
  benchmark:
    "{basedir}/output/figures/scrnaseq_arrowplot_degs_benchmark.txt",
  resources:
    mem_mb=47000,
  shell:
    """
    Rscript workflow/scripts/figures/fig_scrnaseq_arrowplot_degs.R "{input.degs_rds}" "{input.manual_l3_order_xlsx}" "{output.arrowplot_degs_pdf}" &> "{log}"
    """
