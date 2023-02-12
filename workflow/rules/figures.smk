rule fig_scrnaseq_umap_pbmc:
  input:
    seurat_rds="{basedir}/output/scrnaseq/clustered/clustered_SeuratObject.Rds",
    manual_l3_order_xlsx=config['manual_l3_order'],
    renv_library="renv",
  output:
    umap_pbmc_pdf="{basedir}/output/figures/scrnaseq_umap_pbmc.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml"
  log:
    "{basedir}/output/figures/scrnaseq_umap_pbmc.log",
  benchmark:
    "{basedir}/output/figures/scrnaseq_umap_pbmc_benchmark.txt",
  resources:
    mem_mb=47000,
  shell:
    """
    Rscript workflow/scripts/figures/fig_scrnaseq_umap_pbmc.R "{input.seurat_rds}" "{input.manual_l3_order_xlsx}" "{output.umap_pbmc_pdf}" &> "{log}"
    """

rule fig_scrnaseq_dotplot_markers_l3:
  input:
    seurat_rds="{basedir}/output/scrnaseq/clustered/clustered_SeuratObject.Rds",
    manual_l3_order_xlsx=config['manual_l3_order'],
    renv_library="renv",
  output:
    dotplot_markers_pdf="{basedir}/output/figures/scrnaseq_dotplot_markers_l3.pdf",
  threads:
    1
  conda:
    "../envs/r.yaml"
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
    
rule fig_masscytometry_heatmap_markers_l3:
  input:
    masscytometry_annotated_rds="{basedir}/output/masscytometry/cell_metadata/masscytometry_annotated.Rds",
    manual_l3_order_xlsx=config['manual_l3_order'],
    renv_library="renv",
  output:
    dotplot_markers_pdf="{basedir}/output/figures/masscytometry_heatmap_markers_l3.pdf",
  threads:
    1
  conda:
    "../envs/masscytometry/r.yaml",
  log:
    "{basedir}/output/figures/masscytometry_heatmap_markers_l3.log",
  benchmark:
    "{basedir}/output/figures/masscytometry_heatmap_markers_l3_benchmark.txt",
  resources:
    mem_mb=47000,
  shell:
    """
    Rscript workflow/scripts/figures/fig_scrnaseq_dotplot_expression_markers_l3.R {input.seurat_rds} {input.manual_l3_order_xlsx} {output.dotplot_markers_pdf} &> {log}
    """
    
rule fig_scrnaseq_masscytometry_scatterplot_abundance_l3rl0_colmanuall1:
  input:
    seurat_rds="{basedir}/output/scrnaseq/clustered/clustered_SeuratObject.Rds",
    masscytometry_annotated_rds="{basedir}/output/masscytometry/cell_metadata/masscytometry_annotated.Rds",
    renv_library="renv",
  output:
    scatterplot_pdf="{basedir}/output/figures/scrnaseq_masscytometry_scatterplot_abundance_l3rl0_colmanuall1.pdf",
  threads:
    1
  conda:
    "../envs/r.yaml"
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
    renv_library="renv",
  output:
    boxplot_pdf="{basedir}/output/figures/scrnaseq_boxplot_abundance_l3rl0_groupresponse.pdf",
  threads:
    1
  conda:
    "../envs/r.yaml"
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
    renv_library="renv",
  output:
    boxplot_pdf="{basedir}/output/figures/masscytometry_boxplot_abundance_l3rl0_groupresponse.pdf",
  threads:
    1
  conda:
    "../envs/masscytometry/r.yaml",
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

rule fig_scrnaseq_arrowplot_degs:
  input:
    degs_rds="{basedir}/output/scrnaseq/de/degs.Rds",
    manual_l3_order_xlsx=config['manual_l3_order'],
    renv_library="renv",
  output:
    arrowplot_degs_pdf="{basedir}/output/figures/scrnaseq_arrowplot_degs.pdf",
  threads:
    1
  conda:
    "../envs/r.yaml"
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
