rule gse134809_celltype_annotation:
  input:
    seurat_rds="{basedir}/resources/gse134809/gse134809.Rds",
    pbmc_rds="{basedir}/resources/scrnaseq/reference_data/pbmc_multimodal.Rds",
  output:
    gse134809_celltype_annotations_csv="{basedir}/output/gse134809/annotated/gse134809_celltype_annotations.csv",
  threads: 
    1
  conda:
    "../envs/r-seurat.yaml",
  message:
    "--- GSE134809: Manual curation ---",
  log:
    "{basedir}/output/gse134809/annotated/gse134809_celltype_annotation.log",
  benchmark:
    "{basedir}/output/gse134809/annotated/gse134809_celltype_annotation_benchmark.txt",
  resources:
    mem_mb=47000,
  shell:
    """
    Rscript -e "rmarkdown::render('workflow/scripts/gse134809/gse134809_celltype_annotation.Rmd', knit_root_dir = "{basedir}")" &> "{log}"
    """
    
rule gse134809_integrate_annotations:
  input:
    seurat_rds="{basedir}/resources/gse134809/gse134809.Rds",
    seurat_annotated_csv="{basedir}/resources/gse134809/gse134809_celltypes.csv",
    sample_metadata_xlsx="{basedir}/resources/gse134809/1-s2.0-S0092867419308967-mmc2.xlsx"
  output:
    seurat_annotated_rds="{basedir}/output/gse134809/annotated/gse134809_annotated_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r-seurat.yaml",
  log:
    "{basedir}/output/gse134809/annotated/gse134809_integrate_annotations.log",
  message:
    "--- GSE134809: Integrating annotations ---",
  benchmark:
    "{basedir}/output/gse134809/annotated/gse134809_integrate_annotations_benchmark.txt",
  resources:
    mem_mb=47000,
  shell:
    """
    Rscript workflow/scripts/gse134809/gse134809_merge_annotations.R "{input.seurat_rds}" "{input.seurat_annotated_csv}" "{input.sample_metadata_xlsx}" "{output.seurat_annotated_rds}" &> "{log}"
    """

# Analyses

rule gse134809_differential_abundance:
  input:
    seurat_annotated_rds="{basedir}/output/gse134809/annotated/gse134809_annotated_SeuratObject.Rds",
  output:
    dacs_involvedvuninvolved_l1rl0_csv="{basedir}/output/gse134809/da/dacs_involvedvuninvolved_l1rl0.csv",
    dacs_involvedvuninvolved_immune_l3rl1_csv="{basedir}/output/gse134809/da/dacs_involvedvuninvolved_immune_l3rl1.csv",
  threads: 
    1
  conda:
    "../envs/r-speckle.yaml",
  log:
    "{basedir}/output/gse134809/da/dacs.log",
  message:
    "--- GSE134809: Differential abundance analyses ---",
  benchmark:
    "{basedir}/output/gse134809/da/dacs_benchmark.txt",
  resources:
    mem_mb=47000,
  shell:
    """
    Rscript workflow/scripts/gse134809/gse134809_differential_abundance.R "{input.seurat_annotated_rds}" "{output.dacs_involvedvuninvolved_l1rl0_csv}" "{output.dacs_involvedvuninvolved_immune_l3rl1_csv}" &> "{log}"
    """
    
rule gse134809_differential_expression:
  input:
    seurat_annotated_rds="{basedir}/output/gse134809/annotated/gse134809_annotated_SeuratObject.Rds",
    functions_r="{basedir}/workflow/scripts/functions.R",
  output:
    degs_list_rds="{basedir}/output/gse134809/de/degs_manual_l3_list.Rds",
  threads: 
    1
  conda:
    "../envs/r-deseq2.yaml",
  log:
    "{basedir}/output/gse134809/de/gse134809_differential_expression.log",
  message:
    "--- GSE134809: Differential expression analyses at l3 level ---",
  benchmark:
    "{basedir}/output/gse134809/de/gse134809_differential_expression_benchmark.txt",
  resources:
    mem_mb=47000,
  shell:
    """
    Rscript workflow/scripts/gse134809/gse134809_pseudobulk_differential_expression.R "{input.seurat_annotated_rds}" "{input.functions_r}" "{output.degs_list_rds}" &> "{log}"
    """

