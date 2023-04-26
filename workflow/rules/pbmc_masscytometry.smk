rule masscytometry_merge_annotations:
  input:
    masscytometry_rds="{basedir}/resources/masscytometry/cytof.Rds",
    files_pbmc_masscytometry_xlsx=config['masscytometry_files'],
    sample_metadata_xlsx=config['sample_metadata'],
    feature_metadata_xlsx=config['masscytometry_markers'],
  output:
    sce_annotated_rds="{basedir}/output/masscytometry/cell_metadata/sce_annotated.Rds",
  threads: 
    1
  conda:
    "../envs/r-singlecellexperiment.yaml"
  log:
    "{basedir}/output/masscytometry/cell_metadata/masscytometry_merge_annotations.log",
  benchmark:
    "{basedir}/output/masscytometry/cell_metadata/masscytometry_merge_annotations_benchmark.txt",
  resources:
    mem_mb=47000,
  shell:
    """
    Rscript "workflow/scripts/masscytometry/merge_annotations.R" "{input.masscytometry_rds}" "{input.files_pbmc_masscytometry_xlsx}" "{input.sample_metadata_xlsx}" "{input.feature_metadata_xlsx}" "{output.sce_annotated_rds}" &> "{log}"
    """

rule masscytometry_dimred:
  input:
    sce_annotated_rds="{basedir}/output/masscytometry/cell_metadata/sce_annotated.Rds",
  output:
    sce_dimred_ss_rds="{basedir}/output/masscytometry/dimred/sce_dimred_ss.Rds",
    umap_ss_csv="{basedir}/output/masscytometry/dimred/umap_ss.csv",
  threads: 
    1
  conda:
    "../envs/r-scater.yaml"
  log:
    "{basedir}/output/masscytometry/dimred/masscytometry_dimred.log",
  benchmark:
    "{basedir}/output/masscytometry/dimred/masscytometry_dimred_benchmark.txt",
  resources:
    mem_mb=47000,
  params:
    nsubsample=16000, #This value represents the number of cells to subsample from the total experiment. 16000 is approximately equal to what we acquired in scRNAseq.
  shell:
    """
    Rscript "workflow/scripts/masscytometry/dimred.R" "{input.sce_annotated_rds}" "{params.nsubsample}" "{output.sce_dimred_ss_rds}" "{output.umap_ss_csv}" &> "{log}"
    """

# Analyses

rule masscytometry_differential_abundance:
  input:
    sce_annotated_rds="{basedir}/output/masscytometry/cell_metadata/sce_annotated.Rds",
  output:
    dacs_csv="{basedir}/output/masscytometry/da/dacs.csv",
  threads: 
    1
  conda:
    "../envs/r-speckle.yaml"
  log:
    "{basedir}/output/masscytometry/da/dacs.log",
  benchmark:
    "{basedir}/output/masscytometry/da/dacs_benchmark.txt",
  resources:
    mem_mb=47000,
  shell:
    """
    Rscript workflow/scripts/masscytometry/differential_abundance.R "{input.sce_annotated_rds}" "{output.dacs_csv}" &> "{log}"
    """
