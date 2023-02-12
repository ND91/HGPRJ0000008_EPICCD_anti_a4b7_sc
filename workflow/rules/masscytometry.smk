rule masscytometry_merge_annotations:
  input:
    masscytometry_rds="{basedir}/resources/masscytometry/cytof.Rds",
    files_pbmc_masscytometry_xlsx=config['masscytometry_files'],
    sample_metadata_xlsx=config['sample_metadata'],
    feature_metadata_xlsx=config['masscytometry_markers'],
    renv_library="renv",
  output:
    sce_annotated_rds="{basedir}/output/masscytometry/cell_metadata/sce_annotated.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml"
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

rule masscytometry_clustering:
  input:
    sce_annotated_rds="{basedir}/output/masscytometry/cell_metadata/sce_annotated.Rds",
    renv_library="renv",
  output:
    sce_clustered_rds="{basedir}/output/masscytometry/clustered/sce_clustered.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml"
  log:
    "{basedir}/output/masscytometry/clustered/masscytometry_clustering.log",
  benchmark:
    "{basedir}/output/masscytometry/clustered/masscytometry_clustering_benchmark.txt",
  resources:
    mem_mb=47000,
  shell:
    """
    Rscript "workflow/scripts/masscytometry/clustering.R" "{input.sce_annotated_rds}" "{output.sce_clustered_rds}" &> "{log}"
    """

# Analyses

rule masscytometry_differential_abundance:
  input:
    sce_annotated_rds="{basedir}/output/masscytometry/cell_metadata/sce_annotated.Rds",
    renv_library="renv",
  output:
    dacs_csv="{basedir}/output/masscytometry/da/dacs.csv",
  threads: 
    1
  conda:
    "../envs/r.yaml"
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
