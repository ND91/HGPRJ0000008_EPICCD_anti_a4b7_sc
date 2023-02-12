rule scrnaseq_copy_fastqs:
  input:
    fastqdir=lambda w: (scrnaseq_fastq_df[scrnaseq_fastq_df.file == w.scrnaseq_fastq_fn].old_path).tolist(),
  output:
    directory("resources/scrnaseq/fastq/{scrnaseq_fastq_fn}"),
  conda:
    "../envs/rsync.yaml"
  message:
    "--- Copying {input.fastqdir} to scratch ---"
  shell:
    """
    if [ ! -d "{input.fastqdir}" ]; then
      echo "{input.fastqdir} does not exist!"
      exit 1
    fi
    
    rsync -ar "{input.fastqdir}/" "{output}"
    """

rule scrnaseq_copy_cellranger_reference:
  input:
    config['scrnaseq_reference_dir']
  output:
    directory("resources/scrnaseq/reference_genome")
  conda:
    "../envs/rsync.yaml"
  message:
    "--- Copying reference genome to scratch ---"
  shell:
    """
    rsync -ar "{input}/" "{output}"
    """

rule scrnaseq_cellranger_count:
  input:
    fastqdir=lambda wildcards: 'resources/scrnaseq/fastq/' + lookup_fbc.loc[lookup_fbc['Run'] == wildcards.run_fbc, "Run_well"].astype(str),
    fastqfbcdir=lambda wildcards: 'resources/scrnaseq/fastq/' + lookup_fbc.loc[lookup_fbc['Run'] == wildcards.run_fbc, "Run_well"].astype(str) + '-HT',
    libraries="config/samples/libraries_pbmc_scrnaseq.csv",
    features="config/samples/features_pbmc_scrnaseq.csv",
    reference_genome="resources/scrnaseq/reference_genome",
  output:
    cellrangerdir=touch(directory("output/scrnaseq/cellranger/{run_fbc}")),
    filtered_matrix=touch(directory("output/scrnaseq/cellranger/{run_fbc}/outs/filtered_feature_bc_matrix")),
    raw_matrix=touch(directory("output/scrnaseq/cellranger/{run_fbc}/outs/raw_feature_bc_matrix")),
  log:
    "output/scrnaseq/cellranger/{run_fbc}_cellranger.log",
  threads: 
    8
  benchmark:
    "output/scrnaseq/cellranger/{run_fbc}_count.benchmark.txt"
  params:
    runid="{run_fbc}",
    cellranger_basepath=config["cellranger_basepath"],
  resources:
    mem_mb=64000,
  shell:
    """
    rm -rf {output.cellrangerdir}
    mem_gb=$(echo {resources.mem_mb}/1024|bc)
    mkdir -p {output.cellrangerdir}
    cd "output/scrnaseq/cellranger"
    #find . -type d -empty -delete
    {params.cellranger_basepath}/cellranger count \
      --id="{params.runid}" \
      --libraries=$OLDPWD/{input.libraries} \
      --feature-ref=$OLDPWD/{input.features} \
      --transcriptome=$OLDPWD/{input.reference_genome} \
      --nosecondary \
      --disable-ui \
      --localcores={threads} \
      --localmem=$mem_gb &> $OLDPWD/{log}
    """

rule scrnaseq_cellranger_count_mv:
  input:
    cellrangertmpdir="output/scrnaseq/cellranger/{run_fbc}",
  output:
    cellrangerfinaldir=directory("{basedir}/output/scrnaseq/cellranger/{run_fbc}"),
  log:
    "{basedir}/output/scrnaseq/{run_fbc}_cellranger_mv.log",
  threads: 
    1
  benchmark:
    "{basedir}/output/scrnaseq/{run_fbc}_cellranger_mv.benchmark.txt"
  resources:
    mem_mb=10000,
  shell:
    """
    mv "{input.cellrangertmpdir}" &> "{output.cellrangerfinaldir}" &> "{log}"
    """

rule scrnaseq_download_pbmc_reference_h5seurat:
  input:
    config['scrnaseq_reference_dir'],
  output:
    pbmc_reference_h5seurat="{basedir}/resources/scrnaseq/reference_data/pbmc_multimodal.h5seurat",
  threads: 
    1
  conda:
    "../envs/wget.yaml",
  message:
    "--- Downloading PBMC reference dataset in h5seurat format from Hao et al. 2020 ---"
  log:
    "{basedir}/resources/scrnaseq/reference_data/download_pbmc_reference_h5seurat.log",
  benchmark:
    "{basedir}/resources/scrnaseq/reference_data/download_pbmc_reference_h5seurat_benchmark.txt",
  shell:
    """
    wget -O {output.pbmc_reference_h5seurat} 'https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat' &> "{log}"
    """

rule scrnaseq_preparation_pbmc_reference:
  input:
    pbmc_reference_h5seurat="{basedir}/resources/scrnaseq/reference_data/pbmc_multimodal.h5seurat",
    renv_library="renv",
  output:
    pbmc_reference_rds="{basedir}/resources/scrnaseq/reference_data/pbmc_multimodal.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml"
  log:
    "{basedir}/resources/scrnaseq/reference_data/preparation_pbmc_reference.log",
  benchmark:
    "{basedir}/resources/scrnaseq/reference_data/preparation_pbmc_reference_benchmark.txt",
  resources:
    mem_mb=40000,
  shell:
    """
    Rscript workflow/scripts/scrnaseq/pbmc_multimodal_preparation.R "{input.pbmc_reference_h5seurat}" "{output.pbmc_reference_rds}" &> "{log}"
    """

rule scrnaseq_normalization:
  input:
    count_mtx="output/scrnaseq/cellranger/{run_fbc}/outs/raw_feature_bc_matrix",
    renv_library="renv",
  output:
    seurat_rds="{basedir}/output/scrnaseq/normalized/{run_fbc}_normalized_SeuratObject.Rds",
  threads:
    1
  conda:
    "../envs/scrnaseq/r.yaml"
  log:
    "{basedir}/output/scrnaseq/normalized/{run_fbc}.log",
  benchmark:
    "{basedir}/output/scrnaseq/normalized/{run_fbc}_benchmark.txt"
  resources:
    mem_mb=16000,
  params:
    runid="{run_fbc}",
    lowerbound_numis=600,
  shell:
    """
    Rscript workflow/scripts/scrnaseq/normalization.R "{input.count_mtx}" "{output.seurat_rds}" "{params.runid}" "{params.lowerbound_numis}" &> "{log}"
    """

rule scrnaseq_celltype_annotation:
  input:
    seurat_rds="{basedir}/output/scrnaseq/normalized/{run_fbc}_normalized_SeuratObject.Rds",
    pbmc_rds="{basedir}/resources/scrnaseq/reference_data/pbmc_multimodal.Rds",
    renv_library="renv",
  output:
    seurat_annotated_csv="{basedir}/output/scrnaseq/cell_metadata/celltype_annotation/{run_fbc}_celltype_annotated.csv",
  threads: 
    1
  conda:
    "../envs/scrnaseq/r.yaml",
  log:
    "{basedir}/output/scrnaseq/cell_metadata/celltype_annotation/{run_fbc}_celltype_annotation.log",
  benchmark:
    "{basedir}/output/scrnaseq/cell_metadata/celltype_annotation/{run_fbc}_celltype_annotation_benchmark.txt",
  resources:
    mem_mb=47000,
  shell:
    """
    Rscript -e "rmarkdown::render('workflow/scripts/scrnaseq/celltype_annotation.Rmd', knit_root_dir = "${{PWD}}")" &> {log}
    """

rule scrnaseq_sample_demultiplex:
  input:
    seurat_rds="{basedir}/output/scrnaseq/normalized/{run_fbc}_normalized_SeuratObject.Rds",
    renv_library="renv",
  output:
    seurat_demultiplex_csv="{basedir}/output/scrnaseq/cell_metadata/sample_demultiplex/{run_fbc}_demultiplex.csv",
  threads: 
    1
  conda:
    "../envs/scrnaseq/r.yaml",
  log:
    "{basedir}/output/scrnaseq/cell_metadata/sample_demultiplex/{run_fbc}_demultiplex_annotation.log",
  benchmark:
    "{basedir}/output/scrnaseq/cell_metadata/sample_demultiplex/{run_fbc}_demultiplex_benchmark.txt",
  resources:
    mem_mb=47000,
  shell:
    """
    Rscript -e "rmarkdown::render('workflow/scripts/scrnaseq/sample_demultiplex.Rmd', knit_root_dir = "${{PWD}}")" &> {log}
    """

rule scrnaseq_integrate_annotations:
  input:
    seurat_rds=expand("{basedir}/output/scrnaseq/normalized/{run_fbc}_normalized_SeuratObject.Rds", basedir=basedir, run_fbc=run_fbcs),
    seurat_demultiplex_csv=expand("config/annotations/demultiplexed.csv", basedir=basedir),
    seurat_annotated_csv=expand("config/annotations/celltypes.csv", basedir=basedir),
    sample_metadata_csv=config['sample_metadata'],
    renv_library="renv",
  output:
    seurat_annotated_rds="{basedir}/output/scrnaseq/cell_metadata/annotated_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/scrnaseq/r.yaml",
  log:
    "{basedir}/output/scrnaseq/cell_metadata/integrate_annotations.log",
  benchmark:
    "{basedir}/output/scrnaseq/cell_metadata/integrate_annotations_benchmark.txt",
  resources:
    mem_mb=47000,
  shell:
    """
    Rscript workflow/scripts/scrnaseq/merge_annotations.R {input.seurat_rds} {input.seurat_demultiplex_csv} {input.seurat_annotated_csv} {input.sample_metadata_csv} {output.seurat_annotated_rds} &> {log}
    """
    
rule scrnaseq_clustering:
  input:
    seurat_annotated_rds="{basedir}/output/scrnaseq/cell_metadata/annotated_SeuratObject.Rds",
    renv_library="renv",
  output:
    seurat_clustered_rds="{basedir}/output/scrnaseq/clustered/clustered_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/scrnaseq/r.yaml",
  log:
    "{basedir}/output/scrnaseq/cell_metadata/clustering.log",
  benchmark:
    "{basedir}/output/scrnaseq/cell_metadata/clustering_benchmark.txt",
  resources:
    mem_mb=47000,
  shell:
    """
    Rscript workflow/scripts/scrnaseq/clustering.R {input.seurat_annotated_rds} {output.seurat_clustered_rds} &> {log}
    """

# Analyses

rule scrnaseq_differential_abundance:
  input:
    seurat_clustered_rds="{basedir}/output/scrnaseq/clustered/clustered_SeuratObject.Rds",
    renv_library="renv",
  output:
    dacs_csv="{basedir}/output/scrnaseq/da/dacs.csv",
  threads: 
    1
  conda:
    "../envs/scrnaseq/r.yaml",
  log:
    "{basedir}/output/scrnaseq/da/dacs.log",
  benchmark:
    "{basedir}/output/scrnaseq/da/dacs_benchmark.txt",
  resources:
    mem_mb=47000,
  shell:
    """
    Rscript workflow/scripts/scrnaseq/differential_abundance.R "{input.seurat_clustered_rds}" "{output.dacs_csv}" &> "{log}"
    """

rule scrnaseq_differential_expression:
  input:
    seurat_annotated_rds="{basedir}/output/scrnaseq/clustered/clustered_SeuratObject.Rds",
    renv_library="renv",
  output:
    degs_rds="{basedir}/output/scrnaseq/de/degs.Rds",
    pb_rds="{basedir}/output/scrnaseq/de/pb.Rds",
  threads: 
    1
  conda:
    "../envs/scrnaseq/r.yaml",
  log:
    "{basedir}/output/scrnaseq/de/degs.log",
  benchmark:
    "{basedir}/output/scrnaseq/de/degs_benchmark.txt",
  resources:
    mem_mb=47000,
  shell:
    """
    Rscript workflow/scripts/scrnaseq/differential_expression.R "{input.seurat_annotated_rds}" "{output.degs_rds}" "{output.pb_rds}" &> "{log}"
    """
