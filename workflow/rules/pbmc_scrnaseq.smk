rule scrnaseq_copy_fastqs:
  input:
    fastqdir=lambda w: (scrnaseq_fastq_df[scrnaseq_fastq_df.file == w.scrnaseq_fastq_fn].old_path).tolist(),
  output:
    directory("resources/scrnaseq/fastq/{scrnaseq_fastq_fn}"),
  log:
    "resources/scrnaseq/fastq/rsync_{scrnaseq_fastq_fn}.log",
  conda:
    "../envs/rsync.yaml"
  message:
    "--- scRNAseq: Copying {input.fastqdir} to scratch ---"
  shell:
    """
    if [ ! -d "{input.fastqdir}" ]; then
      echo "{input.fastqdir} does not exist!"
      exit 1
    else
      rsync -ar "{input.fastqdir}/" "{output}" &> "{log}"
    fi
    """

rule scrnaseq_copy_cellranger_reference:
  input:
    config['scrnaseq_reference_dir']
  output:
    directory("resources/scrnaseq/reference_genome"),
  log:
    "resources/scrnaseq/rsync_scrnaseq_reference.log",
  conda:
    "../envs/rsync.yaml"
  message:
    "--- scRNAseq: Copying reference genome to scratch ---"
  shell:
    """
    rsync -ar "{input}/" "{output}" &> "{log}"
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
    run_fbc="{run_fbc}",
    cellranger_basepath=config["cellranger_basepath"],
  resources:
    mem_mb=64000,
  params:
    run_fbc="{run_fbc}",
  message:
    "--- scRNAseq: CellRanger count {params.run_fbc} ---"
  shell:
    """
    rm -rf {output.cellrangerdir}
    mem_gb=$(echo {resources.mem_mb}/1024|bc)
    mkdir -p {output.cellrangerdir}
    cd "output/scrnaseq/cellranger"
    #find . -type d -empty -delete
    {params.cellranger_basepath}/cellranger count \
      --id="{params.run_fbc}" \
      --libraries=$OLDPWD/{input.libraries} \
      --feature-ref=$OLDPWD/{input.features} \
      --transcriptome=$OLDPWD/{input.reference_genome} \
      --nosecondary \
      --disable-ui \
      --localcores={threads} \
      --localmem=$mem_gb &> $OLDPWD/{log}
    """

rule scrnaseq_cellranger_count_copy:
  input:
    cellrangertmpdir="output/scrnaseq/cellranger/{run_fbc}",
  output:
    cellrangerfinaldir=directory("{basedir}/output/scrnaseq/cellranger/{run_fbc}"),
  conda:
    "../envs/rsync.yaml"
  log:
    "{basedir}/output/scrnaseq/{run_fbc}_cellranger_copy.log",
  threads: 
    1
  benchmark:
    "{basedir}/output/scrnaseq/{run_fbc}_cellranger_mv.benchmark.txt"
  resources:
    mem_mb=10000,
  params:
    run_fbc="{run_fbc}",
  message:
    "--- scRNAseq: Copying CellRanger output {params.run_fbc} to {basedir} ---"
  shell:
    """
    rsync -ar "{input.cellrangertmpdir}" "{output.cellrangerfinaldir}" &> "{log}"
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
    "--- scRNAseq: Downloading PBMC reference dataset in h5seurat format from Hao et al. 2020 ---"
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
  output:
    pbmc_reference_rds="{basedir}/resources/scrnaseq/reference_data/pbmc_multimodal.Rds",
  threads: 
    1
  conda:
    "../envs/r-seuratdisk.yaml",
  message:
    "--- scRNAseq: Converting PBMC reference dataset from h5seurat to rds ---"
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
  output:
    seurat_rds="{basedir}/output/scrnaseq/normalized/{run_fbc}_normalized_SeuratObject.Rds",
  threads:
    1
  conda:
    "../envs/r-seurat.yaml",
  message:
    "--- scRNAseq: Normalizing {params.run_fbc} ---"
  log:
    "{basedir}/output/scrnaseq/normalized/{run_fbc}.log",
  benchmark:
    "{basedir}/output/scrnaseq/normalized/{run_fbc}_benchmark.txt"
  resources:
    mem_mb=16000,
  params:
    run_fbc="{run_fbc}",
    lowerbound_numis=600,
  shell:
    """
    Rscript workflow/scripts/scrnaseq/pbmc_normalization.R "{input.count_mtx}" "{output.seurat_rds}" "{params.run_fbc}" "{params.lowerbound_numis}" &> "{log}"
    """

rule scrnaseq_celltype_annotation:
  input:
    seurat_rds="{basedir}/output/scrnaseq/normalized/{run_fbc}_normalized_SeuratObject.Rds",
    pbmc_rds="{basedir}/resources/scrnaseq/reference_data/pbmc_multimodal.Rds",
  output:
    seurat_annotated_csv="{basedir}/output/scrnaseq/cell_metadata/celltype_annotation/{run_fbc}_celltype_annotated.csv",
  threads: 
    1
  conda:
    "../envs/r-seurat.yaml",
  message:
    "--- scRNAseq: Automatically annotating {params.run_fbc} ---",
  log:
    "{basedir}/output/scrnaseq/cell_metadata/celltype_annotation/{run_fbc}_celltype_annotation.log",
  benchmark:
    "{basedir}/output/scrnaseq/cell_metadata/celltype_annotation/{run_fbc}_celltype_annotation_benchmark.txt",
  resources:
    mem_mb=47000,
  params:
    run_fbc="{run_fbc}",
  shell:
    """
    Rscript -e "rmarkdown::render('workflow/scripts/scrnaseq/pbmc_celltype_annotation.Rmd', knit_root_dir = "{basedir}")" &> "{log}"
    """

rule scrnaseq_sample_demultiplex:
  input:
    seurat_rds="{basedir}/output/scrnaseq/normalized/{run_fbc}_normalized_SeuratObject.Rds",
  output:
    seurat_demultiplex_csv="{basedir}/output/scrnaseq/cell_metadata/sample_demultiplex/{run_fbc}_demultiplex.csv",
  threads: 
    1
  conda:
    "../envs/r-seurat.yaml",
  message:
    "--- scRNAseq: Demultiplexing {params.run_fbc} ---",
  log:
    "{basedir}/output/scrnaseq/cell_metadata/sample_demultiplex/{run_fbc}_demultiplex_annotation.log",
  benchmark:
    "{basedir}/output/scrnaseq/cell_metadata/sample_demultiplex/{run_fbc}_demultiplex_benchmark.txt",
  resources:
    mem_mb=47000,
  params:
    run_fbc="{run_fbc}",
  shell:
    """
    Rscript -e "rmarkdown::render('workflow/scripts/scrnaseq/pbmc_sample_demultiplex.Rmd', knit_root_dir = "{basedir}")" &> "{log}"
    """

rule scrnaseq_integrate_annotations:
  input:
    seurat_rds=expand("{basedir}/output/scrnaseq/normalized/{run_fbc}_normalized_SeuratObject.Rds", basedir=basedir, run_fbc=run_fbcs),
    seurat_demultiplex_csv=expand("config/annotations/demultiplexed.csv", basedir=basedir),
    seurat_annotated_csv=expand("config/annotations/celltypes.csv", basedir=basedir),
    sample_metadata_csv=config['sample_metadata'],
  output:
    seurat_annotated_rds="{basedir}/output/scrnaseq/cell_metadata/annotated_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r-seurat.yaml",
  log:
    "{basedir}/output/scrnaseq/cell_metadata/integrate_annotations.log",
  message:
    "--- scRNAseq: Integrating annotations ---",
  benchmark:
    "{basedir}/output/scrnaseq/cell_metadata/integrate_annotations_benchmark.txt",
  resources:
    mem_mb=47000,
  shell:
    """
    Rscript workflow/scripts/scrnaseq/pbmc_merge_annotations.R "{input.seurat_rds}" "{input.seurat_demultiplex_csv}" "{input.seurat_annotated_csv}" "{input.sample_metadata_csv}" "{output.seurat_annotated_rds}" &> "{log}"
    """
    
rule scrnaseq_clustering:
  input:
    seurat_annotated_rds="{basedir}/output/scrnaseq/cell_metadata/annotated_SeuratObject.Rds",
  output:
    seurat_clustered_rds="{basedir}/output/scrnaseq/clustered/clustered_SeuratObject.Rds",
    umap_csv="{basedir}/output/scrnaseq/clustered/umap.csv",
  threads: 
    1
  conda:
    "../envs/r-seurat.yaml",
  log:
    "{basedir}/output/scrnaseq/clustered/clustering.log",
  message:
    "--- scRNAseq: Dimension reduction and clustering ---",
  benchmark:
    "{basedir}/output/scrnaseq/clustered/clustering_benchmark.txt",
  resources:
    mem_mb=47000,
  shell:
    """
    Rscript workflow/scripts/scrnaseq/pbmc_clustering.R "{input.seurat_annotated_rds}" "{output.seurat_clustered_rds}" "{output.umap_csv}" &> "{log}"
    """

rule monocyte_subsetting:
  input:
    seurat_annotated_rds="{basedir}/output/scrnaseq/cell_metadata/annotated_SeuratObject.Rds",
  output:
    seurat_monocyte_rds="{basedir}/output/scrnaseq/monocyte/monocyte_SeuratObject.Rds",
    umap_csv="{basedir}/output/scrnaseq/monocyte/umap.csv",
  threads: 
    1
  conda:
    "../envs/r-seurat.yaml",
  log:
    "{basedir}/output/scrnaseq/monocyte/monocyte_subsetting.log",
  message:
    "--- scRNAseq: Monocyte subsetting, dimension reduction, and clustering ---",
  benchmark:
    "{basedir}/output/scrnaseq/monocyte/monocyte_subsetting_benchmark.txt",
  resources:
    mem_mb=47000,
  shell:
    """
    Rscript workflow/scripts/scrnaseq/monocyte_subsetting.R "{input.seurat_annotated_rds}" "{output.seurat_monocyte_rds}" "{output.umap_csv}" &> "{log}"
    """

# Analyses

rule scrnaseq_differential_abundance:
  input:
    seurat_clustered_rds="{basedir}/output/scrnaseq/clustered/clustered_SeuratObject.Rds",
  output:
    dacs_csv="{basedir}/output/scrnaseq/da/dacs.csv",
  threads: 
    1
  conda:
    "../envs/r-speckle.yaml",
  log:
    "{basedir}/output/scrnaseq/da/dacs.log",
  message:
    "--- scRNAseq: Differential abundance analyses ---",
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
  output:
    degs_list_rds="{basedir}/output/scrnaseq/de/degs_{level}_list.Rds",
  threads: 
    1
  conda:
    "../envs/r-deseq2.yaml",
  log:
    "{basedir}/output/scrnaseq/de/scrnaseq_differential_expression_{level}.log",
  message:
    "--- scRNAseq: Differential expression analyses per celltype ---",
  benchmark:
    "{basedir}/output/scrnaseq/de/scrnaseq_differential_expression_{level}_benchmark.txt",
  resources:
    mem_mb=47000,
  params:
    level="{level}",
  shell:
    """
    Rscript workflow/scripts/scrnaseq/pbmc_pseudobulk_differential_expression.R "{input.seurat_annotated_rds}" "{params.level}" "{output.degs_list_rds}" &> "{log}"
    """
