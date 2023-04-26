rule brnaseq_copy_fastqs:
  input:
    fastq=lambda w: (brnaseq_fastq_df[brnaseq_fastq_df.file == w.brnaseq_read_fn].old_path + '.fastq.gz').tolist(),
  output:
    "resources/brnaseq/fastq/{brnaseq_read_fn}.fastq.gz",
  conda:
    "../envs/rsync.yaml"
  message:
    "--- CD14+ monocytes bRNAseq: Copying {input.fastq} to scratch ---",
  log:
    "resources/brnaseq/fastq/rsync_{brnaseq_read_fn}.log",
  params:
    brnaseq_read_fn="{brnaseq_read_fn}",
  shell:
    """
    if [ ! -f "{input.fastq}" ]; then
      echo "{input.fastq} does not exist!"
      exit 1
    else
      rsync -ar "{input.fastq}" "{output}" &> "{log}"
    fi
    """

rule brnaseq_copy_star_reference:
  input:
    expand("{brnaseq_reference_dir}", brnaseq_reference_dir=config['brnaseq_reference_dir']),
  output:
    directory("resources/brnaseq/reference_genome"),
  conda:
    "../envs/rsync.yaml",
  message:
    "--- CD14+ monocytes bRNAseq: Copying reference genome to scratch ---",
  shell:
    """
    rsync -ar "{input}/" "{output}"
    """

rule brnaseq_fastqc:
  input:
    "resources/brnaseq/fastq/{brnaseq_fastq_fn}.fastq.gz"
  output:
    fastqcdir=directory("output/brnaseq/fastqc/{brnaseq_fastq_fn}"),
    fastqchtml="output/brnaseq/fastqc/{brnaseq_fastq_fn}/{brnaseq_fastq_fn}_fastqc.html",
    fastqczip="output/brnaseq/fastqc/{brnaseq_fastq_fn}/{brnaseq_fastq_fn}_fastqc.zip",
  conda:
    "../envs/fastqc.yaml"
  log:
    "output/brnaseq/fastqc/{brnaseq_fastq_fn}/{brnaseq_fastq_fn}_fastqc.log",
  message:
    "--- CD14+ monocytes bRNAseq: FastQC {params.brnaseq_fastq_fn} ---"
  threads: 
    8
  benchmark:
    "output/brnaseq/fastqc/{brnaseq_fastq_fn}/{brnaseq_fastq_fn}_benchmark.txt",
  params:
    brnaseq_fastq_fn="{brnaseq_fastq_fn}",
  shell:
    """
    fastqc --outdir="{output.fastqcdir}" \
           --threads "{threads}" \
           "{input}" \
           &> "{log}"
    """

rule brnaseq_umitools_extract:
  input:
    fastq_r1="resources/brnaseq/fastq/{brnaseq_sampleID}_R1.fastq.gz",
    fastq_r2="resources/brnaseq/fastq/{brnaseq_sampleID}_R2.fastq.gz",
    fastq_umi="resources/brnaseq/fastq/{brnaseq_sampleID}_UMI.fastq.gz",
  output:
    fastq_r1_umi="resources/brnaseq/fastq_umi/{brnaseq_sampleID}_umi_R1.fastq.gz",
    fastq_r2_umi="resources/brnaseq/fastq_umi/{brnaseq_sampleID}_umi_R2.fastq.gz",
  conda:
    "../envs/umitools.yaml"
  message:
    "--- CD14+ monocytes bRNAseq: UMI-tools extract and concatenate UMI to {params.brnaseq_sampleID}.fastq.gz ---"
  params:
    brnaseq_sampleID="{brnaseq_sampleID}"
  shell:
    """
    umi_tools extract --bc-pattern={config[brnaseq_umi]} --stdin={input.fastq_umi} --read2-in={input.fastq_r1} --stdout={output.fastq_r1_umi} --read2-stdout
    umi_tools extract --bc-pattern={config[brnaseq_umi]} --stdin={input.fastq_umi} --read2-in={input.fastq_r2} --stdout={output.fastq_r2_umi} --read2-stdout
    """

rule brnaseq_star_genome_load:
  input:
    reference="resources/brnaseq/reference_genome"
  output:
    touch("output/brnaseq/genome_loaded.done")
  conda:
    "../envs/star.yaml"
  message:
    "--- CD14+ monocytes bRNAseq: STAR load genome ---"
  threads: 
    8
  shell:
    """
    STAR --genomeDir "{input.reference}" \
         --genomeLoad LoadAndExit
    """

rule brnaseq_star_pe:
  input:
    reference="resources/brnaseq/reference_genome",
    fastq_r1="resources/brnaseq/fastq_umi/{brnaseq_sampleID}_umi_R1.fastq.gz",
    fastq_r2="resources/brnaseq/fastq_umi/{brnaseq_sampleID}_umi_R2.fastq.gz",
    idx="output/brnaseq/genome_loaded.done",
  output:
    bamfile="output/brnaseq/bam/{brnaseq_sampleID}/{brnaseq_sampleID}_pe_Aligned.sortedByCoord.out.bam",
    sjfile="output/brnaseq/bam/{brnaseq_sampleID}/{brnaseq_sampleID}_pe_SJ.out.tab",
    logfile="output/brnaseq/bam/{brnaseq_sampleID}/{brnaseq_sampleID}_pe_Log.out",
    logprogressfile="output/brnaseq/bam/{brnaseq_sampleID}/{brnaseq_sampleID}_pe_Log.progress.out",
    logfinalfile="output/brnaseq/bam/{brnaseq_sampleID}/{brnaseq_sampleID}_pe_Log.final.out",
  conda:
    "../envs/star.yaml"
  log:
    "output/brnaseq/bam/{brnaseq_sampleID}/{brnaseq_sampleID}_STAR.log",
  params:
    brnaseq_sampleID="{brnaseq_sampleID}"
  message:
    "--- CD14+ monocytes bRNAseq: STAR paired-ended alignment {params.brnaseq_sampleID} ---"
  threads: 
    10
  shell:
    """
    STAR --runThreadN {threads} \
         --genomeDir "{input.reference}" \
         --genomeLoad LoadAndKeep \
         --readFilesIn "{input.fastq_r1}" "{input.fastq_r2}" \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --limitBAMsortRAM 50000000000 \
         --outFileNamePrefix "output/brnaseq/bam/{params.brnaseq_sampleID}/{params.brnaseq_sampleID}_pe_" \
         |& tee -a "{log}"
    """

rule brnaseq_star_genome_unload:
  input:
    reference="resources/brnaseq/reference",
    bam=expand("output/brnaseq/bam/{brnaseq_sampleID}/{brnaseq_sampleID}_pe_Aligned.sortedByCoord.out.bam", brnaseq_sampleID=brnaseq_sampleIDs),
    idx="output/brnaseq/genome_loaded.done",
  output:
    "output/brnaseq/bam/STAR_genome_unload_Log.out"
  conda:
    "../envs/star.yaml"
  log:
    "output/brnaseq/bam/STAR.log"
  message:
    "--- CD14 monocytes bRNAseq: STAR unload genome ---"
  threads: 
    8
  shell:
    """
    STAR --genomeDir "{input.reference}" \
         --genomeLoad Remove \
         --outFileNamePrefix "output/brnaseq/bam/STAR_genome_unload_" \
         |& tee -a "{log}"
    rm "{input.idx}"
    """

rule brnaseq_samtools_filter:
  input:
    bam="output/brnaseq/bam/{brnaseq_sampleID}/{brnaseq_sampleID}_pe_Aligned.sortedByCoord.out.bam",
  output:
    bam_filtered="output/brnaseq/bam/{brnaseq_sampleID}/{brnaseq_sampleID}_filtered.bam",
  conda:
    "../envs/samtools.yaml"
  log:
    "output/brnaseq/bam/{brnaseq_sampleID}/{brnaseq_sampleID}_samtools.log"
  message:
    "--- CD14 monocytes: SAMtools processing {params.brnaseq_sampleID} ---",
  params:
    brnaseq_sampleID="{brnaseq_sampleID}"
  threads: 
    8
  shell:
    """        
    #Remove:
    # - unmapped reads and multiple mappings (4 + 256 = 260)
	  # - reads with mapping score < 10
	  # - mitochondrial sequences
    
		samtools view -@ {threads} -S -h -F 260 -q 10 "{input.bam}" | awk '($1 ~ /^@/) || ($3 != "MT") {{ print $0 }}' | samtools view -@ {threads} -b -o "{output.bam_filtered}" - |& tee -a "{log}"
    """

rule brnaseq_samtools_index:
  input:
    bam="output/brnaseq/bam/{brnaseq_sampleID}/{brnaseq_sampleID}_filtered.bam",
  output:
    bai="output/brnaseq/bam/{brnaseq_sampleID}/{brnaseq_sampleID}_filtered.bam.bai",
  conda:
    "../envs/samtools.yaml"
  log:
    "output/brnaseq/bam/{brnaseq_sampleID}/{brnaseq_sampleID}_samtools.log",
  params:
    brnaseq_sampleID="{brnaseq_sampleID}",
  message:
    "--- CD14 monocytes: SAMtools indexing {params.brnaseq_sampleID} ---"
  threads: 
    8
  shell:
    """        
		samtools index -b -@ {threads} "{input.bam}" |& tee -a "{log}"
    """

rule brnaseq_umitools_dedup:
  input:
    bam="output/brnaseq/bam/{brnaseq_sampleID}/{brnaseq_sampleID}_filtered.bam",
    bai="output/brnaseq/bam/{brnaseq_sampleID}/{brnaseq_sampleID}_filtered.bam.bai",
  output:
    bam_deduplicated="output/brnaseq/bam/{brnaseq_sampleID}/{brnaseq_sampleID}_filtered_deduplicated.bam",
  conda:
    "../envs/umitools.yaml"
  message:
    "--- UMI-tools: deduplicating {params.brnaseq_sampleID} ---"
  params:
    brnaseq_sampleID="{brnaseq_sampleID}"
  shell:
    """
    umi_tools dedup -I "{input.bam}" --paired --buffer-whole-contig --method directional --output-stats="output/brnaseq/bam/{params.brnaseq_sampleID}/{params.brnaseq_sampleID}_umitools_dedup_stats" -S "{output.bam_deduplicated}"
    """

rule brnaseq_featurecounts:
  input:
    bams=expand("output/brnaseq/bam/{brnaseq_sampleID}/{brnaseq_sampleID}_filtered_deduplicated.bam", brnaseq_sampleID=brnaseq_sampleIDs),
    brnaseq_gtf_file=expand("{brnaseq_gtf_file}", brnaseq_gtf_file=config['brnaseq_gtf_file']),
  output:
    "output/brnaseq/counts/counts.txt"
  conda:
    "../envs/featurecounts.yaml"
  log:
    "output/brnaseq/counts/featurecounts.log"
  message:
    "--- CD14 monocytes: Subread::FeatureCounts feature counting ---"
  threads: 
    8
  shell:
    """    
		featureCounts -T {threads} \
                  -a "{input.brnaseq_gtf_file}" \
                  -t exon \
                  -g gene_id \
                  -p \
                  -o "{output}" \
                  "{input.bams}"
    """
    
rule brnaseq_multiqc:
  input:
    fastqc=expand("output/brnaseq/fastqc/{brnaseq_fastq_fn}/{brnaseq_fastq_fn}_fastqc.html", brnaseq_fastq_fn=brnaseq_fastq_fns),
    featurecounts="output/brnaseq/counts/counts.txt",
  output:
    multiqcdir=directory('output/brnaseq/multiqc'),
  conda:
    "../envs/multiqc.yaml"
  log:
    "output/brnaseq/counts/multiqc.log"
  message:
    "--- CD14 monocytes: MultiQC summarizing statistics ---"
  threads: 
    8
  shell:
    """    
		multiqc "output" -o "{output.multiqcdir}"
    """

rule brnaseq_copy_output:
  input:
    # fastqcdir="output/brnaseq/fastqc/{brnaseq_fastq_fn}",
    # bam="output/brnaseq/bam/{brnaseq_sampleID}/{brnaseq_sampleID}_filtered.bam",
    # bai="output/brnaseq/bam/{brnaseq_sampleID}/{brnaseq_sampleID}_filtered.bam.bai",
    # bam_deduplicated="output/brnaseq/bam/{brnaseq_sampleID}/{brnaseq_sampleID}_filtered_deduplicated.bam",
    multiqcdir="output/brnaseq/multiqc",
    featurecounts="output/brnaseq/counts/counts.txt",
  output:
    outputdir=directory("{basedir}/output/brnaseq"),
    featurecounts_outputdir="{basedir}/output/brnaseq/counts/counts.txt",
  conda:
    "../envs/rsync.yaml",
  log:
    "{basedir}/output/brnaseq/copy_brnaseq_output.log",
  message:
    "--- CD14 monocytes: Copying output to {output.outputdir} ---",
  shell:
    """   
		rsync -ar "{input.multiqcdir}/.." -o "{output.outputdir}"
    """

rule summarizedexperiment_preparation:
  input:
    raw_counts_txt="{basedir}/output/brnaseq/counts/counts.txt",
    sample_metadata_xlsx=config['sample_metadata'],
  output:
    cleaned_counts_se_rds="{basedir}/output/brnaseq/counts/counts_cleaned_se.rds",
  conda:
    "../envs/r-deseq2.yaml",
  log:
    "{basedir}/output/brnaseq/counts/summarizedexperiment_preparation.log",
  message:
    "--- CD14 monocytes: Preparing summarizedExperiment ---"
  shell:
    """   
		Rscript workflow/scripts/brnaseq/summarizedexperiment_preparation.R "{input.raw_counts_txt}" "{input.sample_metadata_xlsx}" "{output.cleaned_counts_se_rds}" &> "{log}"
    """

rule differential_expression:
  input:
    cleaned_counts_se_rds="{basedir}/output/brnaseq/counts/counts_cleaned_se.rds",
  output:
    dds_rds="{basedir}/output/brnaseq/de/dds.rds",
    rld_rds="{basedir}/output/brnaseq/de/rld.rds",
    degs_csv="{basedir}/output/brnaseq/de/degs_rvnr.csv",
  conda:
    "../envs/r-deseq2.yaml",
  log:
    "{basedir}/output/differential_expression.log",
  message:
    "--- CD14 monocytes: Performing differential expression ---"
  shell:
    """   
		Rscript workflow/scripts/brnaseq/differential_expression.R "{input.cleaned_counts_se_rds}" "{output.dds_rds}" "{output.rld_rds}" "{output.degs_csv}" &> "{log}"
    """
