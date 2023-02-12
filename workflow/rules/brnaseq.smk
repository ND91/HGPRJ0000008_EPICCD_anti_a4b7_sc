# Preparation
rule copy_brnaseq_fastqs:
  input:
    fastq=lambda w: (filedata[filedata.new_path == w.readfile].old_path + '.fastq.gz').tolist()
  output:
    "resources/fastq/{readfile}.fastq.gz"
  conda:
    "envs/rsync.yaml"
  message:
    "--- Copying {input.fastq} to scratch ---"
  shell:
    """
    rsync -a {input.fastq} {output}
    """

rule copy_reference:
  input:
    config['reference_dir']
  output:
    directory("resources/reference")
  conda:
    "envs/rsync.yaml"
  message:
    "--- Copying reference index genome to scratch ---"
  shell:
    """
    rsync -ar {input}/ {output}
    """

# QC

rule fastqc:
  input:
    "resources/fastq/{readfile}.fastq.gz"
  output:
    fastqcdir=directory("results/fastqc/{readfile}"),
    fastqchtml="results/fastqc/{readfile}/{readfile}_fastqc.html",
    fastqczip="results/fastqc/{readfile}/{readfile}_fastqc.zip",
  conda:
    "envs/fastqc.yaml"
  log:
    "results/fastqc/{readfile}/{readfile}_fastqc.log",
  message:
    "--- FastQC {wildcards.readfile} ---"
  threads: 
    8
  benchmark:
    "results/fastqc/{readfile}/{readfile}.txt",
  shell:
    """
    fastqc --outdir={output.fastqcdir} \
           --threads {threads} \
           {input} \
           &> {log}
    """

# Alignment, deduplicating, and mapping

rule umitools_extract:
  input:
    fastq_r1="resources/fastq/{sampleID}_R1.fastq.gz",
    fastq_r2="resources/fastq/{sampleID}_R2.fastq.gz",
    fastq_umi="resources/fastq/{sampleID}_UMI.fastq.gz",
  output:
    fastq_r1_umi="resources/fastq_umi/{sampleID}_umi_R1.fastq.gz",
    fastq_r2_umi="resources/fastq_umi/{sampleID}_umi_R2.fastq.gz",
  conda:
    "envs/umitools.yaml"
  message:
    "--- UMI-tools: extracting and concatenating UMI to {wildcards.sampleID}.fastq.gz ---"
  shell:
    """
    umi_tools extract --bc-pattern={config[umi]} --stdin={input.fastq_umi} --read2-in={input.fastq_r1} --stdout={output.fastq_r1_umi} --read2-stdout
    umi_tools extract --bc-pattern={config[umi]} --stdin={input.fastq_umi} --read2-in={input.fastq_r2} --stdout={output.fastq_r2_umi} --read2-stdout
    """

rule star_genome_load:
  input:
    reference="resources/reference"
  output:
    touch("results/genome_loaded.done")
  conda:
    "envs/star.yaml"
  message:
    "--- STAR: load genome ---"
  threads: 
    8
  shell:
    """
    STAR --genomeDir {input.reference} \
         --genomeLoad LoadAndExit
    """

rule star_pe:
  input:
    reference="resources/reference",
    fastq_r1="resources/fastq_umi/{sampleID}_umi_R1.fastq.gz",
    fastq_r2="resources/fastq_umi/{sampleID}_umi_R2.fastq.gz",
    idx="results/genome_loaded.done",
  output:
    bamfile="results/bam/{sampleID}/{sampleID}_pe_Aligned.sortedByCoord.out.bam",
    sjfile="results/bam/{sampleID}/{sampleID}_pe_SJ.out.tab",
    logfile="results/bam/{sampleID}/{sampleID}_pe_Log.out",
    logprogressfile="results/bam/{sampleID}/{sampleID}_pe_Log.progress.out",
    logfinalfile="results/bam/{sampleID}/{sampleID}_pe_Log.final.out",
  conda:
    "envs/star.yaml"
  log:
    "results/bam/{sampleID}/{sampleID}_STAR.log"
  message:
    "--- STAR: paired-ended alignment {wildcards.sampleID} ---"
  threads: 
    10
  shell:
    """
    STAR --runThreadN {threads} \
         --genomeDir {input.reference} \
         --genomeLoad LoadAndKeep \
         --readFilesIn {input.fastq_r1} {input.fastq_r2} \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --limitBAMsortRAM 50000000000 \
         --outFileNamePrefix results/bam/{wildcards.sampleID}/{wildcards.sampleID}_pe_ \
         |& tee -a {log}
    """

rule star_genome_unload:
  input:
    reference="resources/reference",
    bam=expand("results/bam/{sampleID}/{sampleID}_pe_Aligned.sortedByCoord.out.bam", sampleID=sampleIDs),
    idx="results/genome_loaded.done",
  output:
    "results/bam/STAR_genome_unload_Log.out"
  conda:
    "envs/star.yaml"
  log:
    "results/bam/STAR.log"
  message:
    "--- STAR: unload genome ---"
  threads: 
    8
  shell:
    """
    STAR --genomeDir {input.reference} \
         --genomeLoad Remove \
         --outFileNamePrefix results/bam/STAR_genome_unload_ \
         |& tee -a {log}
    rm {input.idx}
    """

rule samtools_filter:
  input:
    bam="results/bam/{sampleID}/{sampleID}_pe_Aligned.sortedByCoord.out.bam",
  output:
    bam_filtered="results/bam/{sampleID}/{sampleID}_filtered.bam",
  conda:
    "envs/samtools.yaml"
  log:
    "results/bam/{sampleID}/{sampleID}_samtools.log"
  message:
    "--- SAMtools: post-alignment processing {wildcards.sampleID} ---"
  threads: 
    8
  shell:
    """        
    #Remove:
    # - unmapped reads and multiple mappings (4 + 256 = 260)
	  # - reads with mapping score < 10
	  # - mitochondrial sequences
    
		samtools view -@ {threads} -S -h -F 260 -q 10 {input.bam} | awk '($1 ~ /^@/) || ($3 != "MT") {{ print $0 }}' | samtools view -@ {threads} -b -o {output.bam_filtered} - |& tee -a {log}
    """

rule samtools_index:
  input:
    bam="results/bam/{sampleID}/{sampleID}_filtered.bam",
  output:
    bai="results/bam/{sampleID}/{sampleID}_filtered.bam.bai",
  conda:
    "envs/samtools.yaml"
  log:
    "results/bam/{sampleID}/{sampleID}_samtools.log"
  message:
    "--- SAMtools: indexing {wildcards.sampleID} ---"
  threads: 
    8
  shell:
    """        
		samtools index -b -@ {threads} {input.bam} |& tee -a {log}
    """

rule umitools_dedup:
  input:
    bam="results/bam/{sampleID}/{sampleID}_filtered.bam",
    bai="results/bam/{sampleID}/{sampleID}_filtered.bam.bai",
  output:
    bam_deduplicated="results/bam/{sampleID}/{sampleID}_filtered_deduplicated.bam",
  conda:
    "envs/umitools.yaml"
  message:
    "--- UMI-tools: deduplicating {wildcards.sampleID}.bam ---"
  shell:
    """
    umi_tools dedup -I {input.bam} --paired --buffer-whole-contig --method directional --output-stats=results/bam/{wildcards.sampleID}/{wildcards.sampleID}_umitools_dedup_stats -S {output.bam_deduplicated}
    """

rule featurecounts:
  input:
    bams=expand("results/bam/{sampleID}/{sampleID}_filtered_deduplicated.bam", sampleID=sampleIDs),
    annotationfile=config['gtf_file'],
  output:
    "results/counts/counts.txt"
  conda:
    "envs/featurecounts.yaml"
  log:
    "results/counts/featurecounts.log"
  message:
    "--- Subread::FeatureCounts: Feature counting ---"
  threads: 
    8
  shell:
    """    
		featureCounts -T {threads} \
                  -a {input.annotationfile} \
                  -t exon \
                  -g gene_id \
                  -p \
                  -o {output} \
                  {input.bams}
    """
    
rule multiqc:
  input:
    fastqc=expand("results/fastqc/{readfile}/{readfile}_fastqc.html", readfile=readfiles),
    featurecounts="results/counts/counts.txt",
  output:
    multiqcdir=directory('results/multiqc'),
  conda:
    "envs/multiqc.yaml"
  log:
    "results/counts/multiqc.log"
  message:
    "--- MultiQC: Summarizing statistics ---"
  threads: 
    8
  shell:
    """    
		multiqc . -o {output.multiqcdir}
    """