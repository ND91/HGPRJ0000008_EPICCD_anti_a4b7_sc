import pandas
import re

sample_metadata = config["sample_metadata"]
sample_metadata_df = pandas.read_excel(sample_metadata, header=0)
fastq_r1 = sample_metadata_df.loc[:, 'R1_path'].tolist()
fastq_r2 = sample_metadata_df.loc[:, 'R2_path'].tolist()
sampleIDs = sample_metadata_df.loc[:, 'Sample_ID'].tolist()

rule all:
  input:
    expand("mapped_reads/{sample_r1}.bam", sample=fastq_r1),
    expand("mapped_reads/{sample_r2}.bam", sample=fastq_r2)

rule fastqc:
	input:
 		"data/samples/{fastq_r1}.fastq",
		"data/samples/{fastq_r2}.fastq"
	output:
		fastqcdir=directory(config["outputdir"]+ "/fastqc/{sampleID}/"),
		fastqchtml=config["outputdir"] + "/fastqc/{sampleID}/{sampleID}.html",
		fastqczip=config["outputdir"] + "/fastqc/{sampleID}/{sampleID}.zip",
	log:
		"output/fastqc/{sampleID}/{sampleID}_fastqc.log",
	threads: 
		8
	benchmark:
		"output/fastqc/{sampleID}/{sampleID}_benchmark.txt",
	params:
		sampleID="{sampleID}"
	shell:
		"""
		#mkdir -p {output.fastqcdir}
		
		#fastqc --outdir={output.fastqcdir} --threads {threads} {input} &> {log}
		
		echo {output.fastqcdir}
		echo {threads}
		echo {input}
		echo {log}
		"""