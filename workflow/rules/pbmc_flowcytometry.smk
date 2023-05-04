rule flowcytometry_import:
  input:
    flowjo_wsp_xlsx="{basedir}/resources/flowcytometry/211209 analysis VEDO .wsp FlowJo table.xls",
    files_pbmc_flowcytometry_xlsx=config['flowcytometry_files'],
    sample_metadata_xlsx=config['sample_metadata'],
  output:
    celltype_percentage_csv="{basedir}/output/flowcytometry/celltype_percentages.csv",
  threads: 
    1
  conda:
    "../envs/r-tidyverse.yaml"
  log:
    "{basedir}/output/flowcytometry/flowcytometry_import.log",
  benchmark:
    "{basedir}/output/flowcytometry/flowcytometry_import_benchmark.txt",
  resources:
    mem_mb=47000,
  shell:
    """
    Rscript "workflow/scripts/flowcytometry/flowcytometry_import.R" "{input.flowjo_wsp_xlsx}" "{input.files_pbmc_flowcytometry_xlsx}" "{input.sample_metadata_xlsx}" "{output.celltype_percentage_csv}" &> "{log}"
    """

# Analyses

rule flowcytometry_differential_abundance:
  input:
    celltype_percentage_csv="{basedir}/output/flowcytometry/celltype_percentages.csv",
  output:
    dacs_csv="{basedir}/output/flowcytometry/da/dacs.csv",
  threads: 
    1
  conda:
    "../envs/r-tidyverse.yaml"
  log:
    "{basedir}/output/flowcytometry/da/dacs.log",
  benchmark:
    "{basedir}/output/flowcytometry/da/dacs_benchmark.txt",
  resources:
    mem_mb=2000,
  shell:
    """
    Rscript workflow/scripts/flowcytometry/flowcytometry_differential_abundance.R "{input.celltype_percentage_csv}" "{output.dacs_csv}" &> "{log}"
    """
