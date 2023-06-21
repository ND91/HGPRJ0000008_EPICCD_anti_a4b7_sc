# Response-associated single-cell characterization of Crohn's disease patients treated with vedolizumab 

Here we present the pipeline used to process the data. The data encompasses:
- Single-cell RNA-sequencing data on PBMCs isolated from Crohn's disease (CD) patients and prepared using BioLegend hashtag multiplexed single cell suspensions prepared using the 10X Genomics (chemistry v3.1) whereafter the libraries were sequenced on the HiSeq4000.
- Mass cytometry time-of-flight (CyTOF) on PBMCs isolated from CD patients as generated using the Helios.
- Flow cytometry on the PBMCs isolated from CD patients as generated using XXX and subsequently processed using analyzed and quantified using FlowJo. 
- Bulk RNA-sequencing data on CD14+ monocytes isolated through flow cytometry from CD patients.

Note that some steps were performed on non-R/Python platforms (flow cytometry was performed using FlowJo). Exact reproducibility is therefore limited to what I got from this software. Similarly, mass cytometry was preprocessed by one of the co-authors.

For this pipeline to work, you will need to have [snakemake installed](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). Furthermore, this pipeline makes use of conda environments. Use this pipeline at your own risk!

## Quickstart

You should first setup the conda environment as defined in the `environment.yml` file, which will install `snakemake`, `Python`, as well as the Python packages `pandas` and `openpyxl`.

```
conda env create -f environment.yml
```

To **test** this pipeline, use.

```
snakemake --cores [cores] --use-conda -np
```

This pipeline was run using the following command. Amend to your own computational facilities as necessary.

```
snakemake --cores 20 -p --use-conda --resources mem_mb=189000 --rerun-triggers mtime
```