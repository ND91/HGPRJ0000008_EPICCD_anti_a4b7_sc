rule renv_preparation:
  input:
    renv_lock=config["renv_lockfile"],
    renv_rprofile=config["renv_rprofile"],
    renv_activate=config["renv_activate"],
  output:
    renv_library=directory("renv"),
  threads:
    1
  conda:
    "../envs/r.yaml"
  log:
    "renv_preparation.log",
  benchmark:
    "renv_preparation_benchmark.log"
  message:
	  "Install R packages in renv"
  shell:
    """
    mkdir "renv"
    cp "{input.renv_activate}" renv/
    cp "{input.renv_rprofile}" ./
    cp "{input.renv_lock}" ./
    Rscript "workflow/scripts/renv_setup.R" &> "{log}"
    """