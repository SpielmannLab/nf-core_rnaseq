suppressPackageStartupMessages(library(rmarkdown))
suppressPackageStartupMessages(library(yaml))

yaml_file <- "params_deseq2.yaml"
rmd_file <- "deseq2.rmd"

params <- read_yaml(yaml_file)
list2env(params, envir = environment())
outdir_folder <- sprintf("%s_vs_%s",alternate_condition,reference_condition)
dir.create(outdir_folder)
outdir <- sprintf("%s/%s",out_dir,outdir_folder)
render(rmd_file, output_file=sprintf("report_%s_vs_%s_deseq2",alternate_condition,reference_condition), output_dir = outdir)
