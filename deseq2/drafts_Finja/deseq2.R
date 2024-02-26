suppressPackageStartupMessages(library(rmarkdown))
suppressPackageStartupMessages(library(yaml))

params <- read_yaml("params_deseq2.yaml")
list2env(params, envir = environment())
render("deseq2.rmd", output_file=sprintf("%s_vs_%s_deseq2",alternate_condition,reference_condition), output_dir = outdir)
