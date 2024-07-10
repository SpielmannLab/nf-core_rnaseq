# providing files to use in which working directory
yaml_file <- "params_limma.yaml"
rmd_file <- "limma.rmd"
workDir <- getwd()

setwd(workDir)

# loading libraries
suppressPackageStartupMessages(library(rmarkdown))
suppressPackageStartupMessages(library(yaml))

# setup output directory
params <- read_yaml(yaml_file)
list2env(params, envir = environment())

compared_A <- sprintf("%s-%s", reference_condition_A, alternate_condition_A) # the minus - is needed for makeContrasts function
compared_B <- sprintf("%s-%s", reference_condition_B, alternate_condition_B)
outdir_folder <- sprintf("%s_vs_%s",compared_A,compared_B)
dir.create(outdir_folder)
outdir <- sprintf("%s/%s",out_dir,outdir_folder)
render(rmd_file, output_file=sprintf("report_limma"), output_dir = outdir)
