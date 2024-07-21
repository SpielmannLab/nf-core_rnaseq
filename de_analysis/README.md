# DE analysis 
To perform a DE analysis the tools DESeq2 and Limma are implemented. Custom R scripts were developed to perform the analyzes.

>To apply the DE analysis to a specific project, only the [parameters file](./limma/params_limma.yaml) file needs to be adapted. More detailed descriptions for the individual parameters can be found in the file.

## DESeq2
To run the DE analysis with DESeq2, the three files [deseq2.R](./deseq2/deseq2.R), [deseq2.rmd](./deseq2/deseq2.rmd) and [params_deseq2.yml](./deseq2/params_deseq2.yaml) must be available. In an environment that contains the packages as in [de_analysis.yml](./de_analysis.yml), the script deseq2.R can be executed with the following oneliner:
```
Rscript deseq2.R
```

## Limma
The procedure for using Limma is similar to that for DESeq2. For the execution, the three files [limma.R](./limma/limma.R), [limma.rmd](./limma/limma.rmd) and [params_limma.yml](./limma/params_limma.yaml) must be available. With the following oneliner, the script can be executed in an environment that contains the packages as in [de_analysis.yml](./de_analysis.yml).
```
Rscript limma.R
```
