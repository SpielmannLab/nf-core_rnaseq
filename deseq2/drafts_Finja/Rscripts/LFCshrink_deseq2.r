"Generate a DDS object from count matrix and sample metadata (coldata) files. Also perform outlier removal based on RLE using the RUV package.

Usage:     Rscript LFCshrink_deseq2.r --reference_condition=${reference_condition} --conditions_compared=${conditions_compared} --MA=${MA} --shrink_LFC_using=${shrink_LFC_using}

Options:
    -h --help               	Show this screen.
    --conditions_compared=<value>     
    --MA=<value>                      if a MA-plot should be made <true|false>
    --shrink_LFC_using=<value>        
" -> doc

# loading libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(apeglm))

## Getting in arguments
arguments <- docopt(doc, quoted_args = TRUE)
compared <- arguments$conditions_compared
MA <- arguments$MA
shrink_LFC_using <- arguments$shrink_LFC_using

message("conditions_compared: ", compared)
message("MA: ", MA)
message("shrink_LFC_using: ", shrink_LFC_using)


###############
# Main script
###############
dds <- readRDS("dds_obj.rds")
res <- readRDS("dds_results.rds")

#the function plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet
pdf_MA <- sprintf("MA_plot_%s.pdf", compared)
pdf(pdf_MA)

if(MA){
    xlim <- c(1,1e5); ylim <- c(-3,3)
    DESeq2::plotMA(res,ylim=c(-5,5), main=compared)
}

# --- shrinkage
if(shrink_LFC_using!=FALSE){
    xlim <- c(1,1e5); ylim <- c(-3,3)

    if(shrink_LFC_using == "apeglm"){
        resLFC <- lfcShrink(dds, type="apeglm")
        res_shrink <- resLFC
        if(MA){
            DESeq2::plotMA(resLFC,ylim=c(-5,5), main=compared)
        }
    }

    if(shrink_LFC_using == "ashr"){
        resAsh <- lfcShrink(dds, type="ashr")
        res_shrink <- resAsh
        if(MA){
            DESeq2::plotMA(resAsh,ylim=c(-5,5), main=compared)
        }
    }

    if(shrink_LFC_using == "normal"){
        resNorm <- lfcShrink(dds, type="normal")
        res_shrink <- resNorm
        if(MA){
            DESeq2::plotMA(resNorm,ylim=c(-5,5), main=compared)
        }
    }
    saveRDS(res_shrink, filename = "dds_results_shrink.rds")
}

dev.off()