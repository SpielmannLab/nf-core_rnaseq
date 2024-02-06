"Generate a DDS object from count matrix and sample metadata (coldata) files. Also perform outlier removal based on RLE using the RUV package.

Usage: build_deseq2_object.R --count_data=<value>

Options:
    -h --help               	Show this screen.
    --conditions_compared=<value>  
    --pValue=<value>              Threshold for p value to determine DEG
" -> doc

# loading libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(writexl))

## Getting in arguments
arguments <- docopt(doc, quoted_args = TRUE)
compared <- arguments$conditions_compared
pVal <- arguments$pValue

message("conditions_compared: ", compared)
message("pValue: ", pVal)


###############
# Main script
###############
dds <- readRDS("dds_obj.rds")

# naming pdf file
pdf_pValue_plot <- sprintf("pVal_plot_%s.pdf", compared)
pdf(pdf_pValue_plot)

# actual comparison
res <- results(dds)
pVal_res <- print(sprintf("DEGs with p-value < %g of %s: %d", pVal,compared,sum(res$padj < pVal, na.rm=TRUE)))
write(pVal_res, append =TRUE, file = "countDEGs.txt")
hist(res$padj, main=conditions_compared)

dev.off()

saveRDS(res, filename "dds_result.rds")
