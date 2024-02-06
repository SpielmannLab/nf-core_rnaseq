"Generate a DDS object from count matrix and sample metadata (coldata) files. Also perform outlier removal based on RLE using the RUV package.

Usage:     Rscript DEG_detection_deseq2.r --reference_condition=${reference_condition} --alternate_condition=${alternate_condition} --pValue=${pValue} --log2fc_threshold_DEG=${log2fc_threshold_DEG} --list_all_genes=${list_all_genes} --list_DEG=${list_DEG}

Options:
    -h --help               	Show this screen.
    --condition_compared<value>     
    --pValue=<value>              Treshold for p value to determine DEG
    --log2fc_threshold_DEG=<value>              Treshold for log2FoldChange to determine DEG
    --list_all_genes=<value>            
    --list_DEG=<value>              
" -> doc

# loading libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(writexl))

## Getting in arguments
arguments <- docopt(doc, quoted_args = TRUE)
compared <- arguments$condition_compared
pVal <- arguments$pValue
list_all_genes <- arguments$list_all_genes
list_DEG <- arguments$list_DEG
log2fc_threshold_DEG <- arguments$log2fc_threshold_DEG

message("condition_compared: ", compared)
message("pValue: ", pVal)
message("list_all_genes: ", list_all_genes)
message("list_DEG: ", list_DEG)
message("log2fc_threshold_DEG: ", log2fc_threshold_DEG)

###############
# Main script
###############
res <- readRDS("dds_results.rds")

#excel table with all genes
if(list_all_genes){
  res_df <- as.data.frame(res)
  res_df <- res_df %>%
      mutate(gene_id = c(rownames(res_df))) %>%
      relocate(gene_id, .before=1)
  xlsx_allGenes <- sprintf("geneList_%s.xlsx", compared)
  write_xlsx(res_df, xlsx_allGenes, col_names = TRUE)
}

#excel table with DEGS
if(list_DEG){
  DEGs_res <- subset(res, res$padj<pVal & abs(log2FoldChange) > log2fc_threshold_DEG)
  DEGs_res <- as.data.frame(DEGs_res)
  DEGs_res <- DEGs_res %>%
      mutate(gene_id = c(rownames(DEGs_res))) %>%
      relocate(gene_id, .before=1)
  xlsx_DEG <- sprintf("DEG_list_%s.xlsx", compared)
  write_xlsx(DEGs_res, xlsx_DEG)
}
