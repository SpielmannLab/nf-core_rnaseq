# --- group comparison
do_blind_stabilization <- TRUE
mean_vs_sd_before_stabilisation <- TRUE
mean_vs_sd_after_stabilisation <- TRUE
plot_sample_clustering <- TRUE
plot_volcano_before_LFCshrink <- TRUE
plot_volcano_after_LFCshrink <- TRUE
shrink_LFC_using <- "apeglm"
plot_p_values <- TRUE
list_all_genes<- TRUE
list_DEG <- TRUE
counts <- TRUE


pVal <- 0.1
genes_of_interest <- c("Slc16a2","Klf9","Hr","Dio2")
# nameing pdf file

pdf_pValue_plot <- sprintf("pVal_plot_%s.pdf", reference_condition)
pdf(pdf_pValue_plot)

## try to convert into functions
group_comparison <- function(dds = dds, reference_condition = reference_condition, alternate_condition = alternate_condition){
  par(mfrow=c(1,3), mar=c(4,4,2,1))
  ref <- reference_condition
  alt <- alternate_condition
  contrast <- c(comparison_key ,ref, alt)
  res <- results(dds,contrast=contrast) 
  pVal_res <- print(sprintf("DEGs with p-value < %g of %s vs %s: %d", pVal ,contrast[2],contrast[3],sum(res$padj < pVal, na.rm=TRUE)))
  write(pVal_res, append =TRUE, file = "countDEGs.txt")
  if(plot_p_values){
    hist(res$padj, main=c(contrast[3],"_vs_",ref))
  }
  return(list(res = res, contrast = contrast, pVal_res = pVal_res))
}

comparison <- group_comparison(dds, reference_condition, alternate_condition)
res <- comparison$res
contrast <- comparison$contrast
pVal_res <- comparison$pVal_res

compared <- sprintf("%s_%s_vs_%s", contrast[1], contrast[3], contrast[2])

dev.off()


