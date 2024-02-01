
# --- Batch correction
if(perform_batch_correction == "RUV"){
  #nameing pdf file
  pdf_ruv_plot <- sprintf("ruv_plot_%s.pdf", reference_condition)
  pdf(pdf_ruv_plot)

  #RUV analysis
  ruv_analysis <- function(dds = dds, reference_condition = reference_condition){
    par(mfrow=c(1,1), mar=c(4,4,2,1))
    set <- newSeqExpressionSet(counts(dds))
    idx <- rowSums(counts(set) > 5) >= 2
    set <- set[idx, ]
    set <- betweenLaneNormalization(set, which = "upper")
    ref <- reference_condition
    not_sig <- rownames(res)[which(res$pvalue > RUV_threshold_not_sig)]
    empirical <- rownames(set)[rownames(set) %in% not_sig]
    sets <- RUVg(set, empirical, k = 2)
    plotPCA(sets, col = colors[dds$group])
    legend("top", legend = col_df$group, pch=16, col= col_df$color, cex=.8, ncol=2, title=comparison_key)
    pData(sets)
    ddsruv <- dds
    ddsruv$W1 <- sets$W_1
    ddsruv$W2 <- sets$W_2
    design(ddsruv) <- ~ W1 + W2 + sex + group
    ddsruv <- DESeq(ddsruv)
    res_ruv <- results(ddsruv, contrast = contrast)
    hist(res_ruv$padj, main=c(contrast[3],"_vs_",ref))
    return(list(ddsruv = ddsruv, res_ruv = res_ruv))
  }

  ruv <- ruv_analysis(dds, reference_condition)
  ddsruv <- ruv$ddsruv
  res_ruv <- ruv$res_ruv

  dev.off()

} else if(perform_batch_correction == "SUV"){

}

