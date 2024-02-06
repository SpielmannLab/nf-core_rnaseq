"Generate a DDS object from count matrix and sample metadata (coldata) files. Also perform outlier removal based on RLE using the RUV package.

Usage:     Rscript batch_correction_deseq2.r --comparison_key=${comparison_key} --RUV_threshold_not_sig=${RUV_threshold_not_sig} --conditions_compared=${conditions_compared}

Options:
    -h --help               	Show this screen.
    --conditions_compared=<value>     
    --comparison_key=<value>      The metadata column that will be used to perform DE analysis
    --RUV_threshold_not_sig=<value>    
" -> doc

# loading libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(RUVSeq))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))


## Getting in arguments
arguments <- docopt(doc, quoted_args = TRUE)
group <- arguments$comparison_key
RUV_threshold_not_sig <- arguments$RUV_threshold_not_sig
compared <- arguments$conditions_compared


message("comparison_key: ", group)
message("RUV_threshold_not_sig: ", RUV_threshold_not_sig)
message("conditions_compared: ", compared)



###############
# Define functions
###############

#RUV analysis
ruv_analysis <- function(dds = dds, compared = compared){
  par(mfrow=c(1,1), mar=c(4,4,2,1))
  set <- newSeqExpressionSet(counts(dds))
  idx <- rowSums(counts(set) > 5) >= 2
  set <- set[idx, ]
  set <- betweenLaneNormalization(set, which = "upper")
  not_sig <- rownames(res)[which(res$pvalue > RUV_threshold_not_sig)]
  empirical <- rownames(set)[rownames(set) %in% not_sig]
  set <- RUVg(set, empirical, k = 2)
  plotPCA(set, col = colors[dds$group])
  legend("top", legend = col_df$group, pch=16, col= col_df$color, cex=.8, ncol=2, title=group)
  pData(set)
  dds_ruv <- dds
  dds_ruv$W1 <- set$W_1
  dds_ruv$W2 <- set$W_2
  design(dds_ruv) <- ~ W1 + W2 + design
  dds_ruv <- DESeq(dds_ruv)
  res_ruv <- results(dds_ruv)
  hist(res_ruv$padj, main=compared)
  return(list(dds_ruv = dds_ruv, res_ruv = res_ruv))
}

sva_analysis <- function(dds = dds, compared = compared)){
  dds_sva <- estimateSizeFactors(dds)
  norm.cts <- counts(dds_sva, normalized=TRUE)
  mm <- model.matrix(~ other, colData(dds_sva))
  mm0 <- model.matrix(~ 1, colData(dds_sva))
  norm.cts <- norm.cts[rowSums(norm.cts) > 0,]
  fit <- svaseq(norm.cts, mod=mm, mod0=mm0)
  library(rafalib)
  bigpar()
  dds_sva$group.int <- as.integer(dds_sva$group) + 15
  df <- as.data.frame(fit$sv[,1:2])
  ggplot(df, aes(x=V1, y=V2, color = colors[dds$group], shape = other)) +
    geom_point(size = 3, color = colors[dds$group]) +
    coord_fixed() +
    labs(color = group) +
    ggtitle("glmpca - sva PCA")
  dds_sva$V1 <- df$V1
  dds_sva$V2 <- df$V2
  design(dds_sva) <- ~ V1 + V2 + design
  dds_sva <- DESeq(dds_sva)
  res_sva <- results(dds_sva)
  hist(res_ruv$padj, main=compared)
  return(list(dds_ruv = dds_ruv, res_ruv = res_ruv))
}


###############
# Main script
###############
dds <- readRDS("dds_obj.rds")
design <- readRDS("dds_design.rds")
res <- readRDS("dds_results.rds")

# --- Batch correction
if(perform_batch_correction == "RUV"){
  #nameing pdf file
  pdf_ruv_plot <- sprintf("ruv_plot_%s.pdf", compared)
  pdf(pdf_ruv_plot)

  ruv <- ruv_analysis(dds, conditions_compared)
  dds_ruv <- ruv$dds_ruv
  res_ruv <- ruv$res_ruv

  dev.off()

  dds_objects <- list(dds = dds, dds_batch = dds_ruv)
  results <- list(res = res, res_batch = res_ruv)

  saveRDS(dds_objects, filename = "dds_obj.rds")
  saveRDS(results, filename = "dds_results.rds")

} else if(perform_batch_correction == "SUV"){
  # performing of sva batch correction still needs to be implemented
  pdf_sva_plot <- sprintf("sva_plot_%s.pdf", compared)
  pdf(pdf_sva_plot)

  sva <- sva_analysis(dds, compared)
  dds_sva <- sva$dds_sva
  res_sva <- sva$res_sva


  #plot(fit$sv[,1:2], col=other, pch=dds_sva$group.int, cex=2, xlab="SV1", ylab="SV2")
  #legend("top", levels(dds_sva$group), pch=16,
  #  col=1:3, cex=.8, ncol=3, title="other")

  dev.off()

  dds_objects <- list(dds = dds, dds_batch = dds_sva)
  results <- list(res = res, res_batch = res_sva)

  saveRDS(dds_objects, filename = "dds_obj.rds")
  saveRDS(results, filename = "dds_results.rds")
}

