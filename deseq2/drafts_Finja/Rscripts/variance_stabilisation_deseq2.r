"Generate a DDS object from count matrix and sample metadata (coldata) files. Also perform outlier removal based on RLE using the RUV package.

Usage:     Rscript variance_stabilisation_deseq2.r --other_keys=${other_keys} --perform_variance_stabilisation=${perform_variance_stabilisation} --do_blind_stabilization=${do_blind_stabilization} --mean_vs_sd_before_stabilisation=${mean_vs_sd_before_stabilisation} --mean_vs_sd_after_stabilisation=${mean_vs_sd_after_stabilisation}

Options:
    -h --help               	Show this screen.
    --perform_variance_stabilisation=<value>     condition other conditions defined in <alternate_condition> should be compared to
    --do_blind_stabilization=<value>    
    --mean_vs_sd_before_stabilisation=<value>
    --mean_vs_sd_after_stabilisation=<value>
    --plot_sample_clustering=<value>
    --other_keys=<value>      All other batch metadata columns

" -> doc

# loading libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(vsn))

## Getting in arguments
arguments <- docopt(doc, quoted_args = TRUE)
perform_variance_stabilisation <- arguments$perform_variance_stabilisation
do_blind_stabilization <- arguments$do_blind_stabilization
mean_vs_sd_before_stabilisation <- arguments$mean_vs_sd_before_stabilisation
mean_vs_sd_after_stabilisation <- arguments$mean_vs_sd_after_stabilisation
plot_sample_clustering <- arguments$plot_sample_clustering
other_keys <- arguments$other_keys %>%
    strsplit(split = ",") %>%
    unlist()

message("perform_variance_stabilisation: ", perform_variance_stabilisation)
message("do_blind_stabilization: ", do_blind_stabilization)
message("mean_vs_sd_before_stabilisation: ", mean_vs_sd_before_stabilisation)
message("mean_vs_sd_after_stabilisation: ", mean_vs_sd_after_stabilisation)
message("plot_sample_clustering: ", plot_sample_clustering)
message("other_key: ", other_key)

###############
# Main script
###############
dds <- readRDS("dds_obj.rds")

# --- stabilization 
if(perform_variance_stabilisation != FALSE){
  #naming pdf file
  pdf_transform <- sprintf("tansform_plot.pdf")
  pdf(pdf_transform)
  #transformations vst, rlog, norm
  #sd plotting

  select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]
  other <- c()
  dds_group <- as.character(dds$group)
  other <- coldata[[other_keys[1]]] %>%
      as.character()
  df <- as.data.frame(dds_group,other)

  sampleDists <- dist(t(assay(dds)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(dds$group)
  colnames(sampleDistMatrix) <- paste(dds$group)

#vst
  if(perform_variance_stabilisation == "vst"){
    vsd <- vst(dds, blind=do_blind_stabilization)
    if(mean_vs_sd_after_stabilisation){
      meanSdPlot(assay(vsd))
    }
    if(plot_sample_clustering){
      pdf_heatmap <- sprintf("heatmap_plot.pdf")
      pdf(pdf_heatmap)
      pheatmap(assay(vsd)[select,],cluster_rows=FALSE, show_rownames=TRUE,cluster_cols=TRUE, annotation_col=df)
      dev.off()
    }
  }

# rlog
  if(perform_variance_stabilisation == "rlog"){
    rld <- rlog(dds, blind=do_blind_stabilization)
    if(mean_vs_sd_after_stabilisation){
      meanSdPlot(assay(rld))
    }
    if(plot_sample_clustering){
      pdf_heatmap <- sprintf("heatmap_plot.pdf")
      pdf(pdf_heatmap)
      pheatmap(assay(rld)[select,],cluster_rows=FALSE, show_rownames=TRUE,cluster_cols=TRUE, annotation_col=df)
      dev.off()
    }
  }

#norm
  if(perform_variance_stabilisation == "norm"){
    ntd <- normTransform(dds)
    if(mean_vs_sd_after_stabilisation){
      meanSdPlot(assay(ntd))
    }
    if(plot_sample_clustering){
      pdf_heatmap <- sprintf("heatmap_plot.pdf")
      pdf(pdf_heatmap)
      pheatmap(assay(ntd)[select,],cluster_rows=FALSE, show_rownames=TRUE,cluster_cols=TRUE, annotation_col=df)
      dev.off()
    }
  }
  dev.off()
} else if(!perform_variance_stabilisation){
    mean_vs_sd_after_stabilisation <- FALSE
    if(mean_vs_sd_before_stabilisation){
      meanSdPlot(assay(dds))
    }
  dev.off()
}
if(plot_sample_clustering){
  pdf_heatmap <- sprintf("heatmap_plot_%s.pdf")
  pdf(pdf_heatmap)
  pheatmap(assay(dds)[select,],cluster_rows=FALSE, show_rownames=TRUE,cluster_cols=TRUE, annotation_col=df)
  pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists)
  dev.off()
}

