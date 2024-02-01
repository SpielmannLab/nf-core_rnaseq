
if(shrink_LFC_using || MA){
#the function plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet
pdf_MA <- sprintf("MA_plot_%s.pdf", reference_condition)
pdf(pdf_MA)

if(MA){
    xlim <- c(1,1e5); ylim <- c(-3,3)
    DESeq2::plotMA(res,ylim=c(-5,5), main=compared)
}

# --- shrinkage
xlim <- c(1,1e5); ylim <- c(-3,3)

if(shrink_LFC_using == "apeglm"){
    resLFC <- lfcShrink(dds, coef=compared, type="apeglm")
    if(MA){
        DESeq2::plotMA(resLFC,ylim=c(-5,5), main=compared)
    }
}

if(shrink_LFC_using == "ashr"){
    resAsh <- lfcShrink(dds, coef=compared, type="ashr")
    if(MA){
        DESeq2::plotMA(resAsh,ylim=c(-5,5), main=compared)
    }
}

if(shrink_LFC_using == "normal"){
    resNorm <- lfcShrink(dds, coef=compared, type="normal")
    if(MA){
        DESeq2::plotMA(resNorm,ylim=c(-5,5), main=compared)
    }
}

dev.off()
}


# --- DEG
#excel table with all genes
if(list_all_genes){
  res_df <- as.data.frame(res)
  res_df <- res_df %>%
      mutate(gene_id = c(rownames(res_df))) %>%
      relocate(gene_id, .before=1)
  xlsx_allGenes <- sprintf("geneList_%s_%s_%s.xlsx", contrast[2],"vs", contrast[3])
  write_xlsx(res_df, xlsx_allGenes, col_names = TRUE)
}

#excel table with DEGS
if(list_DEG){
  DEGs_res <- subset(res, res$padj<pVal & abs(log2FoldChange) > 1)
  DEGs_res <- as.data.frame(DEGs_res)
  DEGs_res <- DEGs_res %>%
      mutate(gene_id = c(rownames(DEGs_res))) %>%
      relocate(gene_id, .before=1)
  xlsx_DEG <- sprintf("DEG_list_%s_%s_%s.xlsx", contrast[2],"vs", contrast[3])
  write_xlsx(DEGs_res, xlsx_DEG)
}




# --- stabilization 
if(perform_variance_stabilisation != FALSE){
  #naming pdf file
  pdf_transform <- sprintf("tansform_plot_%s.pdf", reference_condition)
  pdf(pdf_transform)
  blind <- do_blind_stabilization 
  #transformations vst, rlog, norm
  #sd plotting

  select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]
  df <- as.data.frame(colData(dds)[,c("group","sex")])
  sampleDists <- dist(t(assay(dds)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(dds$group)
  colnames(sampleDistMatrix) <- paste(dds$group)

#vst
  if(perform_variance_stabilisation == "vst"){
    vsd <- vst(dds, blind=blind)
    if(mean_vs_sd_after_stabilisation){
      meanSdPlot(assay(vsd))
    }
    if(plot_sample_clustering){
      pdf_heatmap <- sprintf("heatmap_plot_%s.pdf", reference_condition)
      pdf(pdf_heatmap)
      pheatmap(assay(vsd)[select,],cluster_rows=FALSE, show_rownames=TRUE,cluster_cols=TRUE, annotation_col=df)
      dev.off()
    }
  }

# rlog
  if(perform_variance_stabilisation == "rlog"){
    rld <- rlog(dds, blind=blind)
    if(mean_vs_sd_after_stabilisation){
      meanSdPlot(assay(rld))
    }
    if(plot_sample_clustering){
      pdf_heatmap <- sprintf("heatmap_plot_%s.pdf", reference_condition)
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
      pdf_heatmap <- sprintf("heatmap_plot_%s.pdf", reference_condition)
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
  pdf_heatmap <- sprintf("heatmap_plot_%s.pdf", reference_condition)
  pdf(pdf_heatmap)
  pheatmap(assay(dds)[select,],cluster_rows=FALSE, show_rownames=TRUE,cluster_cols=TRUE, annotation_col=df)
  pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists)
  dev.off()
}



# --- PCA
if(PCA){
  #nameing pdf file
  pdf_gpca <- sprintf("gpca_plot_%s.pdf", reference_condition)
  pdf(pdf_gpca)

  smallestGroupSize <- 10
  keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
  dds_filt <- dds[keep,]

  #run DESeq2 without SVA
  dds_filt <- DESeq(dds_filt)

  gpca <- glmpca(counts(dds_filt), L = 2) # nolint: error.
  gpca.dat <- gpca$factors
  gpca.dat$sex <- dds_filt$sex
  gpca.dat$group <- dds_filt$group
  ggplot(gpca.dat, aes(x = dim1, y = dim2, color = group, shape = sex)) +
    geom_point(size = 3, color = colors[dds_filt$group]) +
    coord_fixed() +
    ggtitle("glmpca - Generalized PCA")

  dev.off()
}


if(counts){
# --- counts plotif()
#nameing pdf file
pdf_geneCount <- sprintf("geneCount_plot_%s.pdf", reference_condition)
pdf(pdf_geneCount)

#specific genes
genes_length <- length(genes_of_interest)
for(i in 1:genes_length){
  gene_counts <- plotCounts(dds, gene=genes_of_interest[i], intgroup="group", returnData = TRUE)
  genes_df <- gene_counts
  mean <- mean(genes_df$count)
  mean_groups <- genes_df %>%
    group_by(group) %>%
    summarize(counts = mean(count))
  sd <- sd(genes_df$count)
  plot <- ggplot(data = mean_groups, aes(x=group, y=counts)) + 
    geom_segment(aes(xend = group,y = 0, yend = counts, color = group), size = 10, arrow.fill = NULL) +
    geom_linerange(aes(
      x = group,
      ymin = mean - sd,
      ymax = mean + sd)) +
    geom_point(data=genes_df, aes(x=group, y=count, color = group), position=position_jitter(w=0.1,h=0)) +
    ggtitle(genes_of_interest[i]) +
    scale_y_continuous(expand = c(0,0)) +
    scale_color_manual(values = colors)
  print(plot)
}

dev.off()
}




