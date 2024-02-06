

#loading libraries
library(DESeq2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(glmpca)
library(apeglm)
library(sva)
library(RUVSeq)
library(ggplot2)
library(vsn)
library(writexl)


geneCts <- "/data/humangen_external/test_area/sygo/project37/salmon.merged.gene_counts_P037_2023.tsv"
sampleAnno <- "/data/humangen_external/test_area/sygo/project37/sheets/annotations_sex_P037_2023.csv"
reference_condition <- "Mct8_Oatp1c1_DKO_AAV_Mct8"
alternate_condition <- "Oatp1c1_KO"
conditions_compared <- "Mct8_Oatp1c1_DKO_AAV_Mct8_vs_Oatp1c1_KO"
#"C57BL_6J" | "Oatp1c1_KO" | "Mct8_Oatp1c1_DKO" | "Mct8_Oatp1c1_DKO_AAV_Mct8"

design <- ~ sex + group
comparison_key <- "group"

detect_sample_outliers <- TRUE # options: <true|false>
RLE_threshold_max <- 1 
RUV_threshold_not_sig <- 0.5
perform_batch_correction <- "RUV"
PCA <- TRUE
MA <- TRUE
dispersion <- TRUE
perform_variance_stabilisation <- "vst"
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

#Pre-filtering
smallestGroupSize <- 3     # -------------- !!!!!!
minCounts <- 10     # -------------- !!!!!!




# --- forming count and coldata 
group <- factor(comparison_key)
cts <- read.csv(geneCts,sep="\t",row.names=1) 
coldata <- read.csv(sampleAnno, row.names=1) %>%
  mutate(group = sub(group, pattern = ":", replacement = "_")) %>%
  mutate(group = relevel(factor(group), ref = reference_condition)) 
if(other_keys){
  for(i in 1:length(other_keys)){
    coldata <- mutate(other_keys[i] <- factor(coldata$other_keys[i]))
  }
  coldata %>% select(group, other_keys[i])
}
cts <- cts[,-1]
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts)) 

design <- ~ group + other_keys[i]

# --- setting the colors
colors <- brewer.pal(4, "Set2")
colors_df <- data.frame(group = coldata$group, color = colors[coldata$group])
col_df <- colors_df %>% distinct(group, color)


# --- run DESeq2 without anything
dds <- DESeqDataSetFromMatrix(countData = round(cts), colData = coldata, design = design)

# --- Pre-filtering
keep <- rowSums(counts(dds) >= minCounts) >= smallestGroupSize
dds <- dds[keep,]

dds <- DESeq(dds)

# --- detection and removing of outliers
if(detect_sample_outliers){
  # nameing pdf file
  pdf_dispersion_plot <- sprintf("dispersion_plot_%s.pdf", reference_condition)
  pdf(pdf_dispersion_plot)

  # Figuring out sample-level outliers using plotRLE
  check_n_remove_outliers <- function(cts = cts, coldata = coldata, color = colors_df) {
    SeqES <- newSeqExpressionSet(round(as.matrix(cts), 0), phenoData = coldata)
    RLE <- plotRLE(SeqES, outline=FALSE, ylim=c(-4, 4), col = colors_df$color)
    legend("top", legend = col_df$group, pch=16, col= col_df$color, cex=.8, ncol=2, title="group")
    rle_stddevs <- apply(RLE, FUN = sd, MARGIN = 2) 
    outliers <- rle_stddevs[rle_stddevs > RLE_threshold_max] %>%
      names()
    to_keep <- rle_stddevs[rle_stddevs < RLE_threshold_max] %>%
      names()
    return(SeqES[, to_keep])
  }

  # using the created function check_n_remove_outliers
  SeqES <- check_n_remove_outliers(cts, coldata)
  dds <- DESeqDataSetFromMatrix(countData = counts(SeqES), colData = data.frame(phenoData(SeqES)@data), design = ~ sex + group)

  # run DESeq2 without SVA
  dds <- DESeq(dds)

  # check dispersion estimate
  if(dispersion){
    plotDispEsts(dds)
  }

  # checking outlier removal
  set <- newSeqExpressionSet(counts(dds), phenoData = data.frame(colData(dds)))
  RLE <- plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors_df$color)
  legend("top", legend = col_df$group, pch=16, col = col_df$color, cex=.8, ncol=2, title="group")

  dev.off()
}

#dds$group <- droplevels(dds$group)

# --- group comparison
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

