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

#providing varables
geneCts <- "/Users/sreenivasan/Documents/Works/scripts/nf-core_rnaseq/deseq2/example_data/salmon_count_data.tsv"
sampleAnno <- "/Users/sreenivasan/Documents/Works/scripts/nf-core_rnaseq/deseq2/example_data/metadata.csv"
cts <- read.csv(geneCts,sep="\t",row.names=1) 
coldata <- read.csv(sampleAnno, row.names=1) %>%
  mutate(group = sub(group, pattern = ":", replacement = "_")) %>%
  mutate(group = relevel(factor(group), ref = "Oatp1c1_KO")) %>%
  mutate(sex = factor(sex)) %>%
  select(sex, group)
cts <- cts[,-1]
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))

#run DESeq2 without anything
dds <- DESeqDataSetFromMatrix(countData = round(cts),colData = coldata,  design = ~ sex + group)
dds <- DESeq(dds)

#Pre-filtering
smallestGroupSize <- 10
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

#----------------------------------------------------------
#   file1 - OUTLIER REMOVAL - disp_n_rle_plots_before_vs_after_removal.pdf
#----------------------------------------------------------


#check dispersion estimate
pdf("disp_n_rle_plots_before_vs_after_removal.pdf")
plotDispEsts(dds)

# Figuring out sample-level outliers using plotRLE
check_n_remove_outliers <- function(cts = cts, coldata = coldata, color_df = color_df) {
  SeqES <- newSeqExpressionSet(round(as.matrix(cts), 0), phenoData = coldata)
  distinct_colors_df <- color_df %>% distinct(group, colors)
  RLE <- plotRLE(SeqES, outline = FALSE, ylim = c(-4, 4), col = color_df$colors)
  legend("top", legend = distinct_colors_df$group, pch=16,
    col=distinct_colors_df$colors, cex=.8, ncol=2, title="group")
  rle_stddevs <- apply(RLE, FUN = sd, MARGIN = 2) #by columns (MARGIN = 2)
  outliers <- rle_stddevs[rle_stddevs > 1] %>%
    names()
  to_keep <- rle_stddevs[rle_stddevs < 1] %>%
    names()
  return(SeqES[, to_keep])
}

colors <- brewer.pal(4, "Set2")
color_df = data.frame(group = coldata$group, colors = colors[coldata$group])
SeqES <- check_n_remove_outliers(cts, coldata, color_df)
dds <- DESeqDataSetFromMatrix(countData = counts(SeqES), colData = data.frame(phenoData(SeqES)@data), design = ~ sex + group)

#Pre-filtering
smallestGroupSize <- 10
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

#run DESeq2 without SVA
dds <- DESeq(dds)

#check dispersion estimate
plotDispEsts(dds)

# Checking outlier removal
set <- newSeqExpressionSet(counts(dds), phenoData = data.frame(colData(dds)))
color_df = data.frame(group = data.frame(colData(dds))$group, colors = colors[data.frame(colData(dds))$group])
distinct_colors_df <- color_df %>% distinct(group, colors)
plotRLE(set, outline = FALSE, ylim = c(-4, 4), col = color_df$colors)
legend("top", legend = distinct_colors_df$group, pch=16, col=distinct_colors_df$colors, cex=.8, ncol=2, title="group")

dev.off()



#----------------------------------------------------------
#   file2 - plots_gpca_pvalHist_ruv_outliersRemoved.pdf
#----------------------------------------------------------

# To not relying on such variance stabilising transfromations, e.g., VST or rlog, use methods that work directly on count data
pdf("plots_gpca_pvalHist_ruv_outliersRemoved.pdf")
#, width= 10, heigth)
gpca <- glmpca(counts(dds), L = 2) # nolint: error.
gpca.dat <- gpca$factors
gpca.dat$sex <- dds$sex
gpca.dat$group <- dds$group
#dev.new()
ggplot(gpca.dat, aes(x = dim1, y = dim2, color = group, shape = sex)) +
  geom_point(size = 3) +
  coord_fixed() +
  ggtitle("glmpca - Generalized PCA")

# Finally run the differential expression analysis

par(mfrow=c(1,3), mar=c(4,4,2,1))

res1 <- results(dds,contrast = c("group", "Oatp1c1_KO", "C57BL_6J"))
hist(res1$padj, main = "C57BL_6J_vs_Oatp1c1_KO")

res2 <- results(dds,contrast=c("group", "Oatp1c1_KO", "Mct8_Oatp1c1_DKO"))
hist(res2$padj, main="Mct8_Oatp1c1_DKO_vs_Oatp1c1_KO")

res3 <- results(dds,contrast=c("group", "Oatp1c1_KO", "Mct8_Oatp1c1_DKO_AAV_Mct8"))
hist(res3$padj, main="Mct8_Oatp1c1_DKO_AAV_Mct8_vs_Oatp1c1_KO")

sum(res1$padj < 0.1, na.rm=TRUE)
sum(res2$padj < 0.1, na.rm=TRUE)
sum(res3$padj < 0.1, na.rm=TRUE)

# RUV analysis
par(mfrow=c(1,1), mar=c(4,4,2,1))
set <- newSeqExpressionSet(counts(dds))
idx <- rowSums(counts(set) > 5) >= 2
set <- set[idx, ]
set <- betweenLaneNormalization(set, which = "upper")
not_sig <- intersect(rownames(res2)[which(res2$pvalue > .9)],
  rownames(res1)[which(res1$pvalue > .5)]) %>%
  intersect(., rownames(res3)[which(res3$pvalue > .9)])
empirical <- rownames(set)[rownames(set) %in% not_sig]
set1 <- RUVg(set, empirical, k = 2)
plotPCA(set, col=color_df$colors)
legend("top", legend = distinct_colors_df$group, pch=16, col=distinct_colors_df$colors, cex=.8, ncol=2, title="group")
pData(set)

ddsruv1 <- dds
ddsruv1$W1 <- set$W_1
ddsruv1$W2 <- set$W_2
design(ddsruv1) <- ~ W1 + W2 + sex + group
ddsruv1 <- DESeq(ddsruv1)
res_ruv1 <- results(ddsruv1, contrast = c("group", "Oatp1c1_KO", "Mct8_Oatp1c1_DKO"))
hist(res_ruv1$padj, main="Mct8_Oatp1c1_DKO_vs_Oatp1c1_KO")

dev.off()



#----------------------------------------------------------
#   file3 - excel tables
#----------------------------------------------------------

#excel table with all genes
res1 <- as.data.frame(res1)
res1 <- res1 %>%
    mutate(gene_id = c(rownames(res1))) %>%
    relocate(gene_id, .before=1)
write_xlsx(res1, "geneList_KO_vs_WT.xlsx", col_names = TRUE)

res2 <- as.data.frame(res2)
res2 <- res2 %>%
    mutate(gene_id = c(rownames(res2))) %>%
    relocate(gene_id, .before=1)
write_xlsx(res2, "geneList_DKO_vs_WT.xlsx", col_names = TRUE)

res3 <- as.data.frame(res3)
res3 <- res3 %>%
    mutate(gene_id = c(rownames(res3))) %>%
    relocate(gene_id, .before=1)
write_xlsx(res3, "geneList_AAV_vs_WT.xlsx", col_names = TRUE)

res4 <- as.data.frame(res4)
res4 <- res4 %>%
    mutate(gene_id = c(rownames(res4))) %>%
    relocate(gene_id, .before=1)
write_xlsx(res4, "geneList_DKO_vs_AAV.xlsx", col_names = TRUE)

#excel table with DEGS
DEGs_res1 <- subset(res1, res1$padj<0.05 & abs(log2FoldChange) > 1)
DEGs_res1 <- as.data.frame(DEGs_res1)
DEGs_res1 <- DEGs_res1 %>%
    mutate(gene_id = c(rownames(DEGs_res1))) %>%
    relocate(gene_id, .before=1)
write_xlsx(DEGs_res1, "DEGsList_KO_vs_WT.xlsx")

DEGs_res2 <- subset(res2, res2$padj<0.05 & abs(log2FoldChange) > 1)
DEGs_res2 <- as.data.frame(DEGs_res2)
DEGs_res2 <- DEGs_res2 %>%
    mutate(gene_id = c(rownames(DEGs_res2))) %>%
    relocate(gene_id, .before=1)
write_xlsx(DEGs_res2, "DEGsList_DKO_vs_WT.xlsx")

DEGs_res3 <- subset(res3, res3$padj<0.05 & abs(log2FoldChange) > 1)
DEGs_res3 <- as.data.frame(DEGs_res3)
DEGs_res3 <- DEGs_res3 %>%
    mutate(gene_id = c(rownames(DEGs_res3))) %>%
    relocate(gene_id, .before=1)
write_xlsx(DEGs_res3, "DEGsList_AAV_vs_WT.xlsx")

DEGs_res4 <- subset(res4, res4$padj<0.05 & abs(log2FoldChange) > 1)
DEGs_res4 <- as.data.frame(DEGs_res4)
DEGs_res4 <- DEGs_res4 %>%
    mutate(gene_id = c(rownames(DEGs_res4))) %>%
    relocate(gene_id, .before=1)
write_xlsx(DEGs_res4, "DEGsList_DKO_vs_AAV.xlsx")

dev.off()




#----------------------------------------------------------
#   file4 - other_plots_outliersRemoved.pdf
#----------------------------------------------------------

pdf("other_plots_outliersRemoved.pdf")
#ploting counts for gene with min p-value
plotCounts(dds, gene=which.min(res1$padj), intgroup="group")
plotCounts(dds, gene=which.min(res2$padj), intgroup="group")
plotCounts(dds, gene=which.min(res3$padj), intgroup="group")

#transformations vst, rlog, norm
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
ntd <- normTransform(dds)

#sd plotting
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

#heatmap
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("group","sex")])
pheatmap(assay(ntd)[select,],cluster_rows=FALSE, show_rownames=TRUE,cluster_cols=FALSE, annotation_col=df)
pheatmap(assay(vsd)[select,],cluster_rows=FALSE, show_rownames=TRUE,cluster_cols=FALSE, annotation_col=df)
pheatmap(assay(rld)[select,],cluster_rows=FALSE, show_rownames=TRUE,cluster_cols=FALSE, annotation_col=df)

sampleDists <- dist(t(assay(vsd)))

#pca plot
pcaData <- plotPCA(vsd, intgroup=c("group", "sex"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=dds$group, shape=sex)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()



dev.off()


#----------------------------------------------------------
#   STOP - does thing don't work
#----------------------------------------------------------

#----------------------------------------------------------
#   file3 - LFC shrinkage - lfcShrink_outliersRemoved.pdf
#----------------------------------------------------------
pdf("lfcShrink_outliersRemoved.pdf")
#shrinkage
res1LFC <- lfcShrink(dds, coef="group_Oatp1c1_KO_vs_C57BL_6J" , type="apeglm")
res1Norm <- lfcShrink(dds, coef="group_Oatp1c1_KO_vs_C57BL_6J" , type="normal")
res1Ash <- lfcShrink(dds, coef="group_Oatp1c1_KO_vs_C57BL_6J" , type="ashr")
res2LFC <- lfcShrink(dds, coef="group_Mct8_Oatp1c1_DKO_vs_C57BL_6J", type="apeglm")
res2Norm <- lfcShrink(dds, coef="group_Mct8_Oatp1c1_DKO_vs_C57BL_6J", type="normal")
res2Ash <- lfcShrink(dds, coef="group_Mct8_Oatp1c1_DKO_vs_C57BL_6J", type="ashr")
res3LFC <- lfcShrink(dds, coef="group_Mct8_Oatp1c1_DKO_AAV_Mct8_vs_C57BL_6J", type="apeglm")
res3Norm <- lfcShrink(dds, coef="group_Mct8_Oatp1c1_DKO_AAV_Mct8_vs_C57BL_6J", type="normal")
res3Ash <- lfcShrink(dds, coef="group_Mct8_Oatp1c1_DKO_AAV_Mct8_vs_C57BL_6J", type="ashr")


#the function plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet
par(mfrow=c(1,2), mar=c(4,4,2,1))
plotMA(res1, ylim=c(-5,5), main="Oatp1c1_KO_vs_C57BL_6J")
plotMA(res1LFC, ylim=c(-5,5), main="Oatp1c1_KO_vs_C57BL_6J")
plotMA(res2, ylim=c(-5,5), main="Mct8_Oatp1c1_DKO_vs_C57BL_6J")
plotMA(res2LFC, ylim=c(-5,5), main="Mct8_Oatp1c1_DKO_vs_C57BL_6J")
plotMA(res3, ylim=c(-5,5), main="Mct8_Oatp1c1_DKO_AAV_Mct8_vs_C57BL_6J")
plotMA(res3LFC, ylim=c(-5,5), main="Mct8_Oatp1c1_DKO_AAV_Mct8_vs_C57BL_6J")



#----------------------------------------------------------
#   doesn't work
#----------------------------------------------------------

#heatmap
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rownames(coldata))
colnames(sampleDistMatrix) <- paste(rownames(coldata))
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)



