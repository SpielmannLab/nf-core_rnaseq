library(DESeq2)
suppressPackageStartupMessages(library(dplyr))
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu)
library(glmpca)
library(apeglm)
library(sva)
library(RUVSeq)
library(EDASeq)
library(KernSmooth)

# Load data
geneCts <- "/Users/sreenivasan/Documents/Works/scripts/nf-core_rnaseq/deseq2/example_data/salmon_count_data.tsv"
sampleAnno <- "/Users/sreenivasan/Documents/Works/scripts/nf-core_rnaseq/deseq2/example_data/metadata.csv" # samplesheet Annotations
cts <- read.csv(geneCts, sep = "\t", row.names = 1)
coldata <- read.csv(sampleAnno, row.names = 1) %>%
  mutate(group = sub(group, pattern = ":", replacement = "_")) %>%
  mutate(group = factor(group)) %>%
  mutate(sex = factor(sex)) %>%
  select(-condition)
cts <- cts[, -1]
cts <- cts[, rownames(coldata)]
if (!identical(rownames(coldata), colnames(cts))) {
  stop(message = "The sample names are not identical across the metadata and count matrix")
}

# Figuring out sample-level outliers using plotRLEO
check_n_remove_outliers <- function(cts = cts, coldata = coldata) {
  SeqES <- newSeqExpressionSet(round(as.matrix(cts), 0), phenoData = coldata)
  colors <- brewer.pal(4, "Set2")
  RLE <- plotRLE(SeqES, outline = FALSE, ylim = c(-4, 4), col = colors[coldata$group])

  rle_stddevs <- apply(RLE, FUN = sd, MARGIN = 2)
  outliers <- rle_stddevs[rle_stddevs > 1] %>%
    names()
  to_keep <- rle_stddevs[rle_stddevs < 1] %>%
    names()
  return(SeqES[, to_keep])
}
SeqES <- check_n_remove_outliers(cts, coldata)
dds <- DESeqDataSetFromMatrix(countData = counts(SeqES), colData = data.frame(phenoData(SeqES)@data), design = ~ sex + group)

# Remove lowly sequenced genes
smallestGroupSize <- 10
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep, ]

# Run first pass DESeq without SVA
dds <- DESeq(dds)
# Check dispersion estimates
plotDispEsts(dds)

# Checking outlier removal
set <- newSeqExpressionSet(counts(dds), phenoData = data.frame(colData(dds)))
colors <- brewer.pal(4, "Set2")
plotRLE(set, outline = FALSE, ylim = c(-4, 4), col = colors[dds$group])

# Make PCA plot using variance stabilising transformation
# First do variance stabilising transformation for PCA plots and other exploratory plots. This will not be used for differential analysis
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = c("sex", "group"))
# To not relying on such variance stabilising transfromations, e.g., VST or rlog, use methods that work directly on count data
gpca <- glmpca(counts(dds), L = 2) # nolint: error.
gpca.dat <- gpca$factors
gpca.dat$sex <- dds$sex
gpca.dat$group <- dds$group
dev.new()
ggplot(gpca.dat, aes(x = dim1, y = dim2, color = group, shape = sex)) +
  geom_point(size = 3) +
  coord_fixed() +
  ggtitle("glmpca - Generalized PCA")

# Finally run the differential expression analysis
res <- results(dds, contrast = c("group", "Oatp1c1_KO", "C57BL_6J"))
hist(res_ruv$padj)


# RUV analysis
set <- newSeqExpressionSet(counts(dds))
idx <- rowSums(counts(set) > 5) >= 2
set <- set[idx, ]
set <- betweenLaneNormalization(set, which = "upper")
not.sig <- rownames(res)[which(res$pvalue > .5)]
empirical <- rownames(set)[rownames(set) %in% not.sig]
set <- RUVg(set, empirical, k = 2)
plotPCA(set, col = colors[dds$group])
pData(set)

ddsruv <- dds
ddsruv$W1 <- set$W_1
ddsruv$W2 <- set$W_2
design(ddsruv) <- ~ W1 + W2 + sex + group
ddsruv <- DESeq(ddsruv)
res_ruv <- results(ddsruv, contrast = c("group", "Oatp1c1_KO", "C57BL_6J"))
hist(res_ruv$padj)

# ---------------------------------------------------
# Run SVA --- BUT FOR SOME reason this does not seem to improve the results in this particular dataset
# ---------------------------------------------------
dat <- counts(dds, normalized = TRUE)
idx <- rowMeans(dat) > 1
dat <- dat[idx, ]
mm <- model.matrix(~group, colData(dds))
mm0 <- model.matrix(~1, colData(dds))
svseq <- svaseq(dat, mod = mm, mod0 = mm0, n.sv = 2)

# Change the design to incorporate surrogate variables from SVA analysis
ddssva <- dds
ddssva$SV1 <- svseq$sv[, 1]
ddssva$SV2 <- svseq$sv[, 2]
design(ddssva) <- ~ sex + SV1 + SV2 + group
ddssva <- DESeq(ddssva)

# Check the coefficients
resultsNames(ddssva)
res1 <- results(ddssva, contrast = c("group", "Oatp1c1_KO", "C57BL_6J"))
res1LFC <- lfcShrink(ddssva, coef = "group_Oatp1c1_KO_vs_C57BL_6J", type = "apeglm")
res3 <- results(ddssva, contrast = c("group", "C57BL_6J", "Mct8_Oatp1c1_DKO"))
res2 <- results(ddssva, contrast = c("group", "C57BL_6J", "Mct8_Oatp1c1_DKO_AAV_Mct8"))

plotDispEsts(ddssva)
