"Generate a DDS object from count matrix and sample metadata (coldata) files. Also perform outlier removal based on RLE using the RUV package.

Usage:     Rscript PCAplot_deseq2.r --other_keys=${other_keys} --PCA=${PCA} --minCounts=${minCounts} --smallestGroupSize=${smallestGroupSize}

Options:
    -h --help               	Show this screen.
    --other_keys=<value>      All other batch metadata columns
    --PCA=<value>      All other batch metadata columns
    --minCounts=<value>        Minimum number of total counts that should be present in at least <smallestGroupSize> number of samples
    --smallestGroupSize=<value>        Minimum number of sammples which should have a the <minCounts> number of total counts
" -> doc

# loading libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(glmpca))

## Getting in arguments
arguments <- docopt(doc, quoted_args = TRUE)
other_keys <- arguments$other_keys %>%
    strsplit(split = ",") %>%
    unlist()
PCA <- arguments$PCA
minCounts <- arguments$minCounts
smallestGroupSize <- arguments$smallestGroupSize

message("other_key: ", other_key)
message("PCA: ", PCA)
message("minCounts: ", minCounts)
message("smallestGroupSize: ", smallestGroupSize)


###############
# Main script
###############
dds <- readRDS("dds_obj.rds")

# --- PCA
if(PCA){
  #nameing pdf file
  pdf_gpca <- sprintf("gpca_plot.pdf")
  pdf(pdf_gpca)

  # --- Pre-filtering
  keep <- rowSums(counts(dds) >= minCounts) >= smallestGroupSize
  dds <- dds[keep, ]

  gpca <- glmpca(counts(dds), L = 2) # nolint: error.
  gpca.dat <- gpca$factors
  gpca.dat$group <- dds$group


# doen't work
if(exists("other_keys")){
  other <- coldata[[other_keys[1]]]
  ggplot(gpca.dat, aes(x = dim1, y = dim2, color = group, shape = other)) +
    geom_point(size = 3, color = colors[dds$group]) +
    coord_fixed() +
    ggtitle("glmpca - Generalized PCA")

  dev.off()
}
} else{
  ggplot(gpca.dat, aes(x = dim1, y = dim2, color = group)) +
    geom_point(size = 3, color = colors[dds$group]) +
    coord_fixed() +
    ggtitle("glmpca - Generalized PCA")

  dev.off()
}

