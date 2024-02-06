"Generate a DDS object from count matrix and sample metadata (coldata) files. Also perform outlier removal based on RLE using the RUV package.

Usage:      Rscript building_deseq2_object.r --count_data=${count_data} --metadata=${metadata} --comparison_key=${comparison_key} --other_keys=${other_keys} --minCounts=${minCounts} --smallestGroupSize=${smallestGroupSize} --detect_sample_outliers=${detect_sample_outliers} --RLE_threshold_max=${RLE_threshold_max}

Options:
    -h --help               	Show this screen.
    --count_data=<value>        TSV file containing the raw count matrix
    --metadata=<value>      TSV file containing the sample metadata. This sheet also needs to have a column with 'color'.
    --comparison_key=<value>      The metadata column that will be used to perform DE analysis
    --other_keys=<value>      All other batch metadata columns
    --minCounts=<value>        Minimum number of total counts that should be present in at least <smallestGroupSize> number of samples
    --smallestGroupSize=<value>        Minimum number of sammples which should have a the <minCounts> number of total counts
    --detect_sample_outliers=<value>         Whether or not to detect and remove sample outliers
    --RLE_threshold_max=<value>      The maximum RLE threshold for outlier removal. Default is 1
" -> doc


# loading libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(RUVSeq))
suppressPackageStartupMessages(library(ggplot2))

## Getting in arguments
arguments <- docopt(doc, quoted_args = TRUE)
count_data <- arguments$count_data
metadata <- arguments$metadata
comparison_key <- arguments$comparison_key
other_keys <- arguments$other_keys %>%
    strsplit(split = ",") %>%
    unlist()
minCounts <- arguments$minCounts
smallestGroupSize <- arguments$smallestGroupSize
detect_sample_outliers <- arguments$detect_sample_outliers
RLE_threshold_max <- arguments$RLE_threshold_max

message("count_data: ", count_data)
message("metadata: ", metadata)
message("comparison_key: ", comparison_key)
message("other_keys: ", other_keys)
message("minCounts: ", minCounts)
message("smallestGroupSize: ", smallestGroupSize)
message("detect_sample_outliers: ", detect_sample_outliers)
message("RLE_threshold_max: ", RLE_threshold_max)

###############
# Define functions
###############
# Figuring out sample-level outliers using plotRLE
check_n_remove_outliers <- function(cts = cts, coldata = coldata, color = colors_df) {
  SeqES <- newSeqExpressionSet(round(as.matrix(cts), 0), phenoData = coldata)
  RLE <- plotRLE(SeqES, outline=FALSE, ylim=c(-4, 4), col = colors_df$color)
  legend("top", legend = comp_col_df, pch=16, col= col_df$color, cex=.8, ncol=2, title=comparison_key)
  rle_stddevs <- apply(RLE, FUN = sd, MARGIN = 2) #by columns (MARGIN = 2)
  outliers <- rle_stddevs[rle_stddevs > 1] %>%
    names()
  to_keep <- rle_stddevs[rle_stddevs < 1] %>%
    names()

  coldata <- coldata[!(rownames(coldata) %in% outliers),]
  #colors_df <- colors_df[!(rownames(coldata) %in% outliers),]
  comp_col <- as.factor(coldata[[comparison_key]])
  #comp_col_df <- col_df[[comparison_key]]

  return(c(SeqES[, to_keep], comp_col, comp_col_df, colors_df, col_df))
}

###############
# Main script
###############
# --- forming count and coldata
cts <- read.csv(count_data, sep = "\t", row.names = 1)
coldata <- read.csv(metadata, row.names = 1) %>%
    select(any_of(c(comparison_key, other_keys)))
# Reordering so that the entries in the metadata and the countdata are in the same order
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))

#following can replace coldata$sex_again
#coldata[[other_keys[2]]]

# --- setting colors
comp_col <- as.factor(coldata[[comparison_key]])

colors <- brewer.pal(length(unique(comp_col)), "Set2")
colors_df <- data.frame(condition = comp_col, color = colors[comp_col])
col_df <- colors_df %>% distinct(condition, color)

comp_col_df <- as.factor(col_df[[comparison_key]])


# --- building DESEq2 DDS object
design <- paste0("~ ", paste(c(comparison_key, other_keys), collapse = " + ")) %>% formula()
dds <- DESeqDataSetFromMatrix(countData = round(cts), colData = coldata, design = design)

# --- Pre-filtering
keep <- rowSums(counts(dds) >= minCounts) >= smallestGroupSize
dds <- dds[keep, ]

dds <- DESeq(dds)

# --- detection and removing of outliers
if (detect_sample_outliers) {
    # nameing pdf file
    pdf("dispersion_plot_before_and_after_outlier_removal.pdf")

    # using the created function check_n_remove_outliers
    SeqES <- check_n_remove_outliers(cts, coldata)
    dds <- DESeqDataSetFromMatrix(countData = counts(SeqES), colData = data.frame(phenoData(SeqES)@data), design = design)

    # run DESeq2 without SVA
    dds <- DESeq(dds)

    # check dispersion estimate
    if (dispersion) {
        plotDispEsts(dds)
    }

    # checking outlier removal
    set <- newSeqExpressionSet(counts(dds), phenoData = data.frame(colData(dds)))
    RLE <- plotRLE(SeqES, outline=FALSE, ylim=c(-4, 4), col = colors_df$color)
    legend("top", legend = comp_col_df, pch=16, col= col_df$color, cex=.8, ncol=2, title=comparison_key)
    
    dev.off()
}

saveRDS(dds, filename = "dds_obj.rds")
saveRDS(design, filename = "dds_design.rds")