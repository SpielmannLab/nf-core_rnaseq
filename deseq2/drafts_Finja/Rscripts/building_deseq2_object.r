"Generate a DDS object from count matrix and sample metadata (coldata) files. Also perform outlier removal based on RLE using the RUV package.

Usage: build_deseq2_object.R --count_data=<value>

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
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(sva))
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
    RLE <- plotRLE(SeqES, outline = FALSE, ylim = c(-4, 4), col = colors_df$color)
    legend("top", legend = col_df$group, pch = 16, col = col_df$color, cex = .8, ncol = 2, title = "group")
    rle_stddevs <- apply(RLE, FUN = sd, MARGIN = 2)
    outliers <- rle_stddevs[rle_stddevs > RLE_threshold_max] %>%
        names()
    to_keep <- rle_stddevs[rle_stddevs < RLE_threshold_max] %>%
        names()
    return(SeqES[, to_keep])
}

###############
# Main script
###############
# --- forming count and coldata
group <- comparison_key
cts <- read.csv(count_data, sep = "\t", row.names = 1)
coldata <- read.csv(metadata, row.names = 1) %>%
    select(any_of(c(group, other_keys)))
# Reordering so that the entries in the metadata and the countdata are in the same order
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))


# --- building DESEq2 DDS object
design <- paste0("~ ", paste(c(other_keys, group), collapse = " + ")) %>% formula()
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
    dds <- DESeqDataSetFromMatrix(countData = counts(SeqES), colData = data.frame(phenoData(SeqES)@data), design = ~ sex + group)

    # run DESeq2 without SVA
    dds <- DESeq(dds)

    # check dispersion estimate
    if (dispersion) {
        plotDispEsts(dds)
    }

    # checking outlier removal
    set <- newSeqExpressionSet(counts(dds), phenoData = data.frame(colData(dds)))
    RLE <- plotRLE(set, outline = FALSE, ylim = c(-4, 4), col = colors_df$color)
    legend("top", legend = col_df$group, pch = 16, col = col_df$color, cex = .8, ncol = 2, title = "group")

    dev.off()
}

saveRDS(dds, filename = "DDSobject.rds")
