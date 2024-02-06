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
suppressPackageStartupMessages(library(RColorBrewer))

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
check_n_remove_outliers <- function(cts = cts, coldata = coldata) {
    SeqES <- newSeqExpressionSet(round(as.matrix(cts), 0), phenoData = coldata)

    RLE_before <- plotRLE(SeqES, outline = FALSE, ylim = c(-4, 4))

    rle_stddevs <- apply(RLE_before, FUN = sd, MARGIN = 2) # by columns (MARGIN = 2)
    outliers <- rle_stddevs[rle_stddevs > 1] %>%
        names()
    to_keep <- rle_stddevs[rle_stddevs < 1] %>%
        names()

    coldata <- coldata[to_keep, ]
    SeqES <- SeqES[, to_keep]

    RLE_after <- plotRLE(SeqES, outline = FALSE, ylim = c(-4, 4))

    return(list(SeqES, RLE_before, RLE_after))
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

# following can replace coldata$sex_again
# coldata[[other_keys[2]]]

# --- adding colors into the dataframe
comp_col <- unique(coldata[[comparison_key]])
colors <- brewer.pal(length(unique(comp_col)), "Set2")
names(colors) <- comp_col
coldata <- coldata %>%
    mutate(color = colors[.data[[comparison_key]]])

# --- building DESEq2 DDS object
design <- paste0("~ ", paste(c(comparison_key, other_keys), collapse = " + ")) %>% formula()
dds <- DESeqDataSetFromMatrix(countData = round(cts), colData = coldata, design = design)

# --- Pre-filtering
keep <- rowSums(counts(dds) >= minCounts) >= smallestGroupSize
dds <- dds[keep, ]

dds <- DESeq(dds)

# --- detection and removing of outliers
if (detect_sample_outliers) {

    # using the created function check_n_remove_outliers
    result <- check_n_remove_outliers(cts, coldata)
    SeqES = result[[1]]
    RLE_before <- result[[2]]
    RLE_after <- result[[3]]

    # naming pdf file
    pdf("dispersion_plot_before_and_after_outlier_removal.pdf")
    # Plot RLE distribution using ggplot2, before outlier removal
    df_rle_scores_before <- RLE_before %>%
        data.frame() %>%
        tidyr::pivot_longer(cols = where(is.numeric)) %>%
        mutate(group = coldata[.data$name, comparison_key])
    plot_before <- ggplot(df_rle_scores_before, aes(x = name, y = value, color = group)) +
        geom_boxplot(outlier.shape = NA) +
        scale_color_manual(breaks = coldata[[comparison_key]], values = coldata$color) +
        ylim(c(-3,3)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        ggtitle("RLE scores before outlier removal")
    df_rle_scores_after <- RLE_after %>%
        data.frame() %>%
        tidyr::pivot_longer(cols = where(is.numeric)) %>%
        mutate(group = coldata[.data$name, comparison_key])
    plot_after <- ggplot(df_rle_scores_after, aes(x = name, y = value, color = group)) +
        geom_boxplot(outlier.shape = NA) +
        scale_color_manual(breaks = coldata[[comparison_key]], values = coldata$color) +
        ylim(c(-3,3)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        ggtitle("RLE scores after outlier removal")
    print(plot_grid(plot_before, plot_after, nrow = 2))

    # Convert the SeqES object back to dds object and perform DESeq2
    dds <- DESeqDataSetFromMatrix(countData = counts(SeqES), colData = data.frame(phenoData(SeqES)@data), design = design)
    dds <- DESeq(dds)

    if (dispersion) {
        plotDispEsts(dds)
    }

    dev.off()
}

saveRDS(dds, filename = "dds_obj.rds")
saveRDS(design, filename = "dds_design.rds")

