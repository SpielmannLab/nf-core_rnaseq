"Generate a DDS object from count matrix and sample metadata (coldata) files. Also perform outlier removal based on RLE using the RUV package.

Usage: build_deseq2_object.R --count_data=<value>

Options:
    -h --help               	Show this screen.
    --comparison_key=<value>      The metadata column that will be used to perform DE analysis
    --genes_of_interest=<value>      genes which counts should be plotted
    --counts=<value>
" -> doc

# loading libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))

## Getting in arguments
arguments <- docopt(doc, quoted_args = TRUE)
genes_of_interest <- arguments$genes_of_interest %>%
    strsplit(split = ",") %>%
    unlist()
group <- arguments$comparison_key
counts <- arguments$counts

message("genes_of_interest: ", genes_of_interest)
message("comparison_key: ", group)
message("genes_of_interest: ", counts)


###############
# Main script
###############
dds <- readRDS("dds_obj.rds")

if(counts){
# --- counts plotif()
#nameing pdf file
pdf_geneCount <- sprintf("geneCount_plot.pdf")
pdf(pdf_geneCount)

#specific genes
genes_length <- length(genes_of_interest)
for(i in 1:genes_length){
  gene_counts <- plotCounts(dds, gene=genes_of_interest[i], intgroup=group, returnData = TRUE)
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

