#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(getopt))

version <- "0.0.3"

# Define options
spec <- matrix(c(
  'help'       , 'h', 0, "logical",
  'version'    , 'v', 0, "logical",
  'input_file' , NA , 1, "character",
  'output_file', NA , 1, "character"
), byrow=TRUE, ncol=4)

# Parse
opt <- getopt(spec)

# help
if (!is.null(opt$help)) {
  cat("
Usage:
    Rscript assignment_heatmap.R [options] --input_file <likelihood tsv file from EMU> --output_file <file.png>

Options:
    -h, --help            Show this help message
    -v, --version         Print version and exit
        --input_file      Input file (required)
        --output_file     Output file (required)
\n")
  quit(status=0)
}

# version
if (!is.null(opt$version)) {
  cat("assignment_heatmap.R version", version, "\n")
  quit(status=0)
}

# Validate required arguments
if (is.null(opt$input_file) || is.null(opt$output_file)) {
  cat("Error: --input_file and --output_file are required\n\n")
  cat(getopt(spec, usage=TRUE))
  quit(status=1)
}

input_file  <- opt$input_file
output_file <- opt$output_file


# Read the data
likelihood_df <- read.table(input_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
likelihood_data <- as.data.table(likelihood_df, keep.rownames = "ReadID")

# Melt the data
melted_data <- melt(likelihood_data, id.vars = "ReadID")

# Create the plot
p <- ggplot(melted_data, aes(x = ReadID, y = variable, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(x = "Reads", y = "Taxon name", fill = "Likelihood") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_blank()
  )


# Save to PNG
ggsave(output_file, plot = p, width = 10, height = 8, dpi = 300)

sessionInfo()
cat("Plot saved to", output_file, "\n")
