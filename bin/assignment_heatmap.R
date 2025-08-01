#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))

version <- "0.0.1"

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Handle --version
if ("--version" %in% args) {
  cat("assignment_heatmap.R version", version, "\n")
  quit(status = 0)
}

# Print help if too few args or help flag is used
if (length(args) < 2 || "--help" %in% args || "-h" %in% args) {
  cat("
Usage:
  Rscript script.R <input_file.tsv> <output_file.png>

Options:
  --version     Print version and exit
  --help, -h    Show this help message
\n")
  quit(status = 1)
}

# Assign positional arguments
input_file <- args[1]
output_file <- args[2]

# Read the data
likelihood_df <- read.table(input_file, sep = "\t", header = TRUE)
likelihood_data <- as.data.table(likelihood_df)

# Melt the data
melted_data <- melt(likelihood_data, id.vars = "X")

# Create the plot
p <- ggplot(melted_data, aes(x = variable, y = X, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(x = "Taxid", y = "Reads", fill = "Likelihood") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# Save to PNG
ggsave(output_file, plot = p, width = 10, height = 8, dpi = 300)

sessionInfo()
cat("Plot saved to", output_file, "\n")
