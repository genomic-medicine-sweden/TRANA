#!/usr/bin/env Rscript
# Author: fwa93
# This script creates bar plots for the taxonomic information from emu-combine-outputs
# One bar plot is created for each sample. The plot also contain all taxa from
# the controls. The user can choose up too two controls. One controls is required.

# ---- Parse command-line args ----
version <- "0.0.7"
# load libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(readr)
    library(glue)
    library(viridis)
    library(getopt)
})

# options
spec <- matrix(c(
    'help'       , 'h', 0, "logical",
    'version'    , 'v', 0, "logical",
    'counts', NA , 0, "logical",
    'input_file', 'f' , 1, "character",
    'ctrl_1', 'c1', 1, "character",
    'ctrl_2', 'c2', 2, "character"
), byrow=TRUE, ncol=4)

# Parse
opt <- getopt(spec)

# Help
if (!is.null(opt$help)) {
    cat("
Usage:
    Rscript ctrl_comparison.R [options] --input_file <emu combine output file> <--ctrl_1>

Options:
    -h, --help            Show this help message
    -v, --version         Print version and exit
        --input_file      Input file. This is the output file from emu-combine-outputs. (required)
        --ctrl_1          control 1, name must match samplesheet control name (required)
        --ctrl_2          control 2, name must match samplesheet control name
        --counts          Output results with abundance (read counts). Input file must have read counts. default:false
\n")
    quit(status=0)
}

# version
if (!is.null(opt$version)) {
    cat("ctrl_comparison_cli.R version", version, "\n")
    quit(status=0)
}

# Validate required arguments
if (is.null(opt$input_file) || is.null(opt$ctrl_1)) {
    cat("Error: --input_file and --ctrl_1 are required\n\n")
    cat(getopt(spec, usage=TRUE))
    quit(status=1)
}
data_file <- opt$input_file
control_name1 <- opt$ctrl_1
control_name2 <- opt$ctrl_2

if (!is.null(opt$counts)) {
    counts = TRUE
} else {
    counts = FALSE
}


# ---- Output folder ----
output_dir <- "plots_vs_selected_ctrl"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Read and reshape data ----
raw <- read_tsv(data_file, col_types = cols(.default = col_character())) %>%
    mutate(across(where(is.character), ~na_if(., "")))  # Treat empty strings as NA

# patterns to remove from the header
# only filtlong read names
remove_this1 <- "_filtered.fastq$"
colnames(raw) <- gsub(remove_this1, "", colnames(raw))
# downsampled reads
remove_this2 <- "_downsampled.fastq$"
colnames(raw) <- gsub(remove_this2, "", colnames(raw))
# reads with only .fastq suffix
remove_this3 <-  "\\.fastq$"
colnames(raw) <- gsub(remove_this3, "", colnames(raw))


##
# Find the last row index
last_row <- nrow(raw)
# Find the first column name
first_col <- colnames(raw)[1]
# If the last cell is NA, set it to "unassigned"
if (is.na(raw[[first_col]][last_row])) {
    raw[[first_col]][last_row] <- "unassigned"
}
##

# Convert numeric columns to numeric
numeric_cols <- setdiff(colnames(raw), c("species", "genus", "family", "order", "class", "phylum", "superkingdom"))
raw[numeric_cols] <- lapply(raw[numeric_cols], as.numeric)

# Pivot longer
long_df <- raw %>%
    pivot_longer(cols = all_of(numeric_cols), names_to = "sample_name", values_to = "abundance") %>%
    filter(!is.na(abundance) & abundance > 0)


# Create an identifier
long_df <- long_df %>%
    unite("Taxonomic_level", colnames(raw)[1], sep = " | ", remove = FALSE)
# controls
control_names <- if (!is.null(control_name2)) {
    c(control_name1, control_name2)
} else {
    c(control_name1)
}

ctrl_df <- filter(long_df, sample_name %in% control_names)
missing_controls <- setdiff(control_names, unique(ctrl_df$sample_name))

if (length(missing_controls) > 0) {
    stop(glue("❌ Missing control(s): {paste(missing_controls, collapse = ', ')}"))
}

# ---- Define sample names (excluding controls) ----
control_names <- c(control_name1, control_name2)
sample_names <- setdiff(unique(long_df$sample_name), control_names)

cat(glue("✅ Loaded {length(sample_names)} samples; comparing to controls: {paste(control_names, collapse = ', ')}\n\n"))

# ---- Loop over each sample ----
for (samp in sample_names) {
    sample_df <- filter(long_df, sample_name == samp)
    # Get union of taxa from sample + controls
    relevant_taxa <- union(sample_df$Taxonomic_level, ctrl_df$Taxonomic_level)
    # Filter for these taxa only
    combined_df <- long_df %>%
        filter(sample_name %in% c(samp, control_name1, control_name2),
            Taxonomic_level %in% relevant_taxa)
    # Optional: factor ordering by total abundance (can be turned off if you want fixed order)
    taxa_order <- combined_df %>%
        group_by(Taxonomic_level) %>%
        summarise(total_abundance = sum(abundance)) %>%
        arrange(desc(total_abundance)) %>%
        pull(Taxonomic_level)

    combined_df$Taxonomic_level <- factor(combined_df$Taxonomic_level, levels = taxa_order)

    # Plot
    p <- ggplot(combined_df, aes(x = Taxonomic_level, y = abundance, fill = sample_name)) +
        geom_bar(stat = "identity", position = "stack") +
        coord_flip() +
        labs(
            title = glue("{samp} vs controls: "),
            subtitle = glue("Control(s): {paste(control_names, collapse = ', ')} "),
            x = "Taxa", y = "Abundance", fill = "Sample"
        ) +
        theme_bw(base_size = 14) +
        theme(
            plot.title = element_text(hjust = 0.5, size = 16),
            axis.text.y = element_text(size = 10),
            plot.subtitle = element_text(hjust = 0.5)
        ) +
        scale_fill_viridis_d()

    # Save plot
    suffix <- if (counts) {
        paste0(samp, "_counts_vs_controls.png")
    } else {
        paste0(samp, "_vs_controls.png")
    }

    out_file <- file.path(output_dir, suffix)
    ggsave(out_file, p, dpi = 300, width = 10, height = 6)
    cat(glue("✅ Saved plot: {out_file}\n"))
}
print("Plots created show all taxa present in the sample, control 1 and control 2 (if specified)")
