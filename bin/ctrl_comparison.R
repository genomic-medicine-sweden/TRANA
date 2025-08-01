#!/usr/bin/env Rscript
# Author: fwa93
# ---- Load libraries ----
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(glue)
  library(viridis)
})

# ---- Parse command-line args ----

version <- "0.0.1"

args <- commandArgs(trailingOnly = TRUE)

# Handle --version
if ("--version" %in% args) {
  cat("assignment_heatmap.R version", version, "\n")
  quit(status = 0)
}

if (length(args) < 2) {
  cat("Usage:\n")
  cat("  Rscript ctrl_comparison_cli.R <tsv result file from emu_combine> <name of control. Name here must end with .fastq>\n  <taxonomic_level> ")
  quit(status = 1)
}
data_file <- args[1]
## Use this instead if you are not running through command line
#data_file <- "M:\\BMA/Frans/Emu_projekt/emu-combined-family-counts_test_file_ctr_comparison2.tsv" 
#control_name <- "my_neg.fastq"
##
control_name <- args[2]

# ---- Output folder ----
output_dir <- "plots_vs_selected_ctrl"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Read and reshape data ----
raw <- read_tsv(data_file, col_types = cols(.default = col_character())) %>%
  mutate(across(where(is.character), ~na_if(., "")))  # Treat empty strings as NA

# Convert numeric columns to numeric
numeric_cols <- setdiff(colnames(raw), c("species", "genus", "family", "order", "class", "phylum", "superkingdom"))
raw[numeric_cols] <- lapply(raw[numeric_cols], as.numeric)

# Pivot longer
long_df <- raw %>%
  pivot_longer(cols = all_of(numeric_cols), names_to = "sample_name", values_to = "abundance") %>%
  filter(!is.na(abundance) & abundance > 0)

##
# Create a "species" identifier
#long_df <- long_df %>%
#  unite("species", c(family), sep = " | ", remove = FALSE)
##
long_df_original <- long_df
# Create an identifier
long_df <- long_df %>%
  unite("Taxonomic_level", colnames(raw)[1], sep = " | ", remove = FALSE)
##
# ---- Split into ctrl/sample ----
ctrl_df <- filter(long_df, sample_name == control_name)
if (nrow(ctrl_df) == 0) {
  stop(glue("❌ No entries found for control column: {control_name}"))
}

sample_names <- setdiff(unique(long_df$sample_name), control_name)

cat(glue("✅ Loaded {length(sample_names)} samples; comparing to  control: {control_name}\n\n"))

# ---- Loop over each sample ----
for (samp in sample_names) {
  sample_df <- filter(long_df, sample_name == samp)
  
  # Find shared species
  common_species <- intersect(sample_df$Taxonomic_level, ctrl_df$Taxonomic_level)
  if (length(common_species) == 0) {
    cat(glue("⚠️  No shared species between {samp} and {control_name}, skipping.\n"))
    next
  }
  
  combined_df <- bind_rows(
    filter(sample_df, Taxonomic_level %in% common_species),
    filter(ctrl_df, Taxonomic_level %in% common_species)
  )
  
  p <- ggplot(combined_df, aes(x = Taxonomic_level, y = abundance, fill = sample_name)) +
    geom_bar(stat = "identity", position = "stack") +
    coord_flip() +
    labs(
      title = glue("Sample: {samp} vs control: {control_name}"),
      x = "Taxa", y = "Relative abundance", fill = "Sample"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.text.y = element_text(size = 10)
    ) +
    scale_fill_viridis_d()
  
  # Save with 300 DPI
  out_file <- file.path(output_dir, paste0(samp, "_vs_", control_name, ".png"))
  ggsave(out_file, p, dpi = 300, width = 10, height = 6)

  cat(glue("✅ Saved plot: {out_file}\n"))
}
