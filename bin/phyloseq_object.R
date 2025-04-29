#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(stringr))


# # # # # # # # # # # # # # # # #
# # # FUNCTIONS USED BELOW    # # #
# # # # # # # # # # # # # # # # #
gms16s_to_phyloseq <- function(counts_file, taxonomy_file,
                        sep_file_sample = "_", n_parts_filename_to_sample = 1) {
    if (endsWith(counts_file, ".xlsx")) {
        gms_counts <- openxlsx::read.xlsx(counts_file)
    } else if (endsWith(counts_file, ".tsv")) {
        gms_counts <- read.csv(counts_file, sep = "\t")
        gms_counts$Sample <- sapply(gms_counts$Source_File,
                                    str_first_n_elements,
                                    n = n_parts_filename_to_sample,
                                    sep = sep_file_sample)
        if (length(unique(gms_counts$Source_File))
                != length(unique(gms_counts$Sample))) {
            warning(
                "Source_File prefix specified by sep_file_sample = \"", sep_file_sample,
                "\" and n_parts_filename_to_sample = ", n_parts_filename_to_sample,
                " not unique. Using the Source_File as sample name instead."
            )
            gms_counts$Sample <- gms_counts$Source_File
        }
    }

    ## first: "unclassified should be made unique to avoid that all are pooled together
    otu_df <- gms_counts[, c("tax_id", "estimated.counts", "Sample")]

    otu_tab_df <- tidyr::pivot_wider(otu_df, id_cols = tax_id,
                                        names_from = Sample,
                                        values_from = estimated.counts)
    otu_tab <- as.matrix(otu_tab_df[,-1])
    rownames(otu_tab) <- otu_tab_df$tax_id
    which(apply(otu_tab, c(1,2), is_numeric), arr.ind = T)

    ## NAs could not be converted to num...
    otu_tab[which(!apply(otu_tab, c(1,2), is_numeric), arr.ind = T)]

    otu_tab_num <- apply(otu_tab, c(1,2), as.numeric)
    # image(otu_tab_num)
    which(otu_tab_num == 0) # no zeroes
    otu_tab_num[which(is.na(otu_tab_num))] <- 0
    # image(otu_tab_num)

    # I decided to sort the samples because their order is randomly determined
    # by Nextflow when the reports are collected in the previous step.
    # In my experience, this is usually not helpful.
    otu_tab_num <- otu_tab_num[,order(colnames(otu_tab_num))]

    gms_taxonomy_general <- read.csv(taxonomy_file, sep = "\t")
    rownames(gms_taxonomy_general) <- gms_taxonomy_general$tax_id
    gms_taxonomy_general <- gms_taxonomy_general[, c("superkingdom", "clade",
        "phylum", "class", "order", "family", "genus", "species", "subspecies",
        "species.subgroup", "species.group")]

    gms_taxonomy <- as.matrix(gms_taxonomy_general)
    ## code adapted from gtools::capwords
    word_caps <- function(s, strict = TRUE) {
        return(paste(toupper(substring(s, 1, 1)), {
            s <- substring(s, 2)
            if (strict)
                tolower(s)
            else s
        }, sep = ""))
    }
    colnames(gms_taxonomy) <- sapply(colnames(gms_taxonomy), word_caps)

    empty_sample_data_table <- as.data.frame(matrix(colnames(otu_tab_num), nrow = ncol(otu_tab_num), dimnames = list(colnames(otu_tab_num), "SampleID")))

    ph <- phyloseq(otu_table = otu_table(otu_tab_num, taxa_are_rows = TRUE),
        tax_table = tax_table(gms_taxonomy))

    sample_data(ph) <- empty_sample_data_table

    return(ph)
}

#' Remove the "x_" prefix from taxonomic terms, where x is an abbreviation
#' for the taxonomic level, e.g. "g" for "genus"
remove_genus_prefix_from_species <- function(tax_tab, i_genus, i_species) {
    rank_nam <- rank_names(tax_tab)

    remove_gen_pref_fom_spec_row <- function(tax_row, i_genus, i_species) {
        genus <- as.character(tax_row[,i_genus])
        species <- as.character(tax_row[,i_species])
        if (!is.na(genus) && nchar(genus) > 0 &&
                !is.na(species) && nchar(species) > 0 &&
                startsWith(species, genus)) {
            genus_pattern = genus
            # genus is followed by blank, _, - or .
            ext_pattern <- paste0("^", genus_pattern, "[ _\\-\\.]")
            if (grepl(ext_pattern, species)) {
                genus_pattern <- ext_pattern
            }
            tax_row[i_species] <- stringr::str_replace(species, pattern = genus_pattern, replacement = "")
        }
        return(tax_row)
    }

    new_tax_table_t <- sapply(rownames(tax_tab),
        function(rn, tax_tab, i_genus, i_species) {
            remove_gen_pref_fom_spec_row(as.matrix(tax_tab)[rn,],
                i_genus = i_genus, i_species = i_species)
        }, tax_tab, i_genus, i_species)

    new_tax_table <- t(new_tax_table_t)
    colnames(new_tax_table) <- rank_nam

    return(new_tax_table)
}


## AUXILIARY FUNCTIONS
str_first_n_elements <- function(string, n = 2, sep = "-") {
    return(paste(stringr::str_split(string, sep)[[1]][1:n], collapse = sep))
}

is_numeric <- function(n) {
    is_num <- suppressWarnings(!is.na(as.numeric(n)))
    return(is_num)
}


# # # # # # # # # #
###### MAIN #######
# # # # # # # # # #
args <- commandArgs(trailingOnly = TRUE)

phyloseq_object <- gms16s_to_phyloseq(counts_file = args[1], # e.g. gms_16s_Samples.xlsx or tsv
    taxonomy_file = args[2], n_parts_filename_to_sample = 3)

tax_table(phyloseq_object) <- remove_genus_prefix_from_species(
    tax_table(phyloseq_object), i_genus = 7, i_species = 8)


save(phyloseq_object, file = "./phyloseq_object.RData")
