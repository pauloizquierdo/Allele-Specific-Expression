# Load necessary libraries
library(tidyverse)

# Clear R environment
rm(list = ls())

# Set working directory
setwd("ASE/ASE_tables")

# Load all ASE tables by day using filename patterns
d1_files <- list.files(pattern = "D1", full.names = TRUE)
d2_files <- list.files(pattern = "D2", full.names = TRUE)
d3_files <- list.files(pattern = "D3", full.names = TRUE)

# Read tables for each day into lists
d1 <- lapply(d1_files, read.table, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
d2 <- lapply(d2_files, read.table, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
d3 <- lapply(d3_files, read.table, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Assign sample names to each table by cleaning up file names
names(d1) <- gsub("\\./|_heteF_D1_table", "", d1_files)
names(d2) <- gsub("\\./|_heteF_D2_table", "", d2_files)
names(d3) <- gsub("\\./|_heteF_D3_table", "", d3_files)

# Select only relevant columns: contig, position, ref/alt alleles, and allele counts
d1 <- lapply(d1, function(x) x[, c(1,2,4:7)])
d2 <- lapply(d2, function(x) x[, c(1,2,4:7)])
d3 <- lapply(d3, function(x) x[, c(1,2,4:7)])

# Rename the 5th and 6th columns to include sample name
d1 <- Map(function(x, nm) {
  colnames(x)[5] <- paste0(nm, "_", colnames(x)[5])
  colnames(x)[6] <- paste0(nm, "_", colnames(x)[6])
  return(x)
}, d1, names(d1))

d2 <- Map(function(x, nm) {
  colnames(x)[5] <- paste0(nm, "_", colnames(x)[5])
  colnames(x)[6] <- paste0(nm, "_", colnames(x)[6])
  return(x)
}, d2, names(d2))

d3 <- Map(function(x, nm) {
  colnames(x)[5] <- paste0(nm, "_", colnames(x)[5])
  colnames(x)[6] <- paste0(nm, "_", colnames(x)[6])
  return(x)
}, d3, names(d3))

# Merge all tables by contig and position across samples per day
d1 <- Reduce(function(x, y) {
  merge(x, y, by = c("contig", "position", "refAllele", "altAllele"), all = TRUE)
}, d1)

d2 <- Reduce(function(x, y) {
  merge(x, y, by = c("contig", "position", "refAllele", "altAllele"), all = TRUE)
}, d2)

d3 <- Reduce(function(x, y) {
  merge(x, y, by = c("contig", "position", "refAllele", "altAllele"), all = TRUE)
}, d3)

# Rename contig and position columns to Chromosome and Position
colnames(d1)[1:2] <- c("Chromosome", "Position")
colnames(d2)[1:2] <- c("Chromosome", "Position")
colnames(d3)[1:2] <- c("Chromosome", "Position")

# Load gene annotation info (from VCF INFO field)
gene_info <- read.csv("gene_info.csv")
gene_info <- gene_info[, c(1, 2, 7)]  # Keep Chromosome, Position, TGN

# Merge ASE tables with gene annotation
d1_gene <- merge(d1, gene_info, by = c("Chromosome", "Position"))
d2_gene <- merge(d2, gene_info, by = c("Chromosome", "Position"))
d3_gene <- merge(d3, gene_info, by = c("Chromosome", "Position"))

# Save merged tables per SNP
write.csv(d1_gene, "day1_F1_hetero_ASE_gene.csv", row.names = FALSE)
write.csv(d2_gene, "day2_F1_hetero_ASE_gene.csv", row.names = FALSE)
write.csv(d3_gene, "day3_F1_hetero_ASE_gene.csv", row.names = FALSE)

# Summarize mean, median, standard deviation, and count per gene across all columns
d1_sum <- d1_gene %>%
  group_by(TGN) %>%
  summarise(across(c(5:34), 
                   list(
                     mean = ~mean(.x, na.rm = TRUE),
                     median = ~median(.x, na.rm = TRUE),
                     sd = ~sd(.x, na.rm = TRUE),
                     n = ~sum(!is.na(.x))
                   ), .names = "{.col}_{.fn}")) %>%
  ungroup()

d2_sum <- d2_gene %>%
  group_by(TGN) %>%
  summarise(across(c(5:34), 
                   list(
                     mean = ~mean(.x, na.rm = TRUE),
                     median = ~median(.x, na.rm = TRUE),
                     sd = ~sd(.x, na.rm = TRUE),
                     n = ~sum(!is.na(.x))
                   ), .names = "{.col}_{.fn}")) %>%
  ungroup()

d3_sum <- d3_gene %>%
  group_by(TGN) %>%
  summarise(across(c(5:32),  # NOTE: fewer columns in d3
                   list(
                     mean = ~mean(.x, na.rm = TRUE),
                     median = ~median(.x, na.rm = TRUE),
                     sd = ~sd(.x, na.rm = TRUE),
                     n = ~sum(!is.na(.x))
                   ), .names = "{.col}_{.fn}")) %>%
  ungroup()

# Save the summarized tables per gene
write.csv(d1_sum, "day1_F1_hetero_ASE_gene_summary.csv", row.names = FALSE)
write.csv(d2_sum, "day2_F1_hetero_ASE_gene_summary.csv", row.names = FALSE)
write.csv(d3_sum, "day3_F1_hetero_ASE_gene_summary.csv", row.names = FALSE)
