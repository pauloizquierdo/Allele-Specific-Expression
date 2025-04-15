# Set working directory and clear the environment
setwd("ASE")
rm(list = ls())

# Load required libraries
library(tidyverse)
library(data.table)

# Load VCF files (no headers) for each genotype across 3 days
# Day 1
AP13_D1 <- read.table("AP13_D1_snps_...vcf", header = F, stringsAsFactors = F)
DAC6_D1 <- read.table("DAC6_D1_snps_...vcf", header = F, stringsAsFactors = F)
F1_D1   <- read.table("F1_D1_snps_...vcf",   header = F, stringsAsFactors = F)

# Day 2
AP13_D2 <- read.table("AP13_D2_snps_Q40_DP20_FM01_annotated_within_genes_biallelic.vcf", header = F, stringsAsFactors = F)
DAC6_D2 <- read.table("DAC6_D2_snps_Q40_DP20_FM01_annotated_within_genes_biallelic.vcf", header = F, stringsAsFactors = F)
F1_D2 <- read.table("F1_D2_snps_Q40_DP20_FM01_annotated_within_genes_biallelic.vcf", header = F, stringsAsFactors = F)

# Day 3
AP13_D3 <- read.table("AP13_D3_snps_Q40_DP20_FM01_annotated_within_genes_biallelic.vcf", header = F, stringsAsFactors = F)
DAC6_D3 <- read.table("DAC6_D3_snps_Q40_DP20_FM01_annotated_within_genes_biallelic.vcf", header = F, stringsAsFactors = F)
F1_D3 <- read.table("F1_D3_snps_Q40_DP20_FM01_annotated_within_genes_biallelic.vcf", header = F, stringsAsFactors = F)

# Function to extract and format relevant VCF annotation fields (AC, AF, AN, TGN, TA)
extract_vcf_info <- function(df) {
  df[, c(1:2, 8)] %>%
    rename(Chromosome = V1, Position = V2) %>%
    mutate(
      AC  = str_extract(V8, "AC=[^;]+") %>% str_replace("AC=", ""),
      AF  = str_extract(V8, "AF=[^;]+") %>% str_replace("AF=", ""),
      AN  = str_extract(V8, "AN=[^;]+") %>% str_replace("AN=", ""),
      TA  = str_extract(V8, "TA=[^;]+") %>% str_replace("TA=", ""),
      TGN = str_extract(V8, "TGN=[^;]+") %>% str_replace("TGN=", "")
    ) %>%
    mutate(across(c(AC, AF, AN), as.numeric)) %>%
    select(-V8)
}

# Extract fields from each VCF
AP13_D1_info <- extract_vcf_info(AP13_D1)
write.csv(AP13_D1_info, "gene_info.csv", row.names = F) # write gene information for all SNPs

DAC6_D1_info <- extract_vcf_info(DAC6_D1)
F1_D1_info <- extract_vcf_info(F1_D1)

AP13_D2_info <- extract_vcf_info(AP13_D2)
DAC6_D2_info <- extract_vcf_info(DAC6_D2) 
F1_D2_info <- extract_vcf_info(F1_D2)

AP13_D3_info <- extract_vcf_info(AP13_D3)
DAC6_D3_info <- extract_vcf_info(DAC6_D3)
F1_D3_info <- extract_vcf_info(F1_D3)

# Assign genotype identity to each dataset
AP13_D1_info$parent <- "AP13"
DAC6_D1_info$parent <- "DAC6"
F1_D1_info$parent   <- "F1"

# Combine day 1 data for plotting allele frequency distribution
combined_info <- rbind(AP13_D1_info, DAC6_D1_info, F1_D1_info)

# Plot histogram of AF for Day 1
ggplot(combined_info, aes(x = AF, fill = parent)) +
  geom_histogram(color = "#e9ecef", bins = 20) +
  facet_wrap(~parent) +
  theme(text = element_text(size = 6), legend.position = "none") +
  scale_fill_manual(values = c("AP13" = "#E69F00", "DAC6" = "#56B4E9", "F1" = "salmon")) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2))


# Filter SNPs that are homozygous (AF < 0.1 or > 0.9) or heterozygous (AF between 0.3–0.7)
filter_AF_homo <- function(df) {
  df %>% filter(AF < 0.1 | AF > 0.9) %>% mutate(Genotype = "Homozygous")
}
filter_AF_hete <- function(df) {
  df %>% filter(AF > 0.3 & AF < 0.7) %>% mutate(Genotype = "Heterozygous")
}

# Apply filters to each sample
AP13_D1_filtered <- filter_AF_homo(AP13_D1_info)
DAC6_D1_filtered <- filter_AF_homo(DAC6_D1_info)
F1_D1_filtered <- filter_AF_hete(F1_D1_info)

AP13_D2_filtered <- filter_AF_homo(AP13_D2_info)
DAC6_D2_filtered <- filter_AF_homo(DAC6_D2_info)
F1_D2_filtered <- filter_AF_hete(F1_D2_info) 

AP13_D3_filtered <- filter_AF_homo(AP13_D3_info)
DAC6_D3_filtered <- filter_AF_homo(DAC6_D3_info)
F1_D3_filtered <- filter_AF_hete(F1_D3_info) 

# Set rownames to "Chr_Position" for matching SNPs
rownames(AP13_D1_filtered) <- paste0(AP13_D1_filtered$Chromosome, "_", AP13_D1_filtered$Position)
rownames(DAC6_D1_filtered) <- paste0(DAC6_D1_filtered$Chromosome,"_", DAC6_D1_filtered$Position)
rownames(F1_D1_filtered) <- paste0(F1_D1_filtered$Chromosome,"_",F1_D1_filtered$Position)

rownames(AP13_D2_filtered) <- paste0(AP13_D2_filtered$Chromosome,"_",AP13_D2_filtered$Position)
rownames(DAC6_D2_filtered) <- paste0(DAC6_D2_filtered$Chromosome,"_", DAC6_D2_filtered$Position)
rownames(F1_D2_filtered) <- paste0(F1_D2_filtered$Chromosome,"_",F1_D2_filtered$Position)

rownames(AP13_D3_filtered) <- paste0(AP13_D3_filtered$Chromosome,"_",AP13_D3_filtered$Position)
rownames(DAC6_D3_filtered) <- paste0(DAC6_D3_filtered$Chromosome,"_", DAC6_D3_filtered$Position)
rownames(F1_D3_filtered) <- paste0(F1_D3_filtered$Chromosome,"_",F1_D3_filtered$Position)

# Identify SNPs common across all 3 genotypes for each day
common_snps_D1 <- intersect(rownames(AP13_D1_filtered), intersect(rownames(DAC6_D1_filtered), rownames(F1_D1_filtered)))
# [repeat for D2 and D3]

# Subset only common SNPs
AP13_D1_filtered <- AP13_D1_filtered[common_snps_D1,]
DAC6_D1_filtered <- DAC6_D1_filtered[common_snps_D1,]
F1_D1_filtered <- F1_D1_filtered[common_snps_D1,]

AP13_D2_filtered <- AP13_D2_filtered[common_snps_D2,]
DAC6_D2_filtered <- DAC6_D2_filtered[common_snps_D2,]
F1_D2_filtered <- F1_D2_filtered[common_snps_D2,]

AP13_D3_filtered <- AP13_D3_filtered[common_snps_D3,]
DAC6_D3_filtered <- DAC6_D3_filtered[common_snps_D3,]
F1_D3_filtered <- F1_D3_filtered[common_snps_D3,]

# Further filter for SNPs where AP13 and DAC6 are opposite homozygotes
index_d1 <- which(AP13_D1_filtered$AF <= 0.1 & DAC6_D1_filtered$AF >= 0.9 |
                  AP13_D1_filtered$AF >= 0.9 & DAC6_D1_filtered$AF <= 0.1)
index_d2 <- which(AP13_D2_filtered$AF <= 0.1 & DAC6_D2_filtered$AF >= 0.9 | 
                  AP13_D2_filtered$AF >= 0.9 & DAC6_D2_filtered$AF <= 0.1)
index_d3 <- which(AP13_D3_filtered$AF <= 0.1 & DAC6_D3_filtered$AF >= 0.9 | 
                  AP13_D3_filtered$AF >= 0.9 & DAC6_D3_filtered$AF <= 0.1)

# Apply above index filters
AP13_D1_filtered <- AP13_D1_filtered[index_d1, ]
DAC6_D1_filtered <- DAC6_D1_filtered[index_d1,]
F1_D1_filtered <- F1_D1_filtered[index_d1,]

AP13_D2_filtered <- AP13_D2_filtered[index_d2,]
DAC6_D2_filtered <- DAC6_D2_filtered[index_d2,]
F1_D2_filtered <- F1_D2_filtered[index_d2,]

AP13_D3_filtered <- AP13_D3_filtered[index_d3,]
DAC6_D3_filtered <- DAC6_D3_filtered[index_d3,]
F1_D3_filtered <- F1_D3_filtered[index_d3,]

# Extract position information
pos_day1 <- AP13_D1_filtered[,1:2]
pos_day2 <- AP13_D2_filtered[,1:2]
pos_day3 <- AP13_D3_filtered[,1:2]

# Save position files for downstream ASEReadCounter
write.table(pos_day1, "pos_homoParents_heteF1_day1.txt", sep="\t", row.names=F, col.names=F, quote=F)
write.table(pos_day2, "pos_homoParents_heteF1_day2.txt", sep="\t", row.names=F, col.names=F, quote=F)
write.table(pos_day3, "pos_homoParents_heteF1_day3.txt", sep="\t", row.names=F, col.names=F, quote=F)












