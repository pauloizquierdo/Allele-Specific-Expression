# ASE Analysis Pipeline

This repository contains a pipeline for identifying and summarizing **allele-specific expression (ASE)** across genotypes using SNPs that are homozygous in parents and heterozygous in F1 hybrids.

---

## 📁 Directory Overview

| Script | Description |
|--------|-------------|
| `01_annotate_and_filter_snps_within_genes.sh` | Annotates SNPs and filters to retain only those within gene regions and that are biallelic. |
| `02_filter_homo_hetero_snps_for_ASE.r` | Filters SNPs based on allele frequency to identify sites that are homozygous in parents (AF < 0.1 or > 0.9) and heterozygous in F1 (AF 0.3–0.7). |
| `03_filter_population_vcf_by_homoParents_heteF1.sh` | Filters a population-level VCF using position lists of selected SNPs from parents and F1. |
| `04_ASEReadCounter_by_sample_array.sh` | SLURM array script to run GATK ASEReadCounter on each F1 sample using filtered VCFs and BAMs. |
| `05_summarize_ASE_gene_level.r` | Combines ASE count tables per day, joins them with gene annotations, and summarizes mean/median/sd/snps expression per gene. |

---

## 🧪 Dependencies

- [`bcftools`](http://samtools.github.io/bcftools/)
- [`GATK 4`](https://gatk.broadinstitute.org/)
- `R` with the following packages:
  - `tidyverse`
  - `data.table`

---

## 🧬 Input Requirements

- VCF files annotated with gene information (e.g., `TGN=...`)
- BAM files for F1 samples, coordinate-sorted and indexed
- Sample lists (`F1_D1_samples.txt`) and filtered SNP position files per day

---

## ⚙️ Output

- Filtered VCFs per sample
- ASEReadCounter tables per sample
- Merged and summarized gene-level ASE tables per day

---

## 📌 Notes

- Ensure all VCFs are compressed and indexed (`.vcf.gz + .tbi`).
- Scripts are meant to be run in order: `01 → 05`.

---

