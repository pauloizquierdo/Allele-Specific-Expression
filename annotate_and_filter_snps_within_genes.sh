# Load bcftools
module purge
module load BCFtools/1.19-GCC-13.2.0

##############################################
# STEP 1: Filter SNPs that are within genes
##############################################
# Keep variants where the INFO/TGN field is NOT just a dot ("."),
# which implies the SNP is annotated as being within a gene.

bcftools view -i 'INFO/TGN~"\."' -o AP13_D1_snps_Q40_DP20_FM01_annotated_within_genes.vcf -O z ../population/AP13_D1_snps_Q40_DP20_FM01_annotated.vcf &
bcftools view -i 'INFO/TGN~"\."' -o DAC6_D1_snps_Q40_DP20_FM01_annotated_within_genes.vcf -O z ../population/DAC6_D1_snps_Q40_DP20_FM01_annotated.vcf &
bcftools view -i 'INFO/TGN~"\."' -o F1_D1_snps_Q40_DP20_FM01_annotated_within_genes.vcf -O z ../population/F1_D1_snps_Q40_DP20_FM01_annotated.vcf &

bcftools view -i 'INFO/TGN~"\."' -o AP13_D2_snps_Q40_DP20_FM01_annotated_within_genes.vcf -O z ../population/AP13_D2_snps_Q40_DP20_FM01_annotated.vcf &
bcftools view -i 'INFO/TGN~"\."' -o DAC6_D2_snps_Q40_DP20_FM01_annotated_within_genes.vcf -O z ../population/DAC6_D2_snps_Q40_DP20_FM01_annotated.vcf &
bcftools view -i 'INFO/TGN~"\."' -o F1_D2_snps_Q40_DP20_FM01_annotated_within_genes.vcf -O z ../population/F1_D2_snps_Q40_DP20_FM01_annotated.vcf &

bcftools view -i 'INFO/TGN~"\."' -o AP13_D3_snps_Q40_DP20_FM01_annotated_within_genes.vcf -O z ../population/AP13_D3_snps_Q40_DP20_FM01_annotated.vcf &
bcftools view -i 'INFO/TGN~"\."' -o DAC6_D3_snps_Q40_DP20_FM01_annotated_within_genes.vcf -O z ../population/DAC6_D3_snps_Q40_DP20_FM01_annotated.vcf &
bcftools view -i 'INFO/TGN~"\."' -o F1_D3_snps_Q40_DP20_FM01_annotated_within_genes.vcf -O z ../population/F1_D3_snps_Q40_DP20_FM01_annotated.vcf &

##############################################
# STEP 2: Add allele frequency annotations
##############################################
# Use `fill-tags` plugin to compute allele frequencies (AF tag).

bcftools +fill-tags -Oz -o AP13_D1_snps_Q40_DP20_FM01_annotated_within_genes_freq.vcf AP13_D1_snps_Q40_DP20_FM01_annotated_within_genes.vcf &
bcftools +fill-tags -Oz -o DAC6_D1_snps_Q40_DP20_FM01_annotated_within_genes_freq.vcf DAC6_D1_snps_Q40_DP20_FM01_annotated_within_genes.vcf &
bcftools +fill-tags -Oz -o F1_D1_snps_Q40_DP20_FM01_annotated_within_genes_freq.vcf F1_D1_snps_Q40_DP20_FM01_annotated_within_genes.vcf &

bcftools +fill-tags -Oz -o AP13_D2_snps_Q40_DP20_FM01_annotated_within_genes_freq.vcf AP13_D2_snps_Q40_DP20_FM01_annotated_within_genes.vcf &
bcftools +fill-tags -Oz -o DAC6_D2_snps_Q40_DP20_FM01_annotated_within_genes_freq.vcf DAC6_D2_snps_Q40_DP20_FM01_annotated_within_genes.vcf &
bcftools +fill-tags -Oz -o F1_D2_snps_Q40_DP20_FM01_annotated_within_genes_freq.vcf F1_D2_snps_Q40_DP20_FM01_annotated_within_genes.vcf &

bcftools +fill-tags -Oz -o AP13_D3_snps_Q40_DP20_FM01_annotated_within_genes_freq.vcf AP13_D3_snps_Q40_DP20_FM01_annotated_within_genes.vcf &
bcftools +fill-tags -Oz -o DAC6_D3_snps_Q40_DP20_FM01_annotated_within_genes_freq.vcf DAC6_D3_snps_Q40_DP20_FM01_annotated_within_genes.vcf &
bcftools +fill-tags -Oz -o F1_D3_snps_Q40_DP20_FM01_annotated_within_genes_freq.vcf F1_D3_snps_Q40_DP20_FM01_annotated_within_genes.vcf &

##############################################
# STEP 3: Keep only biallelic SNPs
##############################################
# Use `-m2 -M2` to retain only biallelic variants
# `-v snps` ensures only SNPs are retained (excludes indels, etc.)

bcftools view -m2 -M2 -v snps -o AP13_D1_snps_Q40_DP20_FM01_annotated_within_genes_biallelic.vcf AP13_D1_snps_Q40_DP20_FM01_annotated_within_genes_freq.vcf &
bcftools view -m2 -M2 -v snps -o DAC6_D1_snps_Q40_DP20_FM01_annotated_within_genes_biallelic.vcf DAC6_D1_snps_Q40_DP20_FM01_annotated_within_genes_freq.vcf &
bcftools view -m2 -M2 -v snps -o F1_D1_snps_Q40_DP20_FM01_annotated_within_genes_biallelic.vcf F1_D1_snps_Q40_DP20_FM01_annotated_within_genes_freq.vcf &

bcftools view -m2 -M2 -v snps -o AP13_D2_snps_Q40_DP20_FM01_annotated_within_genes_biallelic.vcf AP13_D2_snps_Q40_DP20_FM01_annotated_within_genes_freq.vcf &
bcftools view -m2 -M2 -v snps -o DAC6_D2_snps_Q40_DP20_FM01_annotated_within_genes_biallelic.vcf DAC6_D2_snps_Q40_DP20_FM01_annotated_within_genes_freq.vcf &
bcftools view -m2 -M2 -v snps -o F1_D2_snps_Q40_DP20_FM01_annotated_within_genes_biallelic.vcf F1_D2_snps_Q40_DP20_FM01_annotated_within_genes_freq.vcf &

bcftools view -m2 -M2 -v snps -o AP13_D3_snps_Q40_DP20_FM01_annotated_within_genes_biallelic.vcf AP13_D3_snps_Q40_DP20_FM01_annotated_within_genes_freq.vcf &
bcftools view -m2 -M2 -v snps -o DAC6_D3_snps_Q40_DP20_FM01_annotated_within_genes_biallelic.vcf DAC6_D3_snps_Q40_DP20_FM01_annotated_within_genes_freq.vcf &
bcftools view -m2 -M2 -v snps -o F1_D3_snps_Q40_DP20_FM01_annotated_within_genes_biallelic.vcf F1_D3_snps_Q40_DP20_FM01_annotated_within_genes_freq.vcf &
