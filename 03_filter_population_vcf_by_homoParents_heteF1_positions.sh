# Reset environment modules and load BCFtools
module purge
module load BCFtools/1.19-GCC-13.2.0

# Filter population VCF for positions identified as homozygous in parents
# and heterozygous in F1 on Day 1
bcftools view -T ../ASE/pos_homoParents_heteF1_day1.txt \
  population_snps_Q40_DP20_FM01_annotated.vcf \
  -o d1_homoPare_heteF1_population_snps_Q40_DP20_FM01_annotated.vcf &

# Repeat for Day 2
bcftools view -T ../ASE/pos_homoParents_heteF1_day2.txt \
  population_snps_Q40_DP20_FM01_annotated.vcf \
  -o d2_homoPare_heteF1_population_snps_Q40_DP20_FM01_annotated.vcf &

# Repeat for Day 3
bcftools view -T ../ASE/pos_homoParents_heteF1_day3.txt \
  population_snps_Q40_DP20_FM01_annotated.vcf \
  -o d3_homoPare_heteF1_population_snps_Q40_DP20_FM01_annotated.vcf &
