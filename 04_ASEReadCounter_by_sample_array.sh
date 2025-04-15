#!/bin/bash
#SBATCH --job-name=ASE                # Job name
#SBATCH --nodes=1                     # Request 1 node
#SBATCH --ntasks=1                    # Run one task (single process)
#SBATCH --array=1-15                  # SLURM array job from 1 to 15 (for 15 samples)
#SBATCH --cpus-per-task=1             # Use 1 CPU per task
#SBATCH --mem=44GB                    # Request 44 GB of memory
#SBATCH --time=01:00:00               # Max runtime: 1 hour
#SBATCH -A data-machine               # Project account for allocation
#SBATCH --output=F1_D1_%A_%a.out      # Stdout log file per task
#SBATCH --error=F1_D1%A_%a.err        # Stderr log file per task

# set working directory
cd /mnt/research/glbrc_group/lowry_lab/RNASeq_cistrans/ASE

# Load bcftools module
module purge 
module load BCFtools/1.19-GCC-13.2.0

# Get sample name from list using SLURM_ARRAY_TASK_ID
SAMPLE_NAME=$(awk "NR==${SLURM_ARRAY_TASK_ID}" F1_D1_samples.txt | cut -f1)

# Create a filtered VCF for homozygous parents and heterozygous F1 from the prefiltered set (already done, so commented out)
# bcftools view -T pos_homoParents_heteF1_day1.txt F1_D1_snps_Q40_DP20_FM01_annotated_within_genes_biallelic.vcf -o F1_D1_snps_Q40_DP20_FM01_annotated_within_genes_biallelic_heteF_D1.vcf &

# Filter the full multi-sample VCF to retain one sample per run
bcftools view -s ${SAMPLE_NAME} F1_D1_snps_Q40_DP20_FM01_annotated_within_genes_biallelic_heteF_D1.vcf \
  -Oz -o sample_vcfs/${SAMPLE_NAME}_snps_Q40_DP20_FM01_annotated_within_genes_biallelic_heteF_D1.vcf.gz && \
tabix -p vcf sample_vcfs/${SAMPLE_NAME}_snps_Q40_DP20_FM01_annotated_within_genes_biallelic_heteF_D1.vcf.gz &

# Switch to GATK module after using bcftools (there is a module conflict between bcftools and gatk)
module purge
module load GATK/4.5.0.0-GCCcore-12.3.0-Java-17

# Run GATK ASEReadCounter on this sample's BAM and VCF
gatk ASEReadCounter \
   -R /mnt/research/glbrc_group/lowry_lab/WGS/referenceGenome/Pvirgatumvar_AP13HAP1_772_v6.0.hardmasked.fa \
   -I /mnt/research/glbrc_group/lowry_lab/RNASeq_cistrans/STAR_mapping/${SAMPLE_NAME}_2ndPass_Aligned.sortedByCoord.out.dedup.split_rg.bam \
   -V /mnt/research/glbrc_group/lowry_lab/RNASeq_cistrans/ASE/sample_vcfs/${SAMPLE_NAME}_snps_Q40_DP20_FM01_annotated_within_genes_biallelic_heteF_D1.vcf.gz \
   -O /mnt/research/glbrc_group/lowry_lab/RNASeq_cistrans/ASE/ASE_tables/${SAMPLE_NAME}_heteF_D1_table
