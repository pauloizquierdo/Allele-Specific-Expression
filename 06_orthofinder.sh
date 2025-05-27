# This script runs OrthoFinder on the P. virgatum (switchgrass) subgenomes (K and N).
# It identifies homeologous genes by separating subgenomes, extracting primary transcripts, and running orthology analysis.

########################
# Step 1: Prepare the Input Data
########################

# Download the protein (amino acid) sequences for the P. virgatum genome from Phytozome.
# Website: https://phytozome-next.jgi.doe.gov/

# Decompress the downloaded protein file (.gz format) to obtain the raw FASTA file.
gzip -d Pvirgatumvar_AP13HAP1_772_v6.1.protein.fa.gz


########################
# Step 2: Split the Subgenomes (N and K)
########################

# Use AWK to separate protein sequences by subgenome:
# - Headers with 'NG' go to subgenome_N.fa (subgenome N)
# - Headers with 'KG' go to subgenome_K.fa (subgenome K)

awk '
/^>/ {
  if ($0 ~ /NG/) {out="subgenome_N.fa"}
  else if ($0 ~ /KG/) {out="subgenome_K.fa"}
  else {out=""} # Ignore other headers like J or unknowns
}
out {print > out}
' ../../../WGS/referenceGenome/Pvirgatumvar_AP13HAP1_772_v6.1.protein.fa

########################
# Step 3: Prepare Files for OrthoFinder
########################

# Note: The protein files contain multiple transcript isoforms per gene.
# To avoid inflated orthogroups, we keep only the longest (primary) isoform per gene.
# This reduces computation time (~10x faster) and improves orthology accuracy.

# Use OrthoFinder's script (primary_transcript.py) to extract the primary transcript from each file.
# This creates new FASTA files containing only the longest isoform per gene.

for f in *fa ; do python  ../OrthoFinder/tools/primary_transcript.py $f ; done

########################
# Step 4: Run OrthoFinder
########################

# Load OrthoFinder module
module load OrthoFinder/2.5.5-foss-2023a

# Run OrthoFinder on the directory containing the primary transcript FASTA files.
# The '-f' option specifies the input folder.
orthofinder -f primary_transcripts/