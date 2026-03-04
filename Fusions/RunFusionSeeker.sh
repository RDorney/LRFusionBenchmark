#!/bin/bash
set -Eeuo pipefail
trap 'echo "[ERROR] Script failed at line $LINENO"; exit 1' ERR
shopt -s nullglob

############################
# Configuration
############################
#set directory names
START_DIR=$(pwd)
REF_DIR=/bioinformatics/ryley/Gencode44/reference_v44
DATA_DIR=/bioinformatics/ryley/Gencode44/Huh7_Library/FusionSeeker

INPUT_DIR=${DATA_DIR}
ALIGNED_DIR=/bioinformatics/ryley/Gencode44/Huh7_Library/FusionSeeker

#thread
threads=4

#shared default values
min_read_supp=2

############################
# Reference files
############################
reference_genome=${REF_DIR}/GRCh38.primary_assembly.genome.fa
gencode_gtf=${REF_DIR}/gencode.v44.annotation.gtf

############################
# Tool runners
############################
run_fusionseeker () {
  local input_bam="$1"
  local output_prefix="$2"
  local threads="$3"
  
  local stdout_log="fusionseeker.${output_prefix}.${threads}.out.log"
  local stderr_log="fusionseeker.${output_prefix}.${threads}.err.log"
  
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate FusionSeeker 
  
  local OLD_PATH="$PATH"
  export PATH="/opt/FusionSeeker:/opt/FusionSeeker/bsalign:$PATH"
  
    fusionseeker \
    --bam "$input_bam" \
    -o "fusionseeker_${output_prefix}" \
    --datatype nanopore \
    --ref "$reference_genome" \
    --gtf "$gencode_gtf" \
    --geneid \
    --thread "$threads" \
    --minsupp "$min_read_supp" \
    > "$stdout_log" 2> "$stderr_log"
    
    export PATH="$OLD_PATH"
    conda deactivate
}

export -f run_fusionseeker

#############################
#Loop through tools and files
#############################
for input_bam in ${INPUT_DIR}/chopperfiltered_dorado_trimmed_*.sorted.bam; do
    output_prefix=$(basename "${input_bam}" .sorted.bam)
    run_fusionseeker "$input_bam" "$output_prefix" "$threads" 
done

echo ""
echo "============================================"
echo "FusionSeeker complete!"
echo "============================================"

