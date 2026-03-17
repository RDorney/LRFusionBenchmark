#!/bin/bash
set -Eeuo pipefail
trap 'echo "[ERROR] Script failed at line $LINENO"; exit 1' ERR
shopt -s nullglob

############################
# Configuration
############################
#set directory names
WORK_DIR=$(pwd)

# please copy the below directories
REF_DIR=  #/bioinformatics/ryley/Gencode44/reference_v44

INPUT_DIR= #/bioinformatics/ryley/Gencode44/Huh7_Library/FusionSeeker

#thread
threads=12  #8 is default. currently using 5 on Minion. 

#shared default values
min_read_supp=2

############################
# Reference files
############################
reference_genome=${REF_DIR}/GRCh38.primary_assembly.genome.fa
gencode_gtf=${REF_DIR}/gencode.v44.annotation.gtf

############################
# Install FusionSeeker and Set Up environment
############################
source "$(conda info --base)/etc/profile.d/conda.sh"
conda create --name FusionSeeker python=3.11.3 minimap2=2.28 samtools=1.20 pysam=0.22.1 #set up conda environment. this may take a while

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate FusionSeeker 
  
git clone https://github.com/Maggi-Chen/FusionSeeker.git #install FusionSeeker
cd FusionSeeker/
./fusionseeker -h

git clone https://github.com/ruanjue/bsalign.git #install bsalign v1.2.1
cd bsalign
make

local OLD_PATH="$PATH"
export PATH="/${WORK_DIR}/FusionSeeker:/${WORK_DIR}/FusionSeeker/bsalign:$PATH"

cd ${WORK_DIR}

############################
# Tool runners
############################
run_fusionseeker () {
  local input_bam="$1"
  local output_prefix="$2"
  local threads="$3"
  
  local stdout_log="fusionseeker.${output_prefix}.${threads}.out.log"
  local stderr_log="fusionseeker.${output_prefix}.${threads}.err.log"
  
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
    
}

export -f run_fusionseeker

#############################
#Loop through tools and files
#############################
for input_bam in ${INPUT_DIR}/chopperfiltered_dorado_trimmed_*.sorted.bam; do
    output_prefix=$(basename "${input_bam}" .sorted.bam)
    outdir="fusionseeker_${output_prefix}"
    #check if files has been run
    if [[ -d "$outdir" || -d "${INPUT_DIR}/${outdir}" ]]; then
      if ! [[ -d "${outdir}/poa_workspace" || -d "${INPUT_DIR}/${outdir}/poa_workspace" ]]; then
        echo "Skipping FusionSeeker for ${output_prefix} (dir exists)"
        continue
      elif [[ -d "${outdir}/poa_workspace" ]]; then
        rm -r "${outdir}/"
        echo "previous run for ${output_prefix} was incomplete. Redoing"
      fi
    fi
    echo "running FusionSeeker on ${output_prefix}"
    run_fusionseeker "$input_bam" "$output_prefix" "$threads" 
done

export PATH="$OLD_PATH"
conda deactivate
    
echo ""
echo "============================================"
echo "FusionSeeker complete!"
echo "============================================"

