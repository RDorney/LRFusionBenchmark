#!/bin/bash
set -o pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate Chopper_env
Nanopore_folders=/bioinformatics/ryley/Library_Benchmark/raw_files/Huh*/*/2024*/
for trimmed_reads in ${Nanopore_folders}/dorado_trimmed*F*fastq.gz ${Nanopore_folders}/FB*d1*fastq.gz; do
  output_prefix=$(basename ${trimmed_reads} .fastq.gz)
  output_folder=$(dirname ${trimmed_reads})
  output_file=${output_folder}/chopperfiltered_${output_prefix}.fastq.gz

  if [[ -f "$output_file" ]]; then
      echo "$output_file already exists"
  else
      gunzip -c ${trimmed_reads} | chopper --minlength 100 --quality 8 | gzip > ${output_file}
  fi
  
done