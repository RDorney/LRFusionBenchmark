#!/bin/bash
set -Eeuo pipefail
trap 'echo "[ERROR] Script failed at line $LINENO"; exit 1' ERR
shopt -s nullglob

REF_DIR=/bioinformatics/ryley/Gencode44/reference_v44
DATA_DIR=/bioinformatics/ryley/Algorithm_Benchmark
INPUT_DIR=${DATA_DIR}/Adapter_porechop_trimmed
ALIGNED_DIR=/bioinformatics/ryley/Gencode44/Sim_Data/aligned_files
LOGDIR="resource_logs"
mkdir -p "$LOGDIR"

threads=1
METRICS_FILE="${LOGDIR}/benchmark_metrics.csv"
min_read_supp=2
reference_genome=${REF_DIR}/GRCh38.primary_assembly.genome.fa.gz
input_gtf=${REF_DIR}/gencode.v44.annotation.gtf

# CSV functions
init_metrics_csv
init_metrics_csv

run_fusionseeker() {
  conda run -n FusionSeeker /opt/FusionSeeker/fusionseeker \
    --bam "$input_bam" \
    -o "fusionseeker_${output_prefix}" \
    --datatype nanopore \
    --ref "$reference_genome" \
    --gtf "$input_gtf" \
    --geneid \
    --thread "$threads" \
    --minsupp "$min_read_supp" \
    > "$stdout_log" 2> "$stderr_log"
}

for depth in 1G 10G 100G; do
  for seqid in Q85 Q90 Q95; do
    for input_fastq in ${INPUT_DIR}/fastq_files/porechoptrimmed*${depth}*${seqid}.fastq.gz; do
      output_prefix=$(basename "$input_fastq" .fastq.gz)
      input_bam="${ALIGNED_DIR}/${output_prefix}.sorted.bam"
      stdout_log="${LOGDIR}/fusionseeker.${output_prefix}.out.log"
      stderr_log="${LOGDIR}/fusionseeker.${output_prefix}.err.log"
      time_log="${LOGDIR}/fusionseeker.${output_prefix}.time.log"
      output_dir="fusionseeker_${output_prefix}"

      [[ -f "$input_bam" ]] || { echo "Missing BAM: $input_bam"; continue; }

      echo "[$(date)] Running FusionSeeker on $output_prefix"
      /usr/bin/time -v -o "$time_log" run_fusionseeker

      write_metrics "$output_prefix" "fusionseeker" "$input_fastq" "$output_dir" "$time_log" "$threads"
    done
  done
done
