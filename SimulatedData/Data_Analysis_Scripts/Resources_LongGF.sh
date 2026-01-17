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
input_gtf=${REF_DIR}/gencode.v44.annotation.gtf

# CSV functions
init_metrics_csv

run_longgf() {
    conda run -n longgf LongGF "$input_nbam" "$input_gtf" 100 50 200 \
        min_sup_read:"$min_read_supp" \
        > "$stdout_log" 2> "$stderr_log"
}

for depth in 1G 10G 100G; do
  for seqid in Q85 Q90 Q95; do
    for input_fastq in ${INPUT_DIR}/fastq_files/porechoptrimmed*${depth}*${seqid}.fastq.gz; do
      output_prefix=$(basename "$input_fastq" .fastq.gz)
      input_nbam="${ALIGNED_DIR}/${output_prefix}.name.sorted.bam"
      stdout_log="${LOGDIR}/longgf.${output_prefix}.out.log"
      stderr_log="${LOGDIR}/longgf.${output_prefix}.err.log"
      time_log="${LOGDIR}/longgf.${output_prefix}.time.log"
      output_dir="longgf_results_${output_prefix}"

      [[ -f "$input_nbam" ]] || { echo "Missing BAM: $input_nbam"; continue; }

      echo "[$(date)] Running LongGF on $output_prefix"
      /usr/bin/time -v -o "$time_log" run_longgf

      write_metrics "$output_prefix" "longgf" "$input_fastq" "$output_dir" "$time_log" "$threads"
    done
  done
done
