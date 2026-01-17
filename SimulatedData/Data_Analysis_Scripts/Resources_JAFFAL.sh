#!/bin/bash
set -Eeuo pipefail
trap 'echo "[ERROR] Script failed at line $LINENO"; exit 1' ERR
shopt -s nullglob

############################
# Config
############################
REF_DIR=/bioinformatics/ryley/Gencode44/reference_v44
DATA_DIR=/bioinformatics/ryley/Algorithm_Benchmark
INPUT_DIR=${DATA_DIR}/Adapter_porechop_trimmed
ALIGNED_DIR=/bioinformatics/ryley/Gencode44/Sim_Data/aligned_files
LOGDIR="resource_logs"
mkdir -p "$LOGDIR"

threads=1
METRICS_FILE="${LOGDIR}/benchmark_metrics.csv"
min_read_supp=2

############################
# CSV functions
############################
# (copy init_metrics_csv, parse_time_output, write_metrics here)
# initialize CSV
init_metrics_csv
init_metrics_csv

############################
# Reference
############################
input_gtf=${REF_DIR}/gencode.v44.annotation.gtf
reference_genome=${REF_DIR}/GRCh38.primary_assembly.genome.fa.gz

############################
# Jaffa runner
############################
run_jaffal() {
    local OLD_PATH="$PATH"
    export JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64
    export PATH="$JAVA_HOME/bin:$PATH"

    /opt/JAFFA-version-2.3/tools/bin/bpipe run -n "$threads" \
      -p genome=hg38 \
      -p annotation=genCode43 \
      -p exclude=NoSupport \
      -p jaffa_output="jaffal_results_${output_prefix}" \
      /opt/JAFFA-version-2.3/JAFFAL.groovy \
      "$input_fastq" \
      > "$stdout_log" 2> "$stderr_log"

    export PATH="$OLD_PATH"
    unset JAVA_HOME
}

############################
# Loop over files
############################
for depth in 1G 10G 100G; do
  for seqid in Q85 Q90 Q95; do
    for input_fastq in ${INPUT_DIR}/fastq_files/porechoptrimmed*${depth}*${seqid}.fastq.gz; do
      output_prefix=$(basename "$input_fastq" .fastq.gz)
      stdout_log="${LOGDIR}/jaffal.${output_prefix}.out.log"
      stderr_log="${LOGDIR}/jaffal.${output_prefix}.err.log"
      time_log="${LOGDIR}/jaffal.${output_prefix}.time.log"
      output_dir="jaffal_results_${output_prefix}"

      [[ -f "$input_fastq" ]] || { echo "Missing: $input_fastq"; continue; }

      echo "[$(date)] Running Jaffa on $output_prefix"
      /usr/bin/time -v -o "$time_log" run_jaffal

      echo "[$(date)] Recording metrics for Jaffa $output_prefix"
      write_metrics "$output_prefix" "jaffal" "$input_fastq" "$output_dir" "$time_log" "$threads"
    done
  done
done
