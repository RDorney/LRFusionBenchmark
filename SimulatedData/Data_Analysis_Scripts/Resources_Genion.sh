#!/bin/bash
set -Eeuo pipefail
shopt -s nullglob
trap 'echo "[ERROR] Script failed at line $LINENO"; exit 1' ERR

REF_DIR=/bioinformatics/ryley/Gencode44/reference_v44
DATA_DIR=/bioinformatics/ryley/Algorithm_Benchmark
INPUT_DIR=${DATA_DIR}/Adapter_porechop_trimmed
ALIGNED_DIR=/bioinformatics/ryley/Gencode44/Sim_Data/aligned_files
LOGDIR="resource_logs"
mkdir -p "$LOGDIR"

threads=1
METRICS_FILE="${LOGDIR}/benchmark_metrics.csv"
min_read_supp=2

genomicSuperDups=/home/ryleyd/reference_files/genion_reference_files/genomicSuperDups.txt
annotation=${REF_DIR}/Homo_sapiens.GRCh38.110.chr.gtf
cdna_self=${REF_DIR}/cdna.GRCh38v110.selfalign.tsv
input_gtf=${REF_DIR}/gencode.v44.annotation.gtf

# CSV functions
init_metrics_csv

run_genion() {
  conda run -n genion-env /opt/genion/genion -i "$input_fastq" \
         -d "$genomicSuperDups" \
         --gtf "$annotation" \
         -g "$input_paf" \
         -s "$cdna_self" \
         -o "${output_prefix}_genion" \
         --non-coding \
         --min-support "$min_read_supp" \
         > "$stdout_log" 2> "$stderr_log"
}

for depth in 1G 10G 100G; do
  for seqid in Q85 Q90 Q95; do
    for input_fastq in ${INPUT_DIR}/fastq_files/porechoptrimmed*${depth}*${seqid}.fastq.gz; do
      output_prefix=$(basename "$input_fastq" .fastq.gz)
      input_paf="${ALIGNED_DIR}/${output_prefix}.paf"
      stdout_log="${LOGDIR}/genion.${output_prefix}.out.log"
      stderr_log="${LOGDIR}/genion.${output_prefix}.err.log"
      time_log="${LOGDIR}/genion.${output_prefix}.time.log"
      output_dir="${output_prefix}_genion"

      [[ -f "$input_paf" ]] || { echo "Missing PAF: $input_paf"; continue; }

      echo "[$(date)] Running Genion on $output_prefix"
      /usr/bin/time -v -o "$time_log" run_genion

      write_metrics "$output_prefix" "genion" "$input_fastq" "$output_dir" "$time_log" "$threads"
    done
  done
done
