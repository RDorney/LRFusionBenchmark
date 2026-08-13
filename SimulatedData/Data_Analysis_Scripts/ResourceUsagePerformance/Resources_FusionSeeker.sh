#!/bin/bash
set -Eeuo pipefail
trap 'echo "[ERROR] Script failed at line $LINENO"; exit 1' ERR
shopt -s nullglob

REF_DIR=/bioinformatics/ryley/Gencode44/reference_v44
DATA_DIR=/bioinformatics/ryley/Algorithm_Benchmark
INPUT_DIR=${DATA_DIR}/Adapter_porechop_trimmed
ALIGNED_DIR=/bioinformatics/ryley/Gencode44/Sim_Data/aligned_files
LOGDIR="$(pwd)/resource_logs"
mkdir -p "$LOGDIR"

CURRENT_DIR=$(pwd)
SINGLE_DIR="FS_SINGLE"
mkdir -p "$SINGLE_DIR"
MULTI_DIR="FS_MULTI"
mkdir -p "$MULTI_DIR"

min_read_supp=2
reference_genome=${REF_DIR}/GRCh38.primary_assembly.genome.fa.gz
input_gtf=${REF_DIR}/gencode.v44.annotation.gtf

############################
# Benchmark Metrics CSV
############################
METRICS_FILE="${LOGDIR}/benchmark_metrics.csv"

# Initialize CSV if not exists
init_metrics_csv() {
  if [[ ! -f "$METRICS_FILE" ]]; then
    echo "sample_name,tool,dataset_type,input_size_bytes,input_size_gb,size_category,quality,spiked,wall_time_seconds,user_time_seconds,system_time_seconds,cpu_time_seconds,peak_ram_kb,peak_ram_gb,exit_status,output_size_bytes,output_size_mb,timestamp,threads" > "$METRICS_FILE"
    echo "Benchmark CSV initialized: $METRICS_FILE"
  fi
}

# Parse time -v output and extract metrics
parse_time_output() {
  local time_file="$1"

  if [[ ! -f "$time_file" ]]; then
    echo "0,0,0,0,0,0,1"
    return
  fi

  local user_time=$(grep "User time" "$time_file" | awk '{print $NF}')
  local system_time=$(grep "System time" "$time_file" | awk '{print $NF}')
  local wall_time_raw=$(grep "Elapsed (wall clock)" "$time_file" | awk '{print $NF}')
  local peak_ram=$(grep "Maximum resident set size" "$time_file" | awk '{print $NF}')
  local exit_status=$(grep "Exit status" "$time_file" | awk '{print $NF}')

  # Convert wall time to seconds (handles h:mm:ss or m:ss format)
  local wall_seconds=0
  if [[ "$wall_time_raw" =~ ^([0-9]+):([0-9]+):([0-9.]+)$ ]]; then
    wall_seconds=$(echo "${BASH_REMATCH[1]}*3600 + ${BASH_REMATCH[2]}*60 + ${BASH_REMATCH[3]}" | bc)
  elif [[ "$wall_time_raw" =~ ^([0-9]+):([0-9.]+)$ ]]; then
    wall_seconds=$(echo "${BASH_REMATCH[1]}*60 + ${BASH_REMATCH[2]}" | bc)
  fi

  local cpu_time=$(echo "$user_time + $system_time" | bc)
  local peak_ram_gb=$(echo "scale=2; $peak_ram / 1048576" | bc)

  echo "${wall_seconds},${user_time},${system_time},${cpu_time},${peak_ram},${peak_ram_gb},${exit_status}"
}

# Write metrics to CSV
write_metrics() {
  local sample_name="$1"
  local tool="$2"
  local input_file="$3"
  local output_dir="$4"
  local time_file="$5"
  local threading="$6"

  # Get input file size
  local input_size_bytes=$(stat -c%s "$input_file" 2>/dev/null || echo 0)
  local input_size_gb=$(echo "scale=2; $input_size_bytes / 1073741824" | bc)
  
  # Determine size category, quality, and spiked status
  local size_category="NA"
  local quality="NA"
  local spiked="NA"
  local dataset_type="simulated"

  if [[ "$sample_name" =~ 1GB ]]; then
    size_category="1GB"
  elif [[ "$sample_name" =~ 10GB ]]; then
    size_category="10GB"
  elif [[ "$sample_name" =~ 100GB ]]; then
    size_category="100GB"
  fi

  if [[ "$sample_name" =~ Q([0-9]+) ]]; then
    quality="Q${BASH_REMATCH[1]}"
  fi

  if [[ "$sample_name" =~ Spiked ]]; then
    spiked="yes"
  else
    spiked="no"
  fi

  # Get output size
  local output_size_bytes=0
  if [[ -d "$output_dir" ]]; then
    output_size_bytes=$(du -sb "$output_dir" 2>/dev/null | cut -f1 || echo 0)
  fi
  local output_size_mb=$(echo "scale=2; $output_size_bytes / 1048576" | bc)

  # Parse time output
  local time_metrics=$(parse_time_output "$time_file")

  # Get timestamp
  local timestamp=$(date '+%Y-%m-%d %H:%M:%S')

  # Write to CSV
  echo "${sample_name},${tool},${dataset_type},${input_size_bytes},${input_size_gb},${size_category},${quality},${spiked},${time_metrics},${output_size_bytes},${output_size_mb},${timestamp},${threading}" >> "$METRICS_FILE"
}

# Initialize CSV
init_metrics_csv

#######################
# Loop over files
#######################
for threads in 1 10 ; do
  cd "$CURRENT_DIR"
  if [[ "$threads" -eq 10 ]]; then
    cd "$MULTI_DIR"
  else
    cd "$SINGLE_DIR"
  fi
  for depth in 1G 10G 100G; do
    for seqid in Q85 Q90 Q95; do
      for input_fastq in ${INPUT_DIR}/fastq_files/porechoptrimmed*${depth}*${seqid}.fastq.gz; do
        output_prefix=$(basename "$input_fastq" .fastq.gz)
        input_bam="${ALIGNED_DIR}/${output_prefix}.sorted.bam"
        stdout_log="${LOGDIR}/fusionseeker.${output_prefix}.${threads}.out.log"
        stderr_log="${LOGDIR}/fusionseeker.${output_prefix}.${threads}.err.log"
        time_log="${LOGDIR}/fusionseeker.${output_prefix}.${threads}.time.log"
        output_dir="fusionseeker_${output_prefix}"
  
        [[ -f "$input_bam" ]] || { echo "Missing BAM: $input_bam"; continue; }
  
        echo "[$(date)] Running FusionSeeker on $output_prefix"
        /usr/bin/time -v -o "$time_log" conda run -n FusionSeeker /opt/FusionSeeker/fusionseeker \
      --bam "$input_bam" \
      -o "fusionseeker_${output_prefix}" \
      --datatype nanopore \
      --ref "$reference_genome" \
      --gtf "$input_gtf" \
      --geneid \
      --thread "$threads" \
      --minsupp "$min_read_supp" \
      > "$stdout_log" 2> "$stderr_log"
  
        write_metrics "$output_prefix" "fusionseeker" "$input_fastq" "$output_dir" "$time_log" "$threads"
      done
    done
  done
done