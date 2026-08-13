#!/bin/bash
set -e

############################
# GFSeeker Test Pipeline
############################

# Tool paths - use minimap2 v2.28 from JAFFAL
export PATH="/opt/JAFFA-version-2.3/tools/bin:/bioinformatics/siyuan/tools/minigraph-0.21_x64-linux:$PATH"

# Directories
GFSEEKER_DIR="/bioinformatics/siyuan/tools/GFSeeker"
REF_GENOME="/bioinformatics/siyuan/references/GRCh38_gencode_v44_CTAT_lib/ref_genome.fa"
REF_GFF3="/bioinformatics/siyuan/references/gencode.v44.annotation.gff3"
OUTPUT_DIR="/bioinformatics/siyuan/gfseeker_results_v2"
GRAPH_DIR="${OUTPUT_DIR}/ref_graph"
LOGDIR="${OUTPUT_DIR}/logs"

# Threads
THREADS=10

# Data directories (same as CTAT-LR-fusion tests)
CELLLINE_DIR="/bioinformatics/ryley/Algorithm_Benchmark/SGNEx_files/Rasusa_files"
SIMULATED_DIR="/bioinformatics/ryley/Algorithm_Benchmark/Adapter_porechop_trimmed/fastq_files"

# Benchmark metrics CSV (same name as CTAT-LR)
METRICS_FILE="${OUTPUT_DIR}/benchmark_metrics.csv"

############################
# Benchmark Metrics Functions
############################
init_metrics_csv() {
  if [[ ! -f "$METRICS_FILE" ]]; then
    echo "sample_name,tool,dataset_type,input_size_bytes,input_size_gb,size_category,quality,spiked,wall_time_seconds,user_time_seconds,system_time_seconds,cpu_time_seconds,peak_ram_kb,peak_ram_gb,exit_status,output_size_bytes,output_size_mb,timestamp" > "$METRICS_FILE"
    echo "Benchmark CSV initialized: $METRICS_FILE"
  fi
}

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

write_metrics() {
  local sample_name="$1"
  local input_file="$2"
  local output_dir="$3"
  local time_file="$4"

  local input_size_bytes=$(stat -c%s "$input_file" 2>/dev/null || echo 0)
  local input_size_gb=$(echo "scale=2; $input_size_bytes / 1073741824" | bc)

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

  if [[ "$sample_name" =~ SGNex ]]; then
    dataset_type="cellline"
  fi

  local output_size_bytes=0
  if [[ -d "$output_dir" ]]; then
    output_size_bytes=$(du -sb "$output_dir" 2>/dev/null | cut -f1 || echo 0)
  fi
  local output_size_mb=$(echo "scale=2; $output_size_bytes / 1048576" | bc)

  local time_metrics=$(parse_time_output "$time_file")
  local timestamp=$(date '+%Y-%m-%d %H:%M:%S')

  echo "${sample_name},GFSeeker,${dataset_type},${input_size_bytes},${input_size_gb},${size_category},${quality},${spiked},${time_metrics},${output_size_bytes},${output_size_mb},${timestamp}" >> "$METRICS_FILE"
}

############################
# Step 1: Build Reference Graph (only once)
############################
build_reference_graph() {
  if [[ -f "${GRAPH_DIR}/comprehensive.gfa" ]]; then
    echo "Reference graph already exists, skipping build..."
    return 0
  fi

  echo "============================================"
  echo "[$(date)] Building GFSeeker reference graph"
  echo "============================================"

  cd "$GFSEEKER_DIR"

  /usr/bin/time -v -o "${LOGDIR}/graph_build.time" \
    python3 GFSeeker graph \
      "$REF_GFF3" \
      "$REF_GENOME" \
      -o "$GRAPH_DIR" \
      -t "$THREADS" \
      -l "${LOGDIR}/graph_build.log" \
    2>&1 | tee "${LOGDIR}/graph_build.stdout.log"

  # Decompress fusion_annot_lib into graph dir (GFSeeker expects plain text named "fusion_annot_lib")
  if [[ ! -f "${GRAPH_DIR}/fusion_annot_lib" ]]; then
    echo "Decompressing fusion annotation library..."
    gunzip -c /bioinformatics/siyuan/references/GRCh38_gencode_v44_CTAT_lib/fusion_annot_lib.gz > "${GRAPH_DIR}/fusion_annot_lib"
  fi

  echo "[$(date)] Reference graph build complete"
}

############################
# Step 2: Run Fusion Detection
############################
run_gfseeker() {
  local input_file="$1"
  local sample_name=$(basename "$input_file" .fastq.gz)

  # Determine dataset type and set output directory (same structure as CTAT-LR)
  local dataset_subdir
  if [[ "$sample_name" =~ SGNex ]]; then
    dataset_subdir="cellline"
  else
    dataset_subdir="simulated"
  fi

  local result_base="${OUTPUT_DIR}/${dataset_subdir}"
  local result_dir="${result_base}/${sample_name}"
  local log_dir="${result_base}/logs"
  local time_file="${log_dir}/${sample_name}.time"
  local log_file="${log_dir}/${sample_name}.log"

  echo "============================================"
  echo "[$(date)] Running GFSeeker on: $sample_name"
  echo "Input:  $input_file"
  echo "Output: $result_dir"
  echo "============================================"

  mkdir -p "$result_dir" "$log_dir"
  cd "$GFSEEKER_DIR"

  # Decompress fastq.gz (pyfaidx cannot index gzipped files)
  local temp_fastq="${result_dir}/${sample_name}.fastq"
  echo "[$(date)] Decompressing input file..."
  gunzip -c "$input_file" > "$temp_fastq"

  /usr/bin/time -v -o "$time_file" \
    python3 GFSeeker detect \
      "$temp_fastq" \
      "$REF_GENOME" \
      "$GRAPH_DIR" \
      -o "$result_dir" \
      -t "$THREADS" \
      -s 2 \
      -P False \
    > "$log_file" 2>&1

  local exit_code=$?

  # Clean up decompressed file to save disk space
  rm -f "$temp_fastq"

  echo "[$(date)] Finished: $sample_name (exit code: $exit_code)"

  # Write metrics
  write_metrics "$sample_name" "$input_file" "$result_dir" "$time_file"

  return $exit_code
}

############################
# Main
############################
echo "GFSeeker Test Pipeline"
echo "minimap2 version: $(minimap2 --version)"
echo "minigraph version: $(minigraph --version)"
echo ""

init_metrics_csv

# Create necessary directories
mkdir -p "$LOGDIR" "$GRAPH_DIR"

# Build reference graph first
build_reference_graph

echo ""
echo "============================================"
echo "Starting sample processing..."
echo "============================================"

# Process cellline data first (29 samples)
echo ""
echo "=== Step 1: Processing Cellline Data (29 samples) ==="
for input_file in ${CELLLINE_DIR}/*.fastq.gz; do
  if [[ -f "$input_file" ]]; then
    run_gfseeker "$input_file"
  fi
done

# Process simulated data (1GB -> 10GB -> 100GB)
echo ""
echo "=== Step 2: Processing Simulated Data (18 samples) ==="
for size in 1GB 10GB 100GB; do
  for quality in Q85 Q90 Q95; do
    # Normal samples
    for input_file in ${SIMULATED_DIR}/porechoptrimmed_Adapter_Simulated_${size}_nanopore2023_${quality}.fastq.gz; do
      if [[ -f "$input_file" ]]; then
        run_gfseeker "$input_file"
      fi
    done
    # Spiked samples
    for input_file in ${SIMULATED_DIR}/porechoptrimmed_Adapter_Spiked_Simulated_${size}_nanopore2023_${quality}.fastq.gz; do
      if [[ -f "$input_file" ]]; then
        run_gfseeker "$input_file"
      fi
    done
  done
done

echo ""
echo "============================================"
echo "All samples complete!"
echo "Metrics saved to: $METRICS_FILE"
echo "============================================"
