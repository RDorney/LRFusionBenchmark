#!/bin/bash
set -Eeuo pipefail
trap 'echo "[ERROR] Script failed at line $LINENO"; exit 1' ERR
shopt -s nullglob

############################
# Configuration
############################
#set directory names
REF_DIR=/bioinformatics/ryley/Gencode44/reference_v44
DATA_DIR=/bioinformatics/ryley/Algorithm_Benchmark

INPUT_DIR=${DATA_DIR}/Adapter_porechop_trimmed
ALIGNED_DIR=/bioinformatics/ryley/Gencode44/Sim_Data/aligned_files

LOGDIR="resource_logs"
mkdir -p "$LOGDIR"

#thread
threads=1

#shared default values
min_read_supp=2

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

############################
# Reference files
############################
reference_genome=${REF_DIR}/GRCh38.primary_assembly.genome.fa.gz 
input_gtf=${REF_DIR}/gencode.v44.annotation.gtf
#specific to Genion
genomicSuperDups=/home/ryleyd/reference_files/genion_reference_files/genomicSuperDups.txt
annotation=${REF_DIR}/Homo_sapiens.GRCh38.110.chr.gtf
cdna_self=${REF_DIR}/cdna.GRCh38v110.selfalign.tsv

############################
# Tool runners
############################
run_longgf () {
    conda run -n longgf \
    LongGF "$input_nbam" "$input_gtf" \
    100 50 200 \
    min_sup_read:"$min_read_supp" \
    > "$stdout_log" 2> "$stderr_log"

}

run_genion () {
  conda run -n genion-env \
  /opt/genion/genion -i "$input_fastq" \
         -d "$genomicSuperDups" \
         --gtf "$annotation" \
         -g "$input_paf" \
         -s "$cdna_self" \
         -o "${output_prefix}_genion" \
         --non-coding \
         --min-support "$min_read_supp" \
         > "$stdout_log" 2> "$stderr_log"

}

run_fusionseeker () {
  conda run -n FusionSeeker \
  /opt/FusionSeeker/fusionseeker \
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

run_jaffal () {
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

export -f run_jaffal
export -f run_genion
export -f run_longgf
export -f run_fusionseeker

Bench_Tool() {
  local tool="$1"
  local output_prefix="$2"
  local input_fastq="$3"
  local threads="$4"

  local stdout_log="${LOGDIR}/${tool}.${output_prefix}.out.log"
  local stderr_log="${LOGDIR}/${tool}.${output_prefix}.err.log"
  local time_log="${LOGDIR}/${tool}.${output_prefix}.time.log"
  local output_dir="${tool}_results_${output_prefix}"

  echo "[$(date)] Running ${tool} on ${output_prefix}"
  export output_prefix input_fastq input_bam input_nbam input_paf threads
  /usr/bin/time -v -o "$time_log" bash -c "
$(declare -f run_${tool})

input_fastq='$input_fastq'
input_bam='$input_bam'
input_nbam='$input_nbam'
input_paf='$input_paf'
stdout_log='$stdout_log'
stderr_log='$stderr_log'
threads='$threads'
min_read_supp='$min_read_supp'
reference_genome='$reference_genome'
input_gtf='$input_gtf'
genomicSuperDups='$genomicSuperDups'
annotation='$annotation'
cdna_self='$cdna_self'
output_prefix='$output_prefix'

run_${tool}
"

  echo "[$(date)] Recording metrics for ${tool} on ${output_prefix}"
  write_metrics "$output_prefix" "$tool" "$input_fastq" "$output_dir" "$time_log" "$threads"
}


#############################
#Loop through tools and files
#############################
for depth in 1G 10G 100G; do #
  for seqid in Q85 Q90 Q95; do #
    for input_fastq in ${INPUT_DIR}/fastq_files/porechoptrimmed*${depth}*${seqid}.fastq.gz; do
        
        output_prefix=$(basename "${input_fastq}" .fastq.gz)
        input_bam="${ALIGNED_DIR}/${output_prefix}.sorted.bam"
        input_nbam="${ALIGNED_DIR}/${output_prefix}.name.sorted.bam"
        input_paf="${ALIGNED_DIR}/${output_prefix}.paf"
        
        [[ -f "$input_fastq" ]] || { echo "Missing fastq: $input_fastq"; continue; }
        [[ -f "$input_bam" ]] || { echo "Missing BAM: $input_bam"; continue; }
        [[ -f "$input_paf" ]] || { echo "Missing PAF: $input_paf"; continue; }
        [[ -f "$input_nbam" ]] || { echo "Missing name-sorted BAM: $input_nbam"; continue; }
        
        
          for tool in jaffal genion fusionseeker longgf; do

            Bench_Tool "$tool" "$output_prefix" "$input_fastq" "$threads" &
          
          done
      wait
      
    done
  done
done

echo ""
echo "============================================"
echo "Single-Threading Benchmark complete!"
echo "Metrics saved to: $METRICS_FILE"
echo "Continue on to threads = 10 Benchmark "
echo "============================================"

############################
# Multi-threading enabled  #
############################
#thread
threads=10
for depth in 1G 10G ; do #
  for seqid in 85 90 95; do #
    for input_fastq in ${INPUT_DIR}/fastq_files/porechoptrimmed*${depth}*Q${seqid}.fastq.gz; do
      
        output_prefix=$(basename "${input_fastq}" .fastq.gz)
        input_bam="${ALIGNED_DIR}/${output_prefix}.sorted.bam"
        input_nbam="${ALIGNED_DIR}/${output_prefix}.name.sorted.bam"
        input_paf="${ALIGNED_DIR}/${output_prefix}.paf"
        
        [[ -f "$input_bam" ]] || { echo "Missing BAM: $input_bam"; continue; }
        [[ -f "$input_paf" ]] || { echo "Missing PAF: $input_paf"; continue; }
        [[ -f "$input_nbam" ]] || { echo "Missing name-sorted BAM: $input_nbam"; continue; }

      
        for tool in jaffal fusionseeker; do
          Bench_Tool "$tool" "$output_prefix" "$input_fastq" "$threads" &
        done
    done
  done
done

echo ""
echo "============================================"
echo "Multi-threading Benchmark complete!"
echo "Metrics saved to: $METRICS_FILE"
echo "============================================"