#!/bin/bash
#
# CTAT-LR-fusion Run Script
# LongReadFusionCallerBenchmark Project
#
# Configuration (as specified by Ryley):
#   - Gencode v44 (v43 not available as pre-built CTAT library)
#   - minimap2 v2.28 (using /opt/JAFFA-version-2.3/tools/bin/minimap2)
#   - Minimum read support: 2
#   - Annotation filtering disabled (--no_annot_filter)
#

# ============ CONFIGURATION ============

# Tool paths
CTAT_DIR="/home/swu/tools/CTAT-LR-fusion"
MINIMAP2_PATH="/opt/JAFFA-version-2.3/tools/bin"

# Data paths
GENOME_LIB="/bioinformatics/siyuan/references/GRCh38_gencode_v44_CTAT_lib"
OUTPUT_BASE="/bioinformatics/siyuan/ctat_lr_results"
SIMULATED_DATA="/bioinformatics/ryley/Algorithm_Benchmark/Adapter_porechop_trimmed/fastq_files"
CELLLINE_DATA="/bioinformatics/ryley/Algorithm_Benchmark/SGNEx_files"
RASUSA_DATA="/bioinformatics/ryley/Algorithm_Benchmark/SGNEx_files/Rasusa_files"

# Benchmark output
BENCHMARK_CSV="${OUTPUT_BASE}/benchmark_metrics.csv"

# Run parameters
CPU_THREADS=16
MIN_NUM_LR=2  # Minimum read support 

# ============ FUNCTIONS ============

# Initialize benchmark CSV with header
init_benchmark_csv() {
    if [ ! -f "${BENCHMARK_CSV}" ]; then
        echo "sample_name,dataset_type,input_size_bytes,input_size_gb,size_category,quality,spiked,wall_time_seconds,user_time_seconds,system_time_seconds,cpu_time_seconds,peak_ram_kb,peak_ram_gb,output_size_bytes,output_size_mb,num_fusions,exit_status,timestamp" > "${BENCHMARK_CSV}"
        echo "Benchmark CSV initialized: ${BENCHMARK_CSV}"
    fi
}

# Run a single sample with benchmarking
run_sample() {
    local input_fastq=$1
    local output_dir=$2
    local sample_name=$(basename ${input_fastq} .fastq.gz)
    local time_log="${output_dir}/logs/${sample_name}.time"

    echo "============================================"
    echo "[$(date)] Starting analysis: ${sample_name}"
    echo "Input:  ${input_fastq}"
    echo "Output: ${output_dir}/${sample_name}"
    echo "============================================"

    # Get input file size
    local input_size_bytes=$(stat -c%s "${input_fastq}")
    local input_size_gb=$(echo "scale=2; ${input_size_bytes}/1024/1024/1024" | bc)

    # Determine dataset type and metadata
    local dataset_type="unknown"
    local size_category="NA"
    local quality="NA"
    local spiked="NA"

    # Detect metadata from sample_name (path-independent — works for files relocated outside their original dirs)
    if [[ "${sample_name}" == *"SGNEx"* ]] || [[ "${sample_name}" == *"SGNex"* ]]; then
        dataset_type="cellline"
        if [[ "${sample_name}" == *"_1G_"* ]]; then
            size_category="1G"
        elif [[ "${sample_name}" == *"_2.5G_"* ]]; then
            size_category="2.5G"
        elif [[ "${sample_name}" == *"_5G_"* ]]; then
            size_category="5G"
        elif [[ "${sample_name}" == *"_7.5G_"* ]]; then
            size_category="7.5G"
        elif [[ "${sample_name}" == *"_10G_"* ]]; then
            size_category="10G"
        elif [[ "${sample_name}" == *"combined"* ]] || [[ "${sample_name}" == *"concatenated"* ]]; then
            size_category="combined"
        fi
    elif [[ "${sample_name}" == *"Simulated"* ]] || [[ "${sample_name}" == *"Spiked"* ]]; then
        dataset_type="simulated"
        if [[ "${sample_name}" == *"_1GB_"* ]]; then
            size_category="1GB"
        elif [[ "${sample_name}" == *"_10GB_"* ]]; then
            size_category="10GB"
        elif [[ "${sample_name}" == *"_100GB_"* ]]; then
            size_category="100GB"
        fi
        if [[ "${sample_name}" == *"Q85"* ]]; then
            quality="Q85"
        elif [[ "${sample_name}" == *"Q90"* ]]; then
            quality="Q90"
        elif [[ "${sample_name}" == *"Q95"* ]]; then
            quality="Q95"
        fi
        if [[ "${sample_name}" == *"Spiked"* ]]; then
            spiked="yes"
        else
            spiked="no"
        fi
    fi

    # Set PATH to use JAFFA's minimap2 (v2.28)
    export PATH="${MINIMAP2_PATH}:$PATH"

    # Run with /usr/bin/time to capture metrics
    /usr/bin/time -v -o "${time_log}" \
        perl ${CTAT_DIR}/ctat-LR-fusion \
        -T ${input_fastq} \
        --genome_lib_dir ${GENOME_LIB} \
        --output ${output_dir}/${sample_name} \
        --CPU ${CPU_THREADS} \
        --min_num_LR ${MIN_NUM_LR} \
        --no_annot_filter \
        --no_ctat_mm2 \
        2>&1 | tee ${output_dir}/logs/${sample_name}.log

    # Capture perl's exit code, not tee's
    local exit_status=${PIPESTATUS[0]}

    # Parse time output
    local wall_time=$(grep "Elapsed (wall clock)" "${time_log}" | sed 's/.*: //')
    local user_time=$(grep "User time" "${time_log}" | awk '{print $NF}')
    local system_time=$(grep "System time" "${time_log}" | awk '{print $NF}')
    local peak_ram_kb=$(grep "Maximum resident set size" "${time_log}" | awk '{print $NF}')

    # Convert wall time to seconds (format: h:mm:ss or m:ss)
    local wall_seconds=0
    if [[ "${wall_time}" == *:*:* ]]; then
        # h:mm:ss format
        wall_seconds=$(echo "${wall_time}" | awk -F: '{print ($1*3600)+($2*60)+$3}')
    else
        # m:ss format
        wall_seconds=$(echo "${wall_time}" | awk -F: '{print ($1*60)+$2}')
    fi

    # Calculate CPU time and peak RAM in GB
    local cpu_time=$(echo "scale=2; ${user_time} + ${system_time}" | bc)
    local peak_ram_gb=$(echo "scale=2; ${peak_ram_kb}/1024/1024" | bc)

    # Get output size
    local output_size_bytes=0
    local output_size_mb=0
    if [ -d "${output_dir}/${sample_name}" ]; then
        output_size_bytes=$(du -sb "${output_dir}/${sample_name}" | cut -f1)
        output_size_mb=$(echo "scale=2; ${output_size_bytes}/1024/1024" | bc)
    fi

    # Count number of fusions detected
    local num_fusions=0
    local fusion_file="${output_dir}/${sample_name}/ctat-LR-fusion.fusion_predictions.tsv"
    if [ -f "${fusion_file}" ]; then
        num_fusions=$(tail -n +2 "${fusion_file}" | wc -l)
    elif [ "${exit_status}" -eq 0 ]; then
        # CTAT pipeline returned 0 but produced no predictions file — silent failure
        echo "WARNING: ${sample_name} has no ctat-LR-fusion.fusion_predictions.tsv — marking as failed"
        exit_status=2
    fi

    # Get timestamp
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')

    # Append to benchmark CSV
    echo "${sample_name},${dataset_type},${input_size_bytes},${input_size_gb},${size_category},${quality},${spiked},${wall_seconds},${user_time},${system_time},${cpu_time},${peak_ram_kb},${peak_ram_gb},${output_size_bytes},${output_size_mb},${num_fusions},${exit_status},${timestamp}" >> "${BENCHMARK_CSV}"

    echo "[$(date)] Finished: ${sample_name}"
    echo "  Wall time: ${wall_time} (${wall_seconds}s)"
    echo "  CPU time: ${cpu_time}s (user: ${user_time}s, system: ${system_time}s)"
    echo "  Peak RAM: ${peak_ram_kb} KB (${peak_ram_gb} GB)"
    echo "  Output size: ${output_size_mb} MB"
    echo "  Fusions detected: ${num_fusions}"
    echo ""
}

# Prepare reference index (run once before first analysis)
prep_reference() {
    echo "[$(date)] Preparing CTAT-LR-fusion reference index..."

    export PATH="${MINIMAP2_PATH}:$PATH"

    # Create temporary input file (script requires input file parameter)
    echo -e "@dummy\nACGT\n+\nIIII" > /tmp/dummy.fastq

    perl ${CTAT_DIR}/ctat-LR-fusion \
        -T /tmp/dummy.fastq \
        --genome_lib_dir ${GENOME_LIB} \
        --prep_reference_only \
        --CPU ${CPU_THREADS}

    rm -f /tmp/dummy.fastq
    echo "[$(date)] Reference index preparation complete"
}

# Display configuration
show_config() {
    echo ""
    echo "====== CTAT-LR-fusion Configuration ======"
    echo ""
    echo "Tool Paths:"
    echo "  CTAT-LR-fusion: ${CTAT_DIR}"
    echo "  minimap2:       ${MINIMAP2_PATH}/minimap2"
    echo ""
    echo "Data Paths:"
    echo "  Genome library:   ${GENOME_LIB}"
    echo "  Output directory: ${OUTPUT_BASE}"
    echo "  Simulated data:   ${SIMULATED_DATA}"
    echo "  Cell line data:   ${CELLLINE_DATA}"
    echo ""
    echo "Run Parameters:"
    echo "  CPU threads:            ${CPU_THREADS}"
    echo "  Min LR support:         ${MIN_NUM_LR}"
    echo "  Annotation filtering:   Disabled (--no_annot_filter)"
    echo "  External minimap2:      Enabled (--no_ctat_mm2)"
    echo ""
}

# Display help
show_help() {
    echo ""
    echo "CTAT-LR-fusion Run Script"
    echo ""
    echo "Usage: $0 <command> [arguments]"
    echo ""
    echo "Commands:"
    echo "  prep           - Prepare reference index (run once before first analysis)"
    echo "  simulated      - Run on all simulated data files (18 files)"
    echo "  cellline       - Run on all cell line data files (25 files)"
    echo "  single <fastq> - Run on a single sample (output to single/)"
    echo "  rerun <fastq> <subdir> - Run a single sample into OUTPUT_BASE/<subdir>/"
    echo "  all            - Run on all data"
    echo "  config         - Display current configuration"
    echo "  help           - Display this help message"
    echo ""
}

# ============ MAIN ============

# Create output directories
mkdir -p ${OUTPUT_BASE}/simulated/logs
mkdir -p ${OUTPUT_BASE}/cellline/logs
mkdir -p ${OUTPUT_BASE}/single/logs

case "$1" in
    prep)
        prep_reference
        ;;
    simulated)
        init_benchmark_csv
        echo "Processing simulated data (18 files)..."
        for fastq in ${SIMULATED_DATA}/*.fastq.gz; do
            run_sample "${fastq}" "${OUTPUT_BASE}/simulated"
        done
        ;;
    cellline)
        init_benchmark_csv
        echo "Processing cell line data (25 files)..."
        for fastq in ${CELLLINE_DATA}/*.fastq.gz; do
            run_sample "${fastq}" "${OUTPUT_BASE}/cellline"
        done
        ;;
    rasusa)
        mkdir -p ${OUTPUT_BASE}/cellline_rasusa/logs
        init_benchmark_csv
        echo "Processing Rasusa cell line data (29 files)..."
        for fastq in ${RASUSA_DATA}/*.fastq.gz; do
            run_sample "${fastq}" "${OUTPUT_BASE}/cellline_rasusa"
        done
        ;;
    single)
        if [ -z "$2" ]; then
            echo "Error: Please provide a FASTQ file path"
            echo "Usage: $0 single <fastq_file>"
            exit 1
        fi
        init_benchmark_csv
        run_sample "$2" "${OUTPUT_BASE}/single"
        ;;
    rerun)
        if [ -z "$2" ] || [ -z "$3" ]; then
            echo "Error: Please provide FASTQ path and output subdirectory"
            echo "Usage: $0 rerun <fastq_file> <output_subdir>"
            echo "Example: $0 rerun /path/to/cleaned.fastq.gz cellline_rasusa"
            exit 1
        fi
        mkdir -p ${OUTPUT_BASE}/$3/logs
        init_benchmark_csv
        run_sample "$2" "${OUTPUT_BASE}/$3"
        ;;
    all)
        echo "Processing all data..."
        $0 simulated
        $0 cellline
        ;;
    smart)
        init_benchmark_csv
        echo "Processing in smart order: cellline -> 1GB -> 10GB -> 100GB"
        echo ""

        # 1. Cell line data first
        echo "=== Step 1/4: Cell line data (25 files) ==="
        $0 cellline

        # 2. Simulated 1GB files
        echo "=== Step 2/4: Simulated 1GB files (6 files) ==="
        for fastq in ${SIMULATED_DATA}/*_1GB_*.fastq.gz; do
            run_sample "${fastq}" "${OUTPUT_BASE}/simulated"
        done

        # 3. Simulated 10GB files
        echo "=== Step 3/4: Simulated 10GB files (6 files) ==="
        for fastq in ${SIMULATED_DATA}/*_10GB_*.fastq.gz; do
            run_sample "${fastq}" "${OUTPUT_BASE}/simulated"
        done

        # 4. Simulated 100GB files
        echo "=== Step 4/4: Simulated 100GB files (6 files) ==="
        for fastq in ${SIMULATED_DATA}/*_100GB_*.fastq.gz; do
            run_sample "${fastq}" "${OUTPUT_BASE}/simulated"
        done

        echo "=== All processing complete! ==="
        ;;
    config)
        show_config
        ;;
    help|--help|-h|"")
        show_help
        show_config
        ;;
    *)
        echo "Unknown command: $1"
        show_help
        ;;
esac
