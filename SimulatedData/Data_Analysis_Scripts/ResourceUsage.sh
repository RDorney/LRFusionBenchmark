#!/bin/bash
set -e
############################
# Configuration
############################
#set directory names
REF_DIR=/bioinformatics/ryley/Gencode44/reference_v44
DATA_DIR=/bioinformatics/ryley/Algorithm_Benchmark
cd ${DATA_DIR}

INPUT_DIR=${DATA_DIR}/Adapter_porechop_trimmed
ALIGNED_DIR=${INPUT_DIR}/alignment_and_samtools

LOGDIR="resource_logs"
mkdir -p "$LOGDIR"

#thread
threads=1

#shared default values
min_read_supp=2

############################
# Reference files
############################
reference_genome=${REF_DIR}/GRCh38.primary_assembly.genome.fa.gz 
input_gtf=${REF_DIR}/gencode.v44.annotation.gtf
#specific to Genion
genomicSuperDups=${REF_DIR}/genomicSuperDups.txt
annotation=${REF_DIR}/Homo_sapiens.GRCh38.110.chr.gtf
cdna_self=${REF_DIR}cdna.GRCh38v110.selfalign.tsv

############################
# Tool runners
############################
run_longgf () {
  LongGF "$input_nbam" "$input_gtf" 100 50 200 \
    min_sup_read:"$min_read_supp" \
    > "$stdout_log" 2> "$stderr_log"
}

run_genion () {
  genion -i "$input_fastq" \
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
  fusionseeker \
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
  /opt/JAFFA-version-2.3/tools/bin/bpipe run -n "$threads" \
    -p genome=hg38 \
    -p annotation=genCode44 \
    -p exclude=NoSupport \
    -p jaffa_output="jaffal_results_${output_prefix}" \
    /opt/JAFFA-version-2.3/JAFFAL.groovy \
    "$input_fastq" \
    > "$stdout_log" 2> "$stderr_log"
}
#############################
#Loop through tools and files
#############################
for depth in 1 10 ; do #
  for seqid in 85 90 95; do #
    for input_fastq in ${INPUT_DIR}/fastq_files/porechoptrimmed*${depth}*Q${seqid}.fastq.gz; do
      
        output_prefix=$(basename "${input_fastq}" .fastq.gz)
        input_bam="${ALIGNED_DIR}/bam_files/${output_prefix}.sorted.bam"
        input_nbam="${ALIGNED_DIR}/bam_files/${output_prefix}.name.sorted.bam"
        input_paf="${ALIGNED_DIR}/paf_files/${output_prefix}.paf"
      
        for tool in longgf genion fusionseeker jaffal; do

          stdout_log="${LOGDIR}/${tool}.${output_prefix}.out.log"
          stderr_log="${LOGDIR}/${tool}.${output_prefix}.err.log"
          time_log="${LOGDIR}/${tool}.${output_prefix}.time.log"

          echo "[$(date)] Running ${tool} on ${output_prefix}"
  
          /usr/bin/time -v -o "$time_log" "run_${tool}" &
  
          PID=$!
  
          wait "$PID"
      done
    done
  done
done
