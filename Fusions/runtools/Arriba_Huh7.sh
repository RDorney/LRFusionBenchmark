#!/bin/bash
#Set up Arriba
##Activate Arriba Conda Environment
set -Eeuo pipefail
source ~/anaconda3/etc/profile.d/conda.sh 
conda activate Arriba251 

CONDA_PREFIX="$HOME/anaconda3/envs/Arriba251"
ARRIBA_FILES=$CONDA_PREFIX/var/lib/arriba

## Set up directories
PROJ_DIR="/bioinformatics/ryley/Gencode44/Huh7_Library"
mkdir -p ${PROJ_DIR}/Arriba_Huh7
cd ${PROJ_DIR}/Arriba_Huh7 || exit 1

WORK_DIR_ARRIBA="$PROJ_DIR/Arriba_Huh7"
cd "$WORK_DIR_ARRIBA"

###Set up ref shortcuts
GENCODE44_REF_DIR="/bioinformatics/ryley/HCC_cohort_Fusions/Arriba/Gencode44_ref"
GENOME_FASTA="$GENCODE44_REF_DIR/GRCh38.primary_assembly.genome.fa"  #"$WORK_DIR_ARRIBA/GRCh38.fa"
GTF_FILE="$GENCODE44_REF_DIR/gencode.v44.annotation.gtf"  #"$WORK_DIR_ARRIBA/GENCODE38.gtf"
READ_DIR="${PROJ_DIR}/Illumina_reads"
THREADS=4
GEN_INDEX_DIR="$GENCODE44_REF_DIR/STAR_index_GRCh38_GENCODE44"  #"$WORK_DIR_ARRIBA/STAR_index_GRCh38_GENCODE38"

## Run Arriba on AITHM HCC cohort
cd $WORK_DIR_ARRIBA
RESULTS_DIR="$WORK_DIR_ARRIBA/results"
mkdir -p "$RESULTS_DIR"
cd "$RESULTS_DIR"

#Align the RNA-Seq Reads
for sample_number in B1 B2; do
        R1_name="Huh7${sample_number}_R1"
        R2_name="Huh7${sample_number}_R2"
        R1_FASTQ="${R1_name}.fastq.gz"
        R2_FASTQ="${R2_name}.fastq.gz"
        OUT_PREFIX="Huh7${sample_number}"
        
        STAR \
          --runThreadN $THREADS \
        	--genomeDir $GEN_INDEX_DIR --genomeLoad NoSharedMemory \
        	--readFilesIn $READ_DIR/$R1_FASTQ $READ_DIR/$R2_FASTQ --readFilesCommand zcat \
        	--outFileNamePrefix ${OUT_PREFIX}. \
        	--outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 \
          --outFilterMultimapNmax 50 \
          --peOverlapNbasesMin 10 \
          --alignSplicedMateMapLminOverLmate 0.5 \
          --alignSJstitchMismatchNmax 5 -1 5 5 \
          --chimSegmentMin 10 \
          --chimOutType WithinBAM HardClip \
          --chimJunctionOverhangMin 10 \
          --chimScoreDropMax 30 \
          --chimScoreJunctionNonGTAG 0 \
          --chimScoreSeparation 1 \
          --chimSegmentReadGapMax 3 \
          --chimMultimapNmax 50 \
          | arriba \
          -x /dev/stdin \
        	-g $GTF_FILE -a $GENOME_FASTA \
          -b $ARRIBA_FILES/blacklist_hg38_GRCh38_v2.5.1.tsv.gz \
          -k $ARRIBA_FILES/known_fusions_hg38_GRCh38_v2.5.1.tsv.gz \
          -t $ARRIBA_FILES/known_fusions_hg38_GRCh38_v2.5.1.tsv.gz \
          -p $ARRIBA_FILES/protein_domains_hg38_GRCh38_v2.5.1.gff3 \
          -o ${OUT_PREFIX}_fusions.tsv -O ${OUT_PREFIX}_fusions.discarded.tsv 
done

