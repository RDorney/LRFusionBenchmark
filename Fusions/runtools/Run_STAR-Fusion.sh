#!/bin/bash
# STAR 2.7.11b; reference is the CTAT genome library built on GENCODE v44.
#Set up STAR-Fusion
source ~/anaconda3/etc/profile.d/conda.sh
conda activate starfusion_env

PROJ_DIR="/bioinformatics/ryley/Gencode44/Huh7_Library"
#make a directory to store results
WORK_DIR_STARF="$PROJ_DIR/STARFusion"
mkdir -p ${WORK_DIR_STARF}
cd ${WORK_DIR_STARF}

READ_DIR="${PROJ_DIR}/Illumina_reads"

#Set up CTAT Resource Library
#wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play.tar.gz
#tar -zxvf GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play.tar.gz
CTAT_GENOME_LIB="/bioinformatics/ryley/HCC_cohort_Fusions/STARFusion/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play/ctat_genome_lib_build_dir"

#Set up sample names
THREADS=10
for sample_number in B1 B2 ; do 
    R1_name="Huh7${sample_number}_R1"
    R2_name="Huh7${sample_number}_R2"
    R1_FASTQ="${R1_name}.fastq.gz"
    R2_FASTQ="${R2_name}.fastq.gz"
    OUT_PREFIX="Huh7${sample_number}"
    
#Run STAR-Fusion
STAR-Fusion --genome_lib_dir ${CTAT_GENOME_LIB} \
            --left_fq ${READ_DIR}/${R1_FASTQ} \
            --right_fq ${READ_DIR}/${R2_FASTQ} \
            --output_dir starfusion_${OUT_PREFIX} \
            --examine_coding_effect \
            --FusionInspector validate \
            --CPU ${THREADS} 
done
