#!/bin/bash
set -e
export JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64
export PATH="$JAVA_HOME/bin:$PATH"

JAFFA_DIR=/opt/JAFFA-version-2.3
DATA_DIR=/bioinformatics/ryley/Gencode44/Huh7_Library

mkdir -p ${DATA_DIR}/JAFFAdirect_gencode44
cd ${DATA_DIR}/JAFFAdirect_gencode44 || exit 1

ln -sf /bioinformatics/ryley/Library_Benchmark/raw_files/AGRF_Sequencing_Data/Illumina_NovaSeq/fastp_trimmed Illumina_reads

##################
# Important note #
##################
# Make sure the fastq file names are in "%_*.fastq.gz" format. 
# Caution if :
#   more than one or no "_" in name

${JAFFA_DIR}/tools/bin/bpipe run -n 4 -p genome=hg38 -p annotation=genCode44 -p exclude="NoSupport" ${JAFFA_DIR}/JAFFA_direct.groovy Illumina_reads/*.fastq.gz 

unlink Illumina_reads
cd ${DATA_DIR}/JAFFAdirect_gencode44 || exit 1

################
# Message
################
echo "it's all ok :)"
