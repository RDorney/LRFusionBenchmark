#!/bin/bash
set -e
export JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64
export PATH="$JAVA_HOME/bin:$PATH"

JAFFA_DIR=/opt/JAFFA-version-2.3
DATA_DIR=/bioinformatics/ryley/Library_Benchmark

mkdir -p ${DATA_DIR}/JAFFA_direct_gencode43
cd ${DATA_DIR}/JAFFA_direct_gencode43 || exit 1

ln -sf ${DATA_DIR}/raw_files/AGRF_Sequencing_Data/Illumina_NovaSeq/fastp_trimmed Illumina_reads

##################
# Important note #
##################
# Make sure the fastq file names are in "%_*.fastq.gz" format. 
# Caution if :
#   more than one or no "_" in name

${JAFFA_DIR}/tools/bin/bpipe run -p genome=hg38 -p annotation=genCode43 -p exclude="NoSupport" ${JAFFA_DIR}/JAFFA_direct.groovy Illumina_reads/*.fastq.gz 

unlink Illumina_reads
cd ${DATA_DIR}/JAFFA_direct_gencode43 || exit 1

################
# Message
################
echo "it's all ok :)"
