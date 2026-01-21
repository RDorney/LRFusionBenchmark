#!/bin/bash
##############################
# short cuts to Huh7 seq data
##############################

HUH7_FILES=/bioinformatics/ryley/Library_Benchmark/raw_files
############################
# Illumina reads
############################
ln -sf ${HUH7_FILES}/AGRF_Sequencing_Data/Illumina_NovaSeq/fastp_trimmed Illumina_reads
# Huh7B1_R1.fastq.gz  Huh7B1_R2.fastq.gz  Huh7B2_R1.fastq.gz  Huh7B2_R2.fastq.gz 

############################
# PacBio Reads
############################
ln -sf ${HUH7_FILES}/AGRF_Sequencing_Data/PacBio_Revio_full_length_cDNA ~/PacBio_reads
# The Pacbio fastqs we are interested in is PB1_FLNCcDNA.fastq.gz and PB2_FLNCcDNA.fastq/gz

############################
# Nano_reads
############################
mkdir -p ~/Nano_reads
ln -sf ${HUH7_FILES}/Huh-7_Library_Comparision_PCR-cDNA/Huh7_p9_8_7_B1-A2/20241009_1331_MN44892_FAZ80247_03069096/chopperfiltered_dorado_trimmed_FAZ80247_pass_03069096_5394ef1c_concatenated.fastq.gz ~/Nano_reads/Nano_PCR_B1.fastq.gz
ln -sf ${HUH7_FILES}/Huh-7_Library_Comparision_PCR-cDNA/Huh7_p9_8-7_B2/20241024_1501_MN44892_FAZ83542_1344a466/chopperfiltered_dorado_trimmed_FAZ83542_pass_1344a466_45031b9d_concatenated.fastq.gz ~/Nano_reads/Nano_PCR_B2.fastq.gz
ln -sf ${HUH7_FILES}/Huh7_Library_Comparison_direct_cDNA/Huh7_p9_8-7_B1/20241204_1527_MN44892_FBA62655_9d20541a/chopperfiltered_dorado_trimmed_FBA62655_pass_9d20541a_79aab80e_concatenated.fastq.gz ~/Nano_reads/Nano_dcDNA_B1.fastq.gz
ln -sf ${HUH7_FILES}/Huh7_Library_Comparison_direct_cDNA/Huh7_p9_8-7_B2/20241209_1659_MN44892_FBA43334_6846ed3b/chopperfiltered_dorado_trimmed_FBA43334_pass_6846ed3b_6b29293f_concatenated.fastq.gz ~/Nano_reads/Nano_dcDNA_B2.fastq.gz
ln -sf ${HUH7_FILES}/Huh-7_Library_Comparison/Huh7_p9_8-7_B1/20240919_1733_MN44892_FBA22660_d1e1158a/chopperfiltered_FBA22660_pass_d1e1158a_77e47f3b_concatenated.fastq.gz ~/Nano_reads/Nano_dRNA_B1.fastq.gz
ln -sf ${HUH7_FILES}/Huh-7_Library_Comparison/Huh7_p9_8-7_B2/20240924_1732_MN44892_FBA22517_d17fa1e0/chopperfiltered_FBA22517_pass_d17fa1e0_a6486b6e_concatenated.fastq.gz ~/Nano_reads/Nano_dRNA_B2.fastq.gz
