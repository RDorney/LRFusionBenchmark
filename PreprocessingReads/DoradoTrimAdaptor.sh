#!/bin/bash
set -e
# direct cDNA nanopore seq uses the SQK-LSK114 kit
direct_cDNA_Directory=/bioinformatics/ryley/Library_Benchmark/raw_files/Huh7_Library_Comparison_direct_cDNA/

for reads in ${direct_cDNA_Directory}*/*/F*fastq.gz; do
output_prefix=$(basename ${reads} .fastq.gz)
output_folder=$(dirname ${reads})
/opt/dorado-0.9.0-linux-x64/bin/dorado trim --sequencing-kit SQK-LSK114 --emit-fastq ${reads} > ${output_folder}/dorado_trimmed_${output_prefix}.fastq
gzip ${output_folder}/dorado_trimmed_${output_prefix}.fastq &
done
wait

# PCR cDNA nanopore seq uses the SQK-PCS114 kit
PCR_cDNA_Directory=/bioinformatics/ryley/Library_Benchmark/raw_files/Huh-7_Library_Comparision_PCR-cDNA/

for reads in ${PCR_cDNA_Directory}*/*/F*fastq.gz; do
output_prefix=$(basename ${reads} .fastq.gz)
output_folder=$(dirname ${reads})
/opt/dorado-0.9.0-linux-x64/bin/dorado trim --sequencing-kit SQK-PCS114 --emit-fastq ${reads} > ${output_folder}/dorado_trimmed_${output_prefix}.fastq
gzip ${output_folder}/dorado_trimmed_${output_prefix}.fastq &
done
wait

#direct RNA seq trims reads and primers while it is sequencing, so we do not need to trim here