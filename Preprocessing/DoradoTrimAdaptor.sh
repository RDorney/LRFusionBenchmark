#!/bin/bash
#direct cDNA nanopore seq use the  LSK-114 kit
direct_cDNA_Directory=$(/bioinformatics/ryley/Library_Benchmark/raw_files/Huh7_Library_Comparison_direct_cDNA/)

for reads in $(${direct_cDNA_Directory}*/*/F*fastq.gz); do
output_prefix=$(basename ${reads} .fastq.gz)
/opt/dorado-0.9.0-linux-x64/bin/dorado trim --sequencing-kit SQK_LSK_114 --emit-fastq ${reads} > ${direct_cDNA_Directory}_dorado_trimmed_${output_prefix}.fastq
gzip ${direct_cDNA_Directory}_dorado_trimmed_${output_prefix}.fastq &
done

#PCR cDNA nanopore seq use the LSK-114 kit
PCR_cDNA_Directory=$(/bioinformatics/ryley/Library_Benchmark/raw_files/Huh-7_Library_Comparision_PCR-cDNA/)

for reads in $(${PCR_cDNA_Directory}*/*/F*fastq.gz); do
output_prefix=$(basename ${reads} .fastq.gz)
/opt/dorado-0.9.0-linux-x64/bin/dorado trim --sequencing-kit SQK_PCS114 --emit-fastq ${reads} > ${PCR_cDNA_Directory}_dorado_trimmed_${output_prefix}.fastq
gzip ${PCR_cDNA_Directory}_dorado_trimmed_${output_prefix}.fastq &
done

#direct RNA seq trims reads and primers while it is sequencing, so we do not need to trim here