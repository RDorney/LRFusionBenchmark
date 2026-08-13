#!/bin/bash
module load fastp/0.23.4

for index in "B2_E1_22T7FFLT3_ACGGTCAGGA-CGGCCTCGTT" "B1_E1_22T7FFLT3_CATTAGTGCG-TGACTACATA"; do 
    for lane in L005 L006 L007 L008; do

        R1_name="Huh-7_p9_8-7-24_${index}_${lane}_R1"
        R2_name="Huh-7_p9_8-7-24_${index}_${lane}_R2"
        R1_FASTQ="../Huh-7_p9_8-7-24_${index}_${lane}_R1.fastq.gz"
        R2_FASTQ="../Huh-7_p9_8-7-24_${index}_${lane}_R2.fastq.gz"

        fastp -i $R1_FASTQ -I $R2_FASTQ \
            --dont_overwrite \
            -o fastp.${R1_name}.fastq.gz -O fastp.${R2_name}.fastq.gz \
            --unpaired1 unpaired.${R1_name}.fastq.gz --unpaired2 unpaired.${R2_name}.fastq.gz \
            -h "Huh-7_p9_8-7-24_${index}_${lane}_fastp.html" \
            -j "Huh-7_p9_8-7-24_${index}_${lane}_fastp.json" \
            -l 36 --trim_poly_g \
            --adapter_sequence=CTGTCTCTTATACACATCT --adapter_sequence_r2=CTGTCTCTTATACACATCT

    done
done

module load fastqc/0.12.1
fastqc *.fastq.gz

module load multiqc/1.18
multiqc ../* ./*
