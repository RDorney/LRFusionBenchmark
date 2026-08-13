#!/bin/bash
WORK_DIR="$(pwd)"
echo "Working directory is $(pwd)"

###############################################################################
# Input files and references
###############################################################################
RAW_FILES=/bioinformatics/ryley/Library_Benchmark/raw_files
mkdir fastq_dir
DATA_DIR="${WORK_DIR}/fastq_dir"
cd fastq_dir
ln -sf ${RAW_FILES}/AGRF_Sequencing_Data/PacBio_Revio_full_length_cDNA ${DATA_DIR}/PacBioReads

mkdir -p ${DATA_DIR}/ONTReads
ln -sf ${RAW_FILES}/Huh-7_Library_Comparision_PCR-cDNA/Huh7_p9_8_7_B1-A2/20241009_1331_MN44892_FAZ80247_03069096/chopperfiltered_dorado_trimmed_FAZ80247_pass_03069096_5394ef1c_concatenated.fastq.gz ${DATA_DIR}/ONTReads/Nano_PCR_B1.fastq.gz
ln -sf ${RAW_FILES}/Huh-7_Library_Comparision_PCR-cDNA/Huh7_p9_8-7_B2/20241024_1501_MN44892_FAZ83542_1344a466/chopperfiltered_dorado_trimmed_FAZ83542_pass_1344a466_45031b9d_concatenated.fastq.gz ${DATA_DIR}/ONTReads/Nano_PCR_B2.fastq.gz
ln -sf ${RAW_FILES}/Huh7_Library_Comparison_direct_cDNA/Huh7_p9_8-7_B1/20241204_1527_MN44892_FBA62655_9d20541a/chopperfiltered_dorado_trimmed_FBA62655_pass_9d20541a_79aab80e_concatenated.fastq.gz ${DATA_DIR}/ONTReads/Nano_dcDNA_B1.fastq.gz
ln -sf ${RAW_FILES}/Huh7_Library_Comparison_direct_cDNA/Huh7_p9_8-7_B2/20241209_1659_MN44892_FBA43334_6846ed3b/chopperfiltered_dorado_trimmed_FBA43334_pass_6846ed3b_6b29293f_concatenated.fastq.gz ${DATA_DIR}/ONTReads/Nano_dcDNA_B2.fastq.gz
ln -sf ${RAW_FILES}/Huh-7_Library_Comparison/Huh7_p9_8-7_B1/20240919_1733_MN44892_FBA22660_d1e1158a/chopperfiltered_FBA22660_pass_d1e1158a_77e47f3b_concatenated.fastq.gz ${DATA_DIR}/ONTReads/Nano_dRNA_B1.fastq.gz
ln -sf ${RAW_FILES}/Huh-7_Library_Comparison/Huh7_p9_8-7_B2/20240924_1732_MN44892_FBA22517_d17fa1e0/chopperfiltered_FBA22517_pass_d17fa1e0_a6486b6e_concatenated.fastq.gz ${DATA_DIR}/ONTReads/Nano_dRNA_B2.fastq.gz

cd "${WORK_DIR}"

ONT_DIR=${DATA_DIR}/ONTReads
PB_DIR=${DATA_DIR}/PacBioReads

REF_DIR=/bioinformatics/ryley/Gencode44/reference_v44
reference_genome=${REF_DIR}/GRCh38.primary_assembly.genome.fa
input_gtf=${REF_DIR}/gencode.v44.annotation.gtf
#specific to genion
genomicSuperDups=/bioinformatics/ryley/reference_files/genion_reference_files/genomicSuperDups.txt
annotation=${REF_DIR}/Homo_sapiens.GRCh38.110.chr.gtf
cdna_self=${REF_DIR}/cdna.GRCh38v110.selfalign.tsv
###############################################################################
# Tool directory
###############################################################################
# Caller versions:
#   FusionSeeker v1.0.1   (/opt/FusionSeeker/fusionseeker)
#   Genion       v1.1.0   (/opt/genion/genion)
#   LongGF       v0.1.2   (bioconda longgf-0.1.2-h28e74a2_4, in the "longgf" environment)
#   genion-env supplies runtime only: python 3.9.7, minimap2 2.23
FS_CMD="/opt/FusionSeeker/fusionseeker"
Gen_CMD="/opt/genion/genion"

###############################################################################
# Threading
###############################################################################

ALIGN_THREADS=3
threads=3

###############################################################################
# Alignment tools
###############################################################################

M2_CMD="/opt/JAFFA-version-2.3/tools/bin/minimap2"   # minimap2 2.28-r1209
#minimap2 should be used as ${M2_CMD}  

###############################################################################
# Alignment step
###############################################################################

# PacBio alignment (backgrounded)
(
  for input_fastq_file in ${PB_DIR}/*fastq.gz; do
  output_prefix=$(basename "${input_fastq_file}" .fastq.gz)
  
  if [[ -f "${output_prefix}.sorted.bam" \
        && -f "${output_prefix}.sorted.bam.bai" \
        && -f "${output_prefix}.name.sorted.bam" ]]; then
  echo "Skipping ${output_prefix} (all BAM outputs exist)"
  continue
  fi
  
 ${M2_CMD} -cx splice:hq -t ${ALIGN_THREADS} \
  "${reference_genome}" "${input_fastq_file}" > "${output_prefix}.paf"
  
 ${M2_CMD} -ax splice:hq -t ${ALIGN_THREADS} \
  "${reference_genome}" "${input_fastq_file}" \
  | samtools sort -@ ${threads} -O BAM --write-index \
  -o "${output_prefix}.sorted.bam"
  
  samtools sort -n -@ ${threads} \
  "${output_prefix}.sorted.bam" -O BAM \
  -o "${output_prefix}.name.sorted.bam"
  done
) &
  
  # Nanopore alignment generation (backgrounded)
for input_fastq_file in ${ONT_DIR}/*fastq.gz; do
    output_prefix=$(basename "${input_fastq_file}" .fastq.gz)
    
    if [[ -f "${output_prefix}.paf" \
          && -f "${output_prefix}.sorted.bam" \
          && -f "${output_prefix}.sorted.bam.bai" \
          && -f "${output_prefix}.name.sorted.bam" ]]; then
      echo "Skipping ${output_prefix} (all alignment outputs exists)"
      continue
    fi
    
   ${M2_CMD} -cx splice -t ${ALIGN_THREADS} \
    "${reference_genome}" "${input_fastq_file}" > "${output_prefix}.paf"

  ${M2_CMD} -ax splice -t ${ALIGN_THREADS} \
  "${reference_genome}" "${input_fastq_file}" \
  | samtools sort -@ ${threads} -O BAM --write-index \
  -o "${output_prefix}.sorted.bam"
  
  samtools sort -n -@ ${threads} \
  "${output_prefix}.sorted.bam" -O BAM \
  -o "${output_prefix}.name.sorted.bam"
done

###############################################################################
# FusionSeeker
###############################################################################
min_read_supp=2

source "$(conda info --base)/etc/profile.d/conda.sh"

conda activate FusionSeeker

# ONT FusionSeeker
for input_fastq in ${ONT_DIR}/*fastq.gz; do
output_prefix=$(basename "${input_fastq}" .fastq.gz)
outdir="fusionseeker_${output_prefix}"

if [[ -d "${outdir}" ]]; then
echo "Skipping FusionSeeker for ${output_prefix}"
continue
fi

${FS_CMD} \
--bam "${output_prefix}.sorted.bam" \
-o "${outdir}/" \
--datatype nanopore \
--ref "${reference_genome}" \
--gtf "${input_gtf}" \
--geneid \
--thread "${threads}" \
--minsupp "${min_read_supp}"
done

wait

# PacBio FusionSeeker
for input_fastq in ${PB_DIR}/*fastq.gz; do
output_prefix=$(basename "${input_fastq}" .fastq.gz)
outdir="fusionseeker_${output_prefix}"

if [[ -d "${outdir}" ]]; then
echo "Skipping FusionSeeker for ${output_prefix}"
continue
fi

${FS_CMD} \
--bam "${output_prefix}.sorted.bam" \
-o "${outdir}/" \
--datatype isoseq \
--ref "${reference_genome}" \
--gtf "${input_gtf}" \
--geneid \
--thread "${threads}" \
--minsupp "${min_read_supp}"
done

conda deactivate

###############################################################################
# Samtools stats (backgrounded)
###############################################################################

(
  for input_fastq_file in ${ONT_DIR}/*fastq.gz ${PB_DIR}/*fastq.gz; do
  output_prefix=$(basename "${input_fastq_file}" .fastq.gz)
  
  samtools stats "${output_prefix}.sorted.bam" \
  > "${output_prefix}.samtools_stats"
  
  samtools flagstat "${output_prefix}.sorted.bam" -O tsv \
  > "${output_prefix}.samtools_flagstat"
  
  samtools idxstats "${output_prefix}.sorted.bam" \
  > "${output_prefix}.samtools_idxstat"
  done
) &
  
###############################################################################
# Genion
###############################################################################

conda activate genion-env

for input_fastq in ${ONT_DIR}/*fastq.gz ${PB_DIR}/*fastq.gz; do
output_prefix=$(basename "${input_fastq}" .fastq.gz)

if [[ -s "${output_prefix}_genion" ]]; then
echo "Skipping Genion for ${output_prefix}"
continue
fi

${Gen_CMD} \
-i "${input_fastq}" \
-d "${genomicSuperDups}" \
--gtf "${annotation}" \
-g "${output_prefix}.paf" \
-s "${cdna_self}" \
-o "${output_prefix}_genion" \
--non-coding \
--min-support "${min_read_supp}"
done

conda deactivate

###############################################################################
# LongGF
###############################################################################

conda activate longgf

for input_fastq in ${ONT_DIR}/*fastq.gz ${PB_DIR}/*fastq.gz; do
output_prefix=$(basename "${input_fastq}" .fastq.gz)

if [[ -s "LongGF.run.${output_prefix}.log" ]]; then
echo "Skipping LongGF for ${output_prefix}"
continue
fi

LongGF \
"${output_prefix}.name.sorted.bam" \
"${input_gtf}" \
100 50 200 \
min_sup_read:"${min_read_supp}" \
> "LongGF.run.${output_prefix}.log"
done

conda deactivate

###############################################################################
# Final wait
###############################################################################

wait
echo "All processing completed"
