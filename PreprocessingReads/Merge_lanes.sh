#!/bin/bash
###############################################################################
# Merge the per-lane fastp outputs into one FASTQ pair per replicate.
#
# This step sits between Fastp_trim_and_filter.sh, which writes four lanes per
# replicate, and Arriba_Huh7.sh / Run_STAR-Fusion.sh, which each expect a single
# R1/R2 pair named Huh7B1_* and Huh7B2_*.
#
#   in   fastp.Huh-7_p9_8-7-24_<index>_L00{5,6,7,8}_R{1,2}.fastq.gz
#   out  Huh7B{1,2}_R{1,2}.fastq.gz
#
# Concatenating gzip files is valid: the result is a multi-member gzip stream
# that zcat, STAR and Arriba all read as one file. There is no need to
# decompress and recompress, so do not simplify this to a zcat pipeline.
#
# Lanes are concatenated in ascending order. Read order does not affect the
# aligners' output, but fixing the order keeps the step reproducible.
###############################################################################
set -euo pipefail

FASTP_DIR="${1:-.}"     # directory holding the fastp.* outputs
OUT_DIR="${2:-.}"       # directory the aligner scripts read from

# Replicate label -> the index string used in the sequencing file names
declare -A INDEX=(
  [B1]="B1_E1_22T7FFLT3_CATTAGTGCG-TGACTACATA"
  [B2]="B2_E1_22T7FFLT3_ACGGTCAGGA-CGGCCTCGTT"
)

for rep in B1 B2; do
  for read in R1 R2; do

    lane_files=()
    for lane in L005 L006 L007 L008; do
      f="${FASTP_DIR}/fastp.Huh-7_p9_8-7-24_${INDEX[$rep]}_${lane}_${read}.fastq.gz"
      if [[ ! -f "$f" ]]; then
        echo "missing lane file: $f" >&2
        exit 1
      fi
      lane_files+=("$f")
    done

    out="${OUT_DIR}/Huh7${rep}_${read}.fastq.gz"
    echo "[$(date '+%F %T')] ${rep} ${read}: merging 4 lanes -> $(basename "$out")"
    cat "${lane_files[@]}" > "$out"

  done
done

echo "[$(date '+%F %T')] done"
