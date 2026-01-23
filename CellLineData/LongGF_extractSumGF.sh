#!/bin/bash
set -euo pipefail

DIR="/bioinformatics/ryley/Gencode44/SGNex_Fusionv44/LongGF"
LONGGF_SGNEX="LongGF.run.SGNex"
for file in "$DIR"/"${LONGGF_SGNEX}"*.log; do
    [[ -f "$file" ]] || continue

    #Output description
    output_prefix=$(basename ${file} .log )
    # Skip if output already exists
    if [[ -f "${DIR}/${output_prefix}.tsv" ]]; then
        echo "[SKIP] ${output_prefix}.tsv already exists"
        continue
    fi
    
    grep "SumGF" ${file} > "${DIR}/${output_prefix}.tsv"
done