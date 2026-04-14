#!/bin/bash
# Local genetic correlation analysis using SUPERGNOVA
# Applied to primary migraine GWAS and 19 WM phenotypes defined by combined conjFDR loci

# ==============================================================================
# Setup paths
# ==============================================================================
DATA_DIR="./data/cleaned_gwas"
OUT_DIR="./results/supergnova"
REF_DIR="./refdata/1000G_EUR_Phase3_plink"
PARTITION_BED="./data/conjFDR_loci.bed" # BED file containing combined conjFDR loci
MIGRAINE_SUMSTATS="${DATA_DIR}/migraine.sumstats.gz"
TRAITS_LIST="./data/WM_phenotypes.list"

mkdir -p ${OUT_DIR}

for trait in $(cat ${TRAITS_LIST}); do
    python software/SUPERGNOVA/supergnova.py \
        --bfile ${REF_DIR}/1000G.EUR.QC \
        --sumstats1 ${MIGRAINE_SUMSTATS} \
        --sumstats2 ${DATA_DIR}/${trait}.sumstats.gz \
        --partition ${PARTITION_BED} \
        --out ${OUT_DIR}/migraine_vs_${trait}_supergnova
done