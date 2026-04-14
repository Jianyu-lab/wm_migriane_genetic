#!/bin/bash
# Assessment of sample overlap and genetic correlation using bivariate LDSC

# ==============================================================================
# Setup paths
# ==============================================================================
DATA_DIR="./data/cleaned_gwas"
OUT_DIR="./results/ldsc"
REF_DIR="./refdata/eur_w_ld_chr"
TRAITS_LIST="./data/WM_phenotypes.list"
MIGRAINE_SUMSTATS="${DATA_DIR}/migraine_cleaned.sumstats.gz"

mkdir -p ${OUT_DIR}

# Loop over the 19 WM phenotypes
# MHC region excluded (chromosome 6, 26-34 Mb)
for trait in $(cat ${TRAITS_LIST}); do
    # Bivariate LDSC to estimate genetic correlation and cross-trait intercept
    python ldsc.py \
        --rg ${MIGRAINE_SUMSTATS},${DATA_DIR}/${trait}_cleaned.sumstats.gz \
        --ref-ld-chr ${REF_DIR}/ \
        --w-ld-chr ${REF_DIR}/ \
        --out ${OUT_DIR}/migraine_vs_${trait}
done