#!/bin/bash
# Transcriptome-wide association analyses using S-PrediXcan (MetaXcan framework)
# GTEx v8 MASHR across 13 brain-related tissues

# ==============================================================================
# Setup paths
# ==============================================================================
DATA_DIR="./data/cleaned_gwas"
OUT_DIR="./results/twas"
GTEX_DIR="./refdata/GTEx-V8/mashr"
SPREDIXCAN="software/MetaXcan/software/SPrediXcan.py"
TRAITS_LIST="./data/all_phenotypes.list" # List of migraine and WM traits

mkdir -p ${OUT_DIR}

tissues="Brain_Amygdala Brain_Anterior_cingulate_cortex_BA24 Brain_Caudate_basal_ganglia Brain_Cerebellar_Hemisphere Brain_Cerebellum Brain_Cortex Brain_Frontal_Cortex_BA9 Brain_Hippocampus Brain_Hypothalamus Brain_Nucleus_accumbens_basal_ganglia Brain_Putamen_basal_ganglia Brain_Spinal_cord_cervical_c-1 Brain_Substantia_nigra"

for trait in $(cat ${TRAITS_LIST}); do
    for tissue in $tissues; do
        python ${SPREDIXCAN} \
            --model_db_path ${GTEX_DIR}/mashr_${tissue}.db \
            --covariance ${GTEX_DIR}/mashr_${tissue}.txt.gz \
            --gwas_folder ${DATA_DIR} \
            --gwas_file_pattern "${trait}.txt.gz" \
            --snp_column SNP \
            --effect_allele_column A1 \
            --non_effect_allele_column A2 \
            --beta_column BETA \
            --pvalue_column P \
            --output_file ${OUT_DIR}/${trait}_${tissue}_SPrediXcan.csv
    done
done