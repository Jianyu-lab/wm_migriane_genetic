#!/bin/bash
# Transcriptome-wide association study (TWAS) fine-mapping using FOCUS
# Prioritize likely causal genes underlying the shared genetic architecture

# ==============================================================================
# Setup paths
# ==============================================================================
DATA_DIR="./data/cleaned_gwas"
OUT_DIR="./results/focus"
REF_DIR="./refdata/1000G_EUR_Phase3_plink"
TRAITS_LIST="./data/all_phenotypes.list"

mkdir -p ${OUT_DIR}

tissues="Brain_Amygdala Brain_Anterior_cingulate_cortex_BA24 Brain_Caudate_basal_ganglia Brain_Cerebellar_Hemisphere Brain_Cerebellum Brain_Cortex Brain_Frontal_Cortex_BA9 Brain_Hippocampus Brain_Hypothalamus Brain_Nucleus_accumbens_basal_ganglia Brain_Putamen_basal_ganglia Brain_Spinal_cord_cervical_c-1 Brain_Substantia_nigra"

for trait in $(cat ${TRAITS_LIST}); do
    for tissue in $tissues; do
        focus finemap \
            --sumstats ${DATA_DIR}/${trait}.sumstats.gz \
            --tissue ${tissue} \
            --reference ${REF_DIR}/1000G.EUR.QC \
            --out ${OUT_DIR}/${trait}_${tissue}_FOCUS \
            --pip-threshold 0.05
    done
done

# Post-processing: extract PIPs > 0.05 and compute combined PIP for overlap
# This step is typically done in R or Python using FOCUS output files