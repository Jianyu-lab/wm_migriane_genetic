#!/bin/bash
# Bivariate MiXeR analyses for characterizing polygenic overlap between migraine and 19 WM phenotypes
# Polygenic overlap quantified using the Dice coefficient

# ==============================================================================
# Setup paths
# ==============================================================================
DATA_DIR="./data/traitfolder"
OUT_DIR="./results/mixer_bivar"
REF_DIR="./refdata/1000G_EUR_Phase3_plink"
MIXER_PY="software/mixer/precimed/mixer.py"
MIXER_FIG="software/mixer/precimed/mixer_figures.py"
LIB_PATH="software/mixer/src/build/lib/libbgmg.so"
WM_LIST="./data/WM_phenotypes.list"
PRIMARY_TRAIT="migraine"

mkdir -p ${OUT_DIR}

for wm_trait in $(cat ${WM_LIST}); do
    for i in $(seq 1 20); do
        # Fit bivariate model
        python3.8 ${MIXER_PY} fit2 \
            --trait1-file ${DATA_DIR}/${wm_trait}.csv.gz \
            --trait2-file ${DATA_DIR}/${PRIMARY_TRAIT}.csv.gz \
            --trait1-params-file ./results/mixer/${wm_trait}.fit.rep${i}.json \
            --trait2-params-file ./results/mixer/${PRIMARY_TRAIT}.fit.rep${i}.json \
            --out ${OUT_DIR}/${PRIMARY_TRAIT}_vs_${wm_trait}.fit.rep${i} \
            --extract ${REF_DIR}/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${i}.snps \
            --bim-file ${REF_DIR}/1000G.EUR.QC.@.bim \
            --ld-file ${REF_DIR}/1000G.EUR.QC.@.run4.ld \
            --lib ${LIB_PATH}
            
        # Test bivariate model
        python3.8 ${MIXER_PY} test2 \
            --trait1-file ${DATA_DIR}/${wm_trait}.csv.gz \
            --trait2-file ${DATA_DIR}/${PRIMARY_TRAIT}.csv.gz \
            --load-params-file ${OUT_DIR}/${PRIMARY_TRAIT}_vs_${wm_trait}.fit.rep${i}.json \
            --out ${OUT_DIR}/${PRIMARY_TRAIT}_vs_${wm_trait}.test2.rep${i} \
            --bim-file ${REF_DIR}/1000G.EUR.QC.@.bim \
            --ld-file ${REF_DIR}/1000G.EUR.QC.@.run4.ld \
            --lib ${LIB_PATH}
    done
    
    # Combine results
    python3.8 ${MIXER_FIG} combine \
        --json ${OUT_DIR}/${PRIMARY_TRAIT}_vs_${wm_trait}.fit.rep@.json \
        --out ${OUT_DIR}/${PRIMARY_TRAIT}_vs_${wm_trait}.fit
        
    python3.8 ${MIXER_FIG} combine \
        --json ${OUT_DIR}/${PRIMARY_TRAIT}_vs_${wm_trait}.test2.rep@.json \
        --out ${OUT_DIR}/${PRIMARY_TRAIT}_vs_${wm_trait}.test2
        
    python3.8 ${MIXER_FIG} two \
        --json-fit ${OUT_DIR}/${PRIMARY_TRAIT}_vs_${wm_trait}.fit.json \
        --json-test ${OUT_DIR}/${PRIMARY_TRAIT}_vs_${wm_trait}.test2.json \
        --out ${OUT_DIR}/${PRIMARY_TRAIT}_vs_${wm_trait}_bivarres \
        --statistic mean std
done