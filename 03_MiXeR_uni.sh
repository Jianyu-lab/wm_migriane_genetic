#!/bin/bash
# Univariate MiXeR analyses for estimating heritability, polygenicity and discoverability
# repeated 20 times on random subsets of approx 600K SNPs

# ==============================================================================
# Setup paths
# ==============================================================================
DATA_DIR="./data/traitfolder"
OUT_DIR="./results/mixer"
REF_DIR="./refdata/1000G_EUR_Phase3_plink"
MIXER_PY="software/mixer/precimed/mixer.py"
MIXER_FIG="software/mixer/precimed/mixer_figures.py"
LIB_PATH="software/mixer/src/build/lib/libbgmg.so"
TRAITS_LIST="./data/disease.list" # Contains migraine and 19 WM phenotypes

mkdir -p ${OUT_DIR}

for trait in $(cat ${TRAITS_LIST}); do
    for i in $(seq 1 20); do
        # Fit univariate model
        python3.8 ${MIXER_PY} fit1 \
            --trait1-file ${DATA_DIR}/${trait}.csv.gz \
            --out ${OUT_DIR}/${trait}.fit.rep${i} \
            --power-curve --qq-plots \
            --extract ${REF_DIR}/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${i}.snps \
            --bim-file ${REF_DIR}/1000G.EUR.QC.@.bim \
            --ld-file ${REF_DIR}/1000G.EUR.QC.@.run4.ld \
            --lib ${LIB_PATH}
            
        # Test univariate model
        python3.8 ${MIXER_PY} test1 \
            --trait1-file ${DATA_DIR}/${trait}.csv.gz \
            --load-params-file ${OUT_DIR}/${trait}.fit.rep${i}.json \
            --power-curve --qq-plots \
            --out ${OUT_DIR}/${trait}.test.rep${i} \
            --bim-file ${REF_DIR}/1000G.EUR.QC.@.bim \
            --ld-file ${REF_DIR}/1000G.EUR.QC.@.run4.ld \
            --lib ${LIB_PATH}
    done
    
    # Combine results
    python3.8 ${MIXER_FIG} combine \
        --json ${OUT_DIR}/${trait}.fit.rep@.json \
        --out ${OUT_DIR}/${trait}.fit.rep.all
        
    python3.8 ${MIXER_FIG} one \
        --json ${OUT_DIR}/${trait}.fit.rep.all.json \
        --out ${OUT_DIR}/${trait}.fit.rep.all \
        --statistic mean std
done