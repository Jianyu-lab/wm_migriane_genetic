#!/bin/bash
# cond/conjFDR analyses using cfdr.pleio R package (via MATLAB implementation commonly used)
# Significance thresholds: condFDR < 0.01 for trait-specific, conjFDR < 0.05 for joint
# MHC region (6:26-34Mb) and 8p23.1 region (8:7.2-12.5Mb) excluded

# ==============================================================================
# Setup paths
# ==============================================================================
DATA_DIR="./data/traitfolder"
OUT_DIR="./results/cFDR"
SOFTWARE_DIR="./software/pleiofdr"
TRAITS_LIST="./data/cFDR_trait_list.txt" # list of pairs: WM_trait migraine

mkdir -p ${OUT_DIR}
line_count=$(awk 'END{print NR}' "${TRAITS_LIST}")

for i in $(seq 1 $line_count); do
    IDP=$(awk "NR==$i {print \$1}" "${TRAITS_LIST}")
    disease=$(awk "NR==$i {print \$2}" "${TRAITS_LIST}")
    
    # Configure parameters for MATLAB script
    new_traitfolder="${DATA_DIR}/"
    new_traitfile1="${IDP}.mat"
    new_traitname1="${IDP}"
    new_traitfiles="${disease}"
    new_traitnames="${disease}"
    new_outputdir="${OUT_DIR}/${IDP}_${disease}"
    new_reffile="${SOFTWARE_DIR}/1kgPhase3eur_LDr2p1.mat"
    new_randprune_n=500
    new_stattype="conjfdr"
    new_fdrthresh="0.05"
    new_exclude="[6 25000000 34000000; 8 7200000 12500000]"
    
    configfile="${OUT_DIR}/config_${IDP}_${disease}_cjfdr.txt"
    runmefile="${OUT_DIR}/runme_${IDP}_${disease}_cjfdr.m"
    cp -f ${SOFTWARE_DIR}/config_default.txt "$configfile"
    cp -f ${SOFTWARE_DIR}/runme.m "$runmefile"
    
    # Replace parameters in config
    sed -i "s|reffile=.*|reffile=$new_reffile|g" "$configfile"
    sed -i "s|randprune_n=.*|randprune_n=$new_randprune_n|g" "$configfile"
    sed -i "s|traitfolder=.*|traitfolder=$new_traitfolder|g" "$configfile"
    sed -i "s|traitfile1=.*|traitfile1=$new_traitfile1|g" "$configfile"
    sed -i "s|traitname1=.*|traitname1=$new_traitname1|g" "$configfile"
    sed -i "s|traitfiles=.*|traitfiles={'${new_traitfiles}'}|g" "$configfile"
    sed -i "s|traitnames=.*|traitnames={'${new_traitnames}'}|g" "$configfile"
    sed -i "s|stattype=.*|stattype=$new_stattype|g" "$configfile"
    sed -i "s|fdrthresh=.*|fdrthresh=$new_fdrthresh|g" "$configfile"
    sed -i "s|exclude_chr_pos=.*|exclude_chr_pos=$new_exclude|g" "$configfile"
    sed -i "s|outputdir=.*|outputdir=$new_outputdir|g" "$configfile"
    
    sed -i "4s|config=.*|config='$configfile'|g" "$runmefile"
    mv "$runmefile" ${SOFTWARE_DIR}/ -f
    mv "$configfile" ${SOFTWARE_DIR}/ -f
    
    # Run MATLAB
    cmd="args=\"addpath(genpath('${SOFTWARE_DIR}/'));runme_${IDP}_${disease}_cjfdr\"\nmatlab -nodesktop -nosplash -r \"\$args\""
    echo -e $cmd > "${OUT_DIR}/run_${IDP}_${disease}_cfdr.sh"
    bash "${OUT_DIR}/run_${IDP}_${disease}_cfdr.sh"
    
    # Clump results to identify independent genomic loci
    python2 ${SOFTWARE_DIR}/python_convert/sumstats.py clump \
        --clump-field FDR \
        --force \
        --plink plink \
        --sumstats ${OUT_DIR}/${IDP}_${disease}/result.mat.csv \
        --bfile-chr ${SOFTWARE_DIR}/chr@ \
        --exclude-ranges '6:25000000-34000000,8:7200000-12500000' \
        --clump-p1 0.05 \
        --out ${OUT_DIR}/${IDP}_${disease}/result_leadsnp.csv
done