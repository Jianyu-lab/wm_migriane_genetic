# Code Availability for Shared genetic architecture between migraine and white matter phenotypes

This repository contains the scripts used for the analyses described in the paper. Shared genetic architecture between migraine and white matter phenotypes

## Contents

- `01_QC.R`: R script wrapping EasyQC for quality control of GWAS summary statistics. Includes mapping to 1000G reference panels and handling multiallelic/indel variants.
- `02_LDSC_bivar.sh`: Bivariate LDSC analysis for assessing sample overlap and genetic correlation using cross-trait intercepts.
- `03_MiXeR_uni.sh`: Univariate MiXeR analyses to estimate heritability, polygenicity, and discoverability. Repeated 20 times for robustness.
- `04_MiXeR_bivar.sh`: Bivariate MiXeR analyses to quantify polygenic overlap (Dice coefficient) between migraine and WM phenotypes.
- `05_cFDR.sh`: Conditional and conjunctional false discovery rate (cond/conjFDR) analyses to identify shared genomic loci, powered by the `pleiofdr` framework.
- `07_FOCUS.sh`: TWAS fine-mapping using FOCUS to prioritize likely causal genes based on GTEx v8 brain tissues.
- `08_SUPERGNOVA.sh`: Local genetic correlation analysis using SUPERGNOVA for validation of shared conjFDR loci.

## Dependencies

The scripts rely on the following software tools:

- **R Packages**: [data.table](https://github.com/Rdatatable/data.table), [EasyQC](https://www.genepi-regensburg.de/easyqc/), [stringr](https://github.com/tidyverse/stringr)
- **Python Packages**: [pandas](https://github.com/pandas-dev/pandas), [numpy](https://github.com/numpy/numpy), [scipy](https://github.com/scipy/scipy)
- **LDSC (v.1.0.1)**: <https://github.com/bulik/ldsc>
- **MiXeR (v.1.3)**: <https://github.com/precimed/mixer>
- **pleiofdr**: <https://github.com/precimed/pleiofdr>
- **FOCUS**: <https://github.com/bogdanlab/focus>
- **SUPERGNOVA**: <https://github.com/qlu-lab/SUPERGNOVA>
- **FUMA**: <https://fuma.ctglab.nl/> (web-based application, used post-analysis)

<br />

