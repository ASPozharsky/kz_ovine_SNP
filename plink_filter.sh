#!/bin/bash
# === Filter data by HWE and MAF
# === Specify run and path to data
$RUN=test
$bed=\path\to\data

plink1.9 --file ${bed}${RUN} \
      --chr-set 26 \
      --chr 1-26 \
      --nonfounders \
      --mind 0.9 \
      --hwe 0.0001 \
      --maf 0.05 \
      --make-bed --out ${bed}${RUN}_SNP_filt

# === Replace NA with 0 for downstream compatibility  
sed -i "s/NA/0/" ${ped}${RUN}_SNP_filt.fam
sed -i "s/NA/0/" ${ped}${RUN}_SNP_filt.fam

