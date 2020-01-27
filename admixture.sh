#!/bin/bash
# === Run ADMIXTURE for series of K
# === Specify run and path to data and results
# === Data are expected to have _SNP_filt suffix after filtration
RUN=test
bed=/path/to/data/
ADM=/path/to/output/
cd $ADM
for k in `seq 1 20`
do
admixture32 --cv=15 ${bed}${RUN}_SNP_filt.bed $k -j9 | tee log${k}.txt
done
grep -h "CV" log*.txt > converge.txt
