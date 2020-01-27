#!/bin/bash
# === Run FastStructure for series of K
# === Specify run and path to data and results
# === Data are expected to have _SNP_filt suffix after filtration
RUN=test
bed=/path/to/data/
struct=/path/to/output/
FS=/path/to/faststructure/exetutable

cd $struct
for k in `seq 1 20`
do
python $FS/structure.py \
     -K $k \
     --input=${bed}${RUN}_INTER.bed/${RUN}_INTER \
     --output=$struct/output_inter \
     --cv=10 \
     --format=bed \
     --full \
     --seed=43 
done
python $FS/chooseK.py --input=$ADM/struct/output_inter
