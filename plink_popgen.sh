#!/bin/bash
# === Computes He, Ho, Fst and LD by group
# === Specify run and path to data and results
# === Data are expected to have _SNP_filt suffix after filtration
RUN=test
bed=/path/to/data/
stat=/path/to/save/result/

# === Presetting PLINK runs as bash functions
plink_popstat () {      
    if [[ ! -e ${stat}pop_stat/${group} ]] ; then mkdir -p ${stat}pop_stat/${group} ; fi

    plink1.9 --bfile ${bed}${RUN}_SNP_filt \
             --chr-set 26 \
             --nonfounders \
             --family \
             --keep-cluster-names $group \
             --freqx \
             --hardy \
             --het \
             --out ${stat}pop_stat/${folder}/stat_results
}

plink_ld () {
	if [[ ! -e ${stat}LD/${group} ]] ; then mkdir -p ${stat}LD/${group} ; fi
	for chrom in `seq 1 26`
    do
		plink1.9 --bfile ${bed}${RUN}_SNP_filt \
				 --chr-set 26 \
				 --nonfounders \
				 --chr $chrom \
				 --family \
				 --keep-cluster-names $group \
				 --r2  square \
				 --out ${stat}LD/${folder}/${chrom}    #$folder for output
	done
	plink1.9 --bfile ${bed}${RUN}_SNP_filt \
			 --chr-set 26 \
			 --nonfounders \
			 --family \
			 --keep-cluster-names $group \
			 --r2 dprime \
			 --out ${stat}LD/${folder}/ld_table
	}

plink_fst () { 		## Pairwise Fst between two groups. Arguments $group1 $group2
    if [[ ! -e ${stat}Fst ]] ; then mkdir -p ${stat}Fst ; fi
    g1=`basename $group1`
    g2=`basename $group2`
    find ${stat}Fst/${g1}_${g2}* || find ${stat}Fst/${g2}_${g1}*
    if [[ $? -eq 1 && ! $group1 == $group2 ]]
    then
    #touch ${bed}Fst/${group1}_${group2}
    plink1.9 --bfile ${bed}${RUN}_SNP_filt \
             --chr-set 26 \
             --nonfounders \
             --family \
             --keep-cluster-names $group1 $group2 \
             --fst \
             --out ${stat}Fst/${g1}_${g2}
    fi
}


### LD and stat by groups
if [[ $# -eq 0 ]]
then
for breed in `cut -f1 -d " " ${bed}${RUN}_SNP_filt.fam | cut -f1 -d "/" | sort -u`
do
    cnt=$(cut -f1 -d " " ${bed}${RUN}_SNP_filt.fam | grep $breed | sort -u | wc -l)
    echo $cnt
    if [[ $cnt -eq 1 ]] 
    then 
        group=$(cut -f1 -d " " ${bed}${RUN}_SNP_filt.fam | grep $breed | sort -u)
        folder=$breed
        echo $group $folder
        plink_popstat
        plink_ld
    elif [[ $cnt -gt 1 ]] ; then
        group=$(cut -f1 -d " " ${bed}${RUN}_SNP_filt.fam | grep $breed | sort -u)
        folder=$breed
        echo $group $folder
        plink_popstat
        plink_ld
        for subgroup in `cut -f1 -d " " ${bed}${RUN}_SNP_filt.fam | grep $breed | sort -u`  
            do
            group=$subgroup
            folder=$group
            echo $group $folder
            plink_popstat
            plink_ld
            done
    fi 
done
fi
#global Fst
if [[ $1 == F || $# -eq 0 ]]
then
plink1.9 --bfile ${bed}${RUN}_SNP_filt \
             --chr-set 26 \
             --nonfounders \
             --family \
             --fst \
             --out ${stat}Fst/global
# pairwise Fst between populations
for group1 in `cut -f1 -d " " ${bed}${RUN}_SNP_filt.fam | sort -u`
do
    for group2 in `cut -f1 -d " " ${bed}${RUN}_SNP_filt.fam | sort -u`
    do
      plink_fst
    done
done
fi

# Global LD
if [[ $1 == L || $# -eq 0 ]]
then
for chrom in `seq 1 $NCHROM`
do
	plink1.9 --bfile ${bed}${RUN}_SNP_filt \
			 --chr-set 26 \
			 --nonfounders \
			 --chr $chrom \
			 --r2  square \
			 --out ${stat}LD/${chrom}    #$folder for output
done
plink1.9 --bfile ${bed}${RUN}_SNP_filt \
		 --chr-set 26 \
		 --nonfounders \
		 --r2 dprime \
		 --out ${stat}LD/ld_table
fi

