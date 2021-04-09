#!/bin/bash
#$ -N reformat_sumstats
#$ -cwd
#$ -q single.q
#$ -S /bin/bash

## Sumstat prep script, from Else - slightly changed to include beta values instead of OR

## Variables

DataDir="/data/clusterfs/lag/users/gokala/dyslexia-evol/sumstats"
SNPinfoDir="/data/clusterfs/lag/users/gokala/dyslexia-evol/snp_list"

## Reformatting

# Filter based on (all_snp_info.txt + im_snp_stat.txt = snp_IDs.txt) SNP list
awk 'BEGIN{OFS="\t"} NR==FNR{a[$1]=$0;next}{$1=a[$1];print $0,a[$1]}' $SNPinfoDir/snp_IDs.txt $DataDir/dyslexia_passQC.dat > $DataDir/dyslexia_passQC_annotated1.dat

echo -e "\tall.data.id\tim.data.id\trsid\tfreq.a\tfreq.b\tavg.rsqr\tmin.rsqr\tp.batch\tqc.mask\tchr\tposition\talleles\tchrpos\talleleA\talleleB\tsrc\tpvalue\teffect\tstderr\tpass\tim.num.0\tdose.b.0\tim.num.1\tdose.b.1\tAA.0\tAB.0\tBB.0\tAA.1\tAB.1\tBB.1" > $DataDir/dyslexia_passQC_annotated2.dat
awk 'BEGIN{OFS="\t"} NR==FNR{a[$2]=$0;next}{$2=a[$2];print a[$2],$0}' $SNPinfoDir/im_snp_stat.txt $DataDir/dyslexia_passQC_annotated1.dat | tail -n +2 >> $DataDir/dyslexia_passQC_annotated2.dat

# Keep SNPs with min.rsq >= 0.8, and usefull columns
awk 'BEGIN{OFS="\t"} $7>=0.8 {print $3,$4,$5,$6,$7,$10,$11,$13,$14,$15,$17,$18,$19,$21,$22,$23,$24}' $DataDir/dyslexia_passQC_annotated2.dat > $DataDir/dyslexia_passQC_annotated.dat
