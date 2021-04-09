#!/bin/bash
#$ -N reformat_sumstats
#$ -cwd
#$ -q single.q
#$ -S /bin/bash

#-----Reformat Sumstats-----

# Reformat your non-GC (genomic control) corrected GWAS summary stats at the beginning of your analysis
# according to the guideline at https://bitbucket.org/jasonlouisstein/enigmaevolma6/src/master/1000Gphase3_PC_cor/

# The new format should be as the following:
# SNP A1 A2 FREQ1 BETA SE P N MARKER CHR BP
# rs2977670 c g 0.9608 758.3807 485.5590 0.1183 10332 1:723891 1 723891
# rs143225517 t c 0.8454 718.8055 232.3162 0.001974 15846 1:751756 1 751756
# rs146277091 a g 0.0344 -1186.6514 501.2970 0.01793 10501 1:752478 1 752478

#-----Variables-----
sumstat="/data/clusterfs/lag/users/gokala/dyslexia-evol/sumstats/dyslexia_passQC_annotated.dat"
new_sumstat="${sumstat%.dat}_formatted.txt"
tmp="/data/clusterfs/lag/users/gokala/dyslexia-evol/sumstats/tmp.txt"

#-----
echo "$sumstat is being reformatted..."
awk '{print $1, $9, $10, $2, $12, $13, $11, "51801", $8, $6, $7}' $sumstat > $new_sumstat

sed -Ei 's/rsid/SNP/;s/alleleA/A1/;s/alleleB/A2/;s/freq.a/FREQ1/;s/effect/BETA/;s/stderr/SE/;s/pvalue/P/;s/chrpos/MARKER/;s/position/BP/' $new_sumstat

sed -Ei '1s/51801/N/1' $new_sumstat
sed -Ei '1s/chr/CHR/1' $new_sumstat
awk '{sub("chr", "", $10); print}' $new_sumstat > $tmp && mv $tmp $new_sumstat

chmod 777 $new_sumstat
echo "Done!"
#-----
