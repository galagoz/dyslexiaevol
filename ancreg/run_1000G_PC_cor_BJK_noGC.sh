#!/bin/bash
#$ -N run_PC_cor_BJK_noGC
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash

#-----run_1000G_PC_cor_BJK_noGC.sh-----

# Run correlation analysis between 1000G phase 3 PC loadings and GenLang sumstats BETAs for all phenotypes

#-----Variables-----
# $rdataDir - Directory containing GWAS summary statistics (in Rdata format)
# $outDir - Directory to write Spearman's correlation test results
# $rdataList - A txt file containing the list of "/path/to/dir/summary_statisctics.Rdata" of all phenotypes
rdataDir="/data/clusterfs/lag/users/gokala/enigma-evol/sumstatsRdata/global"
outDir="/data/clusterfs/lag/users/gokala/enigma-evol/corvals/"
rdataList="/data/clusterfs/lag/users/gokala/enigma-evol/sumstatsRdata/global/sumstats_rdata_list.txt"

#-----
# Read rdataList
cd /data/clusterfs/lag/users/gokala/enigma-evol
mkdir ${outDir}scripts
while read line; do 
   echo $line
   LINE=$line
   tmp_file_name=$(basename "$line")
   echo $tmp_file_name
   pheno_name="$(cut -d'_' -f4,6,8 <<<"$tmp_file_name")"
   echo $pheno_name
   tmp_run_file="${outDir}scripts/${pheno_name}.sh"
   echo '#!/bin/sh
#$ -N PC_cor_BJK_noGC
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash

#-----Variables-----
#source run_1000G_PC_cor_BJK_noGC.sh

Rscript /data/clusterfs/lag/users/gokala/enigma-evol/1000G_PC_cor_BJK_noGC.R' $LINE $pheno_name $outDir > $tmp_run_file
   chmod a+x $tmp_run_file
   echo "Created the script for cluster ->  submitting ${pheno_name} to the Grid"
   qsub -wd "/data/clusterfs/lag/users/gokala/enigma-evol/corvals/scripts" $tmp_run_file
done < $rdataList
