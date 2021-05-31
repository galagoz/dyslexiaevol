#!/bin/bash

# Make l2.ldscore files for all chromosomes/annotations
# Run this on gridmaster as "bash make_annot_ldscore.sh"

# Variables
inDir="/data/clusterfs/lag/users/gokala/enigma-evol/partherit/amanda_annotations"
annotList="/data/clusterfs/lag/users/gokala/enigma-evol/partherit/amanda_annotations/annot_list.txt"

# Run LD Score Estimation

mkdir $inDir/scripts

while read line; do
   echo $line
   annot=$(basename "$line")
   echo $annot
   tmp_run_file="${inDir}/scripts/${annot}.sh"
   echo '
#$ -N make_annot_ldscore.sh
#$ -cwd
#$ -q single.q
#$ -S /binnDepleted.sh/bash

mkdir '$inDir'/'$annot'
python /data/clusterfs/lag/users/gokala/enigma-evol/partherit/amanda_annotations/Step2_make_evo_based_annots_1KGPhase3.py ' $annot > $tmp_run_file

   chmod a+x $tmp_run_file
   echo "Created the script for cluster -> submitting ${annot} to the Grid"
   qsub -wd "/data/clusterfs/lag/users/gokala/enigma-evol/partherit/amanda_annotations/scripts" $tmp_run_file

done < $annotList
