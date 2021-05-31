# Updated file to work with 1000G Phase 3 data!
# Notes: This is designed to be run one annotation at a time
# If you run it like "> python Step2_make_bed_based_annots_Phase3.py HAR"
# it will substitute HAR at the appropriate spot and finish making the 
# annotation.

import os, sys

#----------------------------------------------------------------------
# Set up the list of chromosomes and other folders
#----------------------------------------------------------------------
mainDir="/data/clusterfs/lag/users/gokala/enigma-evol/partherit/amanda_annotations/"
oneKG="/data/clusterfs/lag/users/gokala/enigma-evol/partherit/annotations/1000G_EUR_Phase3_plink/"
hapmap="/data/workspaces/lag/shared_spaces/Resource_DB/LDscores/Phase3/1000G_EUR_Phase3_baseline/print_snps.txt"
LDSCdir="/home/gokala/ldsc/"

#----------------------------------------------------------------------
# Make the rest of the annotation files needed for running LDSC
# partitioned heritability. Command via this Wiki: 
# https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial
#----------------------------------------------------------------------

resultsDir = mainDir+sys.argv[1]+"/"
print(resultsDir)
os.chdir(resultsDir)

for i in range(1, 23):
	print("Working on annotation: "+sys.argv[1])
	print("Beginning LDScore calculation for chromosome "+str(i))
	os.system("python "+LDSCdir+"ldsc.py \
		--l2 --bfile "+oneKG+"1000G.EUR.QC."+str(i)+" \
		--ld-wind-cm 1 \
		--annot "+resultsDir+sys.argv[1]+"."+str(i)+".annot.gz \
		--out "+resultsDir+sys.argv[1]+"."+str(i)+" \
		--print-snps "+hapmap)
	print("Done with LDScore calculation for chromosome "+str(i))

print("Done with annotation: "+sys.argv[1]+"!")
print("***********************************")
