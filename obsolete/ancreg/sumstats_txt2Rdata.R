## This script converts a txt file to Rdata forma

inputDir="/data/clusterfs/lag/users/gokala/dyslexia-evol/sumstats/"
outputDir="/data/clusterfs/lag/users/gokala/dyslexia-evol/sumstatsRdata"

mergedGR=read.table("/data/clusterfs/lag/users/gokala/dyslexia-evol/sumstats/dyslexia_passQC_annotated_formatted.txt", header=T, fill=T, stringsAsFactors=F)
save(mergedGR,file="/data/clusterfs/lag/users/gokala/dyslexia-evol/sumstatsRdata/dyslexia_passQC_annotated_formatted.Rdata")