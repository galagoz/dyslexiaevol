#################################################################################################################
# By Else Eising, 28 May 2018
# Script for the calculation of genetic correlation between GenLang cohorts using LDSC
#################################################################################################################

# Define paths:
LDSCDir="/data/workspaces/lag/shared_spaces/Resource_DB/ldsc/"
LDScoreDir="/data/workspaces/lag/shared_spaces/Resource_DB/LDscores/eur_w_ld_chr/"
DataDir="/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia/GWAS_results_23andMe/"
SNPinfoDir="/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia/GWAS_results_23andMe/7.2_Annotation/v7.2_european_bundle/"
ResultsDir="/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia/LDSC/"

 

#################################################################################################################
## Step 1: Reformatting Summary Statistics
#################################################################################################################

# Keep all SNPs that pass 23andMe filter
gunzip $DataDir/RST-7764_4/dyslexia.dat.gz
head -n 1 $DataDir/RST-7764_4/dyslexia.dat > $DataDir/RST-7764_4/dyslexia_passQC.dat
awk '$6=="Y" {print $0}' $DataDir/RST-7764_4/dyslexia.dat >> $DataDir/RST-7764_4/dyslexia_passQC.dat

# Change log odss to odds ratios in the original file
echo "all.data.id     src     pvalue  OR  stderr  pass    im.num.0        dose.b.0        im.num.1        dose.b.1    AA.0     AB.0    BB.0    AA.1    AB.1    BB.1" > $DataDir/Filtered_sumstats/dyslexia_passQC_OR.dat
tail -n +2 $DataDir/RST-7764_4/dyslexia_passQC.dat | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, exp($4), $5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' >> $DataDir/Filtered_sumstats/dyslexia_passQC_OR.dat
# use nano $ResultsDir/Clapbeat_02222019_forLDSC.txt to fix the header for the odds ratio column

# merge with all_snp_info.txt and im_snp_stat.txt 
echo -e "all.data.id\trsid\tchr\tposition\talleles\tchrpos\talleleA\talleleB" > $SNPinfoDir/snp_IDs.txt
awk 'BEGIN{OFS="\t"} {print $1,$4,$5,$6,$7, substr($5,4) ":" $6, substr($7,1,1), substr($7, 3)}' $SNPinfoDir/all_snp_info.txt | tail -n +2 >> $SNPinfoDir/snp_IDs.txt
awk 'BEGIN{OFS="\t"} NR==FNR{a[$1]=$0;next}{$1=a[$1];print $0,a[$1]}' $SNPinfoDir/snp_IDs.txt $DataDir/Filtered_sumstats/dyslexia_passQC_OR.dat > $DataDir/Filtered_sumstats/dyslexia_passQC_OR_annotated1.dat


echo -e "\tall.data.id\tim.data.id\trsid\tfreq.a\tfreq.b\tavg.rsqr\tmin.rsqr\tp.batch\tqc.mask\tchr\tposition\talleles\tchrpos\talleleA\talleleB\tsrc\tpvalue\tOR\tstderr\tpass\tim.num.0\tdose.b.0\tim.num.1\tdose.b.1\tAA.0\tAB.0\tBB.0\tAA.1\tAB.1\tBB.1" > $DataDir/Filtered_sumstats/dyslexia_passQC_OR_annotated2.dat
awk 'BEGIN{OFS="\t"} NR==FNR{a[$2]=$0;next}{$2=a[$2];print a[$2],$0}' $SNPinfoDir/im_snp_stat.txt $DataDir/Filtered_sumstats/dyslexia_passQC_OR_annotated1.dat | tail -n +2 >> $DataDir/Filtered_sumstats/dyslexia_passQC_OR_annotated2.dat

# Keep SNPs with min.rsq >= 0.8, and usefull columns
awk 'BEGIN{OFS="\t"} $7>=0.8 {print $3,$4,$5,$6,$7,$10,$11,$13,$14,$15,$17,$18,$19,$21,$22,$23,$24}' $DataDir/Filtered_sumstats/dyslexia_passQC_OR_annotated2.dat > $DataDir/Filtered_sumstats/dyslexia_passQC_OR_annotated.dat
rm $DataDir/Filtered_sumstats/dyslexia_passQC_OR_annotated1.dat
rm $DataDir/Filtered_sumstats/dyslexia_passQC_OR_annotated2.dat


#################################################################################################################
## Step 2: Munge sumstats
#################################################################################################################


$LDSCDir/munge_sumstats.py \
--sumstats $DataDir/Filtered_sumstats/dyslexia_passQC_OR_annotated.dat \
--snp rsid \
--a1 alleleB \
--a2 alleleA \
--N-cas-col im.num.1 \
--N-con-col im.num.0 \
--p pvalue \
--signed-sumstats OR,1 \
--out $ResultsDir/dyslexia_passQC_OR_annotated_forLDSC \
--merge-alleles $LDScoreDir/w_hm3.snplist

#################################################################################################################
## Step 3: Calculate heritability as QC
#################################################################################################################
# results should match results of Niarchou and Gordon 2019 Biorxiv

$LDSCDir/ldsc.py \
--h2 $ResultsDir/dyslexia_passQC_OR_annotated_forLDSC.sumstats.gz \
--ref-ld-chr /data/workspaces/lag/shared_spaces/Resource_DB/LDscores/eur_w_ld_chr/ \
--w-ld-chr /data/workspaces/lag/shared_spaces/Resource_DB/LDscores/eur_w_ld_chr/ \
--samp-prev 0.045 \
--pop-prev 0.05 \
--out $ResultsDir/dyslexia_passQC_OR_annotated_forLDSC.Prev0.05.heritability 

$LDSCDir/ldsc.py \
--h2 $ResultsDir/dyslexia_passQC_OR_annotated_forLDSC.sumstats.gz \
--ref-ld-chr /data/workspaces/lag/shared_spaces/Resource_DB/LDscores/eur_w_ld_chr/ \
--w-ld-chr /data/workspaces/lag/shared_spaces/Resource_DB/LDscores/eur_w_ld_chr/ \
--samp-prev 0.045 \
--pop-prev 0.075 \
--out $ResultsDir/dyslexia_passQC_OR_annotated_forLDSC.Prev0.075.heritability 

$LDSCDir/ldsc.py \
--h2 $ResultsDir/dyslexia_passQC_OR_annotated_forLDSC.sumstats.gz \
--ref-ld-chr /data/workspaces/lag/shared_spaces/Resource_DB/LDscores/eur_w_ld_chr/ \
--w-ld-chr /data/workspaces/lag/shared_spaces/Resource_DB/LDscores/eur_w_ld_chr/ \
--samp-prev 0.045 \
--pop-prev 0.1 \
--out $ResultsDir/dyslexia_passQC_OR_annotated_forLDSC.Prev0.1.heritability 



#################################################################################################################
## Step 4: Check for warnings
#################################################################################################################

cd $ResultsDir
grep 'WARNING' *log
grep 'ERROR' *log
grep 'chi^2' *log 
# if the mean chi2 is not above 1.02, there is not much polygenic signal to work with. if mean chi2 is below 1, ldsc will not work

