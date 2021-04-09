# ========================================================
# Make tables of LDSC partitioned heritability results

# The results files (.results) from the LDSC run should be organized
# by annotation, with separate directories for each annotation. 

# Updated for the nonGC version

# ========================================================

library(tidyverse)

options(stringsAsFactors=FALSE)

annots = list.dirs(path = "/data/clusterfs/lag/users/gokala/dyslexia-evol/ancreg_sumstat/results_selected_annots/", full.names = F, recursive = F)



partheritresults = data.frame(Category = character(0),
                              Prop._SNPs= numeric(0),
                              Prop._h2= numeric(0),
                              Prop._h2_std_error= numeric(0),
                              Enrichment= numeric(0),
                              Enrichment_std_error= numeric(0),
                              Enrichment_p= numeric(0),
                              Annotation=character(0),
                              Analysis=character(0)) #Will have a matrix with rows = number of E3MAs and columns = # of annotations


for (i in 1:length(annots)){
    print(annots[i])
    files = Sys.glob(path = paste0("/data/clusterfs/lag/users/gokala/dyslexia-evol/ancreg_sumstat/results_selected_annots/",annots[i],"/*.gz.results"))
    results = read.table(files,header=TRUE);
    info1 = str_split(files, pattern = "/")
    results$Annotation = info1[[1]][10]
    partheritresults = rbind(partheritresults,results[1,])
  }

### Half the p-values to convert them to one-sided test results
#results = read.table("/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/selected_annots_partherit_results_bonf9.txt", header=T, sep="\t")
partheritresults$Enrichment_p = partheritresults$Enrichment_p/2
partheritresults$Enrichment_p[partheritresults$Enrichment<0] = 1

partheritresults$bonf = p.adjust(partheritresults$Enrichment_p, method = "bonferroni") # correcting for 7 tests/annotations
#partheritresults$annot.p <- if_else(partheritresults$fdr < 0.05, as.character(round(partheritresults$fdr, digits = 4)), "")
partheritresults$significant = if_else(partheritresults$bonf < 0.05, "Yes", "")

write.table(partheritresults, 
            paste0("/data/clusterfs/lag/users/gokala/dyslexia-evol/ancreg_sumstat/results_selected_annots/selected_annots_partherit_results_bonf9_onesided.txt"),
            sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)