options(stringsAsFactors=FALSE)
library(GenomicRanges)

#dircorvals = "/data/clusterfs/lag/users/gokala/beta/ancreg/"

##Output file
foutput = "/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/plots/ancreg/corvalues_before_after_BJK.pdf"

fcorvals = "/data/clusterfs/lag/users/gokala/dyslexia-evol/corvals/corvalues_Dyslexia_BJK.csv"
fancreg_corvals = "/data/clusterfs/lag/users/gokala/dyslexia-evol/corvals/corvalues_Dyslexia_ancreg_BJK.csv"

##Read in correlation values
corvals = read.csv(fcorvals)
ancreg_corvals = read.csv(fancreg_corvals)

##Make a pdf file of correlation estimates
pdf(foutput, width = 8, height = 8)
par(mfrow=c(2,1),new = T)

  ##Make barplots of correlation estimates for each
  corind = grep("BJK_cor",colnames(corvals));
  pind = grep("BJK_P",colnames(corvals));
  seind = grep("BJK_SE",colnames(corvals));
  x = barplot(as.matrix(corvals[1,corind]),main="Dyslexia correlation of effect sizes \nprior to ancestry regression with ancestry PCs",ylab="correlation coefficient",xlab="Ancestry PC",names.arg=paste0("PC",seq(1,20)),ylim=c(-0.05,0.05), yaxp = c(-0.05,0.05,2),las=2, xaxt="n");
  y0 = as.numeric(corvals[1,corind]-corvals[1,seind]);
  y1 = as.numeric(corvals[1,corind]+corvals[1,seind]);
  arrows(x,y0,x,y1,angle=90,length=0)
  ##Add asterisk when significant (Bonferroni corrected)
  bonfsigind = which(corvals[1,pind] < 0.05/20);
  ## draw asterisk for significant ones
  ##text(x[sigind]+((x[sigind+1]-x[sigind])/2),0.095,"*");
  text(x[bonfsigind],0.05,"*");
  ##Add o when nominally significant
  nomsigind = which(corvals[1,pind] >= 0.05/20 & corvals[1,pind] < 0.05);
  text(x[nomsigind],0.04,"o");
  
  ##Make barplots of correlation estimates for each
  corind = grep("BJK_cor",colnames(ancreg_corvals));
  pind = grep("BJK_P",colnames(ancreg_corvals));
  seind = grep("BJK_SE",colnames(ancreg_corvals));
  y = barplot(as.matrix(ancreg_corvals[1,corind]),main="Dyslexia correlation of effect sizes \nafter ancestry regression with ancestry PCs",ylab="correlation coefficient",xlab="Ancestry PC",names.arg=paste0("PC",seq(1,20)),ylim=c(-0.05,0.05), yaxp = c(-0.05,0.05,2),las=2, xaxt="n");
  y0 = as.numeric(ancreg_corvals[1,corind]-ancreg_corvals[1,seind]);
  y1 = as.numeric(ancreg_corvals[1,corind]+ancreg_corvals[1,seind]);
  arrows(x,y0,x,y1,angle=90,length=0)
  ##Add asterisk when significant (Bonferroni corrected)
  bonfsigind = which(ancreg_corvals[1,pind] < 0.05/20);
  ## draw asterisk for significant ones
  ##text(x[sigind]+((x[sigind+1]-x[sigind])/2),0.095,"*");
  text(x[bonfsigind],0.05,"*");
  ##Add o when nominally significant
  nomsigind = which(ancreg_corvals[1,pind] >= 0.05/20 & ancreg_corvals[1,pind] < 0.05);
  text(x[nomsigind],0.04,"o");
  
dev.off()