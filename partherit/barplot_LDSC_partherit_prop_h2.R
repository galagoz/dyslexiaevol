# ========================================================
# Making a barplot of the LDSC partitioned heritability
# "proportion of SNP heritability explained" for a single
# annotation. The human gained enhancers data is what
# I'm plotting here, for manuscript Figure 3B.

# updated for the nonGC version

# ========================================================

library(tidyverse)
library(scales)
#library(BSDA)
library(ggplot2)
library(kableExtra)
library(knitr)
library(cowplot)

options(stringsAsFactors=FALSE)

annots = list.files(path = "/data/clusterfs/lag/users/gokala/dyslexia-evol/ancreg_sumstat/results_selected_annots/", full.names = F,recursive = F)

all_data = read.table("/data/clusterfs/lag/users/gokala/dyslexia-evol/ancreg_sumstat/results_selected_annots/selected_annots_partherit_results_bonf9.txt", header=T, sep="\t")

annotations = c("Human Gained Enhancers - Adult human vs. chimpanzee (prefrontal cortex)",
                "Human Gained Promoters - Adult human vs. chimpanzee (prefrontal cortex)",
                "Human Gained Enhancers - Fetal Brain",
                "Neanderthal Lineage Depleted Regions",
                "Human Accelerated Regions",
                "Human Gained Enhancers - Adult human vs. macaque (prefrontal cortex)",
                "Human Gained Promoters - Adult human vs. macaque (prefrontal cortex)",
                "Neanderthal Introgressed SNPs",
                "Selective Sweep Regions")

# comment out the sig. asterisk label df, because there is nothing significant :(
#y_max <- max(all_data[, 3])
#y_axis_max <- y_max + all_data[all_data$Prop._h2 == y_max, 4] + 0.02
#label.df <- data.frame(Annotation = all_data$Annotation[all_data$significant=="Yes"],Prop._h2=all_data$Prop._h2[all_data$significant=="Yes"]+0.01)

# Put everything into a table

#old_annot_results=read.table("/data/clusterfs/lag/users/gokala/dyslexia-evol/partherit/results_tables/old_annots_partherit_results_FDR7.txt",header=T, sep = "\t")
#new_annot_results=read.table("/data/clusterfs/lag/users/gokala/dyslexia-evol/ancreg_sumstat/results_tables/new_annots_partherit_results_FDR7.txt",header=T, sep = "\t")
#save_kable(kable(new_annot_results, "latex", booktabs = T),"/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/plots/partherit/new_heritability_results.pdf")
#save_kable(kable(selected_annot_results, "latex", booktabs = T),"/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/plots/partherit/selectedAnnots_heritability_results.pdf")

selected_annot_results=read.table("/data/clusterfs/lag/users/gokala/dyslexia-evol/ancreg_sumstat/results_selected_annots/selected_annots_partherit_results_bonf9.txt",header=T, sep = "\t")
save_kable(kable(old_annot_results, "latex", booktabs = T),"/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/plots/partherit/old_heritability_results.pdf")

# Make barplots

p = ggplot(data = all_data, aes(Annotation, Prop._h2)) +
  geom_bar(stat = "identity", position = "dodge", fill = "burlywood4") +
  scale_y_continuous(limits = c(0, y_axis_max), oob = squish) +
  geom_linerange(position = position_dodge(width = 0.9), aes(ymin = Prop._h2 - Prop._h2_std_error, ymax = Prop._h2 + Prop._h2_std_error)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "bottom") +
  labs(
    x = "Annotation",
    y = expression(paste("Proportion ", italic("h"^{2}))),
    title = "Partitioned heritability for dyslexia, following ancestry regression",
    subtitle = "Proportion of heritability explained"
  )
p + geom_text(data = label.df, label = "*")
ggsave("/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/plots/partherit/h2_Prop_dyslexia_ancreg_Phase3_squished.pdf", width = 6.5, height = 3.25, units = "in")

p2 = ggplot(data = all_data, mapping = aes(Annotation, Prop._h2)) +
  geom_bar(stat = "identity", position = "dodge", fill = "burlywood4") +
  geom_linerange(position = position_dodge(width = 0.9), aes(ymin = Prop._h2 - Prop._h2_std_error, ymax = Prop._h2 + Prop._h2_std_error)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "bottom") +
  labs(
    x = "Annotation",
    y = expression(paste("Prop. ", italic("h"^{
      2
    }))),
    title = "Partitioned heritability for dyslexia, following ancestry regression",
    subtitle = "Proportion of heritability explained"
  )
p2# + geom_text(data = label.df, label = "*")
ggsave("/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/plots/partherit/h2_Prop_dyslexia_ancreg_Phase3_squished_v2.pdf", width = 6.5, height = 3.25, units = "in")

p3 = ggplot(data = all_data, mapping = aes(Annotation, Enrichment)) +
  geom_bar(stat = "identity", position = "dodge", fill = "burlywood4") +
  geom_linerange(position = position_dodge(width = 0.9), aes(ymin = Enrichment - Enrichment_std_error, ymax = Enrichment + Enrichment_std_error)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "bottom") +
  labs(
    x = "Annotation",
    y = "Enrichment Score",
    title = "Partitioned heritability for dyslexia, following ancestry regression",
    subtitle = "Heritability Enrichment"
  )
p3# + geom_text(data = label.df, label = "*")
ggsave("/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/plots/partherit/enrichment_dyslexia_ancreg_Phase3_squished_v2.pdf", width =6, height = 3.25, units = "in")

############### Final Plot with enrichment and proportion ###############


p4 = ggplot(data = all_data, aes(Annotation,Prop._h2)) +
  geom_bar(stat = "identity", position = "dodge", fill = "blue") +
  geom_linerange(position = position_dodge(width = 0.9), aes(ymin = Prop._h2 - Prop._h2_std_error, ymax = Prop._h2 + Prop._h2_std_error)) +
  theme_classic() +
  scale_x_discrete(labels= annotations) +
  scale_y_continuous(expand = c(0, 0),limits = c(-0.06,0.06)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "bottom") +
  labs(
    x = "Annotation",
    y = expression(paste("Proportion ", italic("h"^{2})))
  ) +
  #geom_text(data = label.df, label = "*") +
  coord_flip()

p5 = ggplot(data = all_data, aes(Annotation, Enrichment)) +
  geom_bar(stat="identity", position = "dodge", fill = "saddlebrown")+
  geom_linerange(position = position_dodge(width = 0.9), aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error))+
  #geom_text(aes(y = Enrichment+(Enrichment_std_error+2), label = bonferroni), position = position_dodge(width = 0.9), size = 3,  fontface=3, hjust = 0, vjust = 0, angle = 45)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  labs(y = "Enrichment") +
  #geom_text(data = label.df2, label = "*") +
  coord_flip()

final_plot=plot_grid(p4, p5, label_size = 12, align = "h", rel_widths = c(2, 1))

ggsave("/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/plots/partherit/selectedAnnots_23andMe_h2_propAndEnrichment.pdf", width = 13, height = 7, units = "in")