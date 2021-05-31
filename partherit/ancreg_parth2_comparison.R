library(ggplot2)
library(cowplot)

nonancreg=read.table("/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/results/partherit/nonancestry_regressed/newannots_partherit_results_bonf9.txt",header=T)
ancreg=read.table("/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/results/partherit/ancestry_regressed/selected_annots_partherit_results_bonf9.txt",header=T)
nonancreg$ancreg="non-ancestry regressed"
ancreg$ancreg="ancestry regressed"
bef_af=rbind(nonancreg,ancreg)

p1=ggplot(bef_af, aes(Annotation,Enrichment)) + 
  geom_point(aes(col=ancreg)) +
  coord_flip() +
  theme(legend.position = "none")
p2=ggplot(bef_af, aes(Annotation,Prop._h2)) + 
  geom_point(aes(col=ancreg)) +
  coord_flip() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none")
p3=ggplot(bef_af, aes(Annotation,Coefficient)) + 
  geom_point(aes(col=ancreg)) +
  coord_flip() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none")
p4=ggplot(bef_af, aes(Annotation,Enrichment_p)) + 
  geom_point(aes(col=ancreg)) +
  coord_flip() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

plot_row=plot_grid(p1,p2,p3,p4,nrow = 1,rel_widths = c(2,1,1,2))

title <- ggdraw() + 
  draw_label(
    "23andMe-Dyslexia partitioned heritability results",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

dyslexia_rescomp = plot_grid(title,plot_row,ncol=1, rel_heights = c(1,9))

######## Genlang

nonancregGen=read.table("/data/workspaces/lag/workspaces/lg-genlang/Working/Evolution/results/partitioned_heritability/mtag/newannots_partherit_results_bonf9_nonancreg.txt",header=T)
ancregGen=read.table("/data/workspaces/lag/workspaces/lg-genlang/Working/Evolution/results/partitioned_heritability/mtag/new_annots_partherit_results_bonf9.txt",header=T)
nonancregGen$ancreg="non-ancestry regressed"
ancregGen$ancreg="ancestry regressed"
bef_afGen=rbind(nonancregGen,ancregGen)

p1Gen=ggplot(bef_afGen, aes(Annotation,Enrichment)) + 
  geom_point(aes(col=ancreg)) +
  coord_flip() +
  theme(legend.position = "none")
p2Gen=ggplot(bef_afGen, aes(Annotation,Prop._h2)) + 
  geom_point(aes(col=ancreg)) +
  coord_flip() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none")
p3Gen=ggplot(bef_afGen, aes(Annotation,Coefficient)) + 
  geom_point(aes(col=ancreg)) +
  coord_flip() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none")
p4Gen=ggplot(bef_afGen, aes(Annotation,Enrichment_p)) + 
  geom_point(aes(col=ancreg)) +
  coord_flip() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

plot_rowGen=plot_grid(p1Gen,p2Gen,p3Gen,p4Gen,nrow = 1,rel_widths = c(2,1,1,2))

titleGen <- ggdraw() + 
  draw_label(
    "GenLang MTAG partitioned heritability results",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

genlang_rescomp = plot_grid(titleGen,plot_rowGen,ncol=1, rel_heights = c(1,9))

##########
#pdf("/data/workspaces/lag/workspaces/lg-genlang/Working/Evolution/results/partitioned_heritability/mtag/ancreg_parth2_comparison.pdf")
plot_grid(dyslexia_rescomp,genlang_rescomp,nrow=2)
#dev.off()