## Run  BC1_hybrid_female_analysis_revised.R
##      BC1_hybrid_male_analysis_revised.R
##      BC1_hybrid_hermaphrodite_analysis_revised.R
library(ggplot2)
library(patchwork)

p1title<-expression(paste(bold("BC1 females [Pexs x [Pexs x Ppac]]")))
yaxis.title<- ~ atop(paste("Mean log"[2]," difference in normalized"),paste("coverage at introgression sites (BC1-F1)"))
vp1<-violin.data.female %>%
  ggplot(aes(x=Chr_new,y=Intcov))+
  geom_violin()+geom_jitter(size=0.2,alpha=0.3)+
  labs(title=p1title,x=NULL,y=yaxis.title)+
  theme(  text=element_text(size=8),
          panel.border=element_rect(linetype="solid",fill=NA),
          strip.background = element_rect(colour = "black"),
          plot.margin=margin(0,0,0,10),
          panel.spacing = unit(0, "lines"))

p2title<-expression(paste(bold("BC1 males [Ppac x [Pexs x Ppac]]")))

vp2<-violin.data.male %>%
  ggplot(aes(x=Chr_new,y=Intcov))+
  geom_violin()+geom_jitter(size=0.2,alpha=0.3)+
  labs(title=p2title,x=NULL,y=yaxis.title)+
  theme(  text=element_text(size=8),
          panel.border=element_rect(linetype="solid",fill=NA),
          strip.background = element_rect(colour = "black"),
          panel.spacing = unit(0, "lines"))

p3title<-expression(paste(bold("BC1 hermaphrodites [Ppac x [Ppac x Pexs]]")))

vp3<-violin.data.hermaphrodite %>%
  ggplot(aes(x=Chr_new,y=Intcov))+
  geom_violin()+geom_jitter(size=0.2,alpha=0.3)+
  labs(title=p3title,x="Chromosomes",y=yaxis.title)+
  theme(  text=element_text(size=8),
          panel.border=element_rect(linetype="solid",fill=NA),
          strip.background = element_rect(colour = "black"),
          panel.spacing = unit(0, "lines"))
pdf("Fig_coverage_difference_revised.pdf",width=4,height=7)
vp1+vp2+vp3+plot_layout(ncol=1,heights=c(1,1,1))
dev.off()