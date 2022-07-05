library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(grid)

## The version 4 is the assembly fasta file one version before the final.
## The version 5 is the final assembly.
###### Summarise Version 4 Coverage data #######

## Load the 100kb sliding window data of coverage depth of 5 female and 5 male w/ ver.4(one) and ver.5(final)
cov.data<-read.table("Dataset/Stacked_SW_Depth_PexsCon_ver4and5_100kb.txt",header=T)
index<-grep("ver4",cov.data$File)
cov.data<-cov.data[index,]
cov.data<-cov.data %>%
  mutate(Position=(Start+End-1)/2000000,ID=substr(File,29,30)) %>%
  mutate(Sex=ifelse(as.numeric(as.character(ID))<16,"Female","Male"))

mean.cov.data<-cov.data %>%
  group_by(ID) %>%
  summarise(Mean=mean(log2(Count/100000)))

cov.data.x<-cov.data %>%
  right_join(mean.cov.data, by=c("ID")) %>%
  mutate(Nor.cov=log2(Count/100000)-Mean) %>%
  filter(Chr=="ChrX")

female.mean.cov.data<-cov.data.x %>%
  filter(Sex=="Female") %>%
  group_by(Position) %>%
  summarise(Female_mean=mean(Nor.cov))

cov.data.x<-cov.data.x %>%
  right_join(female.mean.cov.data,by="Position") %>%
  mutate(Nor2.cov=Nor.cov-Female_mean)

new.cov.data.x %>%
  ggplot() +
  geom_line(mapping=aes(x=Position,y=fit,col=Sex,group=ID))

ytitle<-expression(paste("Normalized ",italic("log")[2],"-coverage depth"))
pc1<-cov.data.x %>%
  ggplot() +
  geom_line(mapping=aes(x=Position,y=Nor.cov,col=Sex,group=ID),alpha=0.5,size=0.3)+
  labs(x="Position in ChrX*",y=ytitle)+
  theme(panel.border=element_rect(size=0.5,fill=NA,colour="black"),
        plot.margin=margin(0,0,0,20),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=9),
        legend.title=element_blank(),
        legend.key.height=unit(8,"pt"),
        legend.position=c(0.825,0.88),
        legend.spacing.y=unit(0,"pt"),
        legend.text=element_text(size=8,margin=margin(0,0,0,0)))

###### Summarise Version 4 SNP data ######

snp.data<-read.table("Dataset/SW_SNP_count_X_Y_specific_PexsCon_ver4_100kb.txt",header=T)
snp.data.x <- snp.data %>%
  mutate(Position=(Start+End-1)/2000000) %>%
  pivot_longer(cols=c("Male_hetero","Xspecific","Yspecific"),names_to="SNP_type",values_to="Count") %>%
  filter(Chromosome=="ChrX")
snp.data.x$SNP_type<-factor(snp.data.x$SNP_type,levels=c("Male_hetero","Xspecific","Yspecific"))
levels(snp.data.x$SNP_type)<-c("Male heterozygous", "X-derived alternative", "Y-derived alternative")
pc2<-snp.data.x %>%
  ggplot()+
  geom_line(mapping=aes(x=Position,y=Count/100000,col=SNP_type),alpha=0.7)+
  scale_color_manual(values=c("grey","dodgerblue","firebrick"))+
  labs(x="Position in ChrX*",y="Proportion of SNPs")+
  coord_cartesian(ylim=c(0,0.001))+
  theme(panel.border=element_rect(size=0.5,fill=NA,colour="black"),
        legend.title=element_blank(),
        legend.position=c(0.72,0.84),
        legend.key.height=unit(8,"pt"),
        legend.spacing.y=unit(0,"pt"),
        legend.text=element_text(size=8,margin=margin(0,0,0,0)))

###### Drawing Figure for Version 4 ######

pdf("Figure_SW_XY_region_before_manual_curation.pdf",width=4,height=5)
pc1+pc2+plot_layout(ncol=1)
grid.text(label="A",x=unit(0.05,"npc"),y=unit(0.975,"npc"),gp=gpar(cex=2))
grid.text(label="B",x=unit(0.05,"npc"),y=unit(0.500,"npc"),gp=gpar(cex=2))
dev.off()

###### Summarise Version 5 Coverage data #######

## Load the 100kb sliding window data of coverage depth of 5 female and 5 male w/ ver.4(one) and ver.5(final)
cov.data<-read.table("Dataset/Stacked_SW_Depth_PexsCon_ver4and5_100kb.txt",header=T)
index<-grep("ver5",cov.data$File)
cov.data<-cov.data[index,]
cov.data<-cov.data %>%
  mutate(Position=(Start+End-1)/2000000,ID=substr(File,29,30)) %>%
  mutate(Sex=ifelse(as.numeric(as.character(ID))<16,"Female","Male"))

mean.cov.data<-cov.data %>%
  group_by(ID) %>%
  summarise(Mean=mean(log2(Count/100000)))

cov.data.x<-cov.data %>%
  right_join(mean.cov.data, by=c("ID")) %>%
  mutate(Nor.cov=log2(Count/100000)-Mean) %>%
  filter(Chr=="ChrX")

female.mean.cov.data<-cov.data.x %>%
  filter(Sex=="Female") %>%
  group_by(Position) %>%
  summarise(Female_mean=mean(Nor.cov))

cov.data.x<-cov.data.x %>%
  right_join(female.mean.cov.data,by="Position") %>%
  mutate(Nor2.cov=Nor.cov-Female_mean)

new.cov.data.x %>%
  ggplot() +
  geom_line(mapping=aes(x=Position,y=fit,col=Sex,group=ID))

ytitle<-expression(paste("Normalized ",italic("log")[2],"-coverage depth"))
pc1<-cov.data.x %>%
  ggplot() +
  geom_line(mapping=aes(x=Position,y=Nor.cov,col=Sex,group=ID),alpha=0.5,size=0.3)+
  labs(x="Position in ChrX*",y=ytitle)+
  theme(panel.border=element_rect(size=0.5,fill=NA,colour="black"),
        plot.margin=margin(0,0,0,20),
        axis.title.y=element_text(size=9),
        axis.title.x=element_blank(),
        legend.title=element_blank(),
        legend.key.height=unit(8,"pt"),
        legend.position=c(0.825,0.88),
        legend.spacing.y=unit(0,"pt"),
        legend.text=element_text(size=7,margin=margin(0,0,0,0)))

###### Summarise Version 5 SNP data ######

snp.data<-read.table("Dataset/SW_SNP_count_X_Y_specific_PexsCon_ver5_100kb.txt",header=T)
snp.data.x <- snp.data %>%
  mutate(Position=(Start+End-1)/2000000) %>%
  pivot_longer(cols=c("Male_hetero","Xspecific","Yspecific"),names_to="SNP_type",values_to="Count") %>%
  filter(Chromosome=="ChrX")
snp.data.x$SNP_type<-factor(snp.data.x$SNP_type,levels=c("Male_hetero","Xspecific","Yspecific"))
levels(snp.data.x$SNP_type)<-c("Male heterozygous", "X-derived alternative", "Y-derived alternative")
pc2<-snp.data.x %>%
  ggplot()+
  geom_line(mapping=aes(x=Position,y=Count/100000,col=SNP_type),alpha=0.7)+
  scale_color_manual(values=c("grey","dodgerblue","firebrick"))+
  labs(x="Position in ChrX*",y="Proportion of SNPs")+
  coord_cartesian(ylim=c(0,0.0014))+
  theme(panel.border=element_rect(size=0.5,fill=NA,colour="black"),
        legend.title=element_blank(),
        legend.position=c(0.72,0.84),
        legend.key.height=unit(8,"pt"),
        legend.spacing.y=unit(0,"pt"),
        legend.text=element_text(size=7,margin=margin(0,0,0,0)))

###### Drawing Figure for Version 5 ######

pdf("Figure_SW_XY_region_after_manual_curation.pdf",width=4,height=5)
pc1+pc2+plot_layout(ncol=1)
grid.text(label="A",x=unit(0.05,"npc"),y=unit(0.975,"npc"),gp=gpar(cex=2))
grid.text(label="B",x=unit(0.05,"npc"),y=unit(0.500,"npc"),gp=gpar(cex=2))
dev.off()
