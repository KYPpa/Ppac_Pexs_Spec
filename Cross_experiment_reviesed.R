library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)
library(grid)
library(coin)

maintheme<-theme(panel.border=element_rect(size=0.5,fill=NA,colour="black"),
      plot.margin=unit(c(0,0,0,0),"in"),
      axis.text.y=element_text(size=6),
      axis.title.y=element_text(size=6,vjust=-0.3),
      axis.ticks.x=element_blank(),
      axis.title.x=element_blank(),
      axis.text.x=element_blank())

tabletheme<-theme_classic()+
  theme(plot.margin=unit(c(0,0,0,0),"in"),
        plot.background = element_rect(fill="transparent",colour=NA),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=5))

###### Statistics ######

## Data preparation

cdata1<-read.table(file="Dataset/PEEP_F1_data.txt",header=T,sep="\t")
cdata1$CrossName<-factor(cdata1$CrossName,levels=c("PxP","PxE","ExP","ExE"))
cdata1<-cdata1 %>% filter(!is.na(CrossName))

cdata2<-read.table(file="Dataset/CSD_F2_data.txt",header=T,sep="\t")
cdata2$CrossName<-factor(cdata2$CrossName,levels=c("EPS","EPxE","EPxP","ExEP","ExE0"))
cdata2 <-cdata2 %>%  filter(!is.na(CrossName))

cdata3<-read.table(file="Dataset/PEEP_F2_data.txt",header=T,sep="\t")
cdata3$CrossName<-factor(cdata3$CrossName,levels=c("PES","PExE","PExP","ExPE","ExE2"))
cdata3<-cdata3 %>% filter(!is.na(CrossName))

## Hybrid cross

testdata<-cdata1 %>% filter(CrossName=="PxP"|CrossName=="PxE")
wilcox_test(Sum~CrossName,data=testdata)
testdata<-cdata1 %>% filter(CrossName=="ExE"|CrossName=="PxE")
wilcox_test(Sum~CrossName,data=testdata)
testdata<-cdata1 %>% filter(CrossName=="PxP"|CrossName=="ExP")
wilcox_test(Sum~CrossName,data=testdata)
testdata<-cdata1 %>% filter(CrossName=="ExE"|CrossName=="ExP")

cdata1 %>%
  filter(!is.na(AbnormalRatio)) %>%
  group_by(CrossName) %>%
  summarise(mean_se(AbnormalRatio)) %>%
  mutate(se=ymax-y)
cdata1 %>%
  filter(!is.na(SexRatio)) %>%
  group_by(CrossName) %>%
  summarise(mean_se(SexRatio)) %>%
  mutate(se=ymax-y)

## Backcross

testdata<-cdata2 %>% filter(CrossName=="ExEP"|CrossName=="ExE0")
wilcox_test(Sum~CrossName,data=testdata)

testdata1<-cdata2[cdata2$CrossName=="EPxP",]
testdata1$CrossName<-factor(testdata1$CrossName,level=c("EPxP","PExP"))
testdata2<-cdata3[cdata3$CrossName=="PExP",]
testdata2$CrossName<-factor(testdata2$CrossName,level=c("EPxP","PExP"))
testdata<-rbind(testdata1,testdata2)
wilcox_test(Sum~CrossName,data=testdata)

testdata1<-cdata2[cdata2$CrossName=="EPxE",]
testdata1$CrossName<-factor(testdata1$CrossName,level=c("EPxE","PExE"))
testdata2<-cdata3[cdata3$CrossName=="PExE",]
testdata2$CrossName<-factor(testdata2$CrossName,level=c("EPxE","PExE"))
testdata<-rbind(testdata1,testdata2)
wilcox_test(Sum~CrossName,data=testdata)

testdata<-cdata3 %>% filter(CrossName=="ExPE"|CrossName=="PExE")
wilcox_test(Sum~CrossName,data=testdata)
wilcox_test(AbnormalRatio~CrossName,data=testdata)

###### Panel preparation ######

cdata<-read.table(file="Dataset/PEEP_F1_data.txt",header=T,sep="\t")
cdata$CrossName<-factor(cdata$CrossName,levels=c("PxP","PxE","ExP","ExE"))
cdata<-cdata %>% filter(!is.na(CrossName))

cg1<-cdata %>%
  group_by(CrossName) %>%
  summarise(Success=sum(ifelse(Sum>0,1,0))/n()) %>%
  ggplot()+
  geom_bar(aes(x=CrossName,y=Success*100),stat="identity")+
  geom_hline(yintercept=0)+
  coord_cartesian(ylim=c(0,100))+
  labs(y="% Cross w/ progeny")+
  maintheme+
  theme(plot.margin=margin(10,0,0,0))

cg2<-cdata %>%
  ggplot(mapping=aes(x=CrossName,y=Sum))+
  geom_boxplot(outlier.shape = NA,col="cornflowerblue")+
  geom_jitter(width=0.2,height=0,size=0.3)+
  labs(y="No. progeny")+
  maintheme

cg3<-cdata %>%
  ggplot(mapping=aes(x=CrossName,y=AbnormalRatio*100))+
  geom_boxplot(outlier.shape = NA,col="cornflowerblue")+
  geom_jitter(width=0.2,height=0,size=0.3)+
  coord_cartesian(ylim=c(0,100))+
  labs(y="% Immature")+
  maintheme

cg4<-cdata %>%
  group_by(CrossName) %>%
  summarise(N=as.character(n())) %>%
  mutate(Dam=factor(c("P","P","E","E")),Sire=factor(c("P","E","P","E"))) %>%
  pivot_longer(cols=c("Dam","Sire")) %>%
  mutate(name=factor(name,levels=c("Sire","Dam"))) %>%
  ggplot()+geom_text(aes(x=CrossName,y=name,label=value),size=2.5)+
  labs(x=NULL,y=NULL)+tabletheme

cdata<-read.table(file="Dataset/CSD_F2_data.txt",header=T,sep="\t")
cdata$CrossName<-factor(cdata$CrossName,levels=c("EPS","EPxP","EPxE","ExEP","ExE0"))
cdata <-cdata %>%  filter(!is.na(CrossName))

cg5<-cdata %>%
  group_by(CrossName) %>%
  summarise(Success=sum(ifelse(Sum>0,1,0))/n()) %>%
  ggplot()+
  geom_bar(aes(x=CrossName,y=Success*100),stat="identity")+
  geom_hline(yintercept=0)+
  coord_cartesian(ylim=c(0,100))+
  labs(y="% Cross w/ progeny")+
  maintheme+
  theme(plot.margin=margin(10,0,0,0))

cg6<-cdata %>%
  ggplot(mapping=aes(x=CrossName,y=Sum))+
  geom_boxplot(outlier.shape = NA,col="cornflowerblue")+
  geom_jitter(width=0.2,height=0,size=0.3)+
  labs(y="No. progeny")+
  maintheme

cg7<-cdata %>%
  mutate(ModAbnormalRatio=ifelse(CrossName=="EPxP",NA,AbnormalRatio)) %>%
  ggplot()+
  geom_boxplot(mapping=aes(x=CrossName,y=ModAbnormalRatio*100),outlier.shape = NA,col="cornflowerblue")+
  geom_jitter(mapping=aes(x=CrossName,y=AbnormalRatio*100),width=0.2,height=0,size=0.3)+
  coord_cartesian(ylim=c(0,100))+
  labs(y="% Immature")+
  maintheme

cg8<-cdata %>%
  group_by(CrossName) %>%
  summarise(N=as.character(n())) %>%
  mutate(Dam=factor(c("F1","F1","F1","E","E")),Sire=factor(c("-","P","E","F1","E"))) %>%
  pivot_longer(cols=c("Dam","Sire")) %>%
  mutate(name=factor(name,levels=c("Sire","Dam"))) %>%
  ggplot()+geom_text(aes(x=CrossName,y=name,label=value),size=2.5)+
  labs(x=NULL,y=NULL)+
  theme_classic()+tabletheme


cdata<-read.table(file="Dataset/PEEP_F2_data.txt",header=T,sep="\t")
cdata$CrossName<-factor(cdata$CrossName,levels=c("PES","PExP","PExE","ExPE","ExE2"))
cdata <-cdata %>%  filter(!is.na(CrossName))

cg9<-cdata %>%
  group_by(CrossName) %>%
  summarise(Success=sum(ifelse(Sum>0,1,0))/n()) %>%
  ggplot()+
  geom_bar(aes(x=CrossName,y=Success*100),stat="identity")+
  geom_hline(yintercept=0)+
  coord_cartesian(ylim=c(0,100))+
  labs(y="% Cross w/ progeny")+
  maintheme+
  theme(plot.margin=margin(10,0,0,10))

cg10<-cdata %>%
  ggplot(mapping=aes(x=CrossName,y=Sum))+
  geom_boxplot(outlier.shape = NA,col="cornflowerblue")+
  geom_jitter(width=0.2,height=0,size=0.3)+
  labs(y="No. progeny")+
  maintheme

cg11<-cdata %>%
  ggplot(mapping=aes(x=CrossName,y=AbnormalRatio*100))+
  geom_boxplot(outlier.shape = NA,col="cornflowerblue")+
  geom_jitter(width=0.2,height=0,size=0.3)+
  coord_cartesian(ylim=c(0,100))+
  labs(y="% Immature")+
  maintheme

cg12<-cdata %>%
  group_by(CrossName) %>%
  summarise(N=as.character(n())) %>%
  mutate(Dam=factor(c("F1","F1","F1","E","E")),Sire=factor(c("-","P","E","F1","E"))) %>%
  pivot_longer(cols=c("Dam","Sire")) %>%
  mutate(name=factor(name,levels=c("Sire","Dam"))) %>%
  ggplot()+geom_text(aes(x=CrossName,y=name,label=value),size=2.5)+
  labs(x=NULL,y=NULL)+
  theme_classic()+tabletheme

###### Figure Drawing ######

indexsize=0.9

pdf(file="Figure_Cross_experiments.pdf",width=4,height=4)
cg1+cg2+cg3+cg4+cg4+cg4+
  cg5+cg6+cg7+cg8+cg8+cg8+
  cg9+cg10+cg11+cg12+cg12+cg12+
  plot_layout(ncol=3,heights=c(3,1,3,1,3,1))
grid.text("C",x=unit(0.065,"npc"),y=unit(0.97,"npc"),gp=gpar(cex=indexsize))
grid.text("D",x=unit(0.382,"npc"),y=unit(0.97,"npc"),gp=gpar(cex=indexsize))
grid.text("E",x=unit(0.695,"npc"),y=unit(0.97,"npc"),gp=gpar(cex=indexsize))
grid.text("F",x=unit(0.065,"npc"),y=unit(0.65,"npc"),gp=gpar(cex=indexsize))
grid.text("G",x=unit(0.382,"npc"),y=unit(0.65,"npc"),gp=gpar(cex=indexsize))
grid.text("H",x=unit(0.695,"npc"),y=unit(0.65,"npc"),gp=gpar(cex=indexsize))
grid.text("I",x=unit(0.065,"npc"),y=unit(0.328,"npc"),gp=gpar(cex=indexsize))
grid.text("J",x=unit(0.382,"npc"),y=unit(0.328,"npc"),gp=gpar(cex=indexsize))
grid.text("K",x=unit(0.695,"npc"),y=unit(0.328,"npc"),gp=gpar(cex=indexsize))
grid.text("F1[E×P]",x=unit(0.03,"npc"),y=unit(0.53,"npc"),rot=90,gp=gpar(cex=0.7,fontface="bold"))
grid.text("F1[P×E]",x=unit(0.03,"npc"),y=unit(0.21,"npc"),rot=90,gp=gpar(cex=0.7,fontface="bold"))
dev.off()

###### Sex ratio ######
sgtheme<- theme(panel.border=element_rect(size=1,fill=NA,colour="black"),
                 plot.margin=unit(c(5,0,0,20),"pt"),
                 axis.text.y=element_text(size=8),
                 axis.title.y=element_text(size=10),
                 axis.ticks.x=element_blank(),
                 axis.title.x=element_blank(),
                 axis.text.x=element_blank())
sgtabletheme<-theme_classic()+
  theme(plot.margin=unit(c(0,0,0,0),"pt"),
        plot.background = element_rect(fill="transparent",colour=NA),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=8))


cdata<-read.table(file="Dataset/PEEP_F1_data.txt",header=T,sep="\t")
cdata$CrossName<-factor(cdata$CrossName,levels=c("PxP","PxE","ExP","ExE"))
cdata<-cdata %>% filter(!is.na(CrossName))


sg1<-cdata %>%
  ggplot(mapping=aes(x=CrossName,y=SexRatio*100))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width=0.2,height=0)+
  coord_cartesian(ylim=c(0,100))+
  labs(y="% Male progeny")+
  sgtheme


sg2<-cdata %>%
  group_by(CrossName) %>%
  summarise(N=as.character(n())) %>%
  mutate(Dam=factor(c("P","P","E","E")),Sire=factor(c("P","E","P","E"))) %>%
  pivot_longer(cols=c("Dam","Sire")) %>%
  mutate(name=factor(name,levels=c("Sire","Dam"))) %>%
  ggplot()+geom_text(aes(x=CrossName,y=name,label=value),size=3.5)+
  labs(x=NULL,y=NULL)+
  sgtabletheme

cdata<-read.table(file="Dataset/CSD_F2_data.txt",header=T,sep="\t")
cdata$CrossName<-factor(cdata$CrossName,levels=c("EPS","EPxE","EPxP","ExEP","ExE0"))
cdata <-cdata %>%  filter(!is.na(CrossName))

sg3<-cdata %>%
  mutate(ModSexRatio=ifelse(CrossName=="EPxP",NA,SexRatio)) %>%
  ggplot()+
  geom_boxplot(mapping=aes(x=CrossName,y=ModSexRatio*100),outlier.shape = NA)+
  geom_jitter(mapping=aes(x=CrossName,y=SexRatio*100),width=0.2,height=0)+
  coord_cartesian(ylim=c(0,100))+
  labs(y="% Male progeny")+
  sgtheme

sg4<-cdata %>%
  group_by(CrossName) %>%
  summarise(N=as.character(n())) %>%
  mutate(Dam=factor(c("F1","F1","F1","E","E")),Sire=factor(c("-","E","P","F1","E"))) %>%
  pivot_longer(cols=c("Dam","Sire")) %>%
  mutate(name=factor(name,levels=c("Sire","Dam"))) %>%
  ggplot()+geom_text(aes(x=CrossName,y=name,label=value),size=3.5)+
  labs(x=NULL,y=NULL)+
  sgtabletheme

cdata<-read.table(file="Dataset/PEEP_F2_data.txt",header=T,sep="\t")
cdata$CrossName<-factor(cdata$CrossName,levels=c("PES","PExE","PExP","ExPE","ExE2"))
cdata <-cdata %>%  filter(!is.na(CrossName))

sg5<-cdata %>%
  ggplot(mapping=aes(x=CrossName,y=SexRatio*100))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width=0.2,height=0)+
  coord_cartesian(ylim=c(0,100))+
  labs(y="% Male progeny")+
  sgtheme

sg6<-cdata %>%
  group_by(CrossName) %>%
  summarise(N=as.character(n())) %>%
  mutate(Dam=factor(c("F1","F1","F1","E","E")),Sire=factor(c("-","E","P","F1","E"))) %>%
  pivot_longer(cols=c("Dam","Sire")) %>%
  mutate(name=factor(name,levels=c("Sire","Dam"))) %>%
  ggplot()+geom_text(aes(x=CrossName,y=name,label=value),size=3.5)+
  labs(x=NULL,y=NULL)+
  sgtabletheme

indexsize<-2
pdf(file="Figure_Sex_ratio_F1_and_F2_crosses.pdf",width=3.5,height=7)
sg1+sg2+sg3+sg4+sg5+sg6+plot_layout(ncol=1,heights=c(3,1,3,1,3,1))
grid.text("A",x=unit(0.05,"npc"),y=unit(0.96,"npc"),gp=gpar(cex=indexsize))
grid.text("B",x=unit(0.05,"npc"),y=unit(0.64,"npc"),gp=gpar(cex=indexsize))
grid.text("C",x=unit(0.05,"npc"),y=unit(0.31,"npc"),gp=gpar(cex=indexsize))
grid.text("F1[E  X P ]",x=unit(0.05,"npc"),y=unit(0.53,"npc"),rot=90,gp=gpar(cex=1.2,fontface="bold"))
grid.text("F1[P  X E ]",x=unit(0.05,"npc"),y=unit(0.21,"npc"),rot=90,gp=gpar(cex=1.2,fontface="bold"))
dev.off()

###### Inter cross ######

cdata<-read.table(file="Dataset/PEEP_F2_data.txt",header=T,sep="\t")
cdata$CrossName<-factor(cdata$CrossName,levels=c("EPIC","PEIC"))
cdata<-cdata %>% filter(!is.na(CrossName))
levels(cdata$CrossName)<-c("F1[E   x P  ]\nintercross","F1[P    x E  ]\nintercross")

igtheme<-theme(panel.border=element_rect(size=0.5,fill=NA,colour="black"),
                 plot.margin=unit(c(0,0,0,5),"pt"),
                 axis.text.y=element_text(size=10),
                 axis.title.y=element_text(size=10),
                 axis.ticks.x=element_blank(),
                 axis.title.x=element_blank(),
                axis.text.x = element_text(angle = 45, hjust = 1))

ig1<-cdata %>%
  group_by(CrossName) %>%
  summarise(Success=sum(ifelse(Sum>0,1,0))/n()) %>%
  ggplot()+
  geom_bar(aes(x=CrossName,y=Success*100),stat="identity")+
  geom_hline(yintercept=0)+
  coord_cartesian(ylim=c(0,100))+
  labs(y="% Cross w/ progeny")+
  igtheme+
  theme(plot.margin=margin(15,0,0,0))

ig2<-cdata %>%
  ggplot(mapping=aes(x=CrossName,y=Sum))+
  geom_boxplot(outlier.shape = NA,col="cornflowerblue")+
  geom_jitter(width=0.2,size=1,height=0)+
  labs(y="No. progeny")+
  igtheme

ig3<-cdata %>%
  ggplot(mapping=aes(x=CrossName,y=AbnormalRatio*100))+
  geom_boxplot(outlier.shape = NA,col="cornflowerblue")+
  geom_jitter(width=0.2,size=1,height=0)+
  coord_cartesian(ylim=c(0,100))+
  labs(y="% Immature")+
  igtheme

pdf("Figure_intercross.pdf",width=6,height=2.2)
ig1+ig2+ig3+plot_layout(ncol=3)
grid.text("A",x=unit(0.022,"npc"),y=unit(0.95,"npc"),gp=gpar(cex=1.6))
grid.text("B",x=unit(0.36,"npc"),y=unit(0.95,"npc"),gp=gpar(cex=1.6))
grid.text("C",x=unit(0.676,"npc"),y=unit(0.95,"npc"),gp=gpar(cex=1.6))
dev.off()

###
