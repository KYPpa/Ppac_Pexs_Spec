library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(grid)

###### Karyotype data #######
kdata<-read.table(file="Dataset/karyotype_data_for_Yoshida_NEE_2022.txt",header=T,sep="\t")
strainlevel<-c("PS312","RS5221","RSB001","RSA635","RS5522B","RS5811B1","RS5915","RS5920","RS5930")
kdata$StrainID<-factor(kdata$StrainID,levels=strainlevel)
kdata <- kdata %>% pivot_longer(cols=c("X5","X6","X7","X8","X9"),names_to="ChroNum",values_to="Pcell")
kdata$ChroNum<-factor(kdata$ChroNum,levels=c("X5","X6","X7","X8","X9"))
levels(kdata$ChroNum)<-as.character(5:9)

kg1<-kdata %>%
  filter(Cell=="Meiosis") %>%
  ggplot()+
  geom_bar(mapping=aes(x=StrainID,y=Pcell*100,fill=ChroNum),
           position="dodge",
           stat="identity")+
  labs(title="Prophase I cells",y="% Cells",fill="Chromosome\nnumber")+
  coord_cartesian(ylim=c(0,100))+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values=c("dodgerblue","firebrick","darkgoldenrod2","deeppink2","black"))+
  scale_x_discrete(guide = guide_axis(angle = 90))+
  theme(panel.border=element_rect(size=0.5,fill=NA,colour="black"),
        plot.margin=margin(0,0,0,10,"pt"),
        axis.text.x=element_text(size=8),
        axis.title.x=element_blank())

kg2<-kdata %>%
  filter(Cell=="Gamete") %>%
  ggplot()+
  geom_bar(mapping=aes(x=StrainID,y=Pcell*100,fill=ChroNum),
           position="dodge",
           stat="identity")+
  labs(title="Gamete",y="% Cells",fill="Chromosome\nnumber")+
  coord_cartesian(ylim=c(0,100))+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values=c("dodgerblue","firebrick","darkgoldenrod2","deeppink2","black"))+
  scale_x_discrete(guide = guide_axis(angle = 90))+
  theme(panel.border=element_rect(size=0.5,fill=NA,colour="black"),
        axis.text.x=element_text(size=8),
        axis.title.x=element_blank())

pdf(file="Figure_karyotype_data.pdf",width=7,height=4)
kg1+kg2+plot_layout(ncol=1)
grid.text(label="A",x=unit(0.02,"npc"),y=unit(0.96,"npc"),gp=gpar(cex=1.8))
grid.text(label="B",x=unit(0.02,"npc"),y=unit(0.48,"npc"),gp=gpar(cex=1.8))
dev.off()
