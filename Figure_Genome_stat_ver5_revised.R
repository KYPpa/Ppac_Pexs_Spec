## First load Custom_scale_For_facet_grid.R

####### Library ######

library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(grid)

###### LOESS function ######

loessfit<-function(chr,pos,value,span){
  chrlist<-unique(chr)
  res<-rep(NA,length(chr))
  for(c in chrlist){
    index<-which(chr==c)
    df<-data.frame(x=pos[index],y=value[index])
    fit<-loess(y~x,data=df,degree=1,span=span)
    loessfit<-fit$fitted
    na.list<-which(is.na(df$value))
    if(length(na.list)>0){
      na.list<-na.list[order(na.list)]
      for(j in na.list){
        loessfit<-append(loessfit,NA,after=j-1)
      }
    }
    res[index]<-loessfit
  }
  return(res)
}

###### P exspectatus data ######

### Gene_density ###

sw.data<-read.table(file="Dataset/SW_gene_density_PexsCon_ver5_100w-100s.txt",header=T)
sw.data$Pos<-(sw.data$Start+sw.data$End)/2000000
sw.data$Count<-sw.data$Count/100000
selection<-c("ChrI","ChrII","ChrIII","ChrIV","ChrV","ChrX")
rename<-c("ChrI*","ChrII*","ChrIII*","ChrIV*","ChrV*","ChrX*")
sw.data<-sw.data[sw.data$Chromosome %in% selection,]
sw.data$Chromosome<-factor(sw.data$Chromosome,level=selection)
levels(sw.data$Chromosome)<-rename
sw.data$fit<-loessfit(sw.data$Chromosome,sw.data$Pos,sw.data$Count,0.4)
sw.data.gd<-sw.data

###### Expressed gene

sw.data<-read.table(file="Dataset/SW_Exp_Gene_PexsCon_ver5_100w-100s.txt",header=T)
sw.data$Pos<-(sw.data$Start+sw.data$End)/2000000
sw.data$Count<-sw.data$Count/100000
selection<-c("ChrI","ChrII","ChrIII","ChrIV","ChrV","ChrX")
rename<-c("ChrI*","ChrII*","ChrIII*","ChrIV*","ChrV*","ChrX*")
sw.data<-sw.data[sw.data$Chromosome %in% selection,]
sw.data$Chromosome<-factor(sw.data$Chromosome,level=selection)
levels(sw.data$Chromosome)<-rename
sw.data$fit<-loessfit(sw.data$Chromosome,sw.data$Pos,sw.data$Count,0.4)
sw.data.ex<-sw.data

### Repeat Masker

sw.data<-read.table(file="Dataset/SW_masked_RM_stg2_PexsCon_ver5_100w-100s.txt",header=T)
sw.data$Pos<-(sw.data$Start+sw.data$End)/2000000
sw.data$Count<-sw.data$Count/100000
selection<-c("ChrI","ChrII","ChrIII","ChrIV","ChrV","ChrX")
sw.data<-sw.data[sw.data$Chromosome %in% selection,]
sw.data$Chromosome<-factor(sw.data$Chromosome,levels=selection)
rename<-c("ChrI*","ChrII*","ChrIII*","ChrIV*","ChrV*","ChrX*")
levels(sw.data$Chromosome)<-rename
sw.data$fit<-loessfit(sw.data$Chromosome,sw.data$Pos,sw.data$Count,0.4)
sw.data.rs<-sw.data

###### GC content

sw.data<-read.table(file="Dataset/SW_GCcontent_PexsCon_ver5_100w-100s.txt",header=T)
sw.data$Pos<-(sw.data$Start+sw.data$End)/2000000
sw.data$GC_count<-sw.data$GC_count/(sw.data$GC_count+sw.data$AT_count)
colnames(sw.data)[5]<-"Count"
sw.data<-sw.data[,-4]
selection<-c("ChrI","ChrII","ChrIII","ChrIV","ChrV","ChrX")
sw.data<-sw.data[sw.data$Chromosome %in% selection,]
sw.data$Chromosome<-factor(sw.data$Chromosome,levels=selection)
rename<-c("ChrI*","ChrII*","ChrIII*","ChrIV*","ChrV*","ChrX*")
levels(sw.data$Chromosome)<-rename
sw.data$fit<-loessfit(sw.data$Chromosome,sw.data$Pos,sw.data$Count,span=0.4)
sw.data.gc<-sw.data

###### Drawing Figure for full Y-axis of P. exspectatus genome ######

commontheme<-  theme(  text=element_text(size=9),
                       axis.text.x=element_text(size=6),
                       plot.margin=unit(c(0,0,0,0),"in"),
                       axis.title.x=element_blank(),
                       axis.title.y=element_text(size=6),
                       panel.background=element_rect(fill=NA),
                       panel.border=element_rect(linetype="solid",colour="grey40",fill=NA),
                       panel.spacing = unit(0, "lines"))

facet_grid_PexsItoX<-facet_grid_custom(.~Chromosome, scales="free_x",space="free_x", scale_overrides = list(
  scale_override(1, scale_x_continuous(breaks=c(0,10,20),labels=c("",10,20))),
  scale_override(2, scale_x_continuous(breaks=c(0,10,20),labels=c("",10,20))),
  scale_override(3, scale_x_continuous(breaks=c(0,10,20),labels=c("",10,20))),
  scale_override(4, scale_x_continuous(breaks=c(0,10,20,30),labels=c("",10,20,30))),
  scale_override(5, scale_x_continuous(breaks=c(0,10,20),labels=c("",10,20))),
  scale_override(6, scale_x_continuous(breaks=c(0,10,20,30),labels=c("",10,20,30)))
))
pointsize<-0.05
pointalpha<-0.2

gp1<-sw.data.gd %>%
  ggplot(aes(x=Pos))+
  geom_point(aes(y=Count),size=pointsize,alpha=pointalpha)+
  geom_line(aes(y=fit),size=0.5,col="firebrick")+
  labs(y="Gene density")+
  scale_y_continuous(expand=c(0,0))+
  facet_grid_PexsItoX+
  commontheme

gp2<-sw.data.ex %>%
  ggplot(aes(x=Pos))+
  geom_point(aes(y=Count),size=pointsize,alpha=pointalpha)+
  geom_line(aes(y=fit),size=0.5,col="darkorange1")+
  labs(y="Expressed gene density")+
  scale_y_continuous(expand=c(0,0))+
  facet_grid_PexsItoX+
  commontheme+
  theme(strip.background=element_blank(),
        strip.text=element_blank())

gp3<-sw.data.rs %>%
  ggplot(aes(x=Pos))+
  geom_point(aes(y=Count),size=pointsize,alpha=pointalpha)+
  geom_line(aes(y=fit),size=0.5,col="limegreen")+
  labs(y="Repeat density")+
  scale_y_continuous(expand=c(0,0))+
  facet_grid_PexsItoX+
  commontheme+
  theme(strip.background=element_blank(),
        strip.text=element_blank())

gp4<-sw.data.gc %>%
  ggplot(aes(x=Pos))+
  geom_point(aes(y=Count),size=pointsize,alpha=pointalpha)+
  geom_line(aes(y=fit),size=0.5,col="dodgerblue")+
  labs(y="GC contents")+
  scale_y_continuous(expand=c(0,0))+
  facet_grid_PexsItoX+
  commontheme+
  theme(strip.background=element_blank(),
        strip.text=element_blank())

###### Combined data with recombination data ######

pdf(file="Fig_Pexs_genome_stat_revised.pdf",width=3.5,height=7)
{gp1+gp2+gp3+gp4+plot_layout(ncol=1)}+
{ep1+ep2+ep3+ep4+plot_layout(ncol=2,widths=c(3.5,2.33),heights=c(8.5,1))}+
{ep5+ep6+ep7+ep8+plot_layout(ncol=2,widths=c(4.2,1.68),heights=c(8.5,1))}+plot_layout(heights=c(0.6,0.6,0.6,0.6,1,1))
grid.text("A",x=unit(0.05,"npc"),y=unit(0.97,"npc"),gp=gpar(cex=1.5))
grid.text("B",x=unit(0.05,"npc"),y=unit(0.425,"npc"),gp=gpar(cex=1.5))
grid.text("C",x=unit(0.05,"npc"),y=unit(0.213,"npc"),gp=gpar(cex=1.5))
dev.off()



###### P. pacificus data ######

## Gene_density

sw.data<-read.table(file="Dataset/SW_gene_density_El_Paco_V3_100w-100s.txt",header=T)
sw.data$Pos<-(sw.data$Start+sw.data$End)/2000000
sw.data$Count<-sw.data$Count/100000
selection<-c("ChrI","ChrII","ChrIII","ChrIV","ChrV","ChrX")
sw.data<-sw.data[sw.data$Chromosome %in% selection,]
sw.data$Chromosome<-factor(sw.data$Chromosome,level=selection)
sw.data$fit<-loessfit(sw.data$Chromosome,sw.data$Pos,sw.data$Count,0.4)
sw.data.gd<-sw.data

## Expressed gene

sw.data<-read.table(file="Dataset/SW_Exp_Gene_PpacEP_V3_100w-100s.txt",header=T)
sw.data$Pos<-(sw.data$Start+sw.data$End)/2000000
sw.data$Count<-sw.data$Count/100000
selection<-c("ChrI","ChrII","ChrIII","ChrIV","ChrV","ChrX")
rename<-c("ChrI*","ChrII*","ChrIII*","ChrIV*","ChrV*","ChrX*")
sw.data<-sw.data[sw.data$Chromosome %in% selection,]
sw.data$Chromosome<-factor(sw.data$Chromosome,level=selection)
levels(sw.data$Chromosome)<-rename
sw.data$fit<-loessfit(sw.data$Chromosome,sw.data$Pos,sw.data$Count,0.4)
sw.data.ex<-sw.data

## Repeat Masker

sw.data<-read.table(file="Dataset/SW_masked_RM_stg2_El_Paco_100w-100s.txt",header=T)
sw.data$Pos<-(sw.data$Start+sw.data$End)/2000000
sw.data$Count<-sw.data$Count/100000
selection<-c("ChrI","ChrII","ChrIII","ChrIV","ChrV","ChrX")
sw.data<-sw.data[sw.data$Chromosome %in% selection,]
sw.data$Chromosome<-factor(sw.data$Chromosome,levels=selection)
sw.data$fit<-loessfit(sw.data$Chromosome,sw.data$Pos,sw.data$Count,0.4)
sw.data.rs<-sw.data

## GC content

sw.data<-read.table(file="Dataset/SW_GCcontent_El_paco_100w-100s.txt",header=T)
sw.data$Pos<-(sw.data$Start+sw.data$End)/2000000
sw.data$GC_count<-sw.data$GC_count/(sw.data$GC_count+sw.data$AT_count)
colnames(sw.data)[5]<-"Count"
sw.data<-sw.data[,-4]
selection<-c("ChrI","ChrII","ChrIII","ChrIV","ChrV","ChrX")
sw.data<-sw.data[sw.data$Chromosome %in% selection,]
sw.data$Chromosome<-factor(sw.data$Chromosome,levels=selection)
sw.data$fit<-loessfit(sw.data$Chromosome,sw.data$Pos,sw.data$Count,span=0.4)
sw.data.gc<-sw.data

###### Drawing Figure for full Y-axis of P. pacificus genome ######

commontheme<-  theme(  text=element_text(size=9),
                       axis.text.x=element_text(size=6),
                       plot.margin=unit(c(0,0,0,0),"in"),
                       axis.title.x=element_blank(),
                       axis.title.y=element_text(size=6),
                       panel.background=element_rect(fill=NA),
                       panel.border=element_rect(linetype="solid",colour="grey40",fill=NA),
                       panel.spacing = unit(0, "lines"))

facet_grid_PexsItoX<-facet_grid_custom(.~Chromosome, scales="free_x",space="free_x", scale_overrides = list(
  scale_override(1, scale_x_continuous(breaks=c(0,10,20),labels=c("",10,20))),
  scale_override(2, scale_x_continuous(breaks=c(0,10,20),labels=c("",10,20))),
  scale_override(3, scale_x_continuous(breaks=c(0,10,20),labels=c("",10,20))),
  scale_override(4, scale_x_continuous(breaks=c(0,10,20,30),labels=c("",10,20,30))),
  scale_override(5, scale_x_continuous(breaks=c(0,10,20),labels=c("",10,20))),
  scale_override(6, scale_x_continuous(breaks=c(0,10,20,30),labels=c("",10,20,30)))
))
pointsize<-0.05
pointalpha<-0.2

p1<-sw.data.gd %>%
  ggplot(aes(x=Pos))+
  geom_point(aes(y=Count),size=pointsize,alpha=pointalpha)+
  geom_line(aes(y=fit),size=0.5,col="firebrick")+
  labs(y="Gene density")+
  facet_grid_PexsItoX+
  commontheme

p2<-sw.data.ex %>%
  ggplot(aes(x=Pos))+
  geom_point(aes(y=Count),size=pointsize,alpha=pointalpha)+
  geom_line(aes(y=fit),size=0.5,col="darkorange1")+
  labs(y="Expressed gene density")+
  facet_grid_PexsItoX+
  commontheme+
  theme(strip.background=element_blank(),
        strip.text=element_blank())

p3<-sw.data.rs %>%
  ggplot(aes(x=Pos))+
  geom_point(aes(y=Count),size=pointsize,alpha=pointalpha)+
  geom_line(aes(y=fit),size=0.5,col="limegreen")+
  labs(y="Repeat density")+
  facet_grid_PexsItoX+
  commontheme+
  theme(strip.background=element_blank(),
        strip.text=element_blank())

p4<-sw.data.gc %>%
  ggplot(aes(x=Pos))+
  geom_point(aes(y=Count),size=pointsize,alpha=pointalpha)+
  geom_line(aes(y=fit),size=0.5,col="dodgerblue")+
  labs(y="GC contents")+
  facet_grid_PexsItoX+
  commontheme+
  theme(strip.background=element_blank(),
        strip.text=element_blank())

pdf(file="Fig_Ppac_genome_stat_revised.pdf",width=3.7,height=5)
p1+p2+p3+p4+plot_layout(ncol=1)
dev.off()


