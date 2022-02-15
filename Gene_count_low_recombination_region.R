## Requres
# 1. R object, mapdata, for P. exspectatus,
# Dataset/mapdata_Pexspectaus_EMS_hybrids.Robj
# 2. 1kb sliding window data of count of expressed genes in P. exspectatus:
# Dataset/SW_Exp_Gene_PexsCon_ver5.txt
# 3. R object, mapdata, for P. pacificus,
# Dataset/mapdata_Ppac_RILs.Robj
# 4, 1kb sliding window data of count of expressed genes in P. pacificus:
# Dataset/SW_Exp_Gene_PpacEP_V3.txt

library(dplyr)
library(tidyr)
library(ggplot2)
library(qtl)
library(patchwork)

###### loess function ######

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

loessfitse<-function(chr,pos,value,span){
  chrlist<-unique(chr)
  res<-rep(NA,length(chr))
  for(c in chrlist){
    index<-which(chr==c)
    df<-data.frame(x=pos[index],y=value[index])
    fit<-loess(y~x,data=df,degree=1,span=span)
    loessfit<-predict(fit,se=T)$se
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

###### P. exspectatus data ######

load("Dataset/mapdata_Pexspectaus_EMS_hybrids.Robj")

## fitting recombination rate

chrlist<-c("ChrI","ChrII","ChrIII","ChrIV","ChrV","ChrX")
chrlist2<-c("I","II","III","IV","V","X")

res2<-c()
for(i in 1:6){
  value<-as.vector(pull.map(mapdata,chrlist2[i])[[1]])
  marker<-names(pull.map(mapdata,chrlist2[i])[[1]])
  res.sub<-data.frame(matrix(unlist(strsplit(marker,"-")),byrow=T,ncol=2))
  res.sub<-cbind(res.sub,value)
  colnames(res.sub)<-c("Chr","Pos","cM")
  res.sub$Chr<-rep(chrlist[i],nrow(res.sub))
  res.sub$Pos<-as.numeric(res.sub$Pos)
  res.sub2<-data.frame(Chr=rep(chrlist[i]),
                       Pos=res.sub$Pos[1:(nrow(res.sub)-1)]+(diff(res.sub$Pos)/2),
                       RR=diff(res.sub$cM)/diff(res.sub$Pos))
  res2<-rbind(res2,res.sub2)
}

res2$fit<-loessfit(res2$Chr,res2$Pos,res2$RR,span=0.2)
res2$se<-loessfitse(res2$Chr,res2$Pos,res2$RR,span=0.2)

RRdata.Pexs<-res2

## Determine the low recombination region

#Visually, see the position

RRdata.Pexs %>%
  ggplot()+geom_line(aes(x=Pos,y=fit))+
  geom_hline(yintercept=mean(RRdata.Ppac$fit),col="gold")+
  geom_line(aes(x=Pos,y=fit+se),linetype=3)+
  geom_line(aes(x=Pos,y=fit-se),linetype=3)+
  facet_wrap(.~Chr,scales="free",ncol=3)+
  labs(x="Physical position(Mb)",y="Recombination rate (cM/Mb)")+
  theme(text=element_text(size=9))

Pexsbond<-RRdata %>% filter(y<mean(RRdata$y)) %>%
  group_by(Chr) %>%
  summarise(min=min(x),max=max(x))

## ChrIII 
except<-RRdata %>% filter(y<mean(RRdata$y),x<20) %>%
  group_by(Chr) %>%
  summarise(min=min(x),max=max(x))
Pexsbond[3,]<-except[3,]

## ChrI and IV
except<-RRdata %>% filter(y<mean(RRdata$y),x>10) %>%
  group_by(Chr) %>%
  summarise(min=min(x),max=max(x))
Pexsbond[c(1,4),]<-except[c(1,4),]

## Gene count

## 1kb non-overlapping sliding window data for expressed gene count
## is loaded.

sw1k<-read.table("Dataset/SW_Exp_Gene_PexsCon_ver5.txt",header=T)
sw1k$Position<-(sw1k$Position+0.5)/1000
sum(sw1k$Count)
sw1k$Count<-sw1k$Count/sum(sw1k$Count)*100

cdata<-c()
for(chr in chrlist){
  min<-pull(Pexsbond[Pexsbond$Chr==chr,"min"])
  max<-pull(Pexsbond[Pexsbond$Chr==chr,"max"])
  genecount<-sw1k %>%
    dplyr::filter(Chromosome==chr,Position>min,Position<max) %>%
    summarise(sum=sum(Count)) %>% pull
  cdata<-rbind(cdata,c(chr,as.vector(genecount)))
}

# ChrX* (ChrIL and ChrIR)

min<-pull(Pexsbond[Pexsbond$Chr=="ChrX","min"])
max<-pull(Pexsbond[Pexsbond$Chr=="ChrX","max"])
irnum<-sw1k %>%
  dplyr::filter(Chromosome=="ChrX",Position>min,Position<20) %>%
  summarise(sum=sum(Count)) %>% pull
xnum<-sw1k %>%
  dplyr::filter(Chromosome=="ChrX",Position>20,Position<max) %>%
  summarise(sum=sum(Count)) %>% pull

## Data frame format for drawing

cdata<-cdata[-6,]
cdata<-data.frame(cdata)
cdata$Type<-rep("Other")
cdata$Type[cdata$X1=="ChrI"]<-"ChrIL"
cdata<-rbind(cdata,c("ChrX",irnum,"ChrIR"))
cdata<-rbind(cdata,c("ChrX",xnum,"ChrX"))
cdata$X2<-as.numeric(cdata$X2)
chrlist2<-c("ChrX","ChrI","ChrII","ChrIII","ChrIV","ChrV")
cdata$X1<-factor(cdata$X1,levels=chrlist2)
levels(cdata$X1)<-c("X*","I*","II*","III*","IV*","V*")
cdata$Type<-factor(cdata$Type,levels=c("ChrX","ChrIL","ChrIR","Other"))
cdata.pexs<-cdata

###### P. pacificus data ######

load(file="Dataset/mapdata_Ppac_RILs.Robj")

## fitting recombination rate

chrlist<-c("ChrI","ChrII","ChrIII","ChrIV","ChrV","ChrX")
rename_chrlist<-c("c1","c2","c3","c4","c5","c6")

res2<-c()
for(i in 1:6){
  value<-as.vector(pull.map(mapdata,rename_chrlist[i])[[1]])
  marker<-names(pull.map(mapdata,rename_chrlist[i])[[1]])
  res.sub<-data.frame(matrix(unlist(strsplit(marker,"_")),byrow=T,ncol=2))
  res.sub<-cbind(res.sub,value)
  colnames(res.sub)<-c("Chr","Pos","cM")
  res.sub$Pos<-as.numeric(res.sub$Pos)
  res.sub$Chr<-rep(chrlist[i],nrow(res.sub))
  res.sub2<-data.frame(Chr=rep(chrlist[i]),
                       Pos=res.sub$Pos[1:(nrow(res.sub)-1)]+(diff(res.sub$Pos)/2),
                       RR=diff(res.sub$cM)/diff(res.sub$Pos))
  res2<-rbind(res2,res.sub2)
}

res2$fit<-loessfit(res2$Chr,res2$Pos,res2$RR,span=0.2)
res2$se<-loessfitse(res2$Chr,res2$Pos,res2$RR,span=0.2)

## Determine the low recombination region

# Visually, see the position

RRdata.Ppac %>%
  ggplot()+geom_line(aes(x=Pos,y=fit))+
  geom_hline(yintercept=mean(RRdata.Ppac$fit),col="gold")+
  geom_line(aes(x=Pos,y=fit+se),linetype=3)+
  geom_line(aes(x=Pos,y=fit-se),linetype=3)+
  facet_wrap(.~Chr,scales="free",ncol=3)+
  labs(x="Physical position(Mb)",y="Recombination rate (cM/Mb)")+
  theme(text=element_text(size=9))

#ChrII, ChrV
Ppacbond<-RRdata.Ppac %>% dplyr::filter(fit<mean(RRdata.Ppac$fit),Pos>5,Pos<20) %>%
  group_by(Chr) %>%
  summarise(min=min(Pos),max=max(Pos))
#ChrI
exceptI<-RRdata.Ppac %>% filter(fit<mean(RRdata.Ppac$fit),Pos>5,Pos<35) %>%
  group_by(Chr) %>%
  summarise(min=min(Pos),max=max(Pos))
#ChrIII
exceptIII<-RRdata.Ppac %>% filter(fit<mean(RRdata.Ppac$fit),Pos>5,Pos<15) %>%
  group_by(Chr) %>%
  summarise(min=min(Pos),max=max(Pos))
#ChrIV
exceptIV<-RRdata.Ppac %>% filter(fit<mean(RRdata.Ppac$fit),Pos>5,Pos<27) %>%
  group_by(Chr) %>%
  summarise(min=min(Pos),max=max(Pos))
#ChrX
exceptX<-RRdata.Ppac %>% filter(fit<mean(RRdata.Ppac$fit),Pos>7.4,Pos<12) %>%
  group_by(Chr) %>%
  summarise(min=min(Pos),max=max(Pos))

Ppacbond[1,]<-exceptI[1,]
Ppacbond[3,]<-exceptIII[3,]
Ppacbond[4,]<-exceptIV[4,]
Ppacbond[6,]<-exceptX[6,]

## Gene count

## 1kb non-overlapping sliding window data for expressed gene count
## is loaded.

sw1k<-read.table("Dataset/SW_Exp_Gene_PpacEP_V3.txt",header=T)
sw1k$Position<-(sw1k$Position+0.5)/1000
sum(sw1k$Count)
sw1k$Count<-sw1k$Count/sum(sw1k$Count)*100

cdata<-c()
for(chr in chrlist){
  min<-pull(Ppacbond[Ppacbond$Chr==chr,"min"])
  max<-pull(Ppacbond[Ppacbond$Chr==chr,"max"])
  genecount<-sw1k %>%
    dplyr::filter(Chromosome==chr,Position>min,Position<max) %>%
    summarise(sum=sum(Count)) %>% pull
  cdata<-rbind(cdata,c(chr,as.vector(genecount)))
}

# ChrI including ChrIL and ChrIR

min<-pull(Ppacbond[Ppacbond$Chr=="ChrI","min"])
max<-pull(Ppacbond[Ppacbond$Chr=="ChrI","max"])
ilnum<-sw1k %>%
  dplyr::filter(Chromosome=="ChrI",Position>min,Position<21.05) %>%
  summarise(sum=sum(Count)) %>% pull
irnum<-sw1k %>%
  dplyr::filter(Chromosome=="ChrI",Position>21.05,Position<max) %>%
  summarise(sum=sum(Count)) %>% pull

# Format data for drawing

cdata<-cdata[-1,]
cdata<-data.frame(cdata)
cdata$Type<-rep("Other")
cdata$Type[cdata$X1=="ChrX"]<-"ChrX"
cdata<-rbind(cdata,c("ChrI",ilnum,"ChrIL"))
cdata<-rbind(cdata,c("ChrI",irnum,"ChrIR"))
cdata$X2<-as.numeric(cdata$X2)
chrlist2<-c("ChrI","ChrX","ChrII","ChrIII","ChrIV","ChrV")
cdata$X1<-factor(cdata$X1,levels=chrlist2)
levels(cdata$X1)<-c("I","X","II","III","IV","V")
cdata$Type<-factor(cdata$Type,levels=c("ChrX","ChrIL","ChrIR","Other"))
cdata.ppac<-cdata


###### Bar plot ######

bartheme<-theme(text=element_text(size=8),
                legend.key.size=unit(2,"pt"),
                legend.box.margin=margin(0,0,0,0,"pt"),
                legend.margin=margin(0,0,0,0,"pt"),
                axis.title.x=element_text(size=7),
                axis.title.y=element_text(size=6),
                plot.margin=unit(c(4,4,4,4),"pt"),
                panel.background=element_rect(fill=NA),
                panel.border=element_rect(linetype="solid",colour="grey40",fill=NA))

p1<-cdata.pexs %>% ggplot()+
  geom_bar(aes(x=X1,y=X2,fill=Type),stat="identity")+
  coord_cartesian(ylim=c(0,20))+
  scale_fill_manual(values=c("limegreen","darkorange","purple","grey40"))+
  labs(x="Low recombination region",y="% genes (expressed)",fill=NULL)+
  guides(fill=guide_legend(nrow=2,byrow=T))+
  bartheme+
  theme(legend.position="none")

p2<-cdata.ppac %>% ggplot()+
  geom_bar(aes(x=X1,y=X2,fill=Type),stat="identity")+
  coord_cartesian(ylim=c(0,20))+
  scale_fill_manual(values=c("limegreen","darkorange","purple","grey40"))+
  labs(x="Low recombination region",y="% genes (expressed)",fill=NULL)+
  guides(fill=guide_legend(nrow=2,byrow=T))+
  bartheme

pdf("Figure_Gene_count.pdf",width=5,height=2)
p1+p2+plot_layout(ncol=2)
dev.off()


