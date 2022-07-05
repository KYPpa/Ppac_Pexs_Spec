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
library(grid)

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

round_in_point_one<-function(value,updown="up"){
  if(updown=="up"){
    ceiling(value*10)/10
  }else{
    floor(value*10)/10
  }
}


loessfit_100kbdf<-function(chr,pos,value,span){
  chrlist<-unique(chr)
  res<-c()
  for(c in chrlist){
    index<-which(chr==c)
    df<-data.frame(x=pos[index],y=value[index])
    fit<-loess(y~x,data=df,degree=1,span=span)
    newdf<-data.frame(Chr=rep(c),
                      x=seq(round_in_point_one(min(pos[index])),
                            round_in_point_one(max(pos[index]),"down"),
                            0.1))
    newdf$y<-predict(fit,newdata=newdf,by=0.1)
    newdf$se<-predict(fit,newdata=newdf,se=T)$se
    colnames(newdf)<-c("Chr","Pos","fit","se")
    res<-rbind(res,newdf)
  }
  return(res)
}


###### Gene counts w/ P. pacificus genome  ######

load("Dataset/mapdata_Pexspectaus_EMS_hybrids.Robj")

## fitting recombination rate

chrlist<-c("ChrI","ChrII","ChrIII","ChrIV","ChrV","ChrX")
chrlist2<-c("I","II","III","IV","V","X")

res<-c()
res2<-c()
for(i in 1:6){
  value<-as.vector(pull.map(mapdata,chrlist2[i])[[1]])
  marker<-names(pull.map(mapdata,chrlist2[i])[[1]])
  res.sub<-data.frame(matrix(unlist(strsplit(marker,"-")),byrow=T,ncol=2))
  res.sub<-cbind(res.sub,value)
  colnames(res.sub)<-c("Chr","Pos","cM")
  res.sub$Chr<-rep(chrlist[i],nrow(res.sub))
  res.sub$Pos<-as.numeric(res.sub$Pos)
  res<-rbind(res,res.sub)
  res.sub2<-data.frame(Chr=rep(chrlist[i]),
                       Pos=res.sub$Pos[1:(nrow(res.sub)-1)]+(diff(res.sub$Pos)/2),
                       RR=diff(res.sub$cM)/diff(res.sub$Pos))
  res2<-rbind(res2,res.sub2)
}

res_rr<-loessfit_100kbdf(res2$Chr,res2$Pos,res2$RR,span=0.2)
res_marey<-loessfit_100kbdf(res$Chr,res$Pos,res$cM,span=0.2)

res_marey%>%
  ggplot()+geom_line(aes(x=Pos,y=fit))+
  facet_wrap(.~Chr,scales="free",ncol=3)

RRdata.Pexs<-res_rr
## Determine the low recombination region

#Visually, see the position

RRdata.Pexs %>%
  ggplot()+geom_line(aes(x=Pos,y=fit))+
  geom_hline(yintercept=mean(RRdata.Pexs$fit),col="gold")+
  geom_line(aes(x=Pos,y=fit+se),linetype=3)+
  geom_line(aes(x=Pos,y=fit-se),linetype=3)+
  facet_wrap(.~Chr,scales="free",ncol=3)+
  labs(x="Physical position(Mb)",y="Recombination rate (cM/Mb)")+
  theme(text=element_text(size=9))

# Define Low recombination regions in X chromosome

lowrec<-RRdata.Pexs %>% filter(Chr=="ChrX",fit<mean(RRdata.Pexs$fit)) %>%
  summarise(start=min(Pos),end=max(Pos))

# Calculate total cM of the region
lowrec_cM<-res_marey$fit[res_marey$Chr=="ChrX"&res_marey$Pos==lowrec$end[1]]-
  res_marey$fit[res_marey$Chr=="ChrX"&res_marey$Pos==lowrec$start[1]]

## 1kb non-overlapping sliding window data for expressed gene count
## is loaded.

sw.eg<-read.table("Dataset/SW_Exp_Gene_PexsCon_ver5.txt",header=T)
sw.eg$Position<-(sw.eg$Position+0.5)/1000
sum(sw.eg$Count)
sw.eg$Count<-sw.eg$Count/sum(sw.eg$Count)*100

sw.pg<-read.table("Dataset/SW_gene_density_PexsCon_ver5.txt",header=T)
sw.pg$Position<-(sw.pg$Position+0.5)/1000
sum(sw.pg$Count)
sw.pg$Count<-sw.pg$Count/sum(sw.pg$Count)*100

# Function of couting genes

Count_N_Gene<-function(sw.data,chr,start,end){
  selected.sw.data<-
    sw.data %>% dplyr::filter(Chromosome==chr,Position>=start,Position<=end)
  return(sum(selected.sw.data$Count))
}

count_lowrec.eg<-Count_N_Gene(sw.eg,"ChrX",lowrec$start,lowrec$end)
count_lowrec.pg<-Count_N_Gene(sw.pg,"ChrX",lowrec$start,lowrec$end)

# Random pick of the given cM 

#"ChrX" is removed because 34cM doesn't exist out of low recombination regions
res_marey_outside<-res_marey %>% filter(Chr!="ChrX")

try<-10000
res.stat.eg<-c()
res.stat.pg<-c()
while(try>0){
  rind<-sample(1:nrow(res_marey_outside),1)
  chr<-res_marey_outside$Chr[rind]
  start<-res_marey_outside$Pos[rind]
  search_cM<-res_marey_outside$fit[rind]+lowrec_cM
  if(search_cM<=max(res_marey_outside$fit[res_marey_outside$Chr==chr])){
    end<-res_marey_outside %>%
      filter(Chr==chr,fit>search_cM) %>%
      summarise(end=min(Pos)) %>%
      pull(end)
    res.stat.eg<-c(res.stat.eg,Count_N_Gene(sw.eg,chr,start,end))
    res.stat.pg<-c(res.stat.pg,Count_N_Gene(sw.pg,chr,start,end))
    try<-try-1
  }
}

Pexs.res.stat.eg<-res.stat.eg
median(Pexs.res.stat.eg)
Pexs.res.stat.pg<-res.stat.pg
median(Pexs.res.stat.pg)
Pexs_count_lowrec.eg<-count_lowrec.eg #17.05
Pexs_count_lowrec.pg<-count_lowrec.pg #14.83189

## Histogram of the bootstrapping test

ggplot() + geom_histogram(aes(x=res.stat.eg),bins=100)+
  annotate(geom="segment",x=count_lowrec.eg,xend=count_lowrec.eg,y=-10,yend=0,arrow=arrow(length=unit(0.02,"npc")))
ggplot() + geom_histogram(aes(x=res.stat.pg),bins=100)+
  annotate(geom="segment",x=count_lowrec.pg,xend=count_lowrec.pg,y=-10,yend=0,arrow=arrow(length=unit(0.02,"npc")))

quantile(res.stat.pg,probs=0.9999)
quantile(res.stat.eg,probs=0.9999)

###### Gene counts w/ P. pacificus genome ######

load("Dataset/mapdata_Ppac_RILs.Robj")

chrlist<-c("ChrI","ChrII","ChrIII","ChrIV","ChrV","ChrX")
rename_chrlist<-c("c1","c2","c3","c4","c5","c6")

res<-c()
res2<-c()
for(i in 1:6){
  value<-as.vector(pull.map(mapdata,rename_chrlist[i])[[1]])
  marker<-names(pull.map(mapdata,rename_chrlist[i])[[1]])
  res.sub<-data.frame(matrix(unlist(strsplit(marker,"_")),byrow=T,ncol=2))
  res.sub<-cbind(res.sub,value)
  colnames(res.sub)<-c("Chr","Pos","cM")
  res.sub$Pos<-as.numeric(res.sub$Pos)
  res.sub$Chr<-rep(chrlist[i],nrow(res.sub))
  res<-rbind(res,res.sub)
  res.sub2<-data.frame(Chr=rep(chrlist[i]),
                       Pos=res.sub$Pos[1:(nrow(res.sub)-1)]+(diff(res.sub$Pos)/2),
                       RR=diff(res.sub$cM)/diff(res.sub$Pos))
  res2<-rbind(res2,res.sub2)
}

res_rr<-loessfit_100kbdf(res2$Chr,res2$Pos,res2$RR,span=0.2)
res_rr<-res_rr[-nrow(res_marey),]
res_marey<-loessfit_100kbdf(res$Chr,res$Pos,res$cM,span=0.2)

res_marey%>%
  ggplot()+geom_line(aes(x=Pos,y=fit))+
  facet_wrap(.~Chr,scales="free",ncol=3)

RRdata.Ppac<-res_rr

## Determine the low recombination region

#Visually, see the position

RRdata.Ppac %>%
  ggplot()+geom_line(aes(x=Pos,y=fit))+
  geom_hline(yintercept=mean(RRdata.Ppac$fit),col="gold")+
  geom_line(aes(x=Pos,y=fit+se),linetype=3)+
  geom_line(aes(x=Pos,y=fit-se),linetype=3)+
  facet_wrap(.~Chr,scales="free",ncol=3)+
  labs(x="Physical position(Mb)",y="Recombination rate (cM/Mb)")+
  theme(text=element_text(size=9))

# Define Low recombination regions in chromosomeI

lowrec<-RRdata.Ppac %>% filter(Chr=="ChrI",Pos>5,Pos<35,fit<mean(RRdata.Ppac$fit)) %>%
  summarise(start=min(Pos),end=max(Pos))

# Calculate total cM of the region
lowrec_cM<-res_marey$fit[res_marey$Chr=="ChrI"&res_marey$Pos==lowrec$end[1]]-
  res_marey$fit[res_marey$Chr=="ChrI"&res_marey$Pos==lowrec$start[1]]

## 1kb non-overlapping sliding window data for expressed gene count
## is loaded.

sw.eg<-read.table("Dataset/SW_Exp_Gene_PpacEP_V3.txt",header=T)
sw.eg$Position<-(sw.eg$Position+0.5)/1000
sum(sw.eg$Count)
sw.eg$Count<-sw.eg$Count/sum(sw.eg$Count)*100

sw.pg<-read.table("Dataset/SW_gene_density_El_Paco_V3.txt",header=T)
sw.pg$Position<-(sw.pg$Position+0.5)/1000
sum(sw.pg$Count)
sw.pg$Count<-sw.pg$Count/sum(sw.pg$Count)*100

# Function of couting genes

Count_N_Gene<-function(sw.data,chr,start,end){
  selected.sw.data<-
    sw.data %>% dplyr::filter(Chromosome==chr,Position>=start,Position<=end)
  return(sum(selected.sw.data$Count))
}

count_lowrec.eg<-Count_N_Gene(sw.eg,"ChrI",lowrec$start,lowrec$end)
count_lowrec.pg<-Count_N_Gene(sw.pg,"ChrI",lowrec$start,lowrec$end)

# Random pick of the given cM 

#"ChrI" is removed because 23cM doesn't exist out of low recombination regions
res_marey_outside<-res_marey %>% filter(Chr!="ChrI")

try<-10000
res.stat.eg<-c()
res.stat.pg<-c()
while(try>0){
  rind<-sample(1:nrow(res_marey_outside),1)
  chr<-res_marey_outside$Chr[rind]
  start<-res_marey_outside$Pos[rind]
  search_cM<-res_marey_outside$fit[rind]+lowrec_cM
  if(search_cM<=max(res_marey_outside$fit[res_marey_outside$Chr==chr])){
    end<-res_marey_outside %>%
      filter(Chr==chr,fit>search_cM) %>%
      summarise(end=min(Pos)) %>%
      pull(end)
    res.stat.eg<-c(res.stat.eg,Count_N_Gene(sw.eg,chr,start,end))
    res.stat.pg<-c(res.stat.pg,Count_N_Gene(sw.pg,chr,start,end))
    try<-try-1
  }
}

Ppac.res.stat.eg<-res.stat.eg
median(Ppac.res.stat.eg)
Ppac.res.stat.pg<-res.stat.pg
median(Ppac.res.stat.pg)
Ppac_count_lowrec.eg<-count_lowrec.eg #17.90
Ppac_count_lowrec.pg<-count_lowrec.pg #16.12

## Histogram of the bootstrapping test

ggplot() + geom_histogram(aes(x=res.stat.eg),bins=100)+
  annotate(geom="segment",x=count_lowrec.eg,xend=count_lowrec.eg,y=-10,yend=0,arrow=arrow(length=unit(0.02,"npc")))
ggplot() + geom_histogram(aes(x=res.stat.pg),bins=100)+
  annotate(geom="segment",x=count_lowrec.pg,xend=count_lowrec.pg,y=-10,yend=0,arrow=arrow(length=unit(0.02,"npc")))


quantile(res.stat.pg,probs=0.95)
quantile(res.stat.eg,probs=0.9999)

###### Drawing the figure of histograms of bootstrapping #####

p1<-ggplot() + geom_histogram(aes(x=Pexs.res.stat.pg),bins=100)+
  annotate(geom="segment",x=Pexs_count_lowrec.pg,xend=Pexs_count_lowrec.pg,y=250,yend=100,col="red",arrow=arrow(length=unit(0.10,"npc")))+
  labs(x=NULL,y="Count")+
  theme(plot.margin=margin(20,0,20,20))
p2<-ggplot() + geom_histogram(aes(x=Pexs.res.stat.eg),bins=100)+
  annotate(geom="segment",x=Pexs_count_lowrec.eg,xend=Pexs_count_lowrec.eg,y=300,yend=100,col="red",arrow=arrow(length=unit(0.10,"npc")))+
  labs(x=NULL,y=NULL)
p3<-ggplot() + geom_histogram(aes(x=Ppac.res.stat.pg),bins=100)+
  annotate(geom="segment",x=Ppac_count_lowrec.pg,xend=Ppac_count_lowrec.pg,y=300,yend=100,col="red",arrow=arrow(length=unit(0.10,"npc")))+
  labs(x="% Genes across genome",y="Count")
p4<-ggplot() + geom_histogram(aes(x=Ppac.res.stat.eg),bins=100)+
  annotate(geom="segment",x=Ppac_count_lowrec.eg,xend=Ppac_count_lowrec.eg,y=300,yend=100,col="red",arrow=arrow(length=unit(0.10,"npc")))+
  labs(x="% Genes across genome",y=NULL)
pdf("Fig_Gene_count_low_recombination_new.pdf",width=5,height=4)
p1+p2+p3+p4+plot_layout(ncol=2)
grid.text(expression(paste(bold("Predicted gene"))),x=unit(0.35,"npc"),y=unit(0.96,"npc"),gp=gpar(cex=0.8))
grid.text(expression(paste(bold("Expressed gene"))),x=unit(0.80,"npc"),y=unit(0.96,"npc"),gp=gpar(cex=0.8))
grid.text(expression(paste(bolditalic("P. exspectatus"))),x=unit(0.03,"npc"),y=unit(0.75,"npc"),rot=90,gp=gpar(cex=0.8))
grid.text(expression(paste(bolditalic("P. pacificus"))),x=unit(0.03,"npc"),y=unit(0.30,"npc"),rot=90,gp=gpar(cex=0.8))
grid.text("A",x=unit(0.13,"npc"),y=unit(0.94,"npc"),gp=gpar(cex=1.5))
grid.text("B",x=unit(0.58,"npc"),y=unit(0.94,"npc"),gp=gpar(cex=1.5))
grid.text("C",x=unit(0.13,"npc"),y=unit(0.49,"npc"),gp=gpar(cex=1.5))
grid.text("D",x=unit(0.58,"npc"),y=unit(0.49,"npc"),gp=gpar(cex=1.5))
dev.off()