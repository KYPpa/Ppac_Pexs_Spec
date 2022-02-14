library(ggplot2)
library(dplyr)
library(tidyr)
library(qtl)

###### Load genotype data and format ######

# Read file of sliding window data of genotypes (100kb window)
# Column
#1. File name (File)
#2. Chromosome name (Chr)
#3. Start position (Start)
#4. End position (End)
#5. Count of the number of RSB001-specific SNPs (Count)

sw.data<-read.table("Dataset/Genotype_data_Ppac_RILs.txt",sep="\t",header=T)
sw.data$Pos<-(sw.data$Start-1+sw.data$End)/2000000
sw.data$SeqID<-as.numeric(substr(sw.data$File,21,22))
chrlist<-c("ChrI","ChrII","ChrIII","ChrIV","ChrV","ChrX")
sw.data<-sw.data[as.character(sw.data$Chr) %in% chrlist,]

sw.data %>% ggplot(aes(x=Count))+geom_histogram()+coord_cartesian(xlim=c(0,50))
table(sw.data$Count)
sw.data$Genotype<-ifelse(sw.data$Count<=10,"A","B")
sw.data$Chr<-factor(sw.data$Chr,level=chrlist)

rename_chrlist<-c("c1","c2","c3","c4","c5","c6")
levels(sw.data$Chr)<-rename_chrlist

## Histogram of allele frequency

sw.data %>% group_by(Chr,Pos) %>%
  summarise(Prop=sum(Genotype=="B")/n()) %>%
  ggplot()+geom_histogram(aes(x=Prop))

## Removed markers with proportion < 0.25
valid<-sw.data %>% group_by(Chr,Pos) %>%
  summarise(Prop=sum(Genotype=="B")/n()) %>%
  filter(Prop>0.25)

sw.data2<-sw.data %>% semi_join(valid,by=c("Chr","Pos"))

## Loess fitting of allele

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

fitvec<-c()
for(i in 1:96){
  fitvec<-c(fitvec,loessfit(sw.data2$Chr[sw.data2$SeqID==i],
                            sw.data2$Pos[sw.data2$SeqID==i],
                            sw.data2$Count[sw.data2$SeqID==i],0.1))
}
sw.data2$fit<-fitvec
sw.data2$Genotype2<-ifelse(sw.data2$fit<=10,"A","B")

## Cross table for R/qtl

gt.data<-sw.data2 %>%
  mutate(MK=paste(as.character(Chr),as.character(Pos),sep="_")) %>%
  pivot_wider(id_cols=c("SeqID"),names_from=c("MK"),values_from=c("Genotype2"))
colnames(gt.data)[1]<-"id"
gt.data$id<-as.character(gt.data$id)
chrarray<-as.character(sw.data2$Chr[sw.data2$SeqID==1])
mbarray<-sw.data2$Pos[sw.data2$SeqID==1]-0.05
gt.data<-rbind(c("",chrarray),c("",mbarray),gt.data)
write.csv(gt.data,file="Dataset/Genotype_Ppac_RILs_RQTL.csv",row.names=F,quote=F)

###### Linkage analysis ######

mapthis<-read.cross("csv",file="Dataset/Genotype_Ppac_RILs_RQTL.csv",estimate.map=F,crosstype="riself")
gddata <- est.map(mapthis, error.prob=0.001, verbose=FALSE)
mapdata<-replace.map(mapthis,gddata)
mapdata<-calc.errorlod(mapdata)
pull.map(mapdata)
plot(mapdata)

save(mapdata,file="mapdata_Ppac_RILs.Robj")

###### Marey map ######

## Data preparation

res<-c()
for(i in 1:6){
  value<-as.vector(pull.map(mapdata,rename_chrlist[i])[[1]])
  marker<-names(pull.map(mapdata,rename_chrlist[i])[[1]])
  res.sub<-data.frame(matrix(unlist(strsplit(marker,"_")),byrow=T,ncol=2))
  res.sub<-cbind(res.sub,value)
  colnames(res.sub)<-c("Chr","Pos","cM")
  res.sub$Pos<-as.numeric(res.sub$Pos)
  res.sub$Chr<-rep(chrlist[i],nrow(res.sub))
  res<-rbind(res,res.sub)
}

res<-res %>%
  mutate(Chr=factor(Chr))

## Marey map

res %>%
  ggplot(aes(x=as.numeric(Pos),y=cM))+
  geom_line()+
  facet_wrap(.~Chr,scales="free",ncol=6)+
  labs(x="Physical position(Mb)",y="Genetic position(cM)")+
  theme(text=element_text(size=9))
