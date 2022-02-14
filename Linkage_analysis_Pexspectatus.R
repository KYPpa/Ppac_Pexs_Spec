library(dplyr)
library(tidyr)
library(ggplot2)
library(qtl)
library(patchwork)

###### Data loading and formating ######

# Read file of sliding window data of genotypes
# Column
#1. File name (File)
#2. Chromosome name (Chr)
#3. Start position (Start)
#4. End position (End)
#5. SNP count of the allele of the reference (Ref)
#6. SNP count of the allele of the grand-dam (D)
#7. SNP count of the allele of the grand-sire (S)

sw.data<-read.table("Dataset/Genotype_data_Pexs_EMS_hybrids.txt",sep="\t",header=T)
sw.data$ID<-as.numeric(substr(sw.data$File,25,27))
sw.data$Cat<-ifelse(sw.data$D>0,ifelse(sw.data$S>0,"-","A"),ifelse(sw.data$S>0,"H","-"))
sw.data$Position<-((sw.data$Start-1)+250000)/1000000

## Format the data as a csv file imported in R/qtl

cdata<-sw.data[,c("ID","Chr","Position","Cat")]
cdata$Position<-as.character(as.numeric(cdata$Position))
gdata<-cdata %>%
  pivot_wider(id_cols=c("Chr","Position"),names_from=ID,values_from=Cat) %>%
  data.frame

chrlist<-c("ChrI","ChrII","ChrIII","ChrIV","ChrV","ChrX")
gdata<-gdata[gdata$Chr %in% chrlist,]
gdata<-t(gdata)
colnames(gdata)<-paste(gdata[1,],gdata[2,],sep="-")
gdata<-cbind(c("","",as.numeric(substr(rownames(gdata)[3:105],2,5))),gdata)
colnames(gdata)[1]<-"id"

## Save in the csv file
write.csv(gdata,file="Dataset/Genotype_Pexs_EMS_hybrid_RQTL.csv",quote=F,row.names=F)

###### Linkage analysis ######

# Load cross data
mapthis<-read.cross("csv",file="Dataset/Genotype_Pexs_EMS_hybrid_RQTL.csv",estimate.map=F,crosstype="bc")

# Drop errorneous positions (binomial test)
gt<-geno.table(mapthis)
geno.table(mapthis)
gt[gt$P.value < 0.05/totmar(mapthis),]
todrop<-rownames(gt[gt$P.value < 0.05/totmar(mapthis),])
mapthis<-drop.markers(mapthis,todrop)

# Individuals with >200 informative loci were used. (N=103 => N=96)
plotMissing(mapthis)
mapthis<-subset(mapthis,ind=(ntyped(mapthis)>200))

# Linkage analysis
gddata<-est.map(mapthis,error.prob=0.001,verbose=F)
mapdata<-replace.map(mapthis,gddata)
mapdata<-calc.errorlod(mapdata)

## Drop markers that increase genetic distance by more than 10cM or reduce LOD
dropone <- droponemarker(mapdata, error.prob=0.005)
dropone[dropone$LOD>0|dropone$Ldiff>10,]
todrop<-rownames(dropone[dropone$LOD>0|dropone$Ldiff>10,])
mapdata<-drop.markers(mapdata,todrop)

## Reanalyze linkage map without the dropped markers
gddata<-est.map(mapdata,error.prob=0.001,verbose=F)
mapdata<-replace.map(mapdata,gddata)
mapdata<-calc.errorlod(mapdata)

## Drop again
dropone <- droponemarker(mapdata, error.prob=0.005)
dropone[dropone$LOD>0|dropone$Ldiff>10,]
todrop<-rownames(dropone[dropone$LOD>0|dropone$Ldiff>10,])
mapdata<-drop.markers(mapdata,todrop)

## Once more, reanalyze linkage map without the dropped markers
gddata<-est.map(mapdata,error.prob=0.001,verbose=F)
mapdata<-replace.map(mapdata,gddata)
mapdata<-calc.errorlod(mapdata)

## Save mapdata
save(mapdata,file="mapdata_Pexspectaus_EMS_hybrids.Robj")

plotMap(mapdata, show.marker.names=TRUE)

###### Marey map ######

## Data preparation

chrlist2<-c("I","II","III","IV","V","X")
res<-c()
for(i in 1:6){
  j<-chrlist2[i]
  value<-as.vector(pull.map(mapdata,j)[[1]])
  marker<-names(pull.map(mapdata,j)[[1]])
  res.sub<-data.frame(matrix(unlist(strsplit(marker,"-")),byrow=T,ncol=2))
  res.sub<-cbind(res.sub,value)
  colnames(res.sub)<-c("Chr","Pos","cM")
  res.sub$Pos<-as.numeric(res.sub$Pos)
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
