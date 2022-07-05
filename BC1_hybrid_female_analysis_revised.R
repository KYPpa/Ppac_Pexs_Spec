library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(qtl)

###### Load and format genotype data #####

# Read file of sliding window data of genotypes
# Column
#1. File name (File)
#2. Chromosome name (Chr)
#3. Start position (Start)
#4. End position (End)
#5. Count of the number of SNPs with Ppac allele (Ppaccount)
#6. Count of the number of SNPs with Ppac allele (Pexscount)

sw.data<-read.table("Dataset/Genotype_data_BC1_hybrid_female.txt",sep="\t",header=T)

sw.data<-sw.data[sw.data$Chr!="ChrX",]
sw.data$Pos<-(sw.data$End+sw.data$Start-1)/2000000
sw.data$SeqID<-factor(as.numeric(substr(sw.data$File,23,25)))
sw.data$Chr<-factor(sw.data$Chr)
sw.data$Ratio<-sw.data$Pexscount/(sw.data$Pexscount+sw.data$Ppaccount)

sw.data$SmoothPred<-rep(NA)
sample<-levels(sw.data$SeqID)
sel.chr<-c("ChrI","ChrII","ChrIII","ChrIV","ChrV")
for(i in sample){
  for(chr in sel.chr){
    sw.data.sub<-
      sw.data %>%
      filter(Chr==chr,SeqID==i)
    span<-40/nrow(sw.data.sub)
    fit<-loess(Ratio~Pos,data=sw.data.sub,degree=1,span=span)
    loessfit<-fit$fitted
    na.list<-which(is.na(sw.data.sub$Ratio))
    if(length(na.list)>0){
      na.list<-na.list[order(na.list)]
      for(j in na.list){
        loessfit<-append(loessfit,NA,after=j-1)
      }
    }
    sw.data$SmoothPred[sw.data$Chr==chr&sw.data$SeqID==i]<-loessfit
  }
}

# Histogram to decide the boudary of ralt (here SmoothPred)
xtitle<-expression(paste(italic("r")[alt]," in backcross with ", italic("P. exspectatus")))
pdf(file="Fig_ralt_histogram_BC_to_Pexs.pdf",width=4,height=2)
ggplot(sw.data)+geom_histogram(aes(x=SmoothPred))+
annotate(geom="segment",x=0.7,xend=0.7,y=10000,yend=2000,arrow=arrow(length=unit(10,"pt")))+
  labs(x=xtitle,y="Count")+
theme(panel.border=element_rect(linetype="solid",fill=NA))
dev.off()


# Determine introgression
sw.data$Int<-ifelse(sw.data$SmoothPred>0.7,0,1)

###### Binding to phenotype data ######

# Read file of phenotype data
# Column
#1. Sequence ID (SeqID)
#2. Total no. of progeny (Total)

pt.data<-read.table("Dataset/Phenotype_data_BC1_hybrid_female.txt",sep="\t",header=T)
pt.data$SeqID<-factor(as.numeric(substr(pt.data$SeqID,6,9)),level=unique(sw.data$SeqID))
pt.data<-pt.data[!is.na(pt.data$SeqID),]
pt.data$PA_Prog<-ifelse(pt.data$Total>0,"P",
                        ifelse(pt.data$Total==0,"A",NA))
ptindex<-order(pt.data$Total)
orderedseqid<-pt.data$SeqID[ptindex]
orderedlevel<-match(orderedseqid,as.numeric(as.character(unique(sw.data$SeqID))))
sw.data$SeqID<-factor(sw.data$SeqID,level=unique(sw.data$SeqID)[orderedlevel])
pt.data$SeqID<-factor(pt.data$SeqID,level=levels(sw.data$SeqID))
sw.data$Name<-sw.data$SeqID
levels(sw.data$Name)<-paste(levels(sw.data$Name),pt.data$Total[ptindex],sep=":  ")


###### Introgression map ######

hlinepos<-sum(pt.data$PA_Prog[ptindex]=="A")+0.5

g1<-sw.data %>%
  mutate(Int=factor(Int)) %>%
  ggplot()+
  geom_tile(aes(x=Pos,y=Name,fill=Int,col=Int))+
  geom_hline(yintercept=hlinepos,colour="yellow3")+
  facet_grid(~Chr,scales="free_x",space="free_y")+
  scale_colour_manual(values=c("lightskyblue1","magenta4"))+
  scale_fill_manual(values=c("lightskyblue1","magenta4"))+
  labs(x="Position in P. pacificus genome",y="F2 female Individuals [Pexs x (Pexs x Ppac)]")+
  theme(  axis.text.y=element_blank(),
          legend.position="none",      
          panel.border=element_rect(linetype="solid",fill=NA),
          strip.background = element_rect(colour = "black"),
          panel.spacing = unit(0, "lines"))
#  ggsave(file="IntMap_F2male_all_210211.pdf",device="pdf",width=5,height=8)
pt.data$Order<-factor(match(1:nrow(pt.data),ptindex),level=1:length(ptindex))
g2<-pt.data %>%
  ggplot()+geom_bar(aes(x=Order,y=log2(Total+1)),stat="identity")+
  theme_classic()+
  scale_y_continuous(expand = c(0,0))+
  labs(y="log2(Count+1)")+
  coord_flip()+
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        plot.margin = unit(c(22,10,5.5,0), "points"))

ggarrange(g1,g2,nrow=1,widths=c(1,0.25))

###### Proportion of individuals w/ introgressions and Fisher's exact test ######

## Proportion

g1<-sw.data %>%
  left_join(pt.data,by=c("SeqID")) %>%
  group_by(Chr,Pos) %>%
  summarise(Presence=mean(ifelse(PA_Prog=="P",Int,NA),na.rm=T),
            Absence=mean(ifelse(PA_Prog=="A",Int,NA),na.rm=T)) %>%
  pivot_longer(cols=c("Presence","Absence"),names_to="Progeny",values_to="PropInt")%>%
  ggplot()+geom_line(aes(x=Pos,y=PropInt,col=Progeny))+
  facet_grid(.~Chr,scales="free",space="free")+
  labs(x="Position in P. pacificus genome",y="Proportion of individuals\nwith introgression")+
  theme(  text=element_text(size=9),
          axis.title.x=element_blank(),
          panel.border=element_rect(linetype="solid",fill=NA),
          strip.background = element_rect(colour = "black"),
          panel.spacing = unit(0, "lines"))

## Fisher's exact test

fisher.table<-
  sw.data %>%
  left_join(pt.data,by=c("SeqID")) %>%
  group_by(Chr,Pos) %>%
  summarise(PreInt=sum(ifelse(PA_Prog=="P",Int,0),na.rm=T),
            PreNo=sum(ifelse(PA_Prog=="P",1-Int,0),na.rm=T),
            AbInt=sum(ifelse(PA_Prog=="A",Int,0),na.rm=T),
            AbNo=sum(ifelse(PA_Prog=="A",1-Int,0),na.rm=T))

na.row<-attr(na.omit(fisher.table),"na.action")
p.fisher<-function(vec){
  fisher.test(matrix(vec,nrow=2,ncol=2))$p.value
}
fisher.table$logp<--log10(apply(data.matrix(fisher.table[,-c(1,2)]),1,p.fisher))

# QTL peak

fisher.table %>%
  filter(logp>3.261548) %>%
  summarise(min(Pos),max(Pos))

g2<-fisher.table %>%
  ggplot()+geom_line(aes(x=Pos,y=logp))+
  geom_hline(yintercept=2.731716,linetype=2)+
  facet_grid(.~Chr,scales="free",space="free")+
  labs(x="Position in P. pacificus genome",y="log10 p-value")+
  theme(  text=element_text(size=9),
          panel.border=element_rect(linetype="solid",fill=NA),
          strip.background = element_rect(colour = "black"),
          axis.title.y=element_text(margin=margin(0,5,0,13)),
          panel.spacing = unit(0, "lines"),
          plot.margin=unit(c(5.5,79,5.5,5.5), "points"))
ggarrange(g1,g2,nrow=2,heights=c(1,1))


###### Permutation test to calculate significant p-value ######

pnum<-nrow(pt.data[pt.data$PA_Prog=="P",])
maxlogp<-c()
for(i in 1:1000){
  print(i)
  selection<-sample(pt.data$SeqID,pnum)
  fisher.table<-
    sw.data %>%
    mutate(PA_Prog=ifelse(SeqID %in% selection,"P","A")) %>%
    na.omit %>%
    group_by(Chr,Pos) %>%
    summarise(PreInt=sum(ifelse(PA_Prog=="P",Int,0),na.rm=T),
              PreNo=sum(ifelse(PA_Prog=="P",1-Int,0),na.rm=T),
              AbInt=sum(ifelse(PA_Prog=="A",Int,0),na.rm=T),
              AbNo=sum(ifelse(PA_Prog=="A",1-Int,0),na.rm=T))
  p.fisher<-function(vec){
    fisher.test(matrix(vec,nrow=2,ncol=2))$p.value
  }
  logp<--log10(apply(data.matrix(fisher.table[,-c(1,2)]),1,p.fisher))
  maxlogp[i]<-max(logp)
}

quantile(maxlogp,0.95) #2.848579
quantile(maxlogp,0.99) #3.261548

###### R/QTL QTL ######

## Create csv file for cross data

chrvec<-c("c1","c2","c3","c4","c5")
names(chrvec)<-(c("ChrI","ChrII","ChrIII","ChrIV","ChrV"))
gt.data<-sw.data %>%
  mutate(MK=paste(chrvec[as.character(Chr)],Pos,sep="_"),Genotype=as.character(ifelse(is.na(Int),"-",ifelse(Int==1,"H","A")))) %>%
  pivot_wider(id_cols=c("SeqID"),names_from=c("MK"),values_from=c("Genotype"))
ptindex<-match(gt.data$SeqID,pt.data$SeqID)
pheno<-ifelse(pt.data$PA_Prog[ptindex]=="A",0,1)
gt.data<-cbind(pheno,gt.data)
colnames(gt.data)[2]<-"id"
gt.data$id<-as.character(gt.data$id)
gt.data$pheno<-as.character(gt.data$pheno)
chrarray<-chrvec[as.character(sw.data$Chr)]
mbarray<-sw.data$Pos-0.05
gt.data<-rbind(c("","",chrarray),c("","",mbarray),gt.data)
write.csv(gt.data,file="Dataset/Genotype_BC1_female_RQTL_1.csv",row.names=F,quote=F)

## Generate cross data

mapthis<-read.cross("csv",file="Dataset/Genotype_BC1_female_RQTL_1.csv",estimate.map=F,crosstype="bc")
gddata <- est.map(mapthis, error.prob=0.001, verbose=FALSE)
mapdata<-replace.map(mapthis,gddata)
mapdata<-calc.errorlod(mapdata)

## Sigle QTL mapping

out.bin<-scanone(mapdata,model="binary")
plot(out.bin)

## Permutation

operm.bin<-scanone(mapdata,model="binary",n.perm=1000,perm.Xsp=TRUE)
summary(operm.bin, alpha=0.05) #2.28
summary(out.bin, perms=operm.bin, alpha=0.05, pvalues=TRUE) #c1 52.9

## Fit to qtl model

mapdata2<-sim.geno(mapdata, step=2, n.draws=128, err=0.001)
effectplot(mapdata2,mname1="c1_28.75")
qtl<-makeqtl(mapdata2,chr=c("c1"),pos=c(52.9))
qtl
plot(qtl)
out.fq<-fitqtl(mapdata2,qtl=qtl,model="binary")
summary(out.fq) # 29.71%


###### Summarise coverage data of F1 ######

# Read file of sliding window data of coverage of F1 female
# Column
#1. File name (File)
#2. Chromosome name (Chr)
#3. Start position (Start)
#4. End position (End)
#5. Base pair count (Count)

f1.data<-read.table(file="Dataset/Coverage_data_F1_hybrid.txt",header=T)
f1.data<-f1.data[f1.data$Chr!="ChrX",]
f1.data$Pos<-(f1.data$End+f1.data$Start-1)/2000000
f1.data$SeqID<-factor(substr(f1.data$File,7,8))
f1.data$Chr<-factor(f1.data$Chr,levels=c("ChrI","ChrII","ChrIII","ChrIV","ChrV"))
f1.data<-na.omit(f1.data)

# Select 71:75 that are females of Pexs female vs. Ppacificus male
f1.data<-f1.data[as.numeric(as.character(f1.data$SeqID)) %in% 71:75,]
f1ac.info<-f1.data %>%
  group_by(SeqID) %>%
  summarise(Avecov=mean(Count/100000))

f1.data<-f1.data %>%
  left_join(f1ac.info,by=c("SeqID")) %>%
  mutate(Norcov=log2(Count/100000/Avecov))

f1.data$DepthPred<-rep(NA)
sample<-71:75
sel.chr<-levels(f1.data$Chr)
for(i in sample){
  for(chr in sel.chr){
    f1.data.sub<-
      f1.data %>%
      filter(Chr==chr,as.character(SeqID)==i)
    span<-40/nrow(f1.data.sub)
    fit<-loess(Norcov~Pos,data=f1.data.sub,degree=1,span=span)
    loessfit<-fit$fitted
    na.list<-which(is.na(f1.data.sub$Norcov))
    if(length(na.list)>0){
      na.list<-na.list[order(na.list)]
      for(j in na.list){
        loessfit<-append(loessfit,NA,after=j-1)
      }
    }
    f1.data$DepthPred[f1.data$Chr==chr&f1.data$SeqID==i]<-loessfit
  }
}

f1.cov.ref<-
  f1.data %>%
  group_by(Chr,Pos) %>%
  summarise(F1cov=mean(DepthPred))

###### Summarise coverage data of BC1 female ######

# Read file of sliding window data of coverage of BC1 female
# Column
#1. File name (File)
#2. Chromosome name (Chr)
#3. Start position (Start)
#4. End position (End)
#5. Base pair count (Count)

dp.data<-read.table(file="Dataset/Coverage_data_BC1_hybrid_female.txt",header=T)
dp.data<-dp.data[dp.data$Chr!="ChrX",]
dp.data$Pos<-(dp.data$End+dp.data$Start-1)/2000000
dp.data$SeqID<-factor(as.numeric(substr(dp.data$File,24,26)),level=levels(sw.data$SeqID))
dp.data<-dp.data[!is.na(dp.data$SeqID),]
dp.data$Chr<-factor(dp.data$Chr)
ac.info<-dp.data %>%
  group_by(SeqID) %>%
  summarise(Avecov=mean(Count/100000))
dp.data<-dp.data %>%
  left_join(ac.info,by=c("SeqID")) %>%
  mutate(Norcov=log2(Count/100000/Avecov))

dp.data$DepthPred<-rep(NA)
sample<-levels(dp.data$SeqID)
sel.chr<-c("ChrI","ChrII","ChrIII","ChrIV","ChrV")
for(i in sample){
  for(chr in sel.chr){
    dp.data.sub<-
      dp.data %>%
      filter(Chr==chr,SeqID==i)
    span<-40/nrow(dp.data.sub)
    fit<-loess(Norcov~Pos,data=dp.data.sub,degree=1,span=span)
    loessfit<-fit$fitted
    na.list<-which(is.na(dp.data.sub$Norcov))
    if(length(na.list)>0){
      na.list<-na.list[order(na.list)]
      for(j in na.list){
        loessfit<-append(loessfit,NA,after=j-1)
      }
    }
    dp.data$DepthPred[dp.data$Chr==chr&dp.data$SeqID==i]<-loessfit
  }
}

dp.data <-
  dp.data %>% left_join(f1.cov.ref,by=c("Chr","Pos")) %>%
  mutate(NorDepthPred=DepthPred-F1cov)

###### Analysis of ploidy by coverage depth ######

## Violin plot of the difference in normalized coverage depth (BC1-F1)

sel.chr.new<-c("ChrIL","ChrIR","ChrII","ChrIII","ChrIV","ChrV")

violin.data<-sw.data %>%
  left_join(dp.data,by=c("SeqID","Chr","Start","End","Pos")) %>%
  dplyr::filter(!(Chr=="ChrI"&Pos>20.05&Pos<22.05)) %>%
  na.omit %>%
  mutate(Intcov=Int*NorDepthPred) %>%
  mutate(Chr_new=factor(ifelse(Chr=="ChrI",ifelse(Pos>=21.55,"ChrIR",ifelse(Pos<20.55,"ChrIL",NA)),as.character(Chr)),level=sel.chr.new)) %>%
  na.omit %>%
  group_by(Chr_new,SeqID) %>%
  summarise(Intcov=sum(Intcov)/sum(Int))

violin.data.female<-violin.data

## See histogram of the difference in normalized coverage depth (BC1-F1)

sw.data %>%
  left_join(dp.data,by=c("SeqID","Chr","Start","End","Pos")) %>%
  dplyr::filter(!(Chr=="ChrI"&Pos>20.05&Pos<22.05)) %>%
  na.omit %>%
  mutate(Intcov=Int*NorDepthPred) %>%
  mutate(Chr_new=factor(ifelse(Chr=="ChrI",ifelse(Pos>=21.55,"ChrIR",ifelse(Pos<20.55,"ChrIL",NA)),as.character(Chr)),level=sel.chr.new)) %>%
  na.omit %>%
  group_by(Chr_new,SeqID) %>%
  summarise(Intcov=sum(Intcov)/sum(Int)) %>%
  ggplot()+geom_histogram(aes(x=Intcov))

## Identify the trisomy animals

trisomy.data<-sw.data %>%
  left_join(dp.data,by=c("SeqID","Chr","Start","End","Pos")) %>%
  dplyr::filter(!(Chr=="ChrI"&Pos>20.05&Pos<22.05)) %>%
  na.omit %>%
  mutate(Intcov=Int*NorDepthPred) %>%
  mutate(Chr_new=factor(ifelse(Chr=="ChrI",ifelse(Pos>=22.05,"ChrIR","ChrIL"),as.character(Chr)),level=sel.chr.new)) %>%
  group_by(Chr_new,SeqID) %>%
  summarise(Intcov=sum(Intcov)/sum(Int))

other.chr<-c("ChrIL","ChrII","ChrIII","ChrIV","ChrV")
plist<-c()
countlist<-c()
for(j in 1:5){
  testc<-c(0,0,0,0)
  testc[1]<-sum(na.omit(trisomy.data$Intcov>0.35&trisomy.data$Chr_new=="ChrIR"))
  testc[2]<-(75-testc[1])
  testc[3]<-sum(na.omit(trisomy.data$Intcov>0.35&trisomy.data$Chr_new==other.chr[j]))
  testc[4]<-(75-testc[3])
  testsummary<-fisher.test(matrix(testc,nrow=2,byrow=T))
  countlist<-c(countlist,testc[3])
  plist<-c(plist,testsummary$p.value)
}
testc[1]
testc[1]/75
median(countlist)
median(countlist)/75
median(plist)

# Test of association of trisomy of ChrIR with fertility

chrir<-sw.data %>%
  left_join(dp.data,by=c("SeqID","Chr","Start","End","Pos")) %>%
  dplyr::filter(!(Chr=="ChrI"&Pos>20.05&Pos<22.05)) %>%
  na.omit %>%
  mutate(Intcov=Int*NorDepthPred) %>%
  mutate(Chr_new=factor(ifelse(Chr=="ChrI",ifelse(Pos>=22.05,"ChrIR","ChrIL"),as.character(Chr)),level=sel.chr.new)) %>%
  group_by(Chr_new,SeqID) %>%
  summarise(Intcov=sum(Intcov)/sum(Int))%>%
  dplyr::filter(Intcov>0.35,Chr_new=="ChrIR")

pt.data$PA_Prog<-factor(pt.data$PA_Prog,level=c("A","P"))
selected<-table(pt.data$PA_Prog[as.numeric(as.character(pt.data$SeqID)) %in% as.numeric(as.character(chrir$SeqID))])
notselected<-table(pt.data$PA_Prog[!(as.numeric(as.character(pt.data$SeqID)) %in% as.numeric(as.character(chrir$SeqID)))])
fisher.test(x=matrix(c(selected,notselected),nrow=2,byrow=T))
c(selected,notselected)
c(selected,notselected)/rep(c(sum(selected),sum(notselected)),each=2)
