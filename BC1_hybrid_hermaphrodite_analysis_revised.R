library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(qtl)

###### Load and format genotype data ######

# Read file of sliding window data of genotypes
# Column
#1. File name (File)
#2. Chromosome name (Chr)
#3. Start position (Start)
#4. End position (End)
#5. Count of the number of SNPs with Ppac allele (Ppaccount)
#6. Count of the number of SNPs with Ppac allele (Pexscount)

sw.data<-read.table("Dataset/Genotype_data_BC1_hybrid_hermaphrodite.txt",sep="\t",header=T)
sw.data<-sw.data[sw.data$Chr!="ChrX",]
sw.data$Pos<-(sw.data$End+sw.data$Start-1)/2000000
sw.data$SeqID<-factor(as.numeric(substr(sw.data$File,23,25)))
sw.data$Chr<-factor(sw.data$Chr)
sw.data$Ratio<-sw.data$Ppaccount/(sw.data$Pexscount+sw.data$Ppaccount)

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

## histogram
ggplot(sw.data)+geom_histogram(aes(x=SmoothPred))


sw.data$Int<-ifelse(sw.data$SmoothPred<0.90,1,0)

int.no.list<-sw.data %>%
  group_by(SeqID) %>%
  summarise(IntNo=sum(Int,na.rm=TRUE))

int.no.list$SeqID[int.no.list$IntNo==max(int.no.list$IntNo)]
omit.SeqID<-c("45","3","109","119","68")#control
omit.SeqID<-unique(c(omit.SeqID,as.character(int.no.list$SeqID[int.no.list$IntNo==0])))
sum((as.character(sw.data$SeqID) %in% omit.SeqID))
sw.data<-sw.data[!(as.character(sw.data$SeqID) %in% omit.SeqID),]

###### Binding to phenotype data ######


## Read file of phenotype data
# Column
#1. Sequence ID (SeqID)
#2. Total no. of progeny (Total)
#3. State of fertility or control/test (State):
#C, Control individuals (Those are already removed from sw.data above)
#P, presence of progeny
#FE, presence of fertilized eggs but not progeny
#UE, presence of unfertilized eggs but neither fertilized eggs or progeny
#N, no progeny and eggs

pt.data<-read.table("Dataset/Phenotype_data_BC1_hybrid_hermaphrodite.txt",sep="\t",header=T)
pt.data$SeqID<-factor(as.numeric(substr(pt.data$SeqID,8,10)),level=unique(sw.data$SeqID))
pt.data<-pt.data[!is.na(pt.data$SeqID),]
pt.data$PA_Egg<-ifelse(pt.data$State=="C"|pt.data$State=="P"|pt.data$State=="UE"|pt.data$State=="FE","P",
                       ifelse(pt.data$State=="N","A",NA))
pt.data$PA_Prog<-ifelse(pt.data$State=="C"|pt.data$State=="P","P",
                        ifelse(pt.data$State=="F"|pt.data$State=="N"|pt.data$State=="UE"|pt.data$State=="FE","A",NA))
ptindex<-order(pt.data$PA_Egg,pt.data$Total)
orderedseqid<-pt.data$SeqID[ptindex]
orderedlevel<-match(as.character(orderedseqid),unique(sw.data$SeqID))
sw.data$SeqID<-factor(sw.data$SeqID,level=unique(sw.data$SeqID)[orderedlevel])
pt.data$SeqID<-factor(pt.data$SeqID,level=levels(sw.data$SeqID))
sw.data<-sw.data[!(is.na(sw.data$SeqID)),]
sw.data$Name<-sw.data$SeqID
levels(sw.data$Name)<-paste(levels(sw.data$Name),pt.data$Total[ptindex],pt.data$PA_Prog[ptindex],sep=":  ")

###### Introgression map ######

hlinepos1<-sum(pt.data$PA_Prog[ptindex]=="A")+0.5
hlinepos2<-sum(pt.data$PA_Egg[ptindex]=="A")+0.5

g1<-sw.data %>%
  mutate(Int=factor(Int)) %>%
  ggplot()+
  geom_tile(aes(x=Pos,y=Name,fill=Int,col=Int))+
  geom_hline(yintercept=hlinepos1,colour="yellow3")+
  geom_hline(yintercept=hlinepos2,colour="green3")+
  facet_grid(~Chr,scales="free_x",space="free_y")+
  scale_colour_manual(values=c("indianred1","magenta4"))+
  scale_fill_manual(values=c("indianred1","magenta4"))+
  labs(x="Position in P. pacificus genome",y="F2 hermaphrodite Individuals [Ppac x (Ppac x Pexs)]")+
  theme(  axis.text.y=element_blank(),
          legend.position="none",      
          panel.border=element_rect(linetype="solid",fill=NA),
          strip.background = element_rect(colour = "black"),
          panel.spacing = unit(0, "lines"))
#  ggsave(file="IntMap_F2male_all_210211.pdf",device="pdf",width=5,height=8)
pt.data$Order<-factor(match(1:nrow(pt.data),ptindex),level=1:length(ptindex))
pt.data<-pt.data[pt.data$SeqID %in% unique(sw.data$SeqID),]
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

## Figure of proportion of introgression

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
  filter(logp>3.429958) %>%
  summarise(min(Pos),max(Pos))

fisher.table %>%
  dplyr::filter(Chr=="ChrI"|Chr=="ChrV",logp==max(logp)) %>%
  data.frame

# Figure of log p-value of Fisher's exact test

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

quantile(maxlogp,0.95) #2.759388
quantile(maxlogp,0.99) #3.429958

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
write.csv(gt.data,file="Dataset/Genotype_data_BC1_hermaphrodite_RQTL_1.txt",row.names=F,quote=F)

## Generate cross data

mapthis<-read.cross("csv",file="Dataset/Genotype_data_BC1_hermaphrodite_RQTL_1.txt",estimate.map=F,crosstype="bc")
gddata <- est.map(mapthis, error.prob=0.001, verbose=FALSE)
mapdata<-replace.map(mapthis,gddata)
mapdata<-calc.errorlod(mapdata)

## Sigle QTL mapping

out.bin<-scanone(mapdata,model="binary")
plot(out.bin)

## Permutation

operm.bin<-scanone(mapdata,model="binary",n.perm=1000,perm.Xsp=TRUE)
summary(operm.bin, alpha=0.05) #2.33
summary(out.bin, perms=operm.bin, alpha=0.05, pvalues=TRUE)

## Fit to qtl model

mapdata2<-sim.geno(mapdata, step=2, n.draws=128, err=0.001)
effectplot(mapdata2,mname1="c1_31.65")
effectplot(mapdata2,mname1="c5_5.75")
qtl<-makeqtl(mapdata2,chr=c("c1","c5"),pos=c(83.554,0.741))
qtl
plot(qtl)
out.fq<-fitqtl(mapdata2,qtl=qtl,model="binary")
summary(out.fq)
# 20.069%
# 9.979%

###### Summarise coverage data of F1 ######

# Read file of sliding window data of coverage of F1 male
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

## Select 81:85 - F1 female of Ppac hermaphrodite vs. Pexs male
f1.data<-f1.data[as.numeric(as.character(f1.data$SeqID)) %in% 81:85,]
f1ac.info<-f1.data %>%
  group_by(SeqID) %>%
  summarise(Avecov=mean(Count/100000))
f1.data<-f1.data %>%
  left_join(f1ac.info,by=c("SeqID")) %>%
  mutate(Norcov=log2(Count/100000/Avecov))

f1.data %>%
  ggplot()+geom_line(aes(x=Pos,y=Norcov,col=SeqID))+
  facet_grid(.~Chr,scales="free_x",space="free_x")+
  scale_fill_continuous(low="gray85",high="black")

f1.data$DepthPred<-rep(NA)
sample<-81:85
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


f1.data %>%
  ggplot()+geom_line(aes(x=Pos,y=DepthPred,col=SeqID))+
  facet_grid(.~Chr,scales="free_x",space="free_x")+
  scale_fill_continuous(low="gray85",high="black")

f1.cov.ref<-
  f1.data %>%
  group_by(Chr,Pos) %>%
  summarise(F1cov=mean(DepthPred))

###### Summarise coverage data of BC1 male  ######

# Read file of sliding window data of coverage of BC1 male
# Column
#1. File name (File)
#2. Chromosome name (Chr)
#3. Start position (Start)
#4. End position (End)
#5. Base pair count (Count)

dp.data<-read.table("Dataset/Coverage_data_BC1_hybrid_hermaphrodite.txt",sep="\t",header=T)
dp.data<-dp.data[dp.data$Chr!="ChrX",]
dp.data$Pos<-(dp.data$End+dp.data$Start-1)/2000000
#levels(factor(substr(dp.data$File,23,25)))
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
  summarise(Intcov=sum(Intcov)/sum(Int)) %>%
  dplyr::filter(!is.na(Intcov))

violin.data.hermaphrodite<-violin.data

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

other.chr<-c("ChrII","ChrIII","ChrIV","ChrV")
plist<-c()
countlist<-c()
for(j in 1:4){
  testc<-c(0,0,0,0)
  testc[1]<-sum(na.omit(trisomy.data$Intcov>0.25&trisomy.data$Chr_new=="ChrIL"))
  testc[2]<-(136-testc[1])
  testc[3]<-sum(na.omit(trisomy.data$Intcov>0.25&trisomy.data$Chr_new==other.chr[j]))
  testc[4]<-(136-testc[3])
  testsummary<-fisher.test(matrix(testc,nrow=2,byrow=T))
  countlist<-c(countlist,testc[3])
  plist<-c(plist,testsummary$p.value)
}
testc[1]
testc[1]/136
median(countlist)
median(countlist)/136
median(plist)

plist<-c()
countlist<-c()
for(j in 1:4){
  testc<-c(0,0,0,0)
  testc[1]<-sum(na.omit(trisomy.data$Intcov>0.25&trisomy.data$Chr_new=="ChrIR"))
  testc[2]<-(136-testc[1])
  testc[3]<-sum(na.omit(trisomy.data$Intcov>0.25&trisomy.data$Chr_new==other.chr[j]))
  testc[4]<-(136-testc[3])
  testsummary<-fisher.test(matrix(testc,nrow=2,byrow=T))
  countlist<-c(countlist,testc[3])
  plist<-c(plist,testsummary$p.value)
}
testc[1]
testc[1]/136
median(countlist)
median(countlist)/136
median(plist)

# Test of association of trisomy of ChrIL with fertility

chril<-sw.data %>%
  left_join(dp.data,by=c("SeqID","Chr","Start","End","Pos")) %>%
  dplyr::filter(!(Chr=="ChrI"&Pos>20.05&Pos<22.05)) %>%
  na.omit %>%
  mutate(Intcov=Int*NorDepthPred) %>%
  mutate(Chr_new=factor(ifelse(Chr=="ChrI",ifelse(Pos>=22.05,"ChrIR","ChrIL"),as.character(Chr)),level=sel.chr.new)) %>%
  group_by(Chr_new,SeqID) %>%
  summarise(Intcov=sum(Intcov)/sum(Int))%>%
  dplyr::filter(Intcov>0.25,Chr_new=="ChrIL")
pt.data$PA_Prog<-factor(pt.data$PA_Prog,level=c("A","P"))
selected<-table(pt.data$PA_Prog[as.numeric(as.character(pt.data$SeqID)) %in% as.numeric(as.character(chril$SeqID))])
notselected<-table(pt.data$PA_Prog[!(as.numeric(as.character(pt.data$SeqID)) %in% as.numeric(as.character(chril$SeqID)))])
fisher.test(x=matrix(c(selected,notselected),nrow=2,byrow=T))
c(selected,notselected)
c(selected,notselected)/rep(c(sum(selected),sum(notselected)),each=2)

# Test of association of trisomy of ChrIR with fertility

chrir<-sw.data %>%
  left_join(dp.data,by=c("SeqID","Chr","Start","End","Pos")) %>%
  dplyr::filter(!(Chr=="ChrI"&Pos>20.05&Pos<22.05)) %>%
  na.omit %>%
  mutate(Intcov=Int*NorDepthPred) %>%
  mutate(Chr_new=factor(ifelse(Chr=="ChrI",ifelse(Pos>=22.05,"ChrIR","ChrIL"),as.character(Chr)),level=sel.chr.new)) %>%
  group_by(Chr_new,SeqID) %>%
  summarise(Intcov=sum(Intcov)/sum(Int))%>%
  dplyr::filter(Intcov>0.25,Chr_new=="ChrIR")
pt.data$PA_Prog<-factor(pt.data$PA_Prog,level=c("A","P"))
selected<-table(pt.data$PA_Prog[as.numeric(as.character(pt.data$SeqID)) %in% as.numeric(as.character(chrir$SeqID))])
notselected<-table(pt.data$PA_Prog[!(as.numeric(as.character(pt.data$SeqID)) %in% as.numeric(as.character(chrir$SeqID)))])
fisher.test(x=matrix(c(selected,notselected),nrow=2,byrow=T))
c(selected,notselected)
c(selected,notselected)/rep(c(sum(selected),sum(notselected)),each=2)

# Individuals not categorized as trisomy but supicious because of the recombination on the break point
# ChrIL yes => ChrIR no, 9 individuals.
BR_Zero<-sw.data %>%
  dplyr::filter(!(SeqID %in% chril$SeqID),Chr=="ChrI",Pos==22.05,Int==0)
additional_Trisomy1<-sw.data %>%
  dplyr::filter(!(SeqID %in% chril$SeqID),Chr=="ChrI",Pos==20.05,Int==1) %>%
  semi_join(BR_Zero, by="SeqID")

# ChrIL no => ChrIR yes, 13 individuals.
BR_One<-sw.data %>%
  dplyr::filter(!(SeqID %in% chril$SeqID),Chr=="ChrI",Pos==22.05,Int==1)
additional_Trisomy2<-sw.data %>%
  dplyr::filter(!(SeqID %in% chril$SeqID),Chr=="ChrI",Pos==20.05,Int==0) %>%
  semi_join(BR_One, by="SeqID")

###### Map data for Recombination analysis ######

chrvec<-c("c1","c2","c3","c4","c5")
names(chrvec)<-(c("ChrI","ChrII","ChrIII","ChrIV","ChrV"))
gt.data<-sw.data %>%
  filter(!(SeqID %in% chril$SeqID),!(SeqID %in% chrir$SeqID)) %>%
  filter(!(SeqID %in% additional_Trisomy1$SeqID),!(SeqID %in% additional_Trisomy2$SeqID)) %>%
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
write.csv(gt.data,file="Dataset/Genotype_BC1_hermaphrodite_RQTL_2.csv",row.names=F,quote=F)

mapthis<-read.cross("csv",file="Dataset/Genotype_BC1_hermaphrodite_RQTL_2.csv",estimate.map=F,crosstype="bc")
gddata <- est.map(mapthis, error.prob=0.001, verbose=FALSE)
mapdata<-replace.map(mapthis,gddata)
save(mapdata,file="mapdata_11RSG_wo_Trisomy_XXXXXX.Robj")

###### Marey map ######

## Data preparation

res<-c()
for(i in 1:5){
  value<-as.vector(pull.map(mapdata,chrvec[i])[[1]])
  marker<-names(pull.map(mapdata,chrvec[i])[[1]])
  res.sub<-data.frame(matrix(unlist(strsplit(marker,"_")),byrow=T,ncol=2))
  res.sub<-cbind(res.sub,value)
  colnames(res.sub)<-c("Chr","Pos","cM")
  res.sub$Pos<-as.numeric(res.sub$Pos)
  res.sub$Chr<-rep(names(chrvec)[i],nrow(res.sub))
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

