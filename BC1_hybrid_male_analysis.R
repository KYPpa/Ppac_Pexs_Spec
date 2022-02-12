library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(qtl)
library(patchwork)
library(grid)

###### SW data ######

# Read file of sliding window data of genotypes
# Column
#1. File name (File)
#2. Chromosome name (Chr)
#3. Start position (Start)
#4. End position (End)
#5. Count of the number of SNPs with Ppac allele (Ppaccount)
#6. Count of the number of SNPs with Ppac allele (Pexscount)
sw.data<-read.table("Dataset/Genotype_data_BC1_hybrid_male.txt",sep="\t",header=T)

sw.data<-sw.data[sw.data$Chr!="ChrX",]
sw.data$Pos<-(sw.data$End+sw.data$Start-1)/2000000
sw.data$SeqID<-factor(as.numeric(substr(sw.data$File,22,24)))
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

# Histogram to decide the boudary
ggplot(sw.data)+geom_histogram(aes(x=SmoothPred))

# Determine introgression
sw.data$Int<-ifelse(sw.data$SmoothPred<0.90,1,0)

int.no.list<-sw.data %>%
  group_by(SeqID) %>%
  summarise(IntNo=sum(Int,na.rm=TRUE))

# Remove error (F1 genotype)
int.no.list$SeqID[int.no.list$IntNo==max(int.no.list$IntNo)] #597
omit.SeqID<-c("597")
omit.SeqID<-c(omit.SeqID,unique(as.character(int.no.list$SeqID[int.no.list$IntNo==0])))
sum((sw.data$SeqID %in% omit.SeqID))
sw.data<-sw.data[!(as.character(sw.data$SeqID) %in% omit.SeqID),]

###### Binding to phenotype data ######

# Read file of sliding window data
# Column
#1. Sequence ID (SeqID)
#2. State of fertility (State):
  #P, presence of progeny
  #N, presence of eggs but not progeny
  #F, no progeny and eggs
#3. Total no. of progeny (Total)

pt.data<-read.table("Dataset/Phenotype_data_BC1_hybrid_male.txt",sep="\t",header=T)
pt.data<-na.omit(pt.data)
pt.data$SeqID<-factor(as.numeric(substr(pt.data$SeqID,5,7)),level=levels(sw.data$SeqID))
pt.data$PA_Egg<-ifelse(pt.data$State=="P"|pt.data$State=="N","P",
                       ifelse(pt.data$State=="F","A",NA))
pt.data$PA_Prog<-ifelse(pt.data$State=="P","P",
                        ifelse(pt.data$State=="F"|pt.data$State=="N","A",NA))
ptindex<-order(pt.data$PA_Egg,pt.data$Total)
orderedseqid<-pt.data$SeqID[ptindex]
orderedlevel<-match(as.character(orderedseqid),levels(sw.data$SeqID))
sw.data$SeqID<-factor(sw.data$SeqID,level=levels(sw.data$SeqID)[orderedlevel])
pt.data$SeqID<-factor(pt.data$SeqID,level=levels(sw.data$SeqID))
sw.data<-sw.data[!(is.na(sw.data$SeqID)),]
sw.data$Name<-sw.data$SeqID
levels(sw.data$Name)<-paste(levels(sw.data$Name),pt.data$Total[ptindex],pt.data$PA_Prog[ptindex],sep=":  ")

hlinepos<-sum(pt.data$PA_Prog[ptindex]=="A")+0.5

###### Introgression map ######

g1<-sw.data %>%
  mutate(Int=factor(Int)) %>%
  ggplot()+
  geom_tile(aes(x=Pos,y=Name,fill=Int,col=Int))+
  geom_hline(yintercept=hlinepos,colour="yellow3")+
  facet_grid(~Chr,scales="free_x",space="free_y")+
  scale_colour_manual(values=c("indianred1","magenta4"))+
  scale_fill_manual(values=c("indianred1","magenta4"))+
  labs(x="Position in P. pacificus genome",y="F2 male Individuals [Ppac x (Pexs x Ppac)]")+
  theme(  axis.text.y=element_blank(),
          legend.position="none",      
          panel.border=element_rect(linetype="solid",fill=NA),
                strip.background = element_rect(colour = "black"),
                panel.spacing = unit(0, "lines"))
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

###### Recombination sites ######

## Violin plot 
sw.data.8r %>%
  left_join(pt.data.8r,by=c("SeqID")) %>%
  na.omit %>%
  dplyr::filter(Chr=="ChrI",Pos>22,Int==1) %>%
  mutate(Progeny=ifelse(PA_Prog=="A","Absence","Presence")) %>%
  group_by(SeqID,Progeny) %>%
  summarise(BP=min(Pos)-0.1) %>%
  ggplot(aes(y=BP,x=Progeny))+geom_violin()+geom_jitter(width=0.3,size=1)+
  geom_hline(yintercept=35.8,col="black",linetype=3)+
  scale_y_continuous(breaks=c(24,28,32,36,40))+
  stat_summary(fun=mean,geom="crossbar",width=0.2,col="red")+
  coord_flip(ylim=c(24,40))+
  labs(y="ChrIR recombination sites (Mb)",x="Progeny of BC1 males")+
  theme(  text=element_text(size=9),
          panel.border=element_rect(linetype="solid",fill=NA),
          panel.spacing = unit(0, "lines"))

## GLM
glm.data<-sw.data.8r %>%
  left_join(pt.data.8r,by=c("SeqID")) %>%
  na.omit %>%
  dplyr::filter(Chr=="ChrI",Pos>22,Int==1) %>%
  mutate(Prog=ifelse(PA_Prog=="P",1,0)) %>%
  group_by(SeqID,Prog) %>%
  summarise(BP=min(Pos)-0.1)

glm1<-glm(glm.data$Prog~glm.data$BP,family="binomial")
#glm1 => intercept, -6.6384; BP, 0.1777
glm0<-glm(glm.data$Prog~1,family="binomial")
anova(glm1,glm0,test="Chisq")

## Logistic line

logis.data<-data.frame(BP=seq(24,40,length.out=1000))
logis.data$Y<-1/(1+exp(6.6384-0.1777*logis.data$BP))
logis.data$upper<-1/(1+exp(6.6384-(0.1777+0.07284)*logis.data$BP))
logis.data$lower<-1/(1+exp(6.6384-(0.1777-0.07284)*logis.data$BP))

logis.data %>%
  ggplot()+
  geom_line(aes(x=BP,y=Y))+
  geom_point(data=glm.data,aes(x=BP,y=Prog))+
  coord_cartesian(xlim=c(24,40),ylim=c(0,1))+
  scale_x_continuous(breaks=c(24,28,32,36,40),position="top")+
  geom_vline(xintercept=35.8,col="black",linetype=3)+
  labs(y="Proportion of fertile BC1")+
  theme(  text=element_text(size=9),
          panel.border=element_rect(linetype="solid",fill=NA),
          panel.spacing = unit(0, "lines"),
          plot.margin=margin(0,0,0,0),
          axis.title.x=element_blank())

###### Proportion of individuals w/ introgressions and Fisher's exact test ######

## Figure of proportion of introgression
g1<-sw.data %>%
  left_join(pt.data,by=c("SeqID")) %>%
  na.omit %>%
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
  na.omit %>%
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

#QTL peak
fisher.table %>%
  filter(logp>3.921132) %>%
  summarise(min(Pos),max(Pos))

#Graph of log p-value w/ significant level calculated bellow
g2<-fisher.table %>%
  ggplot()+geom_line(aes(x=Pos,y=logp))+
  geom_hline(yintercept=3.1188,linetype=2)+
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

pnum<-nrow(na.omit(pt.data[pt.data$PA_Prog=="P",]))
maxlogp<-c()
for(i in 1:1000){
  print(i)
  selection<-sample(na.omit(pt.data)$SeqID,pnum)
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

quantile(maxlogp,0.95) #3.118815
quantile(maxlogp,0.99) #3.921132 

###### R/QTL QTL ######

chrvec<-c("c1","c2","c3","c4","c5")
names(chrvec)<-c("ChrI","ChrII","ChrIII","ChrIV","ChrV")
gt.data<-sw.data %>%
  mutate(MK=paste(chrvec[as.character(Chr)],"_",Pos,sep=""),Genotype=as.character(ifelse(is.na(Int),"-",ifelse(Int==1,"H","A")))) %>%
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
write.csv(gt.data,file="Genotype_BC1_male_RQTL_1.csv",row.names=F,quote=F)

mapthis<-read.cross("csv",file="Genotype_BC1_male_RQTL_1.csv",estimate.map=F,crosstype="bc")
gddata <- est.map(mapthis, error.prob=0.001, verbose=FALSE)
mapdata<-replace.map(mapthis,gddata)
mapdata<-calc.errorlod(mapdata)

pull.map(mapdata,"I")

## Sigle QTL mapping
out.bin<-scanone(mapdata,model="binary")
plot(out.bin)

operm.bin<-scanone(mapdata,model="binary",n.perm=1000,perm.Xsp=TRUE)
summary(operm.bin, alpha=0.05) #2.61
summary(out.bin, perms=operm.bin, alpha=0.05, pvalues=TRUE)

mapdata2<-sim.geno(mapdata,n.draws=16,error.prob=0.001)
effectplot(mapdata,mname1="c1_31.55")

chrnames<-c("ChrI","ChrII","ChrIII","ChrIV","ChrV")
names(chrnames)<-c("c1","c2","c3","c4","c5")
res<-c()
for(i in 1:5){
  value<-as.vector(pull.map(mapdata2,i)[[1]])
  marker<-names(pull.map(mapdata2,i)[[1]])
  res.sub<-data.frame(matrix(unlist(strsplit(marker,"_")),byrow=T,ncol=2))
  res.sub<-cbind(res.sub,value)
  colnames(res.sub)<-c("Chr","Pos","cM")
  res.sub$Pos<-as.numeric(res.sub$Pos)
  res.sub$Chr<-chrnames[i]
  res<-rbind(res,res.sub)
}

res %>%
  ggplot(aes(x=as.numeric(Pos),y=cM))+
  geom_line()+
  facet_wrap(.~Chr,scales="free",ncol=5)+
  labs(x="Physical position(Mb)",y="Genetic position(cM)")+
  theme(text=element_text(size=9))+
  ggsave(file="F2_male_RR_plot.pdf",width=7,height=1.8)

## Fit to qtl model

mapdata2<-sim.geno(mapdata2, step=2, n.draws=128, err=0.001)
qtl<-makeqtl(mapdata2,chr=c(1),pos=c(74.7))
qtl
plot(qtl)
out.fq<-fitqtl(mapdata2,qtl=qtl,model="binary")
summary(out.fq)
# 13.18%

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

## Select 76:80 - F1 males of Pexs female vs. Ppac male
f1.data<-f1.data[as.numeric(as.character(f1.data$SeqID)) %in% 76:80,]
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
sample<-76:80
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

dp.data<-read.table(file="Dataset/Coverage_data_BC1_hybrid_male.txt",sep="\t",header=T)
dp.data<-dp.data[dp.data$Chr!="ChrX",]
dp.data$Pos<-(dp.data$End+dp.data$Start-1)/2000000
dp.data$SeqID<-factor(substr(dp.data$File,23,25))
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

violin.data %>%
  ggplot(aes(x=Chr_new,y=Intcov))+
  geom_violin()+geom_jitter(size=0.2,alpha=0.3)+
  labs(x="Chromosomes",y="Mean normalized extra coverage\nwith introgression (log2 scale)")+
  theme(  text=element_text(size=10),
          panel.border=element_rect(linetype="solid",fill=NA),
          strip.background = element_rect(colour = "black"),
          panel.spacing = unit(0, "lines"))

trisomy.data<-sw.data %>%
  left_join(dp.data,by=c("SeqID","Chr","Start","End","Pos")) %>%
  dplyr::filter(!(Chr=="ChrI"&Pos>20.05&Pos<22.05)) %>%
  na.omit %>%
  mutate(Intcov=Int*NorDepthPred) %>%
  mutate(Chr_new=factor(ifelse(Chr=="ChrI",ifelse(Pos>=22.05,"ChrIR","ChrIL"),as.character(Chr)),level=sel.chr.new)) %>%
  group_by(Chr_new,SeqID) %>%
  summarise(Intcov=sum(Intcov)/sum(Int))

other.chr<-c("ChrIR","ChrII","ChrIII","ChrIV","ChrV")
plist<-c()
countlist<-c()
for(j in 1:5){
  testc<-c(0,0,0,0)
  testc[1]<-sum(na.omit(trisomy.data$Intcov>0&trisomy.data$Chr_new=="ChrIL"))
  testc[2]<-(291-testc[1])
  testc[3]<-sum(na.omit(trisomy.data$Intcov>0&trisomy.data$Chr_new==other.chr[j]))
  testc[4]<-(291-testc[3])
  testsummary<-fisher.test(matrix(testc,nrow=2,byrow=T))
  countlist<-c(countlist,testc[3])
  plist<-c(plist,testsummary$p.value)
}
testc[1]
testc[1]/291
median(countlist)
median(countlist)/291
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
  dplyr::filter(Intcov>0,Chr_new=="ChrIL")

pt.data$PA_Prog<-factor(pt.data$PA_Prog,level=c("A","P"))
selected<-table(pt.data$PA_Prog[as.numeric(as.character(pt.data$SeqID)) %in% as.numeric(as.character(chril$SeqID))])
notselected<-table(pt.data$PA_Prog[!(as.numeric(as.character(pt.data$SeqID)) %in% as.numeric(as.character(chril$SeqID)))])
fisher.test(x=matrix(c(selected,notselected),nrow=2,byrow=T))
c(selected,notselected)
c(selected,notselected)/rep(c(sum(selected),sum(notselected)),each=2)

# Individuals not categorized as trisomy but supicious because of the recombination on the break point
# 7 individuals total in 225

BR_Zero<-sw.data %>%
  dplyr::filter(!(SeqID %in% chril$SeqID),Chr=="ChrI",Pos==22.05,Int==0)

additional_Trisomy<-sw.data %>%
  dplyr::filter(!(SeqID %in% chril$SeqID),Chr=="ChrI",Pos==20.05,Int==1) %>%
  semi_join(BR_Zero, by="SeqID")




###### Map data for Recombination analysis ######
chrvec<-c("c1","c2","c3","c4","c5")
names(chrvec)<-(c("ChrI","ChrII","ChrIII","ChrIV","ChrV"))
gt.data<-sw.data %>%
  filter(!(SeqID %in% chril$SeqID)) %>%
  filter(!(SeqID %in% additional_Trisomy$SeqID)) %>%
  mutate(MK=paste(chrvec[as.character(Chr)],"_",Pos,sep=""),Genotype=as.character(ifelse(is.na(Int),"-",ifelse(Int==1,"H","A")))) %>%
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
write.csv(gt.data,file="Genotype_BC1_male_RQTL_2.csv",row.names=F,quote=F)

mapthis<-read.cross("csv",file="Genotype_BC1_male_RQTL_2.csv",estimate.map=F,crosstype="bc")
gddata <-est.map(mapthis, error.prob=0.001, verbose=FALSE)
mapdata<-replace.map(mapthis,gddata)
save(mapdata,file="mapdata_BC1_hybrid_male.Robj")