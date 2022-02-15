library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(qtl)
library(patchwork)
library(grid)

###### Data preparation ######

# Pexs

load(file="Dataset/mapdata_Pexspectaus_EMS_hybrids.Robj")

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

res.pexs<-res %>%
  mutate(Chr=factor(Chr))

# Ppac

load(file="Dataset/mapdata_Ppac_RILs.Robj")
chrlist<-c("ChrI","ChrII","ChrIII","ChrIV","ChrV","ChrX")
rename_chrlist<-c("c1","c2","c3","c4","c5","c6")
levels(sw.data$Chr)<-rename_chrlist

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

res.ppac<-res %>%
  mutate(Chr=factor(Chr))

# BC1 male - Pexs vs. Ppac

load(file="Dataset/mapdata_BC1_hybrid_male.Robj")
chrvec<-c("c1","c2","c3","c4","c5")
names(chrvec)<-(c("ChrI","ChrII","ChrIII","ChrIV","ChrV"))

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

res.bc1male<-res %>%
  mutate(Chr=factor(Chr))

# BC1 hermaphrodite - Ppac vs. Pexs

load(file="Dataset/mapdata_BC1_hybrid_hermaphrodite.Robj")
chrvec<-c("c1","c2","c3","c4","c5")
names(chrvec)<-(c("ChrI","ChrII","ChrIII","ChrIV","ChrV"))

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

res.bc1hermaphrodite<-res %>%
  mutate(Chr=factor(Chr))

###### Theme ######

commontheme<-theme(text=element_text(size=9),
                   legend.key.height=unit(3,"pt"),
                   legend.justification="left",
                   legend.margin=margin(0,0,0,0,unit="in"),
                   legend.box.spacing=unit(3, "pt"),
                   plot.margin=unit(c(0,0,3,5),"pt"),
                   axis.title.y=element_text(size=7),
                   panel.background=element_rect(fill=NA),
                   panel.border=element_rect(linetype="solid",colour="grey40",fill=NA),
                   panel.spacing = unit(0.2, "lines"),
                   plot.title=element_text(size=8),
                   strip.text=element_text(margin=margin(1,0,1,0,unit="pt")))
res.pexs$Chr<-factor(res.pexs$Chr,levels=c("ChrI","ChrII","ChrIII","ChrIV","ChrV","ChrX"))
levels(res.pexs$Chr)<-c("ChrI*","ChrII*","ChrIII*","ChrIV*","ChrV*","ChrX*")

###### Generate figure objects ######

ptitle<-expression(bolditalic("P. exspectatus"))
g1<-res.pexs %>%
  ggplot()+geom_line(aes(x=as.numeric(Pos),y=cM),col="black")+
  facet_wrap(.~Chr,ncol=6,scale="free")+
  labs(x=NULL,y="Genetic distance, cM",title=ptitle)+
  commontheme

ptitle<-expression(bolditalic("P. pacificus"))
g2<-res.ppac %>%
  ggplot()+geom_line(aes(x=as.numeric(Pos),y=cM),col="black")+
  facet_wrap(.~Chr,ncol=6,scale="free")+
  labs(x=NULL,y="Genetic distance, cM",title=ptitle)+
  commontheme

ptitle<-expression(paste(bold("F1["),bolditalic("P. exspectatus"),bold(" female ×"),
                         bolditalic("P. pacificus"),bold(" male]")))

g3<-res.bc1male %>%
  ggplot()+geom_line(aes(x=as.numeric(Pos),y=cM),col="black")+
  facet_wrap(.~Chr,ncol=6,scale="free")+
  labs(x=NULL,y="Genetic distance, cM",title=ptitle)+
  commontheme

ptitle<-expression(paste(bold("F1["),bolditalic("P. pacificus"),bold(" hermaphrodite ×"),
                         bolditalic("P. exspectatus"),bold(" male]")))

g4<-res.bc1hermaphrodite %>%
  ggplot()+geom_line(aes(x=as.numeric(Pos),y=cM),col="black")+
  facet_wrap(.~Chr,ncol=6,scale="free")+
  labs(x="Physical distance, Mb",y="Genetic distance, cM",title=ptitle)+
  commontheme

###### Draw the figure ######

pdf(file="Fig_Marey_map.pdf",width=5.4,height=7)
g1+g2+{g3+plot_spacer()+g4+plot_spacer()+plot_layout(ncol=2,widths=(c(5.1,1)))}+plot_layout(ncol=1,heights=c(1,1,2.7))
grid.text("A",x=unit(0.03,"npc"),y=unit(0.976,"npc"),gp=gpar(cex=1.7))
grid.text("B",x=unit(0.03,"npc"),y=unit(0.750,"npc"),gp=gpar(cex=1.7))
grid.text("C",x=unit(0.03,"npc"),y=unit(0.510,"npc"),gp=gpar(cex=1.7))
grid.text("D",x=unit(0.03,"npc"),y=unit(0.260,"npc"),gp=gpar(cex=1.7))
dev.off()
