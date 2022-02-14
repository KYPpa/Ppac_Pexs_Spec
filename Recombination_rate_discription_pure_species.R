library(dplyr)
library(tidyr)
library(ggplot2)
library(qtl)
library(patchwork)

## R object, mapdata, is loaded.
## Those were produced in the scripts,
## Linkage_analysis_PpacificusRILs.R
## Linkage_analysis_Pexspectatus.R

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


###### Figure theme ######

linetheme<-theme(text=element_text(size=9),
                 plot.title=element_text(size=8,face="bold.italic"),
                 panel.background=element_rect(fill=NA),
                 panel.border=element_rect(linetype="solid",colour="grey40",fill=NA),
                 plot.margin=unit(c(0,0,0,0),"pt"),
                 axis.title.x=element_text(size=6),
                 axis.title.y=element_text(size=7))
tiletheme<-theme(legend.position="none",
                 axis.title.x=element_text(size=2,margin=margin(0,0,0,0,"pt")),
                 plot.margin=unit(c(0,0,0,0),"pt"))
chrtextsize<-2.5
tiletextsize<-2
###### P. exspectatus data ######

chrlist<-c("ChrI","ChrII","ChrIII","ChrIV","ChrV","ChrX")
chrlist2<-c("I","II","III","IV","V","X")

load("Dataset/mapdata_Pexspectaus_EMS_hybrids.Robj")

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

sp1<-RRdata.Pexs %>%
  mutate(Chr=sapply(Chr,paste,"*",sep="")) %>%
  ggplot(aes(x=as.numeric(Pos)))+geom_line(aes(y=fit))+
  geom_hline(yintercept=mean(RRdata.Pexs$fit),col="gold")+
  geom_ribbon(aes(ymin=fit-se,ymax=fit+se),alpha=0.2)+
  geom_line(aes(y=fit))+
  facet_grid(.~Chr,scales="free",space="free")+
  labs(x="Physical position(Mb)",y="Recombination rate (cM/Mb)",title="P. exspectatus")+
  linetheme+
  theme(text=element_text(size=9))

RRdata2<-RRdata.Pexs[RRdata.Pexs$Chr %in% c("ChrX"),]
RRdata3<-RRdata.Pexs[RRdata.Pexs$Chr %in% c("ChrI"),]
RRdata3$Pos<-24.218945-RRdata3$Pos

###### P. pacificus data ######

load(file="mapdata_Ppac_RILs.Robj")

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

RRdata.Ppac<-res2
sp2<-RRdata.Ppac %>%
  ggplot(aes(x=as.numeric(Pos)))+geom_line(aes(y=fit))+
  geom_hline(yintercept=mean(RRdata.Ppac$fit),col="gold")+
  geom_ribbon(aes(ymin=fit-se,ymax=fit+se),alpha=0.2)+
  geom_line(aes(y=fit))+
  coord_cartesian(ylim=c(-1,12))+
  facet_grid(.~Chr,scales="free",space="free")+
  labs(x="Physical position(Mb)",y="Recombination rate (cM/Mb)",title="P. pacificus")+
  linetheme+
  theme(text=element_text(size=9))

RRdata4<-RRdata.Ppac[RRdata.Ppac$Chr=="ChrI",]
RRdata4$Pos<-39.556110-RRdata4$Pos
RRdata5<-RRdata.Ppac[RRdata.Ppac$Chr=="ChrX",]

###### Figure of recombination rate of all chromosomes ######

sp1/sp2+plot_layout(ncol=1)

###### Figure of recombination rate of three elements ######

## P. exspectatus

ptitle=expression(paste(bolditalic("P. exspectatus")))
p1<-RRdata2 %>%
  ggplot()+
  geom_hline(yintercept=mean(RRdata.Pexs$fit),col="gold")+
  geom_ribbon(aes(x=Pos,ymin=fit-se,ymax=fit+se),alpha=0.2)+
  geom_line(aes(x=Pos,y=fit))+
  annotate("text",15,12,label="ChrX*",size=chrtextsize)+
  coord_cartesian(xlim=c(0,30))+
  labs(title=ptitle,x="Physical position, Mb",
       y="Recombination\nrate, cM/Mb")+
  linetheme

p2<-RRdata3 %>%
  ggplot()+
  geom_hline(yintercept=mean(RRdata.Pexs$fit),col="gold")+
  geom_ribbon(aes(x=Pos,ymin=fit-se,ymax=fit+se),alpha=0.2)+
  geom_line(aes(x=Pos,y=fit))+
  annotate("text",10,18.5,label="inverted ChrI*",size=chrtextsize)+
  coord_cartesian(xlim=c(0,20))+
  labs(title=NULL,
       x="Physical position, Mb",y=NULL)+
  linetheme

## P. pacificus

ptitle=expression(paste(bolditalic("P. pacificus")))
p5<-RRdata4 %>%
  ggplot(aes(x=as.numeric(Pos)))+
  geom_hline(yintercept=mean(RRdata.Ppac$fit),col="gold")+
  geom_ribbon(aes(ymin=fit-se,ymax=fit+se),alpha=0.2)+
  geom_line(aes(y=fit))+
  annotate("text",20,3.4,label="inverted ChrI",size=chrtextsize)+
  coord_cartesian(xlim=c(0,40))+
  labs(title=ptitle,
       x="Physical position, Mb",y="Recombination\nrate, cM/Mb")+
  linetheme+
  theme(plot.margin=unit(c(0,2,2,0),"pt"))

p6<-RRdata5 %>%
  ggplot(aes(x=as.numeric(Pos)))+
  geom_hline(yintercept=mean(RRdata.Ppac$fit),col="gold")+
  geom_ribbon(aes(ymin=fit-se,ymax=fit+se),alpha=0.2)+
  geom_line(aes(y=fit))+
  annotate("text",8,10.5,label="ChrX",size=chrtextsize)+
  coord_cartesian(xlim=c(0,16))+
  labs(title=NULL,
       x="Physical position, Mb",y=NULL)+
  linetheme+
  theme(plot.margin=unit(c(0,4,2,0),"pt"))

## Color labels

tdata1<-data.frame(x=seq(0,30,length.out=301),
                   y=rep(1,301),
                   col=c(rep("a",200),rep("b",101)))
tdata2<-data.frame(x=seq(0,20,length.out=201),
                   y=rep(1,201),
                   col=rep("a",201))
tdata3<-data.frame(x=seq(0,40,length.out=401),
                   y=rep(1,401),
                   col=c(rep("a",191),rep("b",210)))
tdata4<-data.frame(x=seq(0,16,length.out=161),
                   y=rep(1,161),
                   col=rep("a",161))

p3<-tdata1 %>% ggplot()+geom_tile(aes(x=x,y=y,col=col,fill=col))+
  coord_cartesian(xlim=c(0,30))+
  scale_y_discrete(expand=c(0,0))+
  annotate("text",x=10,y=1,col="white",label="inverted ChrIR",size=tiletextsize)+
  annotate("text",x=25,y=1,col="white",label="ChrX",size=tiletextsize)+
  scale_color_manual(values=c("purple","limegreen"))+
  scale_fill_manual(values=c("purple","limegreen"))+
  labs(x="")+
  theme_void()+
  tiletheme
  
p4<-tdata2 %>% ggplot()+geom_tile(aes(x=x,y=y,col=col,fill=col))+
  coord_cartesian(xlim=c(0,20))+
  scale_y_discrete(expand=c(0,0))+
  annotate("text",x=10,y=1,col="white",label="inverted ChrIL",size=tiletextsize)+
  scale_color_manual(values=c("darkorange"))+
  scale_fill_manual(values=c("darkorange"))+
  labs(x="")+
  theme_void()+
  tiletheme
p7<-tdata3 %>% ggplot()+geom_tile(aes(x=x,y=y,col=col,fill=col))+
  coord_cartesian(xlim=c(0,40))+
  scale_y_discrete(expand=c(0,0))+
  annotate("text",x=10,y=1,col="white",label="inverted ChrIR",size=tiletextsize)+
  annotate("text",x=30,y=1,col="white",label="inverted ChrIL",size=tiletextsize)+
  scale_color_manual(values=c("purple","darkorange"))+
  scale_fill_manual(values=c("purple","darkorange"))+
  labs(x="")+
  theme_void()+
  tiletheme+
  theme(axis.title.x=element_text(margin=margin(2,0,5,0,"pt")))

p8<-tdata4 %>% ggplot()+geom_tile(aes(x=x,y=y,col=col,fill=col))+
  coord_cartesian(xlim=c(0,16))+
  scale_y_discrete(expand=c(0,0))+
  annotate("text",x=7.6,y=1,col="white",label="ChrX",size=tiletextsize)+
  scale_color_manual(values=c("limegreen"))+
  scale_fill_manual(values=c("limegreen"))+
  labs(x="")+
  theme_void()+
  tiletheme

## Combined figure

{p1+p2+p3+p4+plot_layout(ncol=2,widths=c(3.5,2.2),heights=c(8.5,1))}/
{p5+p6+p7+p8+plot_layout(ncol=2,widths=c(4.2,1.5),heights=c(8.5,1))}

