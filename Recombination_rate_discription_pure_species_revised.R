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

res2<-loessfit_100kbdf(res2$Chr,res2$Pos,res2$RR,span=0.2)

RRdata.Pexs<-res2

RRdata2<-RRdata.Pexs[RRdata.Pexs$Chr %in% c("ChrX"),]
RRdata3<-RRdata.Pexs[RRdata.Pexs$Chr %in% c("ChrI"),]
RRdata3$Pos<-24.218945-RRdata3$Pos

###### P. pacificus data ######

load(file="Dataset/mapdata_Ppac_RILs.Robj")

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

res2<-loessfit_100kbdf(res2$Chr,res2$Pos,res2$RR,span=0.2)
res2<-na.omit(res2)

RRdata.Ppac<-res2
RRdata4<-RRdata.Ppac[RRdata.Ppac$Chr=="ChrI",]
RRdata4$Pos<-39.556110-RRdata4$Pos
RRdata5<-RRdata.Ppac[RRdata.Ppac$Chr=="ChrX",]

###### Figure of recombination rate of all chromosomes ######

## P. exspectatus
# Rectangular data showing the position of ChrIL, ChrIR and ChrX
pexs_rectdata<-data.frame(xmin=c(0,0,20),xmax=c(24.218945,20,31.480924),
                          ymin=c(23,23,23),ymax=c(25,25,25),col=c("A","B","C"),
                          Chr=c("ChrI*","ChrX*","ChrX*"))

RRdata.Pexs.c<-RRdata.Pexs %>%
  mutate(Chr=sapply(Chr,paste,"*",sep=""))

sp1<-ggplot()+
  geom_hline(yintercept=mean(RRdata.Pexs$fit),col="gold")+
  geom_ribbon(data=RRdata.Pexs.c,aes(x=as.numeric(Pos),ymin=fit-se,ymax=fit+se),alpha=0.2)+
  geom_line(data=RRdata.Pexs.c,aes(x=as.numeric(Pos),y=fit))+
  geom_rect(data=pexs_rectdata,
            mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=col))+
  scale_fill_manual(values=c("darkorange","purple","limegreen"))+
  facet_grid(.~Chr,scales="free",space="free")+
  labs(x="Physical position(Mb)",y="Recombination rate (cM/Mb)",title="P. exspectatus")+
  linetheme+
  theme(text=element_text(size=9),legend.position="none")

## P. pacificus

ppac_rectdata<-data.frame(xmin=c(0,21.05,0),xmax=c(21.05,39.556110,17.019893),
           ymin=c(12,12,12),ymax=c(13,13,13),col=c("A","B","C"),
           Chr=c("ChrI","ChrI","ChrX"))

sp2<-ggplot()+
  geom_hline(yintercept=mean(RRdata.Ppac$fit),col="gold")+
  geom_ribbon(data=RRdata.Ppac,mapping=aes(x=as.numeric(Pos),ymin=fit-se,ymax=fit+se),alpha=0.2)+
  geom_line(data=RRdata.Ppac,mapping=aes(x=as.numeric(Pos),y=fit))+
  geom_rect(data=ppac_rectdata,
            mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=col))+
  scale_fill_manual(values=c("darkorange","purple","limegreen"))+
  coord_cartesian(ylim=c(-1,13))+
  facet_grid(.~Chr,scales="free",space="free")+
  labs(x="Physical position(Mb)",y="Recombination rate (cM/Mb)",title="P. pacificus")+
  linetheme+
  theme(text=element_text(size=9),legend.position="none")

pdf("Fig_Recombination_rate_all_chromosomes_revised.pdf",width=6,height=4)
sp1/sp2+plot_layout(ncol=1)
dev.off()


###### Figure of recombination rate of three elements ######

## P. exspectatus

ptitle=expression(paste(bolditalic("P. exspectatus")))
ep1<-RRdata2 %>%
  ggplot()+
  geom_hline(yintercept=mean(RRdata.Pexs$fit),col="gold")+
  geom_ribbon(aes(x=Pos,ymin=fit-se,ymax=fit+se),alpha=0.2)+
  geom_line(aes(x=Pos,y=fit))+
  annotate("text",15,18.5,label="ChrX*",size=chrtextsize)+
  coord_cartesian(xlim=c(0,30),ylim=c(-2,20))+
  labs(title=ptitle,x="Physical position, Mb",
       y="Recombination\nrate, cM/Mb")+
  linetheme+
  theme(plot.title=element_text(vjust=-2.5),
        axis.title.y=element_text(vjust=-8))

ep2<-RRdata3 %>%
  ggplot()+
  geom_hline(yintercept=mean(RRdata.Pexs$fit),col="gold")+
  geom_ribbon(aes(x=Pos,ymin=fit-se,ymax=fit+se),alpha=0.2)+
  geom_line(aes(x=Pos,y=fit))+
  annotate("text",10,18.5,label="inverted ChrI*",size=chrtextsize)+
  coord_cartesian(xlim=c(0,20),ylim=c(-2,20))+
  labs(title=NULL,
       x="Physical position, Mb",y=NULL)+
  linetheme+
  theme(plot.margin=unit(c(0,4,2,2),"pt"),
        axis.text.y=element_blank())


## P. pacificus

ptitle=expression(paste(bolditalic("P. pacificus")))
ep5<-RRdata4 %>%
  ggplot(aes(x=as.numeric(Pos)))+
  geom_hline(yintercept=mean(RRdata.Ppac$fit),col="gold")+
  geom_ribbon(aes(ymin=fit-se,ymax=fit+se),alpha=0.2)+
  geom_line(aes(y=fit))+
  annotate("text",20,10.5,label="inverted ChrI",size=chrtextsize)+
  coord_cartesian(xlim=c(0,40),ylim=c(-0.5,11.5))+
  labs(title=ptitle,
       x="Physical position, Mb",y="Recombination\nrate, cM/Mb")+
  linetheme+
  theme(plot.margin=unit(c(1,0,2,0),"pt"),
        plot.title=element_text(vjust=-2.5),
        axis.title.y=element_text(vjust=-8))

ep6<-RRdata5 %>%
  ggplot(aes(x=as.numeric(Pos)))+
  geom_hline(yintercept=mean(RRdata.Ppac$fit),col="gold")+
  geom_ribbon(aes(ymin=fit-se,ymax=fit+se),alpha=0.2)+
  geom_line(aes(y=fit))+
  annotate("text",8,10.5,label="ChrX",size=chrtextsize)+
  coord_cartesian(xlim=c(0,16),ylim=c(-0.5,11.5))+
  labs(title=NULL,
       x="Physical position, Mb",y=NULL)+
  linetheme+
  theme(plot.margin=unit(c(0,4,2,2),"pt"),
        axis.text.y=element_blank())

###### Color labels ######

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

ep3<-tdata1 %>% ggplot()+geom_tile(aes(x=x,y=y,col=col,fill=col))+
  coord_cartesian(xlim=c(0,30))+
  scale_y_discrete(expand=c(0,0))+
  annotate("text",x=10,y=1,col="white",label="inverted ChrIR",size=tiletextsize)+
  annotate("text",x=25,y=1,col="white",label="ChrX",size=tiletextsize)+
  scale_color_manual(values=c("purple","limegreen"))+
  scale_fill_manual(values=c("purple","limegreen"))+
  labs(x="")+
  theme_void()+
  tiletheme
  
ep4<-tdata2 %>% ggplot()+geom_tile(aes(x=x,y=y,col=col,fill=col))+
  coord_cartesian(xlim=c(0,20))+
  scale_y_discrete(expand=c(0,0))+
  annotate("text",x=10,y=1,col="white",label="inverted ChrIL",size=tiletextsize)+
  scale_color_manual(values=c("darkorange"))+
  scale_fill_manual(values=c("darkorange"))+
  labs(x="")+
  theme_void()+
  tiletheme
ep7<-tdata3 %>% ggplot()+geom_tile(aes(x=x,y=y,col=col,fill=col))+
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

ep8<-tdata4 %>% ggplot()+geom_tile(aes(x=x,y=y,col=col,fill=col))+
  coord_cartesian(xlim=c(0,16))+
  scale_y_discrete(expand=c(0,0))+
  annotate("text",x=7.6,y=1,col="white",label="ChrX",size=tiletextsize)+
  scale_color_manual(values=c("limegreen"))+
  scale_fill_manual(values=c("limegreen"))+
  labs(x="")+
  theme_void()+
  tiletheme

###### Combined figure ######

{ep1+ep2+ep3+ep4+plot_layout(ncol=2,widths=c(3.5,2.33),heights=c(8.5,1))}/
{ep5+ep6+ep7+ep8+plot_layout(ncol=2,widths=c(4.2,1.68),heights=c(8.5,1))}

