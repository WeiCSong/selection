library(doParallel)
cl <- makeCluster(40)
registerDoParallel(cl)
memory.limit(999999)
library(data.table)
library(R.utils)
sds=fread("/lustre/home/acct-bmelgn/bmelgn-3/ldsc/SDS.txt",data.table=F)
sds_rho=c()
setwd("/lustre/home/acct-bmelgn/bmelgn-3/ldsc/all_gwas/")
write.csv(c("name","rho_true","p"),"sds_rho.csv")
tmp=list.files(pattern="*.gz")
for (name in tmp){
  gwas=fread(name,data.table=FALSE)
  gwas=gwas[which(gwas$sds!=0),]
  x=match(gwas$SNP,sds$ID)
  gwas$freq=sds[x,"DAF"]

######re-normalized per DAF bin
  gwas$SDS=0
  for(i in seq(0.05,0.94,0.01)){
    x=gwas[which(gwas$freq>i&gwas$freq<i+0.01),"sds"]
    x=(x-mean(x))/sd(x)
    gwas[which(gwas$freq>i&gwas$freq<i+0.01),"SDS"]=x
  }

########spearman rho 
  gwas=gwas[order(gwas$p),]
  res=c()
  for(i in seq(1,nrow(gwas),1000)){
    start=i
    end=min(i+1000,nrow(gwas))
    s=median(gwas[start:end,"SDS"])
    res=rbind(res,c(i,s))
  }
  rho_true=cor(res[,1],res[,2],method="spearman")
  l=order(gwas$CHR,gwas$BP)

##########block jackknife for rho
  bj=foreach(i=round(seq(1,nrow(gwas),length.out=1000)),.combine=rbind) %dopar% {
    L=l[-(i:(i+round(nrow(gwas)/1000)+1))]
    LL=L[order(L)]
    res=c()
    for(j in round(seq(1,length(L),length.out=1000))){
    start=j
    end=min(j+round(length(L)/1000)+1,length(L))

    s=median(gwas[LL[start:end],"SDS"])
    res=rbind(res,c(j,s))
    }
  rho=cor(res[,1],res[,2],method="spearman")
  return(c(i,rho))
  }
  rm(L)
  rm(LL)
  rm(gwas)
  gc()

################p value
  z=rho_true/sd(bj[,2])
  p=2*pnorm(-abs(z))
  line=c(name,rho_true,p)
  write(line,file="sds_rho.csv",append=TRUE)
}

####################analysis
head(res_sds)
res_sds=res_sds[which(res_sds[,1] %in% annot[,1]),]
length(which(res_sds[,3]<0.05/875))
res_sds$ukb=annot[rownames(res_sds),"ukb"]
indraw=data.frame(r=abs(res_sds[,2]),ukb=factor(res_sds$ukb))
indraw[which(is.na(indraw$ukb)),"ukb"]=0
oneway_test(r~ukb,data=indraw)

ggplot(indraw,aes(ukb,r))+geom_boxplot()

draw=data.frame(ID=res_sds[,1],rho=abs(res_sds[,2]),group=annot[match(res_sds[,1],annot[,1]),"abr"])
draw=draw[which(!is.na(draw$group)),]
draw=draw[which(!is.na(draw$rho)),]
#dis=annot[which(annot$disease==1),1]
#draw=draw[which(draw$ID %in% dis),]
f=function(x){
  c(x,median(draw[which(draw$group==x),2]))
}
d=data.frame(lapply(unique(draw$group),f))
d=data.frame(t(d))

draw$group=factor(draw$group,levels=d[order(d[,2],decreasing=T),1])
p=ggplot(draw,aes(group,rho,group=group,color=group))
p=p+geom_boxplot(outlier.size =0)+geom_jitter(size=0.5)
p=p+theme_classic()+theme(legend.position="none",axis.ticks=element_blank(),
axis.title.x=element_blank(),axis.text=element_text(color="black",size=10,
family="sans"),axis.title.y=element_text(color="black",size=10,family="sans"))
p=p+ylab("|Spearman rho|")+scale_y_continuous(breaks=c(0,0.5,1))

p_sds_all=p
#disease only
dis=annot[which(annot$disease==1),1]
draw=draw[which(draw$ID %in% dis),]
f=function(x){
  c(x,median(draw[which(draw$group==x),2]))
}
d=data.frame(lapply(unique(draw$group),f))
d=data.frame(t(d))

draw$group=factor(draw$group,levels=d[order(d[,2],decreasing=T),1])
p=ggplot(draw,aes(group,rho,group=group,color=group))
p=p+geom_boxplot(outlier.shape=NA)+geom_jitter()
p=p+theme_classic()+theme(legend.position="none",axis.ticks=element_blank(),
axis.title.x=element_blank(),axis.text=element_text(color="black",size=10,
family="sans"),axis.title.y=element_text(color="black",size=10,family="sans"))
p=p+ylab("|Spearman rho|")
p=p+scale_color_manual(breaks=grouporder,values=gg_color_hue(15))

p_sds_dis=p

p_sds_all/p_sds_dis
#
res_sds$disease=annot[match(res_sds[,1],annot[,1]),"disease"]
res_sds$onset=annot[match(res_sds[,1],annot[,1]),"onset"]
head(res_sds)
res_sds$abs=abs(res_sds[,2])
oneway.test(abs~disease,data=res_sds)
summary(res_sds[which(res_sds$dis==1),])
res_sds$cat=annot[match(res_sds[,1],annot[,1]),"abr"]
write.csv(res_sds,"res_sds.csv")


res_sds[which(res_sds[,1] %in% anno[which(anno[,3]=="reproduction"),1]),]

summary(abs(res_sds[,2]))
dis=res_sds[which(res_sds[,1] %in% anno[which(anno[,2]==1),1]),]
summary(dis[,2])
dis[,1]=anno[match(dis[,1],anno[,1]),4]

rownames(res_sds)=res_sds[,1]
der=annot[which(annot$abr=="DER"),1]
res_sds$z=unlist(lapply(res_sds[,3],qnorm))
res_sds$sd=abs(res_sds[,2]/res_sds[,4])
res_sds[der,]
res_sds$trait=annot[rownames(res_sds),4]
res_sds[which(res_sds[,3]>0.05/890),]

=density(na.omit(res_sds$abs))
plot(sdsden$x,sdsden$y)


rownames(res_sds)
res_sds$sig=""
res_sds[which(abs(res_sds$z)<4|abs(res_sds[,2])<0.1),"sig"]="insig"
res_sds[which(abs(res_sds$z)>4&res_sds[,2]>0.1),"sig"]="negative"
res_sds[which(abs(res_sds$z)>4&res_sds[,2]<(-0.1)),"sig"]="positive"
res_sds[,c(2,4)]=-res_sds[,c(2,4)]


derd=derd[,c(1,2,5)]
colnames(derd)=c("ID","rho","sd")
derd$low=derd$rho-1.96*derd$sd
derd$up=derd$rho+1.96*derd$sd
derd[,c("rho","low","up")]=apply(derd[,c("rho","low","up")],2,scaling)
derd=derd[order(derd$rho,decreasing=T),]
derd$ID=factor(derd$ID,levels=derd$ID[14:1])
derd$y=as.numeric(derd$ID)

derd2=c()
for(i in rownames(derd)){
int=data.frame(x=c(derd[i,"up"],derd[i,"rho"],derd[i,"low"],derd[i,"rho"]),
y=c(derd[i,"y"],derd[i,"y"]-0.5,derd[i,"y"],derd[i,"y"]+0.5),ID=i)
derd2=rbind(derd2,int)
}
derd2$text=annot[match(derd2$ID,annot[,1]),4]

p=ggplot(derd2,aes(x,y),fill="#F8766D")
p=p+geom_polygon(aes(x,y,group=ID),fill="#F8766D")+theme_classic()+geom_vline(xintercept=0)
p=p+theme(axis.ticks=element_blank(),axis.line.y=element_blank(),axis.title.y=
element_blank(),axis.text.y=element_blank(),axis.text.x=element_text(size=10,
color="black",family="sans"),axis.title.x=element_text(size=10,
color="black",family="sans"))+xlab("Spearman rho")
p=p+ annotate(geom="text", x=-0.03, hjust=1,y=1, label="Skin tanning",size=4)
p=p+ annotate(geom="text", x=0.03,hjust=0, y=14, label="Dark brown hair",size=4)
p=p+ annotate(geom="text", x=0.03,hjust=0, y=13, label="Black hair",size=4)
p=p+ annotate(geom="text", x=0.03,hjust=0, y=12, label="Skin color",size=4)
p=p+ annotate(geom="text", x=0.03,hjust=0, y=11, label="Bald Pattern 4",size=4)
p=p+ annotate(geom="text", x=0.03,hjust=0, y=10, label="Bald Pattern 3",size=4)
p=p+ annotate(geom="text", x=0.03,hjust=0, y=9, label="Baldness",size=4)
p=p+ annotate(geom="text", x=0.03,hjust=0, y=8, label="Vitiligo",size=4)
#negative
p=p+ annotate(geom="text", x=-0.03, hjust=1,y=7, label="Bald Pattern 1",size=4)
p=p+ annotate(geom="text", x=-0.03, hjust=1,y=6, label="Basal cell carcinoma",size=4)
p=p+ annotate(geom="text", x=-0.03, hjust=1,y=5, label="Bald Pattern 2",size=4)
p=p+ annotate(geom="text", x=-0.03, hjust=1,y=4, label="Age first facial hair",size=4)
p=p+ annotate(geom="text", x=-0.03, hjust=1,y=3, label="Blonde hair",size=4)
p=p+ annotate(geom="text", x=-0.03, hjust=1,y=2, label="Brown hair",size=4)

p_sds_forest=p

pdf("test.pdf",height=2,width=2)
p
dev.off()

########spearman rho 
gwas=gwas[order(gwas$p),]
res=c()
for(i in seq(1,nrow(gwas),1000)){
  start=i
  end=min(i+1000,nrow(gwas))
  s=median(gwas[start:end,"SDS"])
  res=rbind(res,c(i,s))
}
res[,1]=4149:1
head(res)
colnames(res)=c("bin","SDS")
p=ggplot(data.frame(res),aes(bin,SDS))+geom_point(color="#F8766D",size=0.3)+theme_classic()
p=p+xlab("p value bin (high to low)")+ylab("bin median SDS")
p=p+theme(axis.ticks=element_blank(),axis.text=element_text(color="black",
size=10,family="sans"),axis.title=element_text(size=10,color="black",family="sans"))
p=p+annotate(geom="text", x=1500, y=0.8, 
label="Skin tanning:\nSpearman rho=0.96\np by block jackknife<e-200",size=4)
p_sds_scatter=p
p_sds_scatter$data[,2]=-p_sds_scatter$data[,2]

#xx=c("IMP","DER","REP","COG","self care","PSY",
"NUT","MET","GI","CIRC","immune","MED","MUSC","NEU","RES")
#p_cf_bar=p_cf_bar+scale_x_discrete(breaks=xx,labels=c("body",xx[2:15]))
#xx=c("DER","NUT","REP","body","GI","self care","PSY","RES","MED","COG",
"MUSC","immune","MET","CIRC","NEU")
#p_sds_all=p_sds_all+scale_x_discrete(labels=xx)

p_f2=((p_cm_bar/p_cf_bar)|p_mr_scatter)/p_sds_all/(p_sds_forest|p_sds_scatter)+
plot_layout(height=c(5,3,4))+ plot_annotation(tag_levels = 'A')

pdf("f2_new.pdf",height=10,width=8)
p_f2
dev.off()


