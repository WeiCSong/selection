tmp=list.files(pattern="*txt")
#tmp=gsub("AA","",tmp)

i="AAID183.txt"
Traj=c()
TxP=c()
for(i in tmp){
if(file.info(i)[1,1]<1000){
traj=rep(0,31)
p_tx=1
}else{
dat=read.table(i,head=T)
l=dat[,1]
dat=dat[,-1]
dat=apply(dat,2,as.numeric)
dat[,3:33]=dat[,3:33]/2
rownames(dat)=l
traj=dat[,1] %*% dat[,3:33]

#p_tx=tryCatch(Qx_test(t(dat[,26:30]), dat[,1],perms =10000)[3],error=function(e){1})
}
Traj=rbind(Traj,traj)
#TxP=c(TxP,p_tx)
}
dim(Traj)
tmp=gsub("AA","",tmp)
tmp=gsub(".txt","",tmp)
rownames(Traj)=tmp
names(TxP)=tmp

#
pvec=dat[which(dat[,2]<0&dat[,3]!=1),3]

#


#
install.packages("dtwclust")
library(tidyr)
library(dtwclust)
library(dplyr)
library(ggplot2)
library(reshape)
# wide to long
Traj[which(is.na(Traj))]=0
which(rowSums(Traj)!=0)
df_long <- gather(data.frame(t(Traj[which(rowSums(Traj)!=0),20:31])))
head(df_long)
df_long[which(df_long[,2]==0),2]=0.001
df_long[which(is.na(df_long[,2])),2]=0.001
# make a timepoint column
df_long$time <- rep(1:12,765)
df_list <- as.list(utils::unstack(df_long, value ~ key))

df_list_z <- dtwclust::zscore(df_list)
pc_k=tsclust(df_list_z, type = "h", k = 4L,  distance = "dtw", 
control = hierarchical_control(method = "complete"), seed = 390, 
preproc = NULL, args = tsclust_args(dist = list(window.size = 5L)))

sapply(pc_k, cvi, type = "internal")

p_f5cn=plot(pc_k, type = "sc")
p_f5cn=p_f5cn+ylab("population mean PRS")+xlab("generation to present")+theme(
plot.title =element_blank(),strip.background = element_blank(),
axis.ticks=element_blank(),strip.text= element_blank())+
scale_x_continuous(breaks=c(3,8),labels=c(500,100))+geom_vline(xintercept=
c(3,8))

plot(pc_k, type = "centroids")
plot(pc_k)
summary(pc_k)
pc_k$order
clures=pc_k@cluster
lab=names(clures)[which(clures==4)]
annot[which(rownames(annot) %in% lab),4:5]

hclust(pc_k)

plot(1:31,Traj["ID47",1:31])
d_s12d=data.frame(gen=gen,prs=t(Traj["ID47",1:31]))
colnames(d_s12d)[2]="prs"

p=ggplot(d_s12d,aes(gen,prs))+geom_point(color=gg_color_hue(15)[15])
p=p+theme_classic()+ scale_x_continuous(trans=reverselog_trans(10))
p=p+geom_line(color=gg_color_hue(15)[15])
p_s12d=p+xlab("generation to present")

d_s12c=data.frame(gen=gen,prs=t(Traj["ID313",1:31]))
colnames(d_s12c)[2]="prs"

p=ggplot(d_s12c,aes(gen,prs))+geom_point(color=gg_color_hue(15)[10])
p=p+theme_classic()+ scale_x_continuous(trans=reverselog_trans(10))
p=p+geom_line(color=gg_color_hue(15)[10])
p_s12c=p+xlab("generation to present")

d_s12c=data.frame(gen=gen,prs=t(Traj["ID313",1:31]))
colnames(d_s12c)[2]="prs"

p=ggplot(d_s12c,aes(gen,prs))+geom_point(color=gg_color_hue(15)[10])
p=p+theme_classic()+ scale_x_continuous(trans=reverselog_trans(10))
p=p+geom_line(color=gg_color_hue(15)[10])
p_s12c=p+xlab("generation to present")

d_s12c=data.frame(gen=gen,prs=t(Traj["ID313",1:31]))
colnames(d_s12c)[2]="prs"

p=ggplot(d_s12c,aes(gen,prs))+geom_point(color=gg_color_hue(15)[10])
p=p+theme_classic()+ scale_x_continuous(trans=reverselog_trans(10))
p=p+geom_line(color=gg_color_hue(15)[10])
p_s12c=p+xlab("generation to present")

d_s12b=data.frame(gen=gen,prs=t(Traj["out35",1:31]))
colnames(d_s12b)[2]="prs"

p=ggplot(d_s12b,aes(gen,prs))+geom_point(color=gg_color_hue(15)[12])
p=p+theme_classic()+ scale_x_continuous(trans=reverselog_trans(10))
p=p+geom_line(color=gg_color_hue(15)[12])
p_s12b=p+xlab("generation to present")

d_s12a=data.frame(gen=gen,prs=t(Traj["ID225",1:31]))
colnames(d_s12a)[2]="prs"

p=ggplot(d_s12a,aes(gen,prs))+geom_point(color=gg_color_hue(15)[2])
p=p+theme_classic()+ scale_x_continuous(trans=reverselog_trans(10))
p=p+geom_line(color=gg_color_hue(15)[2])
p_s12a=p+xlab("generation to present")

p_s12a=p_s12a+ylab("Vegetable intake")
p_s12b=p_s12b+ylab("Atopic dematitis")
p_s12c=p_s12c+ylab("Intelligence")
p_s12d=p_s12d+ylab("Insomnia")
p_s12=p_s12a/p_s12b/p_s12c/p_s12d+
plot_annotation(tag_levels = 'A')

tiff("Figure S12.tiff",units="in",height=8,width=5,res=300)
p_f5
dev.off()

head(Traj)

x=(Traj[,31]-Traj[,27])/abs(Traj[,27])
head(res_sds)
y=res_sds[rownames(Traj),"rho"]
x[which(x>0.1)]=0.1
x[which(x<(-0.1))]=-0.1
plot(x,y)
cor.test(x,y,method="spearman")

y=res_sds[which(res_sds[,3]<0.05/870&rownames(res_sds) %in% rownames(Traj)),c(2,3)]
y$ptx=TxP[rownames(y)]
y$traj=Traj[rownames(y),31]-Traj[rownames(y),27]
head(y)
dim(y)
y[,1]=Traj[rownames(y),27]
y[,2]=Traj[rownames(y),31]
length(which(y[,1]*y[,4]>0&y[,3]<0.05))
write.csv(y,"y.csv")
197/755
y["ID225",]

Traj$cluster=clures[rownames(Traj)]
write.csv(Traj,"Traj.csv")

#
tmp=list.files(pattern="*txt")
psel=c()
for(i in tmp){
if(file.info(i)[1,1]<1000){
pinc=1
minc=1
pdec=1
mdec=1
}else{
dat=read.table(i,head=T)
pvec=dat[which(dat[,2]<0&dat[,3]!=1),3]
if(length(pvec)<3){
pdec=1
mdec=1
}else{
pvec=10^pvec
pdec=ks.test(pvec,seq(1/length(pvec),1,1/length(pvec)))$p.value
mdec=min(pvec)*length(pvec)
}
pvec=dat[which(dat[,2]>0&dat[,3]!=1),3]
if(length(pvec)<3){
pinc=1
minc=1
}else{
pvec=10^pvec
pinc=ks.test(pvec,seq(1/length(pvec),1,1/length(pvec)))$p.value
minc=min(pvec)*length(pvec)
}
}
psel=rbind(psel,c(pinc,minc,pdec,mdec))
}
psel=data.frame(psel)
rownames(psel)=tmp
psel$trait=annot[tmp,4]
psel[which(psel[,3]<0.05/870),]






