library(data.table)
library(MESS)

#res_sds_one=c()
#for(i in 1:40){
#prefix=paste("./burnin",i,sep="")
#prefix=paste(prefix,"/sim1/",sep="")
gwas=fread("plink.assoc.linear",data.table=F)
gwas=gwas[which(gwas$TEST=="ADD"),]
gwas=gwas[!duplicated(gwas[,2]),]
rownames(gwas)=gwas[,2]

gwas$BETA=ifelse(gwas$A1=="T",gwas$BETA,-gwas$BETA)

sds=fread("sds1.txt",data.table=F)
sds[,1]=gwas[match(sds$POS,gwas$BP),2]
sds$p=gwas[as.character(sds[,1]),"P"]
sds$SDS=0
sds=sds[which(sds$DAF>0.05&sds$DAF<0.95),]

######re-normalized per DAF bin
for(i in seq(0.05,0.94,0.01)){
x=sds[which(sds$DAF>i&sds$DAF<i+0.01),"rSDS"]
sds[which(sds$DAF>i&sds$DAF<i+0.01),"SDS"]=x
}
sds[which(is.na(sds$SDS)),"SDS"]=0
########spearman rho 
sds=sds[order(sds$p),]
res=c()
for(i in seq(1,nrow(sds),60)){
start=i
end=min(i+60,nrow(sds))
s=median(sds[start:end,"SDS"])
res=rbind(res,c(i,s))
}
rho_true=cor(res[,1],res[,2],method="spearman")
n=nrow(res)
l=order(sds$POS)
##########block jackknife for rho
bj=c()
for (i in round(seq(1,nrow(sds),length.out=60))){
L=l[-(i:(i+round(nrow(sds)/60)+1))]
LL=L[order(L)]
res=c()
for(j in round(seq(1,length(L),length.out=60))){
start=j
end=min(j+round(length(L)/60)+1,length(L))

s=median(sds[LL[start:end],"SDS"])
res=rbind(res,c(j,s))
}
rho=cor(res[,1],res[,2],method="spearman")
bj=rbind(bj,c(i,rho))
}


################p value
z=rho_true/sd(bj[,2])
p=2*pnorm(-abs(z))
line=data.frame(r=rho_true,z=z,n=n)
write.table(line,"rho.txt",quote=F,col.names=F,row.names=F,sep=" ")
summary(res_sds_one[which(res_sds_one[,3]>29),2])


