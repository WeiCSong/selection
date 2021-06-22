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




