library(data.table)
library(MESS)
load("/lustre/home/acct-bmelgn/bmelgn-3/slim/daf.RData")
gwas=fread("plink.assoc.linear",data.table=F)
gwas=gwas[which(gwas$TEST=="ADD"),]

gwas=gwas[!duplicated(gwas[,2]),]
rownames(gwas)=gwas[,2]
gwas=gwas[which(abs(gwas$STAT)>3),]
gwas$BETA=ifelse(gwas$A1=="T",gwas$BETA,-gwas$BETA)
gwas$isaltinc=ifelse(gwas$BETA>0,1,-1)
write.csv(gwas$SNP,"snp.txt",quote=F,row.names=F,col.names=F)
system("/lustre/home/acct-bmelgn/bmelgn-3/Plink --bfile sample --indep 50 5 1.5 --extract snp.txt --noweb")
pruned=fread("plink.prune.in",data.table=F)
pruned=pruned[,1]

gwas=gwas[as.character(pruned),]
#a=table(gwas$isaltinc)[1]
#b=table(gwas$isaltinc)[2]
#model=binom.test(b,a+b,0.5)
#p=model$p.value
#r=data.frame(p=p,a=a,b=b)
#write.table(r,"daf.txt",quote=F,row.names=F,col.names=F)








daf=daf[order(daf)]
allvariant=fread("allvariant.txt",data.table=F)
allvariant=allvariant[!duplicated(allvariant[,2]),]
rownames(allvariant)=allvariant[,2]
gwas$se=gwas$BETA/gwas$STAT
gwas$r2=(2*gwas$BETA^2)/(2*gwas$BETA^2+6000*gwas$se^2)
gwas$DAF=allvariant[rownames(gwas),10]
gwas=gwas[which(!is.na(gwas$DAF)),]
ord=order(gwas$DAF)
ord=quantile(daf,probs=ord/length(ord))
gwas$DAF=ord
gwas=gwas[which(gwas$DAF>0.05),]


up=gwas[which(gwas$BETA>0),]
down=gwas[which(gwas$BETA<0),]
int=data.frame(ch2=up$r2,DAF=up$DAF)
int=int[order(int$DAF),]
int$ch2=cumsum(int$ch2)/sum(int$ch2)
int$index=1:nrow(int)/nrow(int)
S1=auc(int$index,predict(loess(ch2~index,data=int)))

int=data.frame(ch2=down$r2,DAF=down$DAF)
int=int[order(int$DAF),]
int$ch2=cumsum(int$ch2)/sum(int$ch2)
int$index=1:nrow(int)/nrow(int)
S2=auc(int$index,predict(loess(ch2~index,data=int)))

S=S1-S2

jackknife=function(i){
index=setdiff(1:nrow(gwas),(i*20):(i*20+19))
gwas=gwas[index,]
up=gwas[which(gwas$BETA>0),]
down=gwas[which(gwas$BETA<0),]
int=data.frame(ch2=up$r2,DAF=up$DAF)
int=int[order(int$DAF),]
int$ch2=cumsum(int$ch2)/sum(int$ch2)
S1=auc(int$DAF,predict(loess(ch2~DAF,data=int)))

int=data.frame(ch2=down$r2,DAF=down$DAF)
int=int[order(int$DAF),]
int$ch2=cumsum(int$ch2)/sum(int$ch2)

S2=auc(int$DAF,predict(loess(ch2~DAF,data=int)))
return(c(S1-S2,S1,S2))
}
jkindex=1:ceiling(nrow(gwas)/20)
jk=t(data.frame(sapply(jkindex,jackknife)))
SD=apply(jk,2,sd)
res=data.frame(S=c(S1-S2,S1,S2),sd=SD)
write.table(res,"S.txt",quote=F,col.names=F,row.names=F,sep=" ")
