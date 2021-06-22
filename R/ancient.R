param <- commandArgs(trailingOnly=T)
name= eval(paste(text=param[1]))

library(data.table)
library(R.utils)
setwd("/lustre/home/acct-bmelgn/bmelgn-3/ldsc/all_gwas/")
nam=strsplit(name,".",fixed=T)[[1]][1]
SA=fread(name,data.table=FALSE)
SA=SA[which(SA$p<0.01),]
SA$inc=ifelse(SA$b>0,SA$A1,SA$A2)
sasnp=SA$SNP
SA=SA[,c("SNP","inc")]
gc()

#########neareast
load("/lustre/home/acct-bmelgn/bmelgn-3/neareast/neareast_d.RData")
traitsnp=intersect(sasnp,snp[,1])
gwas=SA[which(SA$SNP %in% traitsnp),]
file=paste("list_",name,sep="")
file=paste(file,".txt",sep="")
write.table(traitsnp,file,quote=FALSE,row.names=F,col.names=F)

newwd=paste("/lustre/home/acct-bmelgn/bmelgn-3/ldsc/int/",nam,sep="")
last=paste("--extract ",file,sep="")
last=paste(last,"--make-bed --noweb --out /lustre/home/acct-bmelgn/bmelgn-3/ldsc/int/",sep=" ")
last=paste(last,nam,sep="")
last=paste(last,"/plink",sep="")
command=paste("/lustre/home/acct-bmelgn/bmelgn-3/plink/plink",
"--bed /lustre/home/acct-bmelgn/bmelgn-3/neareast/plink.bed",
sep=" ")
command=paste(command,
"--bim /lustre/home/acct-bmelgn/bmelgn-3/neareast/plink.bim",
sep=" ")
command=paste(command,
"--fam /lustre/home/acct-bmelgn/bmelgn-3/neareast/plink.fam",
sep=" ")
command=paste(command,last,sep=" ")
system(command)

setwd(newwd)
command=paste("/lustre/home/acct-bmelgn/bmelgn-3/plink/plink",
"--bed plink.bed",sep=" ")
command=paste(command,"--bim plink.bim",sep=" ")
command=paste(command,"--fam plink.fam",sep=" ")
command=paste(command,"--indep 50 5 2 --noweb",sep=" ")
system(command)

list=read.table("plink.prune.in")
list=as.character(list[,1])
gwas=gwas[match(list,gwas$SNP),]
l_snp=match(list,snp[,1])
D=d[l_snp,]-1
snplist=snp[l_snp,]
coef=ifelse(snplist[,2]==gwas[,2],1,-1)
f=function(x){
x*coef
}
D=apply(D,2,f)+1
f=function(x){
nas=length(which(is.na(x)))
non=length(x)-nas
if(nas>nas>0.9*nrow(D)){
return(NA)}else{
return(sum(na.omit(x))/(2*non))
}
}
ftrait=apply(D,2,f)
ind$f=ftrait
colnames(ind)=c("ID","lat","long","sex","type","cov","nsnp","date","f")
model=lm(f~lat+long+sex+type+cov+date,data=ind)
INT=data.frame(summary(model)$coefficient[c(2,3,4,6,7),],set="neareast",trait=nam)
INT$coef=rownames(INT)
rownames(INT)=c()
colnames(INT)=c("est","se","t","p","set","trait","coef")
int=rbind(int,INT)

##########proneo

setwd("/lustre/home/acct-bmelgn/bmelgn-3/ldsc/all_gwas/")
load("/lustre/home/acct-bmelgn/bmelgn-3/proneo/proneo_d.RData")
traitsnp=intersect(sasnp,snp[,1])
gwas=SA[which(SA$SNP %in% traitsnp),]
file=paste("list_",name,sep="")
file=paste(file,".txt",sep="")
write.table(traitsnp,file,quote=FALSE,row.names=F,col.names=F)

newwd=paste("/lustre/home/acct-bmelgn/bmelgn-3/ldsc/int/",nam,sep="")
last=paste("--extract ",file,sep="")
last=paste(last,"--make-bed --noweb --out /lustre/home/acct-bmelgn/bmelgn-3/ldsc/int/",sep=" ")
last=paste(last,nam,sep="")
last=paste(last,"/plink",sep="")
command=paste("/lustre/home/acct-bmelgn/bmelgn-3/plink/plink",
"--bed /lustre/home/acct-bmelgn/bmelgn-3/proneo/plink.bed",
sep=" ")
command=paste(command,
"--bim /lustre/home/acct-bmelgn/bmelgn-3/proneo/plink.bim",
sep=" ")
command=paste(command,
"--fam /lustre/home/acct-bmelgn/bmelgn-3/proneo/plink.fam",
sep=" ")
command=paste(command,last,sep=" ")
system(command)

setwd(newwd)
command=paste("/lustre/home/acct-bmelgn/bmelgn-3/plink/plink",
"--bed plink.bed",sep=" ")
command=paste(command,"--bim plink.bim",sep=" ")
command=paste(command,"--fam plink.fam",sep=" ")
command=paste(command,"--indep 50 5 2 --noweb",sep=" ")
system(command)

list=read.table("plink.prune.in")
list=as.character(list[,1])
gwas=gwas[match(list,gwas$SNP),]
l_snp=match(list,snp[,1])
D=d[l_snp,]-1
snplist=snp[l_snp,]
coef=ifelse(snplist[,2]==gwas[,2],1,-1)
f=function(x){
x*coef
}
D=apply(D,2,f)+1
f=function(x){
nas=length(which(is.na(x)))
non=length(x)-nas
if(nas>0.9*nrow(D)){
return(NA)}else{
return(sum(na.omit(x))/(2*non))
}
}
ftrait=apply(D,2,f)
ind$f=ftrait
colnames(ind)=c("ID","lat","long","date","type","sex","damage","cov","nsnp","f")
ind[23,"long"]=3.45
ind[47,"long"]=5.35
ind$long=as.numeric(ind$long)
model=lm(f~lat+long+sex+type+cov+damage+date,data=ind)
INT=data.frame(summary(model)$coefficient[c(2,3,4,8,10),],set="proneo",trait=nam)
INT$coef=rownames(INT)
rownames(INT)=c()
colnames(INT)=c("est","se","t","p","set","trait","coef")
int=rbind(int,INT)

########neolithic
setwd("/lustre/home/acct-bmelgn/bmelgn-3/ldsc/all_gwas/")
load("/lustre/home/acct-bmelgn/bmelgn-3/neolithic/neolithic_d.RData")
traitsnp=intersect(sasnp,snp[,1])
gwas=SA[which(SA$SNP %in% traitsnp),]
file=paste("list_",name,sep="")
file=paste(file,".txt",sep="")
write.table(traitsnp,file,quote=FALSE,row.names=F,col.names=F)

newwd=paste("/lustre/home/acct-bmelgn/bmelgn-3/ldsc/int/",nam,sep="")
last=paste("--extract ",file,sep="")
last=paste(last,"--make-bed --noweb --out /lustre/home/acct-bmelgn/bmelgn-3/ldsc/int/",sep=" ")
last=paste(last,nam,sep="")
last=paste(last,"/plink",sep="")
command=paste("/lustre/home/acct-bmelgn/bmelgn-3/plink/plink",
"--bed /lustre/home/acct-bmelgn/bmelgn-3/neolithic/plink.bed",
sep=" ")
command=paste(command,
"--bim /lustre/home/acct-bmelgn/bmelgn-3/neolithic/plink.bim",
sep=" ")
command=paste(command,
"--fam /lustre/home/acct-bmelgn/bmelgn-3/neolithic/plink.fam",
sep=" ")
command=paste(command,last,sep=" ")
system(command)

setwd(newwd)
command=paste("/lustre/home/acct-bmelgn/bmelgn-3/plink/plink",
"--bed plink.bed",sep=" ")
command=paste(command,"--bim plink.bim",sep=" ")
command=paste(command,"--fam plink.fam",sep=" ")
command=paste(command,"--indep 50 5 2 --noweb",sep=" ")
system(command)

list=read.table("plink.prune.in")
list=as.character(list[,1])
gwas=gwas[match(list,gwas$SNP),]
l_snp=match(list,snp[,1])
D=d[l_snp,]-1
snplist=snp[l_snp,]
coef=ifelse(snplist[,2]==gwas[,2],1,-1)
f=function(x){
x*coef
}
D=apply(D,2,f)+1
f=function(x){
nas=length(which(is.na(x)))
non=length(x)-nas
if(nas>0.9*nrow(D)){
return(NA)}else{
return(sum(na.omit(x))/(2*non))
}
}
ftrait=apply(D,2,f)
ind$f=ftrait
colnames(ind)=c("ID","lat","long","date","type","sex","damage","cov","nsnp","f")
ind[23,"long"]=3.45
ind[47,"long"]=5.35
ind$long=as.numeric(ind$long)
model=lm(f~lat+long+sex+type+cov+damage+date,data=ind)
INT=data.frame(summary(model)$coefficient[c(2,3,4,8,10),],set="neolithic",trait=nam)
INT$coef=rownames(INT)
rownames(INT)=c()
colnames(INT)=c("est","se","t","p","set","trait","coef")
int=rbind(int,INT)

write.csv(int,"res_ancient.csv")


