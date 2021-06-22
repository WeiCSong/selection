#########some of the ldsc annotations are trait-specific and were generated using this file.
library(data.table)
library(R.utils)
library(readr)

param <- commandArgs(trailingOnly=T)
i= eval(paste(text=param[1]))
trait=eval(paste(text=param[2]))
path=paste("//lustre//home//acct-bmelgn//bmelgn-3//ldsc//int//",trait,sep="")
i=normalizePath(i)
print(i)
gwas=fread(i,data.table=FALSE)
load("//lustre//home//acct-bmelgn//bmelgn-3//ldsc//annot.RData") #uploaded to github
hm=fread("w_hm3.snplist",data.table=FALSE)




for(j in 1:22){
h=head[[j]]
cadd=CADD[[j]]
cms=CMS[[j]]
hard=HARD[[j]]
soft=SOFT[[j]]
int=gwas[which(gwas$CHR==j),]
isai=int[match(h$SNP,int$SNP),"isaltinc"] #is the ancestral allele trait-increasing?
isai[which(is.na(isai))]=0
sds=int[match(h$SNP,int$SNP),"sds"]
sds[which(is.na(sds))]=0
AI=ifelse(isai==1,1,0)
AD=ifelse(isai==-1,1,0)
h$cadd_pos=cadd*AI
h$cadd_neg=cadd*AD
h$cms_pos=cms*AI
h$cms_neg=cms*AD
h$hard_pos=hard*AI
h$hard_neg=hard*AD
h$soft_pos=soft*AI
h$soft_neg=soft*AD
h$sds=sds
out=paste("//spe",j,sep=".")
out=paste(out,"annot.gz",sep=".")
out=paste(path,out,sep="")
write.table(h,gzfile(out),row.names=FALSE,quote=FALSE)
}




