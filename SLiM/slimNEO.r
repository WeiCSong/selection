library(data.table)

gwas=fread("plink.assoc.linear",data.table=F)
gwas=gwas[which(gwas$TEST=="ADD"),]
pruned=fread("plink.prune.in",data.table=F)
pruned=pruned[,1]
gwas=gwas[!duplicated(gwas[,2]),]
rownames(gwas)=gwas[,2]
gwas=gwas[intersect(pruned,rownames(gwas)),]
gwas=gwas[which(gwas$P<0.01),]
gwas$BETA=ifelse(gwas$A1=="T",gwas$BETA,-gwas$BETA)
tmp=list.files(pattern="*gen*")
tmp=tmp[grep("txt",tmp)]
pData=t(data.frame(strsplit(tmp,".",fixed=T)))[,1]
pData=t(data.frame(strsplit(pData,"gen",fixed=T)))
pData=data.frame(pData,lat="",long="")
colnames(pData)=c("pop","gen","lat","long")
pData[which(pData[,1]=="p2"),3]=rnorm(length(which(pData[,1]=="p2")),42,0.42)
pData[which(pData[,1]=="p2"),4]=rnorm(length(which(pData[,1]=="p2")),-2.47,1.37)
pData[which(pData[,1]=="p4"),3]=rnorm(length(which(pData[,1]=="p4")),47,0.7)
pData[which(pData[,1]=="p4"),4]=rnorm(length(which(pData[,1]=="p4")),19.5,1.2)
pData[which(pData[,1]=="p5"),3]=rnorm(length(which(pData[,1]=="p5")),52,0.9)
pData[which(pData[,1]=="p5"),4]=rnorm(length(which(pData[,1]=="p5")),10.6,1.2)
pData[,2:4]=apply(pData[,2:4],2,as.numeric)

splitVar=function(x){
return(as.vector(strsplit(x," ",fixed=T)))
}

getburden=function(ind){
variant=fread(ind,data.table=F)

rownames(variant)=variant[,2]
variant=variant[rownames(gwas),9,drop=F]
variant[which(is.na(variant[,1])),1]=0
rownames(variant)=rownames(gwas)
variant=data.frame(gt=variant[,1],coef=ifelse(gwas$BETA>0,1,-1),stringsAsFactors=F)
variant=(as.numeric(variant[,1])-1)*variant[,2]+1
hq=runif(1,0.1,0.9)
hq=ceiling(hq*length(variant))
variant=sample(variant,hq,replace=F)
burden=sum(variant)/(2*length(variant))
return(burden)
}

burden=unlist(sapply(tmp,getburden))
data=data.frame(pData,burden=burden)
model=lm(burden~lat+long+gen,data=data)
write.table(summary(model)$coefficients["gen",],"neo.txt",quote=F,col.names=F,row.names=F,sep=" ")




