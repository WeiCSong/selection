library(data.table)

h2=list.files(pattern="a*")[1]
h2=strsplit(h2,"a",fixed=T)[[1]][2]
h2=as.numeric(h2)
all=readLines("final.txt")

N=6000
variant=all[grep(" m",all)]
individual=all[grep(" H ",all)]
genome=all[grep(" A ",all)][-1]

#variant and individual as data.frame
splitVar=function(x){
return(as.vector(strsplit(x," ",fixed=T)))
}

variant=data.frame(t(data.frame(sapply(variant,splitVar))),stringsAsFactors=F)
colnames(variant)=c("tid","pid","type","pos","s","dom","pop","gen","freq")
variant$freq=as.numeric(variant$freq)/N
variant$pos=as.numeric(variant$pos)+1

individual=data.frame(t(data.frame(sapply(individual,splitVar))),stringsAsFactors=F)
rownames(individual)=individual[,1]


#genotype as list
g=as.vector(sapply(genome,function(x){strsplit(x," ",fixed=T)[[1]][1]}))
genome=sapply(genome,splitVar)
names(genome)=g


#define phenotype
qtl=variant[which(variant$type=="m3"),1]
bqtl=variant[which(variant$type=="m3"),5]
qtlpos=as.numeric(variant[which(variant$type=="m3"),4])

vcf=fread("sample.vcf",data.table=F)
QTL=vcf[match(qtlpos,vcf[,2]),]
ref=data.frame(c("0|0","0|1","1|0","1|1"),c(0,1,1,2))
rownames(ref)=ref[,1]

getgt=function(x){
x=ref[x,2]
return(x)
}

QTL=QTL[,10:3009]

QTL=apply(QTL,2,getgt)
prs=as.numeric(bqtl) %*% QTL
prs=t(prs)[,1]
V_A = sd(prs)^2
V_E = (V_A - h2 * V_A) / h2
env = rnorm(N/2, 0.0, sqrt(V_E))
phenotypes = prs + env
individual$phenotype=phenotypes

###############################################################
pop=t(data.frame(strsplit(rownames(individual),":",fixed=T)))[,1]
ind=paste("i",0:2999,sep="")
individual=data.frame(individual,pop,ind)
rownames(individual)=individual$ind
fam=fread("sample.fam",data.table=F)
#main function
#singleLoci=function(pid){
#pos=variant[which(variant$pid==pid),"pos"]
#rsid=paste(1,pid,sep="_")
#tid=variant[which(variant$pid==pid),"tid"]
#freq=variant[which(variant$pid==pid),"freq"]
#getgt=function(ind){
#g=individual[ind,3:4]
#a=ifelse(tid %in% genome[[g[1,1]]],1,0)
#b=ifelse(tid %in% genome[[g[1,2]]],1,0)
#return(a+b)
#}
#gt=unlist(sapply(individual[,1],getgt))
#d=data.frame(gt=gt,phenotype=individual$phenotype,pop=pop)
#...
#return(c(rsid,pos,b,se,DAF,isdainc))
#}
fam[,6]=phenotypes
write.table(fam,"sample.fam",quote=F,col.names=F,row.names=F,sep="\t")
bim=fread("sample.bim",data.table=F)
bim[,2]=variant[match(bim[,4],variant$pos),2]
write.table(bim,"sample.bim",quote=F,col.names=F,row.names=F,sep=" ")
bim=data.frame(bim,variant[match(bim[,4],variant$pos),c(3,7,8,9)])
rownames(bim)=c()
write.table(bim,"allvariant.txt",quote=F,col.names=F,row.names=F,sep=" ")

################################SDS file###################################
obs=rep(1,1000)
write.table(t(obs),"obs.txt",quote=F,col.names=F,row.names=F,sep=" ")
variant=variant[order(as.numeric(variant$pos)),]
SNP=variant[which(variant$freq>0.05&variant$freq<0.95),]
rownames(SNP)=SNP[,1]
p5=individual[which(individual[,5]=="p5"),]



##################tfile######################################################
vcf=vcf[match(SNP[which(SNP$pos>2000000&SNP$pos<7000000),4],vcf[,2]),]
header=data.frame(SNP[which(SNP$pos>2000000&SNP$pos<7000000),2],"A","T",as.numeric(SNP[which(SNP$pos>2000000&SNP$pos<7000000),4]))

getgt=function(x){
x=ref[x,2]
return(x)
}

vcf=vcf[,2010:3009]
vcf=apply(vcf,2,getgt)
header=data.frame(header,vcf)
write.table(header,"tfile.txt",quote=F,col.names=F,row.names=F,sep=" ")

################s file#####################################################
P5=readLines("p5.txt")

v5=P5[grep(" m",P5)]
v5=data.frame(t(data.frame(sapply(v5,splitVar))),stringsAsFactors=F)
colnames(v5)=c("tid","pid","type","pos","s","dom","pop","gen","freq")
ST=v5[which(v5$freq==1),2]

g5=list()
for(i in 1:1000){
g5[[i]]=unique(c(genome[[3999+2*i]],genome[[4000+2*i]]))
}



singleton=variant[which(variant[,2] %in% ST),]
rownames(singleton)=singleton[,1]

getst=function(x){
int=intersect(x,rownames(singleton))

s=as.numeric(singleton[int,4])
if(length(s)<2){
up=sample(1:100000,1)
down=sample(9899999:9999999,1)
s=c(s,up,down)
}
s=s[order(s)]

return(s)
}
st=lapply(g5,getst)
lapply(st, write, "sfile.txt", append=TRUE, ncolumns=1000)







