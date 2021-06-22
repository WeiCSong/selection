library(data.table)
library(R.utils)
library(readr)

cadd="F:\\selection\\anno_cadd\\cadd."
cms="F:\\selection\\anno_cms\\cms."
hard="F:\\selection\\anno_hard\\hard."
sds="F:\\selection\\anno_sds\\sds."
soft="F:\\selection\\anno_soft\\soft."
base="F:\\selection\\baseline\\baselineLD."
CADD=list()
CMS=list()
HARD=list()
SDS=list()
SOFT=list()
head=list()
for(i in 1:22){
file=paste(hard,i,sep="")
file=paste(file,".annot.gz",sep="")
hardannot=fread(file,data.table=FALSE)
HARD[[i]]=hardannot[,1]
}

file=paste(soft,i,sep="")
file=paste(file,".annot.gz",sep="")
sannot=fread(file,data.table=FALSE)
SOFT[[i]]=sannot[,1]
}

for(i in 1:22){
file=paste(base,i,sep="")
file=paste(file,".annot.gz",sep="")
bannot=fread(file,data.table=FALSE)
head[[i]]=bannot[,1:4]
}

rm(bannot)
rm(cannot)
rm(cmsannot)
rm(hardannot)
rm(prob)
rm(sannot)
save.image("partition.RData")
