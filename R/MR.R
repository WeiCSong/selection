library(MRPRESSO)
library(TwoSampleMR)

library(doParallel)
cl <- makeCluster(14,outfile="log.txt")
registerDoParallel(cl)
#################data preprocess and harmonization#############################
#All gwas in the current study was manually reformatted to include columns: SNP,CHR,BP,A1,A2,b,se,p,N. 

helper1=function(int1,dis){
  exposure=data.frame(SNP=int1$SNP,beta=int1$b,se=int1$se,
  pos=int1$BP,other_allele=as.character(int1$A2),pval=int1$p,
  effect_allele=as.character(int1$A1),samplesize=int1$N,chr=int1$CHR)
  exposure<- format_data(exposure, type="exposure")
  exposure=clump_data(exposure)
  if(nrow(exposure)==0){return(c())}else{
    exposure$exposure=dis
    return(exposure)
  }
}

helper_harmonize=function(exp,expname,out){
  exposure=helper1(exp,expname)
  int2=out[which(out[,1] %in% exposure$SNP),]
  n=nrow(int2)
  if(n==0){d_strong=c()}else{
    outcome_dat <- data.frame(
    SNP=int2[,1],
    effect_allele.outcome=int2[,"A1"],
    other_allele.outcome=int2[,"A2"],
    beta.outcome=int2[,"BETA"],
    se.outcome=int2[,"StdErr"],
    pval.outcome=int2[,"P"],
    outcome=rep(expname,n),
    mr_keep.outcome=rep(TRUE,n),
    pval_origin.outcome=rep("reported",n),
    id.outcome=rep(expname,n),
    eaf.outcome=int2[,"frequency"],
    data_source.outcome=rep("textfile",n))
    d_strong <- harmonise_data(
      exposure_dat = exposure, 
      outcome_dat = outcome_dat)
    }
return(d_strong)
}

harmonized_ncm=list()
for(file in list.files(pattern="*gz")){
  name=sub(".gz","",file)
  gwas=fread(file,data.table=F)
  harmonized_ncm[[name]]=helper_harmonized(gwas,name,ncm) #ncm: formatted gwas summary of number of children for men. Similar gwas file name for ncf,nca,nsm,nsf,nsa.
}


################analysis#########################
head(harmonized_ncm[[1]])

helper_sensitivity=function(dat){
het=mr_heterogeneity(dat, method_list=c("mr_egger_regression", "mr_ivw"))
p_cochrane=het[2,8]
p_rucker=het[1,8]
ple_p=NA
if(nrow(dat)>2){
ple_p=mr_pleiotropy_test(dat)[1,7]
}
return(c(p_cochrane,p_rucker,ple_p))
}

helper_MR=function(dat){
mr=mr(dat, method_list=c("mr_egger_regression","mr_weighted_median", "mr_ivw"))
ivw_b=mr[3,7]
ivw_se=mr[3,8]
ivw_p=mr[3,9]
if(nrow(dat)<3){
	return(c(NA,NA,NA,NA,NA,NA,ivw_b,ivw_se,ivw_p))
	}else{
	egger_b=mr[1,7]
	egger_se=mr[1,8]
	egger_p=mr[1,9]
	wm_b=mr[2,7]
	wm_se=mr[2,8]
	wm_p=mr[2,9]
	return(c(egger_b,
	egger_se,egger_p,wm_b,wm_se,wm_p,ivw_b,ivw_se,ivw_p))
	}
}

helper_PRESSO=function(dat){
outlier=mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
SdOutcome = "se.outcome", SdExposure = "se.exposure", data=dat,
OUTLIERtest = TRUE,NbDistribution =nrow(dat)+1,SignifThreshold =1)
int=data.frame(outlier$`MR-PRESSO results`$`Outlier Test`)
if(nrow(int)==0){
n=nrow(dat)
l=1:n
sen=helper_sensitivity(dat[l[1:n],])}else{
	l=order(int[,1])
	n=length(l)
	repeat{
		sen=helper_sensitivity(dat[l[1:n],])
		if(min(sen[which(!is.na(sen))])>0.05|n==3){
				break
			}
		n=n-1
		}
	}
mrres=helper_MR(dat[l[1:n],])
gc()
return(c(paste(dat[l[1:n],1],collapse="|"),n,sen,mrres))
}

singleTrait_MR=function(dat){
exp=dat[1,"exposure"]
out=dat[1,"outcome"]
n=nrow(dat)
if(n==1|length(which(is.na(dat[,"se.exposure"])))!=0){
return(c(exp,out,n,rep(NA,26)))
}else{
res_presso=rep(NA,14)
if(n>3){
res_presso=helper_PRESSO(dat)
}
sen=helper_sensitivity(dat)
mrres=helper_MR(dat)
return(c(exp,out,n,sen,mrres,res_presso))
}
}

head(ncf_presso)
res_ncm=rbind(res_ncm,res_ncm1)

res_ncf=foreach(i=1:729,.combine=rbind,.packages=c("TwoSampleMR","MRPRESSO")) %dopar% {
int=singleTrait_MR(dat=harmonized_ncf[[i]])
return(int)
}
write.csv(res_ncf,"res_ncf.csv")
write.csv(res_ncm,"res_ncm.csv")

#####################################CAUSE#####################################
tmp=list.files(pattern="*gz")
install.packages("mixsqp", version = "0.1-97",INSTALL_opts="--no-multiarch", repos = "http://cran.us.r-project.org")
install.packages("ashr",INSTALL_opts="--no-multiarch", repos = "http://cran.us.r-project.org")
devtools::install_github("jean997/cause")
library(readr)
library(dplyr)
library(cause)
library(data.table)
hm3=fread("E:\\selection\\w_hm3.snplist",data.table=FALSE)
out=fread("E:\\all_GWAS\\rep\\nsex_female.gz",data.table=FALSE)
out6=out[which(out$SNP %in% hm3$SNP),c(1,4,5,6,7)]
outlist=list(out1,out2,out3,out4,out5,out6)

cause_helper=function(SA,j,expid){
out=outlist[[j]]
outid=c("nca","ncm","ncf","nsa","nsm","nsf")[j]
exp=SA[which(SA$SNP %in% hm3$SNP),c("SNP","A1","A2","b","se")]
X <- gwas_merge(exp,out, snp_name_cols = c("SNP", "SNP"), 
                       beta_hat_cols = c("b", "b"), 
                       se_cols = c("se", "se"), 
                       A1_cols = c("A1", "A1"), 
                       A2_cols = c("A2", "A2"))
params <- est_cause_params(X, X$snp)
rm(exp)
colnames(SA)=c("SNP","chr","pos","effect_allele","other_allele","pval","samplesize",
"beta","se")
exposure=SA[which(SA$pval<1e-5),]
exposure=exposure[which(exposure$SNP %in% X$snp),]
exposure<- format_data(exposure, type="exposure")
exposure=clump_data(exposure)
l=exposure$SNP
res <- cause(X=X, variants = l, param_ests = params,force=TRUE)
rm(SA)
rm(exposure)
gc()
return(c(expid,outid,res$elpd[,5]))
}

CAUSE=c()
for (i in unique(res_presso$exp)){
SA=fread(tmp[i],data.table=FALSE)
SA=SA[,c("SNP","CHR","BP","A1","A2","p","N","b","se")]
expid=strsplit(tmp[i],".gz",fixed=T)[[1]][1]
cause_single=function(j){
cause_helper(SA=SA,j=j,expid=expid)
}
res=data.frame(sapply(1:6,cause_single))

CAUSE_nca=rbind(CAUSE_nca,t(res))
rm(SA)
gc()
}

#################integration########################
library(ggplot2)
MR1$type=factor(MR1$type,levels=c("IMP","DER","REP","COG","self care","PSY",
"NUT","MET","GI","CIRC","immune","MED","MUSC","NEU","RES"))
p=ggplot(MR1[which(!is.na(MR1$type)&!is.na(MR1$ncm)),],aes(type,fill=ncm))+
geom_bar(stat="count",position="fill", width = 1)+ labs(fill = "effect on\nmale's\nN.child")+
theme_classic()+ scale_y_continuous(expand = c(0, 0),breaks=c(0,0.5,1))+ 
theme(axis.text.x = element_blank(),legend.text = element_text(size =10))
p=p+scale_fill_manual(values=brewer.pal(4,"Set3")[c(1,4)])
p=p+theme(axis.title=element_blank(),axis.ticks=element_blank(),legend.title = element_text(size =10),
axis.text.y = element_text(size=10,family="sans",color="black"),
legend.margin=margin(-2,-2,-2,-2),legend.box.margin=margin(-2,-2,-2,-2))
p_cm_bar=p

p=ggplot(MR1[which(!is.na(MR1$type)&!is.na(MR1$ncf)),],aes(type,fill=ncf))+
geom_bar(stat="count",position="fill", width = 1)+ labs(fill = "effect on\nfemale's\nN.child")+
theme_classic()+ scale_y_continuous(expand = c(0, 0),breaks=c(0,0.5,1))+ 
theme(axis.text.x = element_text(angle = -60, vjust = 0.2, hjust=0,size=10,
family="sans",color="black"),legend.text = element_text(size =10))
p=p+scale_fill_manual(values=brewer.pal(4,"Set3")[c(1,4)])
p=p+theme(axis.title=element_blank(),axis.ticks=element_blank(),legend.title = element_text(size =10),
axis.text.y = element_text(size=10,family="sans",color="black"),
legend.margin=margin(-2,-2,-2,-2),legend.box.margin=margin(-2,-2,-2,-2))
p_cf_bar=p

library(patchwork)
MR1[which(MR1$type=="RES"&MR1$nsa=="yes"),]
#################
head(MR)
hdat=data.frame(value=c(-4,4),x="dashed")
p=ggplot(MR,aes(ncm,ncf))+geom_vline(xintercept=c(-4,4),linetype="dashed")+
geom_hline(data=hdat,mapping=aes(yintercept=value,linetype=x),show.legend=TRUE)+geom_point()
p=p+theme_classic()+theme(axis.ticks=element_blank(),
axis.text=element_text(color="black",family="sans",size=10),
axis.title=element_text(color="black",family="sans",size=10))+
scale_linetype_identity()
p=p+theme(axis.title.x= element_blank())
p=p+ylab("Z for female")
p=p+ annotate(geom="text", x=7.5, y=5.6, label="Height",size=4)
p=p+ annotate(geom="text", x=-1.3, y=7.6, label="Body fat percentage",size=4)
p=p+ annotate(geom="text", x=5.5, y=-3.7, label="Skin tanning",size=4)
p=p+ annotate(geom="text", x=5, y=-5.5, label="Sitting height",size=4)
p=p+ annotate(geom="text", x=-6.5, y=-8.1, label="Intelligence",size=4)
p=p+ annotate(geom="text", x=-0.2, y=-6.9, label="Income",size=4)
p=p+ annotate(geom="text", x=-6, y=2.9, label="Left leg",size=4)
p_mr_scatter=p

length(which(!is.na(MR1$nsf)))
table(MR1[which(!is.na(MR1$nsm)),c("nsm","type")])
(p_cm_bar/p_cf_bar)|p_mr_scatter
table(MR1[,c("disease","nsf")])
cor.test(MR$nsf,MR$nsm)$p.value


