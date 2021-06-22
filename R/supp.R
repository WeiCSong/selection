#Generation of supplementary figures. Similar to the main figures.

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
gg_color_hue(15)[7]
#figure sMR
p=ggplot(MR1[which(!is.na(MR1$type)&!is.na(MR1$nsm)),],aes(type,fill=nsm))+
geom_bar(stat="count",position="fill", width = 1)+ labs(fill = "effect on\nmale's\nN.partner")+
theme_classic()+ scale_y_continuous(expand = c(0, 0),breaks=c(0,0.5,1))+ 
theme(axis.text.x = element_blank(),legend.text = element_text(size =10))
p=p+scale_fill_manual(values=brewer.pal(4,"Set3")[c(1,4)])
p=p+theme(axis.title=element_blank(),axis.ticks=element_blank(),legend.title = element_text(size =10),
axis.text.y = element_text(size=10,family="sans",color="black"),
legend.margin=margin(-2,-2,-2,-2),legend.box.margin=margin(-2,-2,-2,-2))
p_sm_bar=p

p=ggplot(MR1[which(!is.na(MR1$type)&!is.na(MR1$nsf)),],aes(type,fill=nsf))+
geom_bar(stat="count",position="fill", width = 1)+ labs(fill = "effect on\nfemale's\nN.child")+
theme_classic()+ scale_y_continuous(expand = c(0, 0),breaks=c(0,0.5,1))+ 
theme(axis.text.x = element_text(angle = -60, vjust = 0.2, hjust=0,size=10,
family="sans",color="black"),legend.text = element_text(size =10))
p=p+scale_fill_manual(values=brewer.pal(4,"Set3")[c(1,4)])
p=p+theme(axis.title=element_blank(),axis.ticks=element_blank(),legend.title = element_text(size =10),
axis.text.y = element_text(size=10,family="sans",color="black"),
legend.margin=margin(-2,-2,-2,-2),legend.box.margin=margin(-2,-2,-2,-2))
p_sf_bar=p

p_sm_bar/p_sf_bar

p=ggplot(MR,aes(nsm,nsf))+geom_vline(xintercept=c(-4,4),linetype="dashed")+
geom_hline(data=hdat,mapping=aes(yintercept=value,linetype=x),show.legend=TRUE)+geom_point()
p=p+theme_classic()+theme(axis.ticks=element_blank(),
axis.text=element_text(color="black",family="sans",size=10),
axis.title=element_text(color="black",family="sans",size=10))+
scale_linetype_identity()
p=p+theme(axis.title.x= element_blank())
p=p+ylab("Z for female")
p=p+ annotate(geom="text", x=17, y=19, label="Ever smoker",size=4)
p=p+ annotate(geom="text", x=8, y=8, label="Risk intolerance",size=4)
p_nsex_scatter=p

p=ggplot(MR,aes(ncm,nsm))+geom_vline(xintercept=c(-4,4),linetype="dashed")+
geom_hline(data=hdat,mapping=aes(yintercept=value,linetype=x),show.legend=TRUE)+geom_point()
p=p+theme_classic()+theme(axis.ticks=element_blank(),
axis.text=element_text(color="black",family="sans",size=10),
axis.title=element_text(color="black",family="sans",size=10))+
scale_linetype_identity()
p=p+ylab("Z for N.partner")+xlab("Z for N.child")
p=p+ annotate(geom="text", x=5.5, y=19, label="Ever smoker",size=4)
p=p+ annotate(geom="text", x=-7.5, y=-6, label="Left arm",size=4)
p=p+ annotate(geom="text", x=-5.5, y=3, label="Loneliness",size=4)
p_m_scatter=p
cor.test(MR$ncm,MR$nsm)

p=ggplot(MR,aes(ncf,nsf))+geom_vline(xintercept=c(-4,4),linetype="dashed")+
geom_hline(data=hdat,mapping=aes(yintercept=value,linetype=x),show.legend=TRUE)+geom_point()
p=p+theme_classic()+theme(axis.ticks=element_blank(),
axis.text=element_text(color="black",family="sans",size=10),
axis.title=element_text(color="black",family="sans",size=10))+
scale_linetype_identity()
p=p+ylab("Z for N.partner")+xlab("Z for N.child")
p=p+ annotate(geom="text", x=8, y=-5, label="Leg fat",size=4)
p=p+ annotate(geom="text", x=-6, y=9, label="Intelligence",size=4)
p_f_scatter=p

cor.test(MR$ncf,MR$nsf)
p_smr=((p_sm_bar/p_sf_bar)|p_nsex_scatter)/(p_m_scatter|p_f_scatter)

pdf("fsmr.pdf",height=8,width=8)
p_smr
dev.off()
annot["ID600","h2"]=0.0136
annot=annot[which(annot$h2>0.01),]

mr=mr[which(mr$expid %in% annot[,1]),]
mr$type=annot[match(mr$expid,annot[,1]),"abr"]
mr$disease=annot[match(mr$expid,annot[,1]),"disease"]
mr$name=annot[match(mr$expid,annot[,1]),"uniqTrait"]
mr$time=annot[match(mr$expid,annot[,1]),"onset"]

mr$rgsig=ifelse(abs(mr$z)>7,"yes","no")
mr$hz=mr$h2/mr$h2se
mr$hsig=ifelse(abs(mr$hz)>7,"yes","no")
head(mr)
mr=mr[which(mr$out.x %in% c("ncm","ncf","nsm","nsf")),]
write.csv(mr,"res_cause.csv")

########cause
mr_ncm=mr[which(mr$out.x=="ncm"),]
length(which(abs(mr_ncm$z)>7))
length(which(abs(mr_ncm$z)>7&mr_ncm$estimation!=0))
length(which(abs(mr_ncm$z)>7&mr_ncm$estimation!=0&mr_ncm$cs<(-1.96)))
length(which(abs(mr_ncm$z)>7&mr_ncm$estimation!=0&mr_ncm$cs<0))
mr_ncm[which(abs(mr_ncm$z)>7&mr_ncm$estimation!=0),"cs"]

##############sds
indraw=data.frame(r=abs(res_sds[,2]),ukb=factor(res_sds$ukb),fil=
ifelse(res_sds$cat %in% c("DER","NUT"),0,1))
indraw[which(is.na(indraw$ukb)),"ukb"]=0
oneway_test(r~ukb,data=indraw[which(indraw$fil==1),])

tiff("fssda.tiff",res=300,units="in",height=4,width=4)
plot(sdsden$x,sdsden$y)
dev.off()
p_fsds_a=ggplot(indraw[which(indraw$fil==1),],aes(ukb,r))+geom_boxplot()+theme_classic()+ylab("|rho|")+
scale_x_discrete(labels=c("no","yes"))
########spearman rho 
library(data.table)
library(R.utils)
library(ggplot2)

gwas=fread("ID572.gz",data.table=F)
gwas=gwas[order(gwas$p),]
gwas$SDS=0
gwas=gwas[which(gwas$sds!=0),]
x=match(gwas$SNP,sds$ID)
gwas$freq=sds[x,"DAF"]
for(i in seq(0.05,0.94,0.01)){
x=gwas[which(gwas$freq>i&gwas$freq<i+0.01),"sds"]
x=(x-mean(x))/sd(x)
gwas[which(gwas$freq>i&gwas$freq<i+0.01),"SDS"]=x
}
res=c()
for(i in seq(1,nrow(gwas),1000)){
start=i
end=min(i+1000,nrow(gwas))
s=median(gwas[start:end,"SDS"])
res=rbind(res,c(i,s))
}
dim(res)
res[,1]=4149:1
head(res)
colnames(res)=c("bin","SDS")
p=ggplot(data.frame(res),aes(bin,SDS))+geom_point(color=gg_color_hue(15)[1],size=0.3)+theme_classic()
p=p+xlab("p value bin (high to low)")+ylab("bin median SDS")
p=p+theme(axis.ticks=element_blank(),axis.text=element_text(color="black",
size=10,family="sans"),axis.title=element_text(size=10,color="black",family="sans"))
p=p+annotate(geom="text", x=2500, y=0.3, 
label="Ulcerative colitis:\nSpearman rho=0.36\np by block jackknife<e-200",size=4)
p_fdis_c=p

sdsden=data.frame(x=sdsden$x,y=sdsden$y)
p=ggplot(sdsden,aes(x,y))+geom_point()+theme_classic()
p=p+xlab("|rho|")+ylab("density")
p_fsds_b=p

p_fsds=(p_fsds_a|p_fsds_b)/(p_fsds_c|p_fsds_d)/(p_fsds_e|p_fsds_f)
pdf("fsds.pdf",height=10,width=8)
p_fsds
dev.off()

library(patchwork)
p_fdis=(p_fsdis_a|p_fdis_b)/(p_fdis_c|p_fdis_d)
pdf("fis.pdf",height=7,width=7)
p_fdis
dev.off()

######sds disease
oneway.test(rho~disease,data=res_sds)
oneway.test(abs~disease,data=res_sds)
res_sds$ons=ifelse(res_sds$onset=="early",1,0)
oneway.test(abs~onset,data=res_sds[which(res_sds$onset %in% c("early","old")),])
median(res_sds[which(res_sds$onset=="early"),"rho"])
res_sds[which(res_sds$onset==""),"onset"]="trait"

p=ggplot(res_sds,aes(onset,rho,color=onset))+geom_boxplot()+theme_classic()
p=p+xlab("disease onset")+theme(legend.position="none")
p_fsdis_a=p

p=ggplot(res_sds,aes(onset,abs,color=onset))+geom_boxplot()+theme_classic()
p=p+xlab("disease onset")+theme(legend.position="none")
p_fsdis_b=p+ylab("|rho|")

median(res_sds[which(res_sds$disease==1),"abs"])

###########sanc
rownames(d_man)=1:3436
d_man$orig=10^(-d_man$p)
d_man=d_man[which(d_man$ID %in% annot[,1]),]
ks.test(res_count[which(res_count$coef=="date"&res_count$set=="proneo"),"p"], "punif", 0, 1)$p.value

d_man[which(duplicated(d_man$p)),]

d_man[which(d_man$p==0.511751504),]
annot=annot[-which(annot[,1]=="ID311"),]
res_prs=read.csv("res_prs.csv")
head(res_prs)
res_count=res_count[order(res_count$trait),]
res_prs=res_prs[which(res_prs$trait %in% annot[,1]),]

dsan=data.frame(count=res_count[which(res_count$coef=="hg"),"t"],
prs=res_prs[which(res_prs$coef=="hg"),"t"])
cor(dsan)
p=ggplot(dsan,aes(count,prs))+geom_point()+theme_classic()
p=p+xlab("t for count")+ylab("t for prs")
p=p+annotate(geom="text", x=-10, y=10, label="%HG: PCC=0.98",size=4)
p_fsanc_a=p

int1=res_count[which(res_count$coef=="date"),]
rownames(int1)=paste(int1$set,int1$trait,sep=".")
int2=res_prs[which(res_prs$coef=="date"),]
rownames(int2)=paste(int2$set,int2$trait,sep=".")
l=intersect(rownames(int1),rownames(int2))
int1=int1[l,]
int2=int2[l,]
dsan=data.frame(count=int1[which(int1$set=="neolithic"),"t"],
prs=int2[which(int2$set=="neolithic"),"t"])
cor(dsan)
p=ggplot(dsan,aes(count,prs))+geom_point()+theme_classic()
p=p+xlab("t for count")+ylab("t for prs")
p=p+annotate(geom="text", x=0, y=4, label="neolithic: PCC=0.96",size=4)
p_fsanc_b=p

dsan=data.frame(count=int1[which(int1$set=="proneo"),"t"],
prs=int2[which(int2$set=="proneo"),"t"])
cor(dsan)
p=ggplot(dsan,aes(count,prs))+geom_point()+theme_classic()
p=p+xlab("t for count")+ylab("t for prs")
p=p+annotate(geom="text", x=0, y=5, label="proneolithic: PCC=0.96",size=4)
p_fsanc_d=p

p_fsan=(p_fsanc_a|p_fsanc_b)/(p_fsanc_c|p_fsanc_d)
pdf("fsanc.pdf",width=6,height=6)
p_fsan
dev.off()

#####sane
head(d_man)
median(d_man[which(d_man$disease==1&d_man$var=="neolithic"),"t"])
cor(d_man[which(d_man$disease==1&d_man$var=="proneo"),"t"],
d_man[which(d_man$disease==1&d_man$var=="hg"),"t"])

dper=d_man[which(d_man$var=="neolithic"),c("t","disease")]
oneway.test(t~disease,data=dper)

p=ggplot(pro,aes(date,f))+geom_point(color=gg_color_hue(15)[1])+theme_classic()
p=p+theme(axis.ticks=element_blank(),
          axis.text= element_text(size=10,family="sans",color="black"),
	    axis.title= element_text(size=10,family="sans",color="black"))
p=p+xlab("year to present (proneolithic)")+ylab("scaled polygenic burden")+scale_x_log10()
p_fsane_a=p+scale_x_reverse()

dneo=neo[,c(14,17)]
colnames(dneo)=c("date","f")
dneo$f=1-dneo$f
p=ggplot(dneo,aes(date,f))+geom_point(color=gg_color_hue(15)[1])+theme_classic()
p=p+theme(axis.ticks=element_blank(),
          axis.text= element_text(size=10,family="sans",color="black"),
	    axis.title= element_text(size=10,family="sans",color="black"))
p=p+xlab("year to present (neolithic dataset)")+ylab("scaled polygenic burden")
p_fsane_b=p+scale_x_reverse()

p=ggplot(pro,aes(date,f))+geom_point(color=gg_color_hue(15)[11])+theme_classic()
p=p+theme(axis.ticks=element_blank(),
          axis.text= element_text(size=10,family="sans",color="black"),
	    axis.title= element_text(size=10,family="sans",color="black"))
p=p+xlab("year to present (proneolithic)")+ylab("scaled polygenic burden")
p_fsane_c=p+scale_x_reverse()+ggtitle("Crohn's disease")

p=ggplot(pro,aes(date,f))+geom_point(color=gg_color_hue(15)[11])+theme_classic()
p=p+theme(axis.ticks=element_blank(),
          axis.text= element_text(size=10,family="sans",color="black"),
	    axis.title= element_text(size=10,family="sans",color="black"))
p=p+xlab("year to present (neolithic dataset)")+ylab("scaled polygenic burden")
p_fsane_d=p+scale_x_reverse()+ggtitle("atopic dermatitis")

p=ggplot(pro,aes(date,f))+geom_point(color=gg_color_hue(15)[11])+theme_classic()
p=p+theme(axis.ticks=element_blank(),
          axis.text= element_text(size=10,family="sans",color="black"),
	    axis.title= element_text(size=10,family="sans",color="black"))
p=p+xlab("year to present (proneolithic)")+ylab("scaled polygenic burden")
p_fsane_e=p+scale_x_reverse()+ggtitle("periodontis")

pro=pro[which(pro$date<12000),]
p=ggplot(pro,aes(date,f))+geom_point(color=gg_color_hue(15)[12])+theme_classic()
p=p+theme(axis.ticks=element_blank(),
          axis.text= element_text(size=10,family="sans",color="black"),
	    axis.title= element_text(size=10,family="sans",color="black"))
p=p+xlab("year to present (neareast)")+ylab("scaled polygenic burden")
p_fsane_f=p+scale_x_reverse()+ggtitle("fracture")

p_fsane_top=p_fsane_a|p_fsane_b
p_fsane_top=p_fsane_top+ggtitle("Ease of skin tanning")
p_fsane=p_fsane_top/(p_fsane_c|p_fsane_d)/(p_fsane_e|p_fsane_f)
pdf("fsane.pdf",height=8,width=6)
p_fsane
dev.off()

##############
library(circlize)
library("ComplexHeatmap")
res_ks[1:3,1:3]"#3288BD", 
col_fun = colorRamp2(c(0.5, 0.55), c("white","#B2182B" ))
p_fsl_hm=Heatmap(res_ks[,l],col = col_fun,name="ks D", cluster_rows = FALSE, 
cluster_columns = F,border = TRUE,row_names_side="left", row_names_gp = gpar(fontsize =10),
column_names_gp = gpar(fontsize =10),heatmap_legend_param = list(direction = "horizontal", 
    title_position = "lefttop"))

pdf("fshm_ks.pdf",height=4,width=4)
draw(p_fsl_hm,heatmap_legend_side ="top")
dev.off()

head(ldsc)
dshmc=ldsc[,c(4,7)]
colnames(dshmc)=c("z","dis")
dshmc$disease=ifelse(dshmc$dis==1,"yes","no")
p=ggplot(dshmc,aes(z,group=disease,color=disease,fill=disease))+geom_density()
p=p+theme_classic()+xlab("LDSC z for Recombination rate")
p


l=column_order(p_fsl_hm)
#################
ldsc=ldsc[which(ldsc$trait %in% annot[,1]),]
ldsc$type=annot[ldsc$trait,"abr"]
Lpli=ldsc[which(ldsc[,1]=="PLIL2_1"),c(4,5,6,7,8)]
colnames(Lpli)[1]="p"
Lpli[,1]=1-pnorm(abs(Lpli[,1]))
Lpli[,1]=-log10(Lpli[,1])
colorannot=data.frame(type=grouporder,color=gg_color_hue(15))
Lpli$label=""
Lpli[which(Lpli$trait=="ID503"),"label"]="Schizophrenia"
Lpli[which(Lpli$trait=="ID313"),"label"]="Intelligence"
Lpli[which(Lpli$trait=="ID55"),"label"]="Ease of getting up"
Lpli[which(Lpli$trait=="ID405"),"label"]="Neutrophil count"
Lpli[which(Lpli$trait=="ID244"),"label"]="Worry"
Lpli$type=factor(Lpli$type,levels=grouporder)
p=ggplot(Lpli,aes(type,p,label=label,color=type))+geom_jitter(aes(group=type),width=0.4)
p=p+theme_classic()+theme(
	axis.text.x = element_text(angle = -60, vjust = 0.2,hjust=0,size=10,
		family="sans",color="black"),
	axis.ticks=element_blank(),
	axis.title.x=element_blank(),
	axis.text.y = element_text(size=10,family="sans",color="black"),
	axis.title.y = element_text(size=10,family="sans",color="black"),
	legend.position="none")
p=p+geom_hline(yintercept=-log10(0.05/873))+ylab("-log10(p for pLi)")
p=p+geom_text_repel(size=3,segment.size = 0.1)
p_pli_mahatten=p

Lcon=ldsc[which(ldsc[,1]=="Conserved_LindbladTohL2_1"),c(4,5,6,7,8)]
colnames(Lcon)[1]="p"

Lcon[,1]=1-pnorm(abs(Lcon[,1]))
Lcon[,1]=-log10(Lcon[,1])
rownames(Lcon)=1:nrow(Lcon)
colorannot=data.frame(type=grouporder,color=gg_color_hue(15))
Lcon$label=""
Lcon[which(Lcon$trait=="ID152"),"label"]="Impedance of left leg"
Lcon[which(Lcon$trait=="ID205"),"label"]="Fresh fruit intake"
Lcon[which(Lcon$trait=="ID229"),"label"]="Snoring"
Lcon[which(Lcon$trait=="ID249"),"label"]="Age first intercourse"
Lcon[which(Lcon$trait=="ID39"),"label"]="Forced vital capacity"
Lcon[which(Lcon$trait=="ID5"),"label"]="BMI"
Lcon[which(Lcon$trait=="ID631"),"label"]="Anorexia Nervosa"
Lcon[which(Lcon$trait=="ID92"),"label"]="Dentures"


Lcon$type=factor(Lcon$type,levels=grouporder)
p=ggplot(Lcon,aes(type,p,label=label,color=type))+geom_jitter(aes(group=type),width=0.4)
p=p+theme_classic()+theme(
	axis.text.x = element_text(angle = -60, vjust = 0.2,hjust=0,size=10,
		family="sans",color="black"),
	axis.ticks=element_blank(),
	axis.title.x=element_blank(),
	axis.text.y = element_text(size=10,family="sans",color="black"),
	axis.title.y = element_text(size=10,family="sans",color="black"),
	legend.position="none")
p=p+geom_hline(yintercept=-log10(0.05/873))+ylab("-log10(p for Con.region)")
p=p+geom_text_repel(size=3,segment.size = 0.1)
p_con_mahatten=p

p_f4_down=p_con_mahatten|p_pli_mahatten

Lcadd=ldsc[which(ldsc[,1]=="cadd_posL2_0"),c(4,5,6,7,8)]
colnames(Lcadd)[1]="p"

Lcadd[,1]=1-pnorm(abs(Lcadd[,1]))
Lcadd[,1]=-log10(Lcadd[,1])
rownames(Lcadd)=1:nrow(Lcadd)
colorannot=data.frame(type=grouporder,color=gg_color_hue(15))
Lcadd$label=""
Lcadd[which(Lcadd$trait=="out63"),"label"]="Large artery stroke"
Lcadd[which(Lcadd$trait=="ID225"),"label"]="Vegetable intake"
Lcadd[which(Lcadd$trait=="ID53"),"label"]="Ever drinker"
Lcadd[which(Lcadd$trait=="out35"),"label"]="Eczema"
Lcadd[which(Lcadd$trait=="ID445"),"label"]=""
Lcadd[which(Lcadd$trait=="ID7"),"label"]="Coronary artery disease"


Lcadd$type=factor(Lcadd$type,levels=grouporder)
p=ggplot(Lcadd,aes(type,p,label=label,color=type))+geom_jitter(aes(group=type),width=0.4)
p=p+theme_classic()+theme(
	axis.text.x = element_text(angle = -60, vjust = 0.2,hjust=0,size=10,
		family="sans",color="black"),
	axis.ticks=element_blank(),
	axis.title.x=element_blank(),
	axis.text.y = element_text(size=10,family="sans",color="black"),
	axis.title.y = element_text(size=10,family="sans",color="black"),
	legend.position="none")
p=p+geom_hline(yintercept=-log10(0.05/873))+ylab("-log10(p for CADD(+))")
p=p+geom_text_repel(size=3,segment.size = 0.1)
p_cadd_mahatten=p

#SLA
head(ldsc)
useanno
asmc=ldsc[which(ldsc[,1]=="PLIL2_1"),4]
pasmc=1-pnorm(abs(asmc))
exp=runif(length(pasmc),0,1)
ks.test(pasmc,"punif",0,1)$p.value
dsmana=data.frame(true=pasmc[order(pasmc)],exp=exp[order(exp)])
dsmana=-log10(dsmana)
p=ggplot(dsmana,aes(exp,true))+geom_point()+theme_classic()
p=p+annotate(geom="text", x=0,hjust=0, y=8.5, label="pLi>0.9\nKolmogorov-Smirnov\nD=0.52,p<e-100",size=3)
p=p+theme(axis.ticks=element_blank(),
		   axis.text=element_text(color="black",family="sans",size=10),
               axis.title=element_text(color="black",family="sans",size=10))
p_smana=p+xlab("-log10 expected p")+ylab("-log10 true p")

asmc=ldsc[which(ldsc[,1]=="cadd_posL2_0"),4]
pasmc=1-pnorm(abs(asmc))
exp=runif(length(pasmc),0,1)
ks.test(pasmc,"punif",0,1)$p.value
dsmanc=data.frame(true=pasmc[order(pasmc)],exp=exp[order(exp)])
dsmanc=-log10(dsmanc)
p=ggplot(dsmanc,aes(exp,true))+geom_point()+theme_classic()
p=p+annotate(geom="text", x=0,hjust=0, y=11, label="CADD(+)\nKolmogorov-Smirnov\nD=0.50,p<e-100",size=3)
p=p+theme(axis.ticks=element_blank(),
		   axis.text=element_text(color="black",family="sans",size=10),
               axis.title=element_text(color="black",family="sans",size=10))
p_smanc=p+xlab("-log10 expected p")+ylab("-log10 true p")

p_sman=(p_smana|p_pli_mahatten)/(p_smanc|p_cadd_mahatten)
pdf("sman.pdf",width=8,height=8)
p_sman
dev.off()


disp=c()
for(i in useanno){
a=length(which(ldsc[,1]==i&ldsc$disease==1&ldsc$sig==0))
b=length(which(ldsc[,1]==i&ldsc$disease==1&ldsc$sig==1))

c=length(which(ldsc[,1]==i&ldsc$disease==0&ldsc$sig==0))
d=length(which(ldsc[,1]==i&ldsc$disease==0&ldsc$sig==1))
e=matrix(c(b,a,d,c),nr=2)
model=fisher.test(e)

disp=rbind(disp,c(i,model$estimate,model$p.value,model$conf.int))
}
disp=data.frame(disp)
disp=disp[which(disp[,4]!=0),]
disp=disp[-9,]
disp[,2:5]=apply(disp[,2:5],2,as.numeric)
disp[,1]=c("Con.region","PLI","Rec.rate","ASMCavg","CADD(+)","LLD","soft(-)",
"soft(+)","CADD(-)")
disp=disp[order(disp[,2],decreasing=T),]
disp$annotation=factor(disp$annotation,levels=disp[9:1,1])
dldisa=disp
colnames(disp)=c("annotation","OR","p","low","up")
p=ggplot(disp,aes(OR,annotation))+geom_point(size=5)+theme_classic()
p=p+geom_errorbarh(aes(xmin=low,xmax=up),height=0)+scale_x_log10()
p=p+theme(axis.ticks=element_blank(),
		   axis.text=element_text(color="black",family="sans",size=10),
               axis.title=element_text(color="black",family="sans",size=10),
               axis.line.y=element_blank())
p_sldisa=p+geom_vline(xintercept=1)

dldisd=ldsc[which(ldsc[,1]=="MAF_Adj_ASMCL2_1"),c(4,7)]
colnames(dldisd)=c("z","dis")
dldisd$disease=ifelse(dldisd$dis==1,"yes","no")
oneway.test(z~dis,data=dldisd)
p=ggplot(dldisd,aes(z,group=disease,color=disease))+stat_ecdf(geom="step")
p=p+theme_classic()+theme(axis.ticks=element_blank(),legend.position="none",
		   axis.text=element_text(color="black",family="sans",size=10),
               axis.title=element_text(color="black",family="sans",size=10))
p=p+xlab("z for ASMCavg")+ylab("cumulative proportion")
p_sldisd=p+geom_vline(xintercept=-3.86)+annotate(geom="text", x=0.5,hjust=0.5, 
y=1, label="ASMCavg\npermutation p=0.07",size=3)+scale_x_reverse()
p_sldisd
p_sldisc=p_sldisc+theme(legend.position="none")
p_sldis=(p_sldisa|p_sldisb)/(p_sldisc|p_sldisd)
pdf("sldis.pdf",width=8,height=8)
p_sldis
dev.off()


##########
dfdaa=res_daf[,c("S","disease","time")]
dfdaa$onset=ifelse(dfdaa$disease==0,"trait",dfdaa$time)
dfdaa$S=-dfdaa$S
p=ggplot(dfdaa,aes(onset,S,color=onset,group=onset))+geom_boxplot()+theme_classic()
p_fdaa=p+theme(legend.position="none")

dfdab=dfdab[order(dfdab$S,decreasing=T),]
dfdab=dfdab[order(dfdab$S,decreasing=T),]
dfdab$trait=factor(dfdab$trait,levels=dfdab$trait[9:1])
p=ggplot(dfdab,aes(S,trait,color=onset))+geom_point(size=3)+theme_classic()
p=p+geom_errorbarh(aes(xmin=low,xmax=up),height=0)
p=p+theme(axis.ticks=element_blank(),
		   axis.text=element_text(color="black",family="sans",size=10),
               axis.title=element_text(color="black",family="sans",size=10),
               axis.line.y=element_blank(),legend.position="none")
p_fdab=p+geom_vline(xintercept=0)+scale_color_manual(values=gg_color_hue(4)[1:3])
p_fdaa|p_fdab

int=data.frame(ch2=up$r,DAF=up$DAF)
int=int[order(int$DAF),]
int$ch2=cumsum(int$ch2)/sum(int$ch2)
int1=int
int=data.frame(ch2=down$r,DAF=down$DAF)
int=int[order(int$DAF),]
int$ch2=cumsum(int$ch2)/sum(int$ch2)
auc=rbind(data.frame(int1,effect="up"),data.frame(int,effect="down"))
p=ggplot(auc,aes(DAF,ch2,color=effect))+geom_point(size=1)+theme_classic()
p=p+theme(axis.ticks=element_blank(),
	axis.text= element_text(size=10,family="sans",color="black"),
	axis.title= element_text(size=10,family="sans",color="black"),
	legend.position=c(0.8,0.2),legend.background=element_blank())
p=p+xlab("DAF")+ylab("cumulative % of variance explained")
p=p+geom_segment(x=0.5,y=0.5,xend=0.35,yend=0.5,color="black")
p=p+annotate("text",x=0.6,y=0.5,label="S=0.11")
p=p+annotate("text",x=0.25,y=0.85,label="Anorexia nervosa")
p=p+scale_y_continuous(expand = c(0, 0))
p=p+scale_color_manual(name = "DA effect", values=c("#BEBADA","#FB8072"),labels = c("trait decresing","trait increasing"))
p_fdac=p

dfdae=res_daf[,c("fold","disease","time")]
dfdae$onset=ifelse(dfdae$disease==0,"trait",dfdae$time)

oneway.test(fold~disease,data=dfdae)

p=ggplot(dfdae,aes(onset,fold,color=onset,group=onset))+geom_boxplot()+theme_classic()
p_fdae=p+theme(legend.position="none")

p_fda=(p_fdaa|p_fdab|p_fdae)/(p_fdac|p_da_vol)
pdf("fsda.pdf",height=6,width=8)
p_fda
dev.off()

y$group=annot[y[,1],"abr"]
y$disease=ifelse(y$disease==1,"yes","no")
y$group=factor(y$group,levels=grouporder)
p=ggplot(y,aes(group,z,group=group,color=group))
p=p+geom_boxplot(outlier.shape =NA)+geom_jitter(size=0.5)
p=p+theme_classic()+theme(legend.position="none",axis.ticks=element_blank(),
axis.title.x=element_blank(),axis.text=element_text(color="black",size=10,
family="sans"),axis.title.y=element_text(color="black",size=10,family="sans"))
p=p+ylab("|S|")
p_s_box=p

summary(neareast_sel_model)
R2(predict(neareast_sel_model,clud_adj),clud_adj$neareast)

x=resid(proneo_sel_model)
x=(x-mean(x))/sd(x)
y=data.frame(rownames(clud),annot[rownames(clud),"uniqTrait"],
predict(proneo_sel_model,clud_adj),clud_adj$proneo,x)
y=y[which(y[,1] %in% annot[,1]),]
colnames(y)=c("id","trait","predict","true","z")
y$disease=annot[y$id,"disease"]
y[order(y$z)[1:13],]
R2= function(preds,vals) {
  1 - (sum((vals - preds)^2) / sum((vals - mean(vals))^2))
}
summary(lm(true~predict,data=y))
cor.test(y$predict,y$true)
bj=c()
for(i in 1:1000){
set.seed(i)
k=sample(1:length(which(y$disease==0)),length(which(y$disease==1)),replace=F)
bj=c(bj,R2(y[k,"predict"],y[k,"true"]))
}

mean(bj)
1-pnorm((mean(bj)-summary(lm(true~predict,data=y[which(y$disease==1),]))$r.squared)/sd(bj))

bj=data.frame(ncf=bj)
p=ggplot(bj,aes(ncf))+geom_density()+theme_classic()
p

colnames(dssel)=c()
dssel=rbind(dssel,c(R2(y[which(y$disease==0),"predict"],y[which(y$disease==0),"true"]),
R2(y[which(y$disease==1),"predict"],y[which(y$disease==1),"true"]),
1-pnorm((mean(bj)-R2(y[which(y$disease==1),"predict"],y[which(y$disease==1),"true"]))/sd(bj))))

####
x=resid(step.model)
x=(x-mean(x))/sd(x)

y=data.frame(rownames(datuse),annot[rownames(datuse),"uniqTrait"],
predict(f_sel_model,datuse),datuse$ncf,x)
y=y[which(y[,1] %in% annot[,1]),]
colnames(y)=c("id","trait","predict","true","z")
y$disease=annot[y$id,"disease"]
write.csv(y,"f_sel_res.csv")

p=ggplot(y,aes(predict,true,color=z))+geom_point()+theme_classic()
p=p+scale_color_gradient2(low="#3288BD",mid="grey",high="#B2182B")
p=p+geom_smooth(method="lm",se=F)+ labs(color = "scaled\nresidual")
p=p+ annotate(geom="text", x=-1.9, hjust=0,y=-8, label="Intelligence",size=3)
p=p+ annotate(geom="text", x=-0.8, hjust=0,y=-6.9, label="Wholemeal or wholegrain bread",size=3)
p=p+ annotate(geom="text", x=1.43, hjust=1,y=7.3, label="Leg fat percentage",size=3)
p=p+ annotate(geom="text", x=-0.84, hjust=0.5,y=4.4, label="Mood swing",size=3)
p=p+theme(axis.ticks=element_blank(),
	axis.text= element_text(size=10,family="sans",color="black"),
	axis.title= element_text(size=10,family="sans",color="black"))
p_f_sel=p+xlab("predict z_ncf")+ylab("true z_ncf")
############
####
x=resid(neareast_sel_model)
x=(x-mean(x))/sd(x)
y=data.frame(rownames(clud),annot[rownames(clud),"uniqTrait"],
predict(neareast_sel_model,clud_adj),clud_adj$neareast,x)
y=y[which(y[,1] %in% annot[,1]),]
colnames(y)=c("id","trait","predict","true","z")
y$disease=annot[y$id,"disease"]
write.csv(y,"neareast_sel_res.csv")

p=ggplot(y,aes(predict,true,color=z))+geom_point()+theme_classic()
p=p+scale_color_gradient2(low="#3288BD",mid="grey",high="#B2182B")
p=p+geom_smooth(method="lm",se=F)+ labs(color = "scaled\nresidual")
p=p+ annotate(geom="text", x=-0.8, hjust=0,y=-5.46, label="Back pain",size=3)
p=p+ annotate(geom="text", x=-2.1, hjust=1,y=-5.8, label="Skin color",size=3)
p=p+ annotate(geom="text", x=0.17, hjust=1,y=3.5, label="urine KCR",size=3)
p=p+theme(axis.ticks=element_blank(),
	axis.text= element_text(size=10,family="sans",color="black"),
	axis.title= element_text(size=10,family="sans",color="black"))
p_near_sel=p+xlab("predict t_neareast")+ylab("true t_neareast")
############
####
x=resid(proneo_sel_model)
x=(x-mean(x))/sd(x)
y=data.frame(rownames(clud),annot[rownames(clud),"uniqTrait"],
predict(proneo_sel_model,clud_adj),clud_adj$proneo,x)
y=y[which(y[,1] %in% annot[,1]),]
colnames(y)=c("id","trait","predict","true","z")
y$disease=annot[y$id,"disease"]
write.csv(y,"proneo_sel_res.csv")

p=ggplot(y,aes(predict,true,color=z))+geom_point()+theme_classic()
p=p+scale_color_gradient2(low="#3288BD",mid="grey",high="#B2182B")
p=p+geom_smooth(method="lm",se=F)+ labs(color = "scaled\nresidual")
p=p+ annotate(geom="text", x=0.42, hjust=0.5,y=-6, label="Voice broke age",size=3)
p=p+ annotate(geom="text", x=-0.5, hjust=0,y=-5.4, label="Back pain",size=3)
p=p+ annotate(geom="text", x=0.5, hjust=0,y=-4.2, label="Intelligence",size=3)
p=p+ annotate(geom="text", x=-0.63, hjust=0.5,y=4.6, label="tiredness",size=3)
p=p+theme(axis.ticks=element_blank(),
	axis.text= element_text(size=10,family="sans",color="black"),
	axis.title= element_text(size=10,family="sans",color="black"))
p_proneo_sel=p+xlab("predict t_proneo")+ylab("true t_proneo")
############
p=ggplot(y,aes(predict,true,color=z))+geom_point()+theme_classic()
p=p+scale_color_gradient2(low="#3288BD",mid="grey",high="#B2182B")
p=p+geom_smooth(method="lm",se=F)+ labs(color = "scaled\nresidual")
p=p+ annotate(geom="text", x=-2.6, hjust=0.5,y=-12, label="Black hair",size=3)
p=p+ annotate(geom="text", x=-2, hjust=0,y=-11, label="Cheese intake",size=3)
p=p+ annotate(geom="text", x=0.78, hjust=0,y=-6.1, label="Intelligence",size=3)
p=p+ annotate(geom="text", x=2.4, hjust=0.5,y=16, label="Skin tanning",size=3)
p=p+theme(axis.ticks=element_blank(),
+ axis.text= element_text(size=10,family="sans",color="black"),
+ axis.title= element_text(size=10,family="sans",color="black"))
p_hg_sel=+xlab("predict t_%HG")+ylab("true t_%HG")
#########################################
####
x=resid(sds_sel_model)
x=(x-mean(x))/sd(x)
y=data.frame(rownames(clud),annot[rownames(clud),"uniqTrait"],
predict(sds_sel_model,clud_adj),clud_adj$sds,x)
y=y[which(y[,1] %in% annot[,1]),]
colnames(y)=c("id","trait","predict","true","z")
y$disease=annot[y$id,"disease"]
write.csv(y,"sds_sel_res.csv")

p=ggplot(y,aes(predict,true,color=z))+geom_point()+theme_classic()
p=p+scale_color_gradient2(low="#3288BD",mid="grey",high="#B2182B")
p=p+geom_smooth(method="lm",se=F)+ labs(color = "scaled\nresidual")
p=p+ annotate(geom="text", x=0.42, hjust=0.5,y=-6, label="Voice broke age",size=3)
p=p+ annotate(geom="text", x=-0.5, hjust=0,y=-5.4, label="Back pain",size=3)
p=p+ annotate(geom="text", x=0.5, hjust=0,y=-4.2, label="Intelligence",size=3)
p=p+ annotate(geom="text", x=-0.63, hjust=0.5,y=4.6, label="tiredness",size=3)
p=p+theme(axis.ticks=element_blank(),
	axis.text= element_text(size=10,family="sans",color="black"),
	axis.title= element_text(size=10,family="sans",color="black"))
p_sds_sel=p+xlab("predict rho_sds")+ylab("true rho_sds")
#########################################
p_f_sel=p_f_sel+scale_color_gradient2(low="#3288BD",mid="grey",high="#B2182B")
p_sds_sel=p_sds_sel+scale_color_gradient2(low="#3288BD",mid="grey",high="#B2182B",limits=c(-4,4))
p_near_sel=p_near_sel+scale_color_gradient2(low="#3288BD",mid="grey",high="#B2182B",limits=c(-4,4))
p_proneo_sel=p_proneo_sel+scale_color_gradient2(low="#3288BD",mid="grey",high="#B2182B",limits=c(-4,4))
p_hg_sel=p_hg_sel+scale_color_gradient2(low="#3288BD",mid="grey",high="#B2182B",limits=c(-4,4))
p_ssel=p_f_sel + p_sds_sel + p_near_sel + p_proneo_sel+p_hg_sel+guide_area() + 
  plot_layout(guides = 'collect')
pdf("fssel.pdf",width=8,height=5)
p_ssel
dev.off()

