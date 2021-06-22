#code for figure 5 and 6, the integration of all selection metrices.
load("annot.RData")
library(circlize)
library("ComplexHeatmap")
library(data.table)
s1=fread("S1.csv",data.table=F) #A manually curated meta file for all GWAS included.
rownames(s1)=s1[,1]
head(s1)
res_daf$fold=res_daf[,3]/(res_daf[,2]+res_daf[,3])
annot$fold=res_daf[rownames(annot),"fold"]
s1$prop=ifelse(is.na(s1$case),NA,s1$case/s1$N)
annot$ukb=s1[rownames(annot),"ukb"]
annot$meta=s1[rownames(annot),"meta"]
annot$prop=s1[rownames(annot),"prop"]

allmet=c("ncm","ncf","sds","neolithic","neareast","proneo","HG","cmsp","cmsn",
"hardp","hardn","softp","softn","b2","asmc","fold")
#quantitative
hm_quan=c()
r2_quan=c()
apply(annot,2,function(x){length(which(x==0))})

annot["ID618","h2"]=0.1465
annot["out43","h2"]=0.4
annot["ID627","h2"]=1
annot["ID631","h2"]=0.243

ggdensity(log10(annot[which(annot$quant==1),"N"]))#4.6
ggdensity(log10(annot[which(annot$quant==0),"N"]))#5.3

ggdensity(annot[which(annot$quant==0),"h2"])#0.1
ggdensity(annot[which(annot$quant==1),"h2"])

ggdensity(annot[which(annot$quant==1),"brain"])#2

ggdensity(annot[which(annot$quant==0),"prop"])#0.2


form="~h2+log10(N)+brain+meta"
for(metric in allmet){
formula=paste(metric,form,sep="")
d=annot[which(annot$quant==1&annot[,metric]!=0),]
d[,metric]=abs(d[,metric])
model=lm(formula,data=d)
r2=summary(model)$r.squared
int=summary(model)$coefficients
tstat=int[match(c("h2","log10(N)","brain","meta","lgc","int","chisq"),rownames(int)),3]
names(tstat)=c("h2","log10(N)","brain","meta","lgc","int","chisq")
hm_quan=rbind(hm_quan,tstat)
r2_quan=c(r2_quan,r2)
}

#qulitative
hm_qual=c()
r2_qual=c()

form="~h2+log10(N)+brain+meta+prop"
for(metric in allmet){
formula=paste(metric,form,sep="")
d=annot[which(annot$quant==0 &annot[,metric]!=0),]
d[,metric]=abs(d[,metric])
model=lm(formula,data=d)
r2=summary(model)$r.squared
tstat=summary(model)$coefficients[c("h2","log10(N)","brain","meta","prop"),3]
hm_qual=rbind(hm_qual,tstat)
r2_qual=c(r2_qual,r2)
}
rownames(hm_qual)=allmet
rownames(hm_quan)=allmet

names(r2_qual)=allmet
names(r2_quan)=allmet

rownames(hm_quan)[16]="%DA"
rownames(hm_qual)[16]="%DA"
###########Figure 5A
library(circlize)
library("ComplexHeatmap")
col_fun = colorRamp2(c(-10, 0, 10), c("#3288BD", "white","#B2182B" ))
ha = rowAnnotation(totalR2= anno_barplot(r2_quan),width = unit(0.7, "in"),
annotation_name_gp=gpar(fontsize=10))
colnames(hm_quan)[2]="lg(N)"
rownames(hm_quan)[7:14]=c("%HG","cms(+)","cms(-)","hardp(+)","hard(-)","soft(+)",
"soft(-)","Bal.sel")

pf5a_left=Heatmap(hm_quan[,1:4],col = col_fun,name="t value", cluster_rows = FALSE, 
cluster_columns = FALSE,heatmap_legend_param = list(at = c(-10,-5,0,5,10)),
row_names_gp = gpar(fontsize = 9),column_names_gp = gpar(fontsize = 9), 
row_names_side = "left",right_annotation=ha)

pdf("pf5a_left.pdf",height=3.8,width=3.1)
draw(pf5a_left)
dev.off()

ha = rowAnnotation(totalR2= anno_barplot(r2_qual),width = unit(0.7, "in"),
annotation_name_gp=gpar(fontsize=10))

rownames(hm_qual)[7:14]=c("%HG","cms(+)","cms(-)","hardp(+)","hard(-)","soft(+)",
"soft(-)","Bal.sel")

pf5a_right=Heatmap(hm_qual,col = col_fun,name="t value", cluster_rows = FALSE, 
cluster_columns = FALSE,heatmap_legend_param = list(at = c(-10,-5,0,5,10)),
row_names_gp = gpar(fontsize = 9),column_names_gp = gpar(fontsize = 9), 
row_names_side = "left",right_annotation=ha)

pdf("pf5a_right.pdf",height=3.8,width=3.6)
draw(pf5a_right)
dev.off()

############################

annot=annot[rownames(s1),]


dsc=dsc[rownames(annot),]
dsc[]=apply(dsc,2,function(x){
x[which(is.na(x))]=0
return(x)})

annot[]=apply(annot,2,function(x){
x[which(is.na(x))]=0
return(x)})

allmetric=colnames(annot)[c(10:13,17:37,41,48:76)]
dsc=dsc[,c(6:7,13,15,16,18,20,21,25,29,30,32,34,37,40,43,49,64,67,69,70,72,
76,78,80,82,84,86,88)]

annot1=cbind(annot,dsc)

d_adj=data.frame(sapply(allmetric,adj))
dim(d_adj)

d_f5b=c()
r_f5b=c()
coef_list=c("sds","neolithic","neareast","proneo","HG","cmsp","cmsn","hardp",
"hardn","softp","softn","asmc","b2","fold")
form="ncm~sds+neolithic+neareast+proneo+HG+cmsp+cmsn+hardp+hardn+softp+softn+asmc+b2+fold"
model=lm(form,data=d_adj[which(d_adj$ncm!=0),])
d_f5b=rbind(d_f5b,summary(model)$coefficients[coef_list,3])
r_f5b=c(r_f5b,summary(model)$r.squared)

form="ncf~sds+neolithic+neareast+proneo+HG+cmsp+cmsn+hardp+hardn+softp+softn+asmc+b2+fold"
model=lm(form,data=d_adj[which(d_adj$ncf!=0),])
d_f5b=rbind(d_f5b,summary(model)$coefficients[coef_list,3])
r_f5b=c(r_f5b,summary(model)$r.squared)

form="sds~neolithic+neareast+proneo+HG+cmsp+cmsn+hardp+hardn+softp+softn+asmc+b2+fold"
model=lm(form,data=d_adj[which(d_adj$sds!=0),])
coef=summary(model)$coefficients
int=coef[match(coef_list,rownames(coef)),3]
names(int)=coef_list
d_f5b=rbind(d_f5b,int)
r_f5b=c(r_f5b,summary(model)$r.squared)

form="neolithic~proneo+HG+cmsp+cmsn+hardp+hardn+softp+softn+asmc+b2+fold"
model=lm(form,data=d_adj[which(d_adj$neolithic!=0),])
coef=summary(model)$coefficients
int=coef[match(coef_list,rownames(coef)),3]
names(int)=coef_list
d_f5b=rbind(d_f5b,int)
r_f5b=c(r_f5b,summary(model)$r.squared)

form="neareast~cmsp+cmsn+hardp+hardn+softp+softn+asmc+b2+fold"
model=lm(form,data=d_adj[which(d_adj$neareast!=0),])
coef=summary(model)$coefficients
int=coef[match(coef_list,rownames(coef)),3]
names(int)=coef_list
d_f5b=rbind(d_f5b,int)
r_f5b=c(r_f5b,summary(model)$r.squared)

form="proneo~cmsp+cmsn+hardp+hardn+softp+softn+asmc+b2+fold"
model=lm(form,data=d_adj[which(d_adj$proneo!=0),])
coef=summary(model)$coefficients
int=coef[match(coef_list,rownames(coef)),3]
names(int)=coef_list
d_f5b=rbind(d_f5b,int)
r_f5b=c(r_f5b,summary(model)$r.squared)

form="HG~cmsp+cmsn+hardp+hardn+softp+softn+asmc+b2+fold"
model=lm(form,data=d_adj[which(d_adj$HG!=0),])
coef=summary(model)$coefficients
int=coef[match(coef_list,rownames(coef)),3]
names(int)=coef_list
d_f5b=rbind(d_f5b,int)
r_f5b=c(r_f5b,summary(model)$r.squared)

rownames(d_f5b)=c("ncm","ncf","sds","neolithic","neareast","proneo","%HG")
colnames(d_f5b)[c(5,13,14)]=c("%HG","Bal.sel","%DA")
d_f5b[]=apply(d_f5b,2,function(x){
x[which(is.na(x))]=0
x[which(x>10)]=10
x[which(x<(-10))]=-10
return(x)})

##################Figure 5B

ha = rowAnnotation(totalR2= anno_barplot(r_f5b),width = unit(1, "in"),
annotation_name_gp=gpar(fontsize=10))
colnames(d_f5b)[6:11]=c("cms(+)","cms(-)","hardp(+)","hard(-)","soft(+)",
"soft(-)")
pf5b=Heatmap(d_f5b,col = col_fun,name="t value", cluster_rows = FALSE, 
cluster_columns = FALSE,heatmap_legend_param = list(at = c(-10,-5,0,5,10)),
row_names_gp = gpar(fontsize = 9),column_names_gp = gpar(fontsize = 9), 
row_names_side = "left",right_annotation=ha)

pdf("pf5b.pdf",height=2,width=6)
draw(pf5b)
dev.off()

##############################
library(MASS)
allform=paste(colnames(d_adj)[c(10,11,17,19:22,27:55)],collapse="+")
conform=paste(colnames(d_adj)[c(10,11,17,19:22,47)],collapse="+")
caddform=paste(colnames(d_adj)[c(10,11)],collapse="+")

d_f6a=c()
r_fg=c()

for(i in allmet){
allmodel=lm(paste(i,allform,sep="~"),data=d_adj[which(d_adj[,50]!=0&d_adj[,i]!=0),])
r_all=summary(allmodel)$r.squared
d_f6a=rbind(d_f6a,summary(allmodel)$coefficients[,3])

conmodel=lm(paste(i,conform,sep="~"),data=d_adj[which(d_adj[,50]!=0&d_adj[,i]!=0),])
r_con=summary(conmodel)$r.squared


caddmodel=lm(paste(i,caddform,sep="~"),data=d_adj[which(d_adj[,50]!=0&d_adj[,i]!=0),])
r_cadd=summary(caddmodel)$r.squared

r_fg=rbind(r_fg,c(r_all,r_con,r_cadd))
}

allannot=sapply(form_fg,function(x){unlist(strsplit(x,"+",fixed=T))})
allannot=unique(unlist(allannot))


############Figure 6A
d_f6a=d_f6a[,-1]
rownames(d_f6a)=allmet
d_f6a=d_f6a[,c(1:7,27,8:26,28:35)]
colnames(d_f6a)
colnames(d_f6a)=c("CADD(+)","CADD(-)","Neanderthal","LLD","alleleAge","Nuc.diversity",
"pLi","Recombination","H3K27acQTL","H3K4me1QTL","Coding","Conserved","CpG",
"CTCFsite","DGF","DHS","Enhancer","fetalDHS","GERP","eQTL","H3K27ac","H3K4me",
"H3K4m3","H3K9ac","Intron","Promoter","flankPromoter","Repressed","superEnhancer",
"TFBS","Transcript","TSS","UTR3","UTR5","weakEnhancer")
d_f6a=d_f6a[,c(1:8,12,9:11,13:35)]

for(i in allmet){
index=length(allannot[[i]])
xx=abs(d_f6a[i,])
xx=t(t(xx))
xx=xx[order(xx[,1],decreasing=T),]
xx=names(xx)[1:index]
xx=setdiff(colnames(d_f6a),xx)
d_f6a[i,xx]=0
}

rownames(d_f6a)[c(7,16)]=c("%HG","%DA")

ha = rowAnnotation(totalR2= anno_barplot(cbind(r_fg[,3],r_fg[,2]-r_fg[,3],r_fg[,1]-r_fg[,2]),
gp = gpar(fill = col, col = col)),width = unit(1, "in"))
rownames(d_f6a)[8:14]=c("cms(+)","cms(-)","hard(+)","hard(-)","soft(+)",
"soft(-)","Bal.sel")
p_f6a=Heatmap(d_f6a,col = col_fun,name="t value", cluster_rows = FALSE, 
cluster_columns = FALSE,heatmap_legend_param = list(at = c(-10,0,10)),
border = TRUE,row_names_gp = gpar(fontsize = 9),column_split=factor(c(
rep("CADD",2),rep("conservation",7),rep("functional genomics",26))),column_names_gp = gpar(fontsize = 9,col=col), 
row_names_side = "left",right_annotation=ha,column_title_gp=gpar(fontsize=10,col=col))
pdf("f6a.pdf",height=3,width=7)
draw(p_f6a)
dev.off()

col=brewer.pal(4,"Set3")[c(1,3,4)]
col=c(col,"black")

library(ggplot2)

iduse=rownames(annot)[which(annot$sds!=0&annot$caddp!=0)]
df6d=data.frame(id=iduse,sds=annot[iduse,"sds"],cadd=annot[iduse,"caddp"])

p=ggplot(df6d,aes(cadd,sds))+geom_point(size=0.8)+theme_classic()
p=p+theme(axis.ticks=element_blank(),
          axis.text=element_text(color="black",family="sans",size=10),
          axis.title=element_text(color="black",family="sans",size=10))
p=p+xlab("ldsc Z_CADD(+)")+ylab("rho_SDS")
p_f6b=p

iduse=rownames(annot)[which(annot$asmc!=0&annot$age!=0)]
df6d=data.frame(id=iduse,sds=annot[iduse,"asmc"],cadd=annot[iduse,"CpG_Content_50kbL2_1"])
p=ggplot(df6d,aes(cadd,sds))+geom_point(size=0.8)+theme_classic()
p=p+theme(axis.ticks=element_blank(),
          axis.text=element_text(color="black",family="sans",size=10),
          axis.title=element_text(color="black",family="sans",size=10))
p=p+xlab("ldsc Z_CpG")+ylab("ldsc Z_ASMCavg")
p_f6c=p

iduse=rownames(annot)[which(annot$asmc!=0&annot$age!=0)]
df6d=data.frame(id=iduse,sds=annot[iduse,"asmc"],cadd=annot[iduse,"div"])
p=ggplot(df6d,aes(cadd,sds))+geom_point(size=0.8)+theme_classic()
p=p+theme(axis.ticks=element_blank(),
          axis.text=element_text(color="black",family="sans",size=10),
          axis.title=element_text(color="black",family="sans",size=10))
p=p+xlab("ldsc Z_Nucleotide.diversity")+ylab("ldsc Z_ASMCavg")
p_f6d=p

library(patchwork)
pdf("Figure 6.pdf",height=5.5,width=7)
plot_spacer()+grid::grid.grabExpr(draw(p_f6a))+(p_f6b|p_f6d|p_f6c)+plot_layout(height=c(0.01,2,1))
dev.off()


###############
iduse=rownames(annot)[which(annot$proneo!=0&annot$Enhancer_AnderssonL2_1!=0)]
df6d=data.frame(id=iduse,sds=annot[iduse,"proneo"],cadd=annot[iduse,"Enhancer_AnderssonL2_1"])
p=ggplot(df6d,aes(cadd,sds))+geom_point(size=0.8)+theme_classic()
p=p+theme(axis.ticks=element_blank(),
          axis.text=element_text(color="black",family="sans",size=10),
          axis.title=element_text(color="black",family="sans",size=10))
p=p+xlab("ldsc Z_Enhancer")+ylab("t for proneolithic")
p_sf13a=p

iduse=rownames(annot)[which(annot$hardp!=0&annot$BLUEPRINT_H3K27acQTL_MaxCPPL2_1!=0)]
df6d=data.frame(id=iduse,sds=annot[iduse,"hardp"],cadd=annot[iduse,"BLUEPRINT_H3K27acQTL_MaxCPPL2_1"])
p=ggplot(df6d,aes(cadd,sds))+geom_point(size=0.8)+theme_classic()
p=p+theme(axis.ticks=element_blank(),
          axis.text=element_text(color="black",family="sans",size=10),
          axis.title=element_text(color="black",family="sans",size=10))
p=p+xlab("ldsc Z_H3K27 QTL")+ylab("ldsc Z_hard(+)")
p_sf13b=p

iduse=rownames(annot)[which(annot$neolithic!=0&annot$Intron_UCSCL2_1!=0)]
df6d=data.frame(id=iduse,sds=annot[iduse,"hardp"],cadd=annot[iduse,"Intron_UCSCL2_1"])
p=ggplot(df6d,aes(cadd,sds))+geom_point(size=0.8)+theme_classic()
p=p+theme(axis.ticks=element_blank(),
          axis.text=element_text(color="black",family="sans",size=10),
          axis.title=element_text(color="black",family="sans",size=10))
p=p+xlab("ldsc Z_intron")+ylab("t for neolithic")
p_sf13c=p

pdf("fs13.pdf",height=3,width=9)
p_sf13a|p_sf13b|p_sf13c
dev.off()




