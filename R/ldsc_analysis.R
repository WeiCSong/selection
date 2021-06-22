head(ldsc) #commands for running ldsc and cleaning of ldsc output are standard and trivial.
ldsc=ldsc[which(ldsc$trait %in% annot[,1]),]
ldsc=read.csv("ldsc.csv")
useanno=unique(ldsc[,1])[c(18,21,23,25,32,43,50,68,70,88,2,3,5,6,10,11,14,16,1)]
useanno=useanno[-1]
library(ComplexHeatmap)
df4a=read.csv("df4a.csv",row.names=1,encoding="UTF-8")
head(df4a)
ldsc$p=1-pnorm(abs(ldsc[,4]))
ldsc$sig=ifelse(ldsc$p<0.05/(871),1,0)
table(ldsc$sig)
countsig=as.data.frame.matrix(table(ldsc[which(ldsc[,1] %in% useanno),c(1,11)]))
rownames(df4a)
rownames(df4a)=c("Con.region","ASMCavg","LLD_AFR","Rec.Rate","B statistic",
"PLI","AlleleAge","Neanderthal","hardsweep(+)","softsweep(+)","Nuc.Diversity",
"hardsweep(-)","cms(-)","softsweep(-)","cms(+)","cadd(+)","cadd(-)","Bal.sel")

match(rownames(countsig),rownames(df4a))
col_fun = colorRamp2(c(-4, 0, 4), c("#3288BD", "white","#B2182B" ))
ha = HeatmapAnnotation(which="row",N.trait.significant= anno_barplot(countsig[rownames(df4a),2]),
annotation_name_gp=gpar(fontsize=10),width = unit(1.5, "in"))
p_f4a=Heatmap(df4a,col = col_fun,name="median z", cluster_rows = FALSE, 
cluster_columns = FALSE,heatmap_legend_param = 
list(at = c(-4,0,4)),column_split = c(rep(1,15),0,0),column_title= NULL,
border = TRUE,right_annotation=ha,row_names_side="left",
row_names_gp = gpar(fontsize = 9),column_names_gp = gpar(fontsize =9))

gb_f4a = grid.grabExpr(draw(p_f4a))




asmc=ldsc[which(ldsc[,1]=="MAF_Adj_ASMCL2_1"),4]
pasmc=1-pnorm(abs(asmc))
exp=runif(length(pasmc),0,1)
ks.test(pasmc,"punif",0,1)$p.value
df4b=data.frame(true=pasmc[order(pasmc)],exp=exp[order(exp)])
df4b=-log10(df4b)
p=ggplot(df4b,aes(exp,true))+geom_point()+theme_classic()
p=p+annotate(geom="text", x=0,hjust=0, y=12, label="ASMCavg\nKolmogorov-Smirnov\nD=0.54,p<e-100",size=3)
p=p+theme(axis.ticks=element_blank(),
		   axis.text=element_text(color="black",family="sans",size=10),
               axis.title=element_text(color="black",family="sans",size=10))
p_f4b=p+xlab("-log10 expected p")+ylab("-log10 true p")

Lcadd=ldsc[which(ldsc[,1]=="MAF_Adj_ASMCL2_1"),c(4,5,6,7,8)]
colnames(Lcadd)[1]="p"

Lcadd[,1]=1-pnorm(abs(Lcadd[,1]))
Lcadd$p=-log10(Lcadd[,1])
rownames(Lcadd)=1:nrow(Lcadd)
colorannot=data.frame(type=grouporder,color=gg_color_hue(15))
Lcadd$label=""
Lcadd[which(Lcadd$trait=="ID5"),"type"]="IMP"
Lcadd[which(Lcadd$trait=="ID148"),"label"]="Body water mass"
Lcadd[which(Lcadd$trait=="ID123"),"label"]="Reaction time"
Lcadd[which(Lcadd$trait=="ID313"),"label"]="Intelligence"
Lcadd[which(Lcadd$trait=="ID503"),"label"]="Schizophrenia"
Lcadd[which(Lcadd$trait=="ID187"),"label"]="Vigorous activity"
Lcadd[which(Lcadd$trait=="ID100"),"label"]="Health rating"
Lcadd[which(Lcadd$trait=="ID60"),"label"]="Processed meat intake"
Lcadd[which(Lcadd$trait=="ID319"),"label"]="Age at menarche"


Lcadd$type=factor(Lcadd$type,levels=grouporder)
pos <- position_jitter(width = 0.4, seed = 1)
p=ggplot(Lcadd,aes(type,p,label=label,color=type))+geom_jitter(aes(group=type),position = pos)
p=p+theme_classic()+theme(
	axis.text.x = element_text(angle = -60, vjust = 0.2,hjust=0,size=10,
		family="sans",color="black"),
	axis.ticks=element_blank(),
	axis.title.x=element_blank(),
	axis.text.y = element_text(size=10,family="sans",color="black"),
	axis.title.y = element_text(size=10,family="sans",color="black"),
	legend.position="none")
p=p+geom_hline(yintercept=-log10(0.05/873))+ylab("-log10(p for ASMCavg)")
p=p+geom_text_repel(size=3,segment.size = 0.1, position = pos)
p_asmc_mahatten=p
p_asmc_mahatten=p_asmc_mahatten+scale_x_discrete(labels=xx)
