#analysis of neolithic prs result from ancient.r

library(ggplot2)
library(ggrepel)
head(d_man)
d_man=int
d_man=d_man[!is.na(d_man$trait),]
d_man$p=-log10(d_man$p)
d_man$type=factor(d_man$type,levels=grouporder)
d_man$label=""
ll=which(d_man$ID %in% id1 & d_man$coef=="hg")
d_man[ll,"label"]=d_man[ll,"trait"]
d_man[679,"label"]="Skin tanning"
d_man[843,"label"]="Time watching TV"
d_man[847,"label"]="raw vegetable intake"
d_man[1031,"label"]=""
d_man[1459,"label"]="N.children"
d_man[1931,"label"]="Heavy drinking"
d_man[1315,"label"]="Shift work"
d_man[1931,"label"]="Heavy drinking"
d_man[2541,"label"]="Back pain"
d_man[2029,"label"]=""
d_man[1459,"label"]=""
d_man[ll,]
d_man=d_man[which(d_man$ID %in% annot[,1]),]
p=ggplot(d_man,aes(type,p,label=label,color=type))+geom_jitter(aes(group=type,shape=coef),width=0.4)
p=p+theme_classic()+theme(
	axis.text.x = element_text(angle = -60, vjust = 0.2,hjust=0,size=10,
		family="sans",color="black"),
	axis.ticks=element_blank(),
	axis.title.x=element_blank(),
	axis.text.y = element_text(size=10,family="sans",color="black"),
	axis.title.y = element_text(size=10,family="sans",color="black"),
	legend.position=c(0.7,0.8))
p=p+scale_shape_manual(values=c(16,4),
	labels=c("year to present (3 datasets)","% of hunter-gatherer ancestry"),
	name="tested coefficient")+guides(color=F)
p=p+geom_hline(yintercept=-log10(0.05/3600))+ylab("-log10(p in linear regression)")
p=p+geom_text_repel(size=3,segment.size = 0.1)
p_anc_mahatten=p

neo=read.csv("183_neolithic.csv")
p=ggplot(N,aes(HG,f))+geom_point(color=gg_color_hue(15)[1],shape=4)+theme_classic()
p=p+theme(axis.ticks=element_blank(),
          axis.text= element_text(size=10,family="sans",color="black"),
	    axis.title= element_text(size=10,family="sans",color="black"))
p=p+xlab("%HG")+ylab("scaled polygenic burden")
p_f3b=p+ annotate(geom="text", x=10,hjust=0, y=0.6, label="Ease of skin tanning",size=3)

head(d_cor)
summary(abs(d_cor[which(d_cor$disease=="dis"),4]))
d_cor$abs=abs(d_cor[,4])
oneway.test(abs~disease,data=d_cor)

near=read.csv("183_neareast.csv")
near$latitude=ifelse(near$lat>50,"high","low")
p=ggplot(near[which(near$date<12000),],aes(date,f,shape=latitude))+geom_point(color=gg_color_hue(15)[1])+theme_classic()
p=p+theme(axis.ticks=element_blank(),
          axis.text= element_text(size=10,family="sans",color="black"),
	    axis.title= element_text(size=10,family="sans",color="black"),
          legend.position=c(0.25,0.6))
p=p+xlab("year to present (neareast dataset)")+ylab("scaled polygenic burden")
p=p+scale_shape_manual(values=c(16,1))
p_f3c=p+scale_x_reverse()

p_anc_mahatten=p_anc_mahatten+scale_x_discrete(labels=xx)
p_f3_left=p_anc_mahatten/(p_f3b|p_f3c)
pdf("f3_left.pdf",height=6,width=6)
p_f3_left
dev.off()
d_man$var=ifelse(d_man$coef=="hg","hg",d_man$set)
d_cor=spread(d_man[,c("ID","disease","type","t","var")],var,t)
d_cor[,5:7]=-d_cor[,5:7]
d_p=spread(d_man[,c("ID","disease","type","p","var")],var,p)

anc_dig=d_p[which(d_p$proneo>1.3|d_p$hg>1.3),1]
anc_dig=intersect(anc_dig,annot[,1])
sig_d_cor=d_cor[anc_dig,]

head(sig_d_cor)
sig_d_cor$sign=apply(sig_d_cor[,4:7],1,median)
d_cor$sig=""
d_cor[-which(rownames(d_cor) %in% rownames(sig_d_cor)),"sig"]="insig"
d_cor[rownames(sig_d_cor),"sig"]="insig"

head(d_cor)

cor.test(d_cor[,4],d_cor$neolithic)$p.value
d_cor=d_cor[which(d_cor$ID %in% annot[,1]),]
library(GGally)
colnames(d_cor)[4]="%HG"
d_cor$disease=ifelse(d_cor$disease=="disease","dis","trait")
d_cor$disease=factor(d_cor$disease,levels=c("dis","trait"))
p_anc_cor=ggpairs(d_cor, 
     columns = 4:7, 
     aes(color=disease,alpha=0.2),upper = list(continuous = wrap("cor", size=3)))+
               theme_classic()+theme(axis.ticks=element_blank(),
		   axis.text=element_text(color="black",family="sans",size=10),
               axis.title=element_text(color="black",family="sans",size=10))+
               xlab("t in linear regression")+ylab("t in linear regression")

pdf("f3_right.pdf",height=5.5,width=5.5)
p_anc_cor
dev.off()
