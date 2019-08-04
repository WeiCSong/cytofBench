#This file contains code for analysis and plotting of clustering resolution, corresponding to Figure 4 & S4.
#colon.fcs is the output of Xshift standalone.bat with whole colon data as input.
library(cytofCore)
olon=read.FCS("colon.fcs")
olon=exprs(olon)
head(olon)
D=data.frame(olon[,20:22],stringsAsFactors = FALSE)
D[,3]=as.numeric(D[,3])+1
m=max(D[,2:3])
pred=factor(D[,2],levels=1:m)
true=factor(D[,3],levels=1:m)
conf=confusionMatrix(pred,true)
d=as.data.frame.matrix(conf$table)
unique(olon[,"cluster_id"])
colSums(d)
d=d[,1:13] #Number of anual gating is 13
f=function(x){
x/sum(x)
}

d=t(apply(d,1,f))
library(pheatmap)

p=pheatmap(d,cluster_row=FALSE,cluster_col=FALSE) #Visual inspection
head(d)

g=function(x){
which(x==max(x))
}
assign=apply(d,1,g)
assign=assign[order(assign)]
d=d[names(assign),]
D[,2]=assign[as.character(D[,2])]
pred=factor(D[,2],levels=1:13)
true=factor(D[,3],levels=1:13)
conf=confusionMatrix(pred,true)
conf$overall
ari=ARI(pred,true)
nmi=NMI(pred,true)
fm=FMeasure(as.integer(pred),as.integer(true),silent=TRUE)
accense 
#ac:0.788,ari:0.655,nmi:0.605,f:0.823
#Column "sample" marks the sample origin of each cell.
olon[which(olon[,"sample"] %in% 1:8),"sample"]="early"
olon[which(olon[,"sample"] %in% 25:32),"sample"]="polyp"
olon[which(olon[,"sample"] %in% 9:16),"sample"]="late"
olon[which(olon[,"sample"] %in% 17:24),"sample"]="normal"
D=data.frame(olon[,20:22])
D[,3]=as.numeric(D[,3])+1
sample_cluster=as.data.frame.matrix(table(D[,c(1,3)]))
sample_cluster=apply(sample_cluster,2,f)

h=function(x){
c(which(x==max(x)),max(x))
}
ord=apply(sample_cluster,2,h)
sample_cluster=t(sample_cluster)
ord=t(ord)
head(ord)
ord=data.frame(type=assign[as.character(1:123)],ord)
ord[,2]=factor(ord[,2],levels=c(3,4,1,2))
colnames(ord)=c("type","sample","max")

d=d[with(ord, order(type,sample,-max)),]
ord_clu=with(ord, order(type,sample,-max))
colnames(d)[1:13]=c("CD4 T naive","CD4 Tcm","CD4 Tem","Th0","Th1","Th2","CD8 T naive","CD8 Tc0","CD8 Tc1","CD8 Tc2"
,"CD8 Tcm","CD8 Tem","Treg")
draw=data.frame(d[as.character(ord_clu),],
sample_cluster[ord_clu,c(3,4,1,2)])
#Figure 4A
tiff("refine_xshift.tiff",units="in",height=8,width=5,res=300)
pheatmap(draw,cluster_row=FALSE,cluster_col=FALSE,gaps_col=13)
dev.off()

#colon.RData is the output of cytofkit GUI and contains result for Phenograph with whole colon data as input.
load("colon.RData")
exp=analysis_results$expressionData
predict=analysis_results$clusterRes[[1]]
x=as.vector(exp[,ncol(exp)])
D=data.frame(cell=1:length(x),cluster=predict,label=x,stringsAsFactors = FALSE)
D[,3]=as.numeric(D[,3])+1
m=length(unique(D[,2]))
pred=factor(D[,2],levels=1:m)
true=factor(D[,3],levels=1:m)
conf=confusionMatrix(pred,true)
d=as.data.frame.matrix(conf$table)
unique(olon[,"cluster_id"])
colSums(d)
d=d[,1:13]
f=function(x){
x/sum(x)
}

d=t(apply(d,1,f))
library(pheatmap)
p=pheatmap(d,cluster_row=FALSE,cluster_col=FALSE)
head(d)

g=function(x){
which(x==max(x))
}
assign=apply(d,1,g)
D[,2]=assign[as.integer(D[,2])]
pred=factor(D[,2],levels=1:13)
true=factor(D[,3],levels=1:13)
conf=confusionMatrix(pred,true)
conf$overall
ari=ARI(pred,true)
nmi=NMI(pred,true)
fm=FMeasure(as.integer(pred),as.integer(true),silent=TRUE)
#ac:0.74,ari:0.59,nmi:0.55,fm:0.78
assign=assign[order(assign)]
d=d[names(assign),]

olon[which(olon[,"sample"] %in% 1:8),"sample"]="early"
olon[which(olon[,"sample"] %in% 25:32),"sample"]="polyp"
olon[which(olon[,"sample"] %in% 9:16),"sample"]="late"
olon[which(olon[,"sample"] %in% 17:24),"sample"]="normal"
D=data.frame(olon[,"sample"],predict)

sample_cluster=as.data.frame.matrix(table(D))
sample_cluster=apply(sample_cluster,2,f)

h=function(x){
c(which(x==max(x)),max(x))
}
ord=apply(sample_cluster,2,h)
sample_cluster=t(sample_cluster)
ord=t(ord)
head(ord)
ord=data.frame(type=assign[as.character(1:31)],ord)
ord[,2]=factor(ord[,2],levels=c(3,4,1,2))
colnames(ord)=c("type","sample","max")
ord_clu=with(ord, order(type,sample,-max))
colnames(d)[1:13]=c("CD4 T naive","CD4 Tcm","CD4 Tem","Th0","Th1","Th2","CD8 T naive","CD8 Tc0","CD8 Tc1","CD8 Tc2"
,"CD8 Tcm","CD8 Tem","Treg")
draw=data.frame(d[ord_clu,],
sample_cluster[ord_clu,c(3,4,1,2)])
#Figure S4A
tiff("refine_phenograph.tiff",units="in",height=8,width=5,res=300)
pheatmap(draw,cluster_row=FALSE,cluster_col=FALSE,gaps_col=13)
dev.off()

#int_l13 recorded clustering result of DEPECHE on whole levine13dim data.
D=int_l13
max(D[,2])
pred=factor(D[,2],levels=1:24)
true=factor(D[,3],levels=1:24)
conf=confusionMatrix(pred,true)
d=as.data.frame.matrix(conf$table)
d=t(d)
d=d[,1:6]
#F:0.56
f=function(x){
x/sum(x)
}
d=t(apply(d,1,f))
g=function(x){
which(x==max(x))
}
assign=apply(d,1,g)
assign=assign[order(assign)]
d=d[names(assign),]
D[,3]=assign[as.character(D[,3])]

library(pheatmap)
p=pheatmap(d,cluster_row=FALSE,cluster_col=FALSE)
rownames(d)=as.character(anno[rownames(d),1])
#Figure S4B
tiff("depeche_l13.tiff",units="in",height=8,width=5,res=300)
p
dev.off()
dev.new()

head(int_sam)
D=int_sam
max(D[,2])
pred=factor(D[,2],levels=1:24)
true=factor(D[,3],levels=1:24)
conf=confusionMatrix(pred,true)
d=as.data.frame.matrix(conf$table)
d=t(d)
d=d[,1:8]
#F:0.88
f=function(x){
x/sum(x)
}
d=t(apply(d,1,f))
g=function(x){
which(x==max(x))
}
d[11,8]=0.333333334
assign=apply(d,1,g)
assign=assign[order(assign)]
d=d[names(assign),]
D[,3]=assign[as.character(D[,3])]

library(pheatmap)
p=pheatmap(d,cluster_row=FALSE,cluster_col=FALSE)
rownames(d)=as.character(anno[rownames(d),1])
anno=read.csv("anno.csv",head=FALSE,row.name=1)
#Figure 4B
tiff("depeche_sam.tiff",units="in",height=8,width=5,res=300)
p
dev.off()
