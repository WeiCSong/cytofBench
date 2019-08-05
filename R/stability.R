#This file contains code for Figure 2. The input data 40k.csv and 20k.csv are the concatenation of the output of evaluation.R. The columns are:
#tool	data	time	runtime accuracy	f	nmi	ari	CH	DB	XB
#Left panel of Figure 2A&B was generated in PowerPoint.
data_40k=read.csv("40k.csv",head=TRUE)
data_20k=read.csv("20k.csv",head=TRUE)
library(tidyr)
#wide to long
data_40k <- gather(data_40k, index, value, runtime:XB, factor_key=TRUE)
data_20k <- gather(data_20k, index, value, runtime:XB, factor_key=TRUE)
library(dplyr)
#calculate sd and mean
Data_40k=data_40k %>% 
   group_by(.dots=c("tool","data","index")) %>% summarize(sd=sd(value),mean=mean(value))
Data_20k=data_20k %>% 
   group_by(.dots=c("tool","data","index")) %>% summarize(sd=sd(value),mean=mean(value))
Data_40k=data.frame(Data_40k)
Data_20k=data.frame(Data_20k)
#CV
draw=data.frame(Data_40k[,c("tool","data","index")],CV=Data_40k[,"sd"]/Data_40k[,"mean"])
library(ggplot2)
library(ggridges)
draw$tool=factor(draw$tool,level=c("kmeans","accense","DEPECHE","xshift",
"phenograph","flowmean","flowsom","ACDC","LDA"))
#figure 2A
p=ggplot(draw[which(draw$index %in% c("accuracy",	"f",	"nmi",	"ari")),], aes(x=CV, y=tool, fill=tool)) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_manual(
    values = c(ggplotColours(n = 7),"grey","black"))
p=p+theme_bw() +theme(panel.border = element_blank(),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
p=p+xlab("CV_for_external_evaluations")
p=p+geom_vline(xintercept=0.001,color="brown")
p
p_ex=p

#figure S2A
p=ggplot(draw[which(draw$index %in% c("CH",	"DB",	"XB")),], aes(x=CV, y=tool, fill=tool)) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_manual(
    values = c(ggplotColours(n = 7),"grey","black"))
p=p+theme_bw() +theme(panel.border = element_blank(),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
p=p+xlab("CV_for_internal_evaluations")
p=p+geom_vline(xintercept=0.001,color="brown")
p
p_in=p

#RD
draw=data.frame(Data_40k[,c("tool","data","index")],RD=(Data_40k[which(Data_40k$tool!="flowmean"),"mean"]-Data_20k[,"mean"])/Data_20k[,"mean"]
#Figure 2B
draw$tool=factor(draw$tool,level=c("kmeans","accense","DEPECHE","xshift",
"phenograph","flowsom","ACDC","LDA"))
p=ggplot(draw[which(draw$index %in% c("accuracy",	"f",	"nmi",	"ari")),], aes(x=RD, y=tool, fill=tool)) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_manual(
    values = c(ggplotColours(n = 7),"grey","black"))
p=p+theme_bw() +theme(panel.border = element_blank(),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
p=p+xlab("RD_for_external_evaluations")
p=p+geom_vline(xintercept=0,color="brown")
p
p_er=p

#Figure S2B

p=ggplot(draw[which(draw$index %in% c("CH",	"DB",	"XB")),], aes(x=RD, y=tool, fill=tool)) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_manual(
    values = c(ggplotColours(n = 7),"grey","black"))
p=p+theme_bw() +theme(panel.border = element_blank(),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
p=p+xlab("RD_for_internal_evaluations")
p=p+geom_vline(xintercept=0,color="brown")
p
p_ir=p
