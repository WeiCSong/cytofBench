#This file contains R code for external & internal evaluations. Clustering results for DEPECHE, FlowSOM and flowMeans were stored in R as list,
#whereas clustering results from other tools were exported as .mat, .fcs or .csv. They were processed separately.

#External evaluations for DEPECHE, FlowSOM, flowMeans
library(doParallel)
cl <- makeCluster(2)
registerDoParallel(cl)
#assign predicted labels to true labels by hungarian algorithm#
assignALL=function(int){
  assigned=foreach(i=1:5,.combine=cbind,.packages=c("caret","clue")) %dopar% {
    n=3*i-1
    int[,n]=int[,n]
    m=max(int[,n:(n+1)])
    pred=factor(int[,n],levels=1:m)
    true=factor(int[,n+1],levels=1:m)
    conf=confusionMatrix(pred,true)
    d=as.data.frame.matrix(conf$table)
    assign=solve_LSAP(as(nrow(int)-d,"matrix"))
    assign=as.vector(assign)
    names(assign)=1:length(assign)
    assign=assign[which(rowSums(d)!=0)]
    d=data.frame(cell=int[,n-1],cluster=assign[as.integer(int[,n])],label=int[,n+1])
    colnames(d)[1]=colnames(int)[n-1]
    colnames(d)[2]=colnames(int)[n]
    colnames(d)[3]=colnames(int)[n+1]
    d
  }
assigned
}
l_depeche_assigned=lapply(depeche_res,assignALL)
l_flowsom_assigned=lapply(flowsom_res,assignALL)
l_flowmeans_assigned=lapply(flowmeans_res,assignALL)
#summary(l_phenograph_assigned)

#external evaluation
#Repeated for each tool
ee_depeche=data.frame(study=c(),nCell=c(),time=c(),runtime=c(),Accuracy=c(),
fmeasure=c(),nmi=c(),ari=c())
for(j in 1:5){
mod=(j-1)%%3
mod=2^mod
int=l_depeche_assigned[[j]]

INT=foreach(i=1:5,.combine=rbind,.packages=c("FlowSOM","aricode","caret","MixGHD")) %dopar% {
  n=3*i-1
  m=max(int[,n:(n+1)])

  pred=factor(int[,n],levels=1:m)
  true=factor(int[,n+1],levels=1:m)
  conf=confusionMatrix(pred,true)
  accuracy=conf$overall["Accuracy"]

  pred=factor(int[,n],levels=1:m)
  true=factor(int[,n+1],levels=1:m)
  ari=ARI(pred,true)
  nmi=NMI(pred,true)
  fm=FMeasure(as.integer(pred),as.integer(true),silent=TRUE)
  runtime=strsplit(colnames(int)[n-1],":",fixed=TRUE)[[1]][2]
  nCell=colnames(int)[n+1]
  study=colnames(int)[n]
  row=c(study=study,nCell=nCell,time=i,runtime=runtime,accuracy,fmeasure=fm,
  nmi=nmi,ari=ari)
  row
  }
ee_depeche=rbind(ee_depeche,INT)
}

#Internal evaluation
exp=list(l13,l32,sam,muscle,CC,colon) #Names of data sets were defined at pre-processing step 
ie_depeche=data.frame(study=c(),nCell=c(),time=c(),runtime=c(),calinski_harabasz=c(),
davies_bouldin=c(),xie_beni=c())
for(j in 1:12){ #Two sample sizes, six data sets, 12 results (length of list) in total.
  ne=ceiling(j/2)
  e=exp[ne][[1]]
  clures=depeche_res[j][[1]]

  int=foreach(nc=1:5,.combine=rbind,.packages="clusterCrit") %dopar% {
    sta=3*nc-2
    cellID=as.integer(clures[,sta])
    rt=strsplit(colnames(clures)[sta],":",fixed=TRUE)[[1]][2]
    nCell=colnames(clures)[sta+2]
    study=colnames(clures)[sta+1]
    lab=as.integer(clures[,sta+1])
    intIdx <- intCriteria(as(e[cellID,],"matrix"),lab,c("Calinski_Harabasz",
    "Davies_Bouldin","Xie_Beni"))
    row=c(study=study,nCell=nCell,time=nc,runtime=rt,unlist(intIdx))
    row
  }
  ie_depeche=rbind(ie_depeche,int)
}

#Evaluations for Accense
#The output of Accense consists of two .mat file: .result and sub. They were renamed as "**_*k_*.mat" and "**_*_*.mat".
library(R.matlab)
library(MixGHD)
library(aricode)
library(FlowSOM)
library(clue)
library(caret)
library(clusterCrit)
ie_accense=c()
ee_accense=c()
for(i in c("l13","l32","sam","muscle","CC","colon")){
  for(j in c(2,4)){
    for(k in 1:5){
      file=paste(i,"_",j,"0k_",k,".mat",sep="")
      file1=paste(i,"_",j,"0_",k,".mat",sep="")
      data=readMat(file)
      data1=readMat(file1)
      clures=foreach(i=1:length(data1$Subpop),.combine=rbind) %dopar% {
        L=unlist(data1$Subpop[[i]][[1]][4])
        int=data.frame(cell=L,cluster=rep(i,length(L)))
        int
      }
      exp=data$X.all
      rownames(clures)=clures[,1]
      clures=clures[order(clures[,1]),]
      clures=data.frame(clures,label=data1$rlabels[as.integer(clures[,1]),1])
      #assignment
      D=clures
      m=max(D[,2:3])
      pred=factor(D[,2],levels=1:m)
      true=factor(D[,3],levels=1:m)
      conf=confusionMatrix(pred,true)
      d=as.data.frame.matrix(conf$table)
      assign=solve_LSAP(as(20000*j-d,"matrix"))
      assign=as.vector(assign)
      names(assign)=1:length(assign)
      D[,2]=assign[as.integer(D[,2])]
      pred=factor(D[,2],levels=1:m)
      true=factor(D[,3],levels=1:m)
      conf=confusionMatrix(pred,true)
      ###ee
      ac=conf$overall[["Accuracy"]]
      ari=ARI(pred,true)
      nmi=NMI(pred,true)
      fm=FMeasure(as.integer(pred),as.integer(true),silent=TRUE)
      row=as.vector(c(study=i,nCell=paste(j,"0k",sep=""),time=k,runtime="runtime",ac,fm, #runtimes were recorded in a separate .txt file
      nmi,ari))
      ee_accense=rbind(ee_accense,row)
      #ie
    intIdx <- intCriteria(as(exp,"matrix"),clures[,2],c("Calinski_Harabasz",
    "Davies_Bouldin","Xie_Beni"))
    row=c(study=i,nCell=paste(j,"0k",sep=""),time=k,runtime=rt,unlist(intIdx))
    ie_accense=rbind(ie_accense,as.vector(row))
    }
  }
}

#Evaluations for ACDC
ie_acdc=c()
ee_acdc=c()
for study in c("l13","l32","sam","muscle","CC","colon")){
    for (i in c(2,4)){
      for(j in 1:5){
        file=paste(study,"_",i,"0k_",j,".mat",sep="")
        md=readMat(file)
        exp=md$X
        #assign
        int=data.frame(as.integer(md$y.pred.oracle),as.integer(md$y.true))
        int=int+1
        m=max(int)
        pred=factor(int[,1],levels=1:m)
        true=factor(int[,2],levels=1:m)
        conf=confusionMatrix(pred,true)
        d=as.data.frame.matrix(conf$table)
        assign=solve_LSAP(as(nrow(int)-d,"matrix"))
        assign=as.vector(assign)
        names(assign)=1:length(assign)
        int[,1]=assign[as.integer(int[,1])]
        #ee
        pred=factor(int[,1],levels=1:m)
        true=factor(int[,2],levels=1:m)
        conf=confusionMatrix(pred,true)
        accuracy=conf$overall["Accuracy"]
        ari=ARI(pred,true)
        nmi=NMI(pred,true)
        fm=FMeasure(as.integer(pred),as.integer(true),silent=TRUE)
        row=as.vector(c(nCell=10000*i,time=j,runtime="runtime",accuracy,fmeasure=fm,
        nmi=nmi,ari=ari))
        ee_acdc=rbind(ee_acdc,row)
        intIdx <- intCriteria(as(exp,"matrix"),int[,2],c("Calinski_Harabasz",
        "Davies_Bouldin","Xie_Beni"))
        row=c(study=study,nCell=paste(i,"0k",sep=""),time=j,runtime=rt,unlist(intIdx))
        ie_acdc=rbind(ie_acdc,as.vector(row))
      }
    }
  }
