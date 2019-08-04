#This file contains R code for DEPECHE, FlowSOM and flowMeans. Phenograph was run on cytofKit GUI without command line.
#All .fcs files for subsampling tests had been moved to one directiory
library(doParallel)
cl <- makeCluster(2)
registerDoParallel(cl)
library("cytofCore")

depeche_res=list()
flowsom_res=list()
flowmeans_res=list()
for(study in c("l13","l32","sam","muscle","CC","colon")){   
  for( nCell in c("20k","40k")){
    depeche_int=foreach(i=1:5,.combine=cbind,.packages=c("cytofCore","DepecheR")) %dopar% {
      file=paste(study,"_",nCell,"_",i,".fcs",sep="")
      data=exprs(read.FCS(file))
      start=Sys.time()
      testDataDepeche <- depeche(data[,1:(ncol(data)-1)])
      end=Sys.time()
      time=difftime(end, start, units = "secs")
      int=data.frame(cell=1:nrow(data),cluster=testDataDepeche$clusterVector,label=as.integer(data[,"label"]))
      colnames(int)[1]=paste("i",as.numeric(time),sep=":")
      colnames(int)[2]=study
      colnames(int)[3]=nCell
      int
      }
      
    depeche_res=c(depeche_res,depeche_int)
    flowsom_int=foreach(i=1:5,.combine=cbind,.packages=c("cytofCore","FlowSOM")) %dopar% {
      file=paste(study,"_",nCell,"_",i,".fcs",sep="")
      data=exprs(read.FCS(file))
      start=Sys.time()
      out <- FlowSOM::ReadInput(data, transform = FALSE, scale = FALSE)
      out <- FlowSOM::BuildSOM(out, colsToUse = 1:(ncol(data)-1))
      out <- FlowSOM::BuildMST(out)
      labels_pre <- out$map$mapping[, 1]
      k <- 14 #Modified for each data set
      out <- FlowSOM::metaClustering_consensus(out$map$codes, k = k, seed = i) 
      #out <- FlowSOM::metaClustering(out$map$codes,"metaClustering_consensus",max=28) 
      labels <- out[labels_pre]
      end=Sys.time()
      time=difftime(end, start, units = "secs")
      int=data.frame(cell=1:nrow(data),cluster=labels,label=data[,"label"])
      colnames(int)[1]=paste("i",as.numeric(time),sep=":")
      colnames(int)[2]=study
      colnames(int)[3]=nCell
      int
      }
      
    flowsom_res=c(flowsom_res,flowsom_int)
    flowmeans_int=foreach(i=1:5,.combine=cbind,.packages=c("cytofCore","flowMeans")) %dopar% {
      file=paste(study,"_",nCell,"_",i,".fcs",sep="")
      data=exprs(read.FCS(file))
      start=Sys.time()
      out=flowMeans(data[,1:(ncol(data)-1)],Standardize=FALSE)
      end=Sys.time()
      time=difftime(end, start, units = "secs")
      int=data.frame(cell=1:nrow(data),cluster=out@Label,label=as.integer(sam[,"label"]))
      colnames(int)[1]=paste("i",as.numeric(time),sep=":")
      colnames(int)[2]=study
      colnames(int)[3]=nCell
      int
      }    
     
     flowmeans_res=c(flowmeans_res,flowmeans_int)
   }
}


 
