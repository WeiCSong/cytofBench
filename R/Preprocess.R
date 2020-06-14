library("cytofCore")
#Original .fcs data were separate files for each manual cluster downloaded from flowRepository (colon data) our cytoBank (muscle & Cell Cycle),
#So the first step is to combine them and record the manual label as the last colomn. The following code was modified from Nolan lab's code (see
#cytoBank suggestion)
#Levine13dim, Levine32dim and Samusik01 data were obtained from flowRepository as processed by Weber et.al.

c=c(16,19:52) #Markers used for clustering, modified for each dataset. 
fcsFiles=list.files(pattern="*.fcs")
	begin=T
	for (file in fcsFiles){
            l=strsplit(file,"_",fixed=TRUE)[[1]]
            cluster=strsplit(l[4],".",fixed=TRUE)[[1]][1]  #extract label names
            cell=l[2]
            time=l[3]           		
            if (begin){
			combinedData=exprs(read.FCS(file))[,c]
                  combinedcluster=rep(cluster,nrow(combinedData))
                  combinedcell=rep(cell,nrow(combinedData))
                  combinedtime=rep(time,nrow(combinedData))
			begin=F
		} else {
			tmpData = exprs(read.FCS(file))[,c]
                  combinedData = rbind(combinedData,tmpData)
                  combinedcluster=c(combinedcluster,rep(cluster,nrow(tmpData)))
                  combinedcell=c(combinedcell,rep(cell,nrow(tmpData)))
                  combinedtime=c(combinedtime,rep(time,nrow(tmpData)))	}
	}
combinedData=data.frame(combinedData,cluster=combinedcluster,cell=combinedcell,time=combinedtime)
write.FCS(flowFrame(exprs=as.matrix(combinedData)),"CC.fcs")  #file name was manually set for each data set


#Data transformation and random subsmpling
CC=read.FCS("CC.fcs")  #Repeated for each data set
CC=exprs(CC)
dim(CC) #81594    38 
CC[,1:35]=asinh(CC[,1:35]/5)

cell_CC_20k=1:20000 #Repeated for 5k, 10k, 20k, 40k, 60k & 80k when required
for (i in 1:5){
set.seed(i)
l=sample(nrow(CC),20000,replace=FALSE)
d=CC[l,]
d=data.frame(d,cell=l) #Record rownames for subsampled cells for internal evaluation
file=paste("CC_","20k_",i,".fcs")
write.FCS(flowFrame(exprs=as.matrix(d)),file)
file=paste("CC_","20k_",i,".csv") # Some tools require csv as input, see correponding files
write.csv(d,file)
cell_CC_20k=cbind(cell_CC_20k,l) # Recorded the selected cells 
}
