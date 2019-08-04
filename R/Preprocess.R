library("cytofCore")
#Original .fcs data were separate files for each manual cluster downloaded from flowRepository (colon data) our cytoBank (muscle & Cell Cycle),
#So the first step is to combine them and record the manual label as the last colomn. The following code was modified from Nolan lab's code (see
#cytoBank suggestion)
#Levine13dim, Levine32dim and Samusik01 data were obtained from flowRepository as processed by Weber et.al.
fcsFiles=list.files(pattern="*.fcs")

	begin=T
	for (file in fcsFiles){
            l=strsplit(file,"_",fixed=TRUE)[[1]]
            cluster=strsplit(l[3],".",fixed=TRUE)[[1]][1]  #extract label names
		if (begin){
			combinedData=read.FCS(file);
                  combinedLabel=rep(cluster,nrow(combinedData))
			begin=F
		} else {
			tmpData = exprs(read.FCS(file));
                  combinedData = rbind(combinedData,tmpData)
                  combinedLabel=c(combinedLabel,rep(cluster,nrow(tmpData)))
		}
	}
combinedData=data.frame(combinedData,	label=combinedLabel)
write.FCS(flowFrame(exprs=as.matrix(combinedData)),"muscle.fcs")  #file name was manually set for each data set
}

#Data transformation and random subsmpling
muscle=read.FCS("muscle.fcs")  #Repeated for each data set
muscle=exprs(muscle)
dim(muscle) #585133   26 
muscle[,1:25]=asinh(muscle[,1:25]/5)

cell_muscle_20k=1:20000 #Repeated for 5k, 10k, 20k, 40k, 60k & 80k when required
for (i in 1:5){
set.seed(i)
l=sample(nrow(muscle),20000,replace=FALSE)
d=muscle[l,]
d=data.frame(d,cell=l) #Record rownames for subsampled cells for internal evaluation
file=paste("muscle_","20k_",i,".fcs")
write.FCS(flowFrame(exprs=as.matrix(d)),file)
file=paste("muscle_","20k_",i,".csv") # Some tools require csv as input, see correponding files
write.csv(d,file)
cell_muscle_20k=cbind(cell_muscle_20k,l) # Recorded the selected cells 
}
