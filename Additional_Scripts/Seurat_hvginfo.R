args = commandArgs(trailingOnly=TRUE)
dataname = args[1]
# filename ='pbmc_large/temp'

library('Matrix')
library('RcppCNPy')

batchid = npyLoad(paste(dataname, '.batch.npy',sep=''), type="integer")
labels = npyLoad(paste(dataname, '.labels.npy',sep=''), type="integer")
data=paste(dataname, '.X.mtx',sep='')
X = readMM(data)
# genenames = read.csv('easy2/Macosko_Regev.genenames.txt',header=F,as.is=T)[,1]
data = lapply(unique(batchid),function(i){
	X=t(X[batchid==i,])
	genenames = paste('gene',c(1:length(X[,1])),sep='_')
	rownames(X) = genenames
	colnames(X) = paste(i,c(1:length(X[1,])),sep='_')
	return(X)
})
labels = lapply(unique(batchid),function(i){
	labels[batchid==i]		# colnames(X) = genenames
})
seurat_data <- lapply(c(1:length(data)),function(i){
	X = SeuratPreproc(data[[i]], label=labels[[i]], batchname=i)
	return(X)
})

genes = lapply(c(1,2), function(i){
	X = seurat_data[[i]]
	write.csv(X@hvg.info,quote=F,file=paste(dataname,'.',i,'.hvg_info.csv',sep=''))
	return(X@hvg.info)
})

