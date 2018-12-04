# Title     : SeuratPCA
# Objective : Running Seurat PCA on a list of dataset (without merging)
# Created by: chenlingantelope
# Created on: 11/30/18

args = commandArgs(trailingOnly=TRUE)
dataname = args[1]
# filename ='pbmc_large/temp'

library('Matrix')
library('RcppCNPy')

batchid = npyLoad(paste(dataname, '.batch.npy',sep=''), type="integer")
labels = npyLoad(paste(dataname, '.labels.npy',sep=''), type="integer")
data=paste(dataname, '.X.mtx',sep='')
X = readMM(data)
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
##############################################
# Seurat
##############################################
library("Seurat")
SeuratPreproc <- function(X,label,batchname,zero_cells,genenames=NA){
	X = as.matrix(X)
	if(is.na(genenames[1])==T){
		genenames = paste('gene',c(1:length(X[,1])),sep='_')
	}
	rownames(X) = genenames
	colnames(X) = paste(batchname,c(1:length(X[1,])),sep='_')
	X <- CreateSeuratObject(raw.data = X,discard=F,min.cells=0,min.genes=0)
	X@meta.data[, "batch"] <- batchname
	X@meta.data[,"labels"] <- label[as.numeric(sapply(strsplit(colnames(X@data),'_'),function(X){X[2]}))]
	X <- NormalizeData(object = X)
	X <- FindVariableGenes(X, do.plot = F, display.progress = F)
	X <- ScaleData(object = X)
	return(X)
}

seurat_data <- lapply(c(1:length(data)),function(i){
	X = SeuratPreproc(data[[i]], label=labels[[i]], batchname=i)
	return(X)
})

pc = lapply(c(1:length(seurat_data)),function(i){
    genes.use <- c()
    for (i in 1:length(seurat_data)) {
      genes.use <- c(genes.use, head(rownames(seurat_data[[i]]@hvg.info), 1000))
    }
    if(length(seurat_data)<=2){n_shared=0}else{n_shared=2}
    genes.use <- names(which(table(genes.use) > n_shared))
	X = seurat_data[[i]]
	X <- RunPCA(X, pc.genes = genes.use, do.print = FALSE)
	PC <- GetDimReduction(object = X, reduction.type = "pca", slot = "cell.embeddings")[,c(1:10)]
	write.table(PC,file=paste(dataname,'.',i,'.CCA.txt',sep=''),quote=F,col.names=F,row.names=F)
	return(PC)
})

genes = lapply(c(1:length(seurat_data)), function(i){
	X = seurat_data[[i]]
	geneinfo = X@hvg.info
	genename = read.table(paste(dataname, '.genenames.txt',sep=''))[,1]
	geneinfo = cbind(genename,geneinfo)
	write.csv(geneinfo,quote=F,file=paste(dataname,'.',i,'.hvg_info.csv',sep=''))
	return(X@hvg.info)
})
