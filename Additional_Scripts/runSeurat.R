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
# sub_data2 = sample(c(1,0),68579,TRUE,c(.1,.9))
# data[[2]] = data[[2]][,sub_data2]
# labels[[2]]= labels[[2]][sub_data2]
# sub_data3 = sample(c(1,0),105868,TRUE,c(.1,.9))
# data[[3]] = data[[3]][,sub_data3]
# labels[[3]]= labels[[3]][sub_data3]
##############################################
# Seurat
##############################################
library("Seurat")

hvg_CCA <- function(data,ndim=10,plotting=F,getlabels=T,filter_genes=TRUE,ngenes=1000){
	if(filter_genes==TRUE){
		genes.use <- c()
		for (i in 1:length(data)) {
		  genes.use <- c(genes.use, head(rownames(data[[i]]@hvg.info), ngenes))
		}
		if(length(data)==2){n_shared=0}else{n_shared=2}
		genes.use <- names(which(table(genes.use) > n_shared))
	}else{
		genes.use = rownames(data[[1]]@data)
	}
	pc = lapply(c(1:length(data)),function(i){
		X = data[[i]]
		X <- RunPCA(X, pc.genes = genes.use, do.print = FALSE)
		PC <- GetDimReduction(object = X, reduction.type = "pca", slot = "cell.embeddings")[,c(1:10)]
		write.table(PC,file=paste(dataname,'.',i,'.CCA.txt',sep=''),quote=F,col.names=F,row.names=F)
		return(PC)
	})
	if(length(data)==2){
		combined <- RunCCA(object = data[[1]], object2 = data[[2]], genes.use=genes.use,num.cc = ndim)
	}else{
		combined <- RunMultiCCA(data, genes.use = genes.use,num.cc = ndim)
	}
	combined <- CalcVarExpRatio(object = combined, reduction.type = "pca", grouping.var = "batch",
    dims.use = 1:ndim)
	combined <- AlignSubspace(object = combined, reduction.type = "cca", grouping.var = "batch",
	    dims.align = 1:ndim)
    latent <- GetDimReduction(object = combined,
        reduction.type = "cca.aligned",
        slot = "cell.embeddings"
    )
    batch <- attributes(combined)$meta.data$batch
    if(getlabels==T){
	    labels <- attributes(combined)$meta.data$label
		return(list(latent,batch,labels,genes.use))
    }else{
    	return(list(latent,batch,genes.use))
    }
}

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

start = proc.time()
seurat_data <- lapply(c(1:length(data)),function(i){
	X = SeuratPreproc(data[[i]], label=labels[[i]], batchname=i)
	return(X)
})

genes = lapply(c(1,2), function(i){
	X = seurat_data[[i]]
	geneinfo = X@hvg.info
	genename = read.table(paste(dataname, '.genenames.txt',sep=''))[,1]
	geneinfo = cbind(genename,geneinfo)
	write.csv(geneinfo,quote=F,file=paste(dataname,'.',i,'.hvg_info.csv',sep=''))
	return(X@hvg.info)
})


combined = hvg_CCA(data=seurat_data,filter_genes=T)
end = proc.time()
print(paste('runtime for CCA =',end-start))

labels = npyLoad(paste(dataname, '.labels.npy',sep=''), type="integer")
print(sum(labels==combined[[3]]))

cells = do.call(c,lapply(seurat_data,function(X){colnames(X@data)}))
cells = cbind(sapply(strsplit(cells,'_'),function(X){X[1]}),sapply(strsplit(cells,'_'),function(X){X[2]}))
write.table(combined[[1]],file=paste(dataname,'.CCA.txt',sep=''),quote=F,col.names=F,row.names=F)
write.table(cells,file=paste(dataname,'.CCA.cells.txt',sep=''),quote=F,col.names=F,row.names=F)
write.table(sapply(strsplit(combined[[4]],'_'),function(X){X[2]}),file=paste(dataname,'.CCA.genes.txt',sep=''),quote=F,col.names=F,row.names=F)

