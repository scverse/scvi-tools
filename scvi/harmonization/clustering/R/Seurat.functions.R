library("Seurat")

hvg_CCA <- function(data,ndim=10,plotting=F,getlabels=T){
	combined <- RunCCA(object = data[[1]], object2 = data[[2]], genes.use = rownames(data[[1]]@data),num.cc = ndim)
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
		return(list(latent,batch,labels))
    }else{
    	return(list(latent,batch))
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
	X <- NormalizeData(object = X)
	X <- ScaleData(object = X)
	X@meta.data[, "batch"] <- batchname
	X@meta.data[,"labels"] <- label
	return(X)
}
