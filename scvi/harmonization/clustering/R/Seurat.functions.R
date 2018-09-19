library("Seurat")

hvg_CCA <- function(data1,data2,ndim=10,plotting=F,discard=T,getlabels=T){
    data = list(data1,data2)
	hvg <- lapply(data,function(X){
		rownames(x = head(x = X@hvg.info, n = 2000))
	})
	hvg <- union(x = hvg[[1]], y = hvg[[2]])
	combined <- RunCCA(object = data[[1]], object2 = data[[2]], genes.use = hvg,num.cc = ndim)
	if(plotting==T){
		p1 <- DimPlot(object = combined, reduction.use = "cca", group.by = "batch", pt.size = 0.5,
		    do.return = TRUE)
		p2 <- VlnPlot(object = combined, features.plot = "CC1", group.by = "batch", do.return = TRUE)
		plot_grid(p1, p2)
	}
	combined <- CalcVarExpRatio(object = combined, reduction.type = "pca", grouping.var = "batch",
    dims.use = 1:ndim)
	if(discard==T){
		combined <- SubsetData(object = combined, subset.name = "var.ratio.pca", accept.low = 0.5)
	}
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

SeuratPreproc <- function(X,label,batchname,genenames=NA){
	X = as.matrix(X)
	if(is.na(genenames[1])==T){
		genenames = paste('gene',c(1:length(X[,1])),sep='_')
	}
	rownames(X) = genenames
	colnames(X) = paste(batchname,c(1:length(X[1,])),sep='_')
	X <- CreateSeuratObject(raw.data = X)
	X <- NormalizeData(object = X)
	X <- ScaleData(object = X)
	X <- FindVariableGenes(object = X, do.plot = FALSE)
	X@meta.data[, "batch"] <- batchname
	if(!is.na(label)){
		X@meta.data$label <- label
	}
	return(X)
}

