library("Seurat")

hvg_CCA <- function(data,ndim=10,plotting=F,filter_genes = FALSE,ngenes=1000){
    if (filter_genes ==TRUE){
        genes.use <- c()
		for (i in 1:length(data)) {
		  genes.use <- c(genes.use, head(rownames(data[[i]]@hvg.info), ngenes))
		}
		if(length(data)==2){n_shared=0}else{n_shared=2}
		genes.use <- names(which(table(genes.use) > n_shared))
		for (i in 1:length(data)) {
		  genes.use <- genes.use[genes.use %in% rownames(data[[i]]@scale.data)]
		}
    }
    else{
        genes.use = rownames(data[[1]]@data)
    }
    combined <- RunCCA(object = data[[1]], object2 = data[[2]], genes.use = genes.use,num.cc = ndim)
	combined <- CalcVarExpRatio(object = combined, reduction.type = "pca", grouping.var = "batch",
    dims.use = 1:ndim)
	combined <- AlignSubspace(object = combined, reduction.type = "cca", grouping.var = "batch",
	    dims.align = 1:ndim)
    latent <- GetDimReduction(object = combined,
        reduction.type = "cca.aligned",
        slot = "cell.embeddings"
    )
    cells = do.call(c,lapply(data,function(X){colnames(X@data)}))
    batch <- sapply(strsplit(cells,'_'),function(X){X[1]})
    cells = sapply(strsplit(cells,'_'),function(X){X[2]})
    return(list(latent,genes.use,batch,cells))
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
    X <- FindVariableGenes(X, do.plot = F, display.progress = F)
	X <- ScaleData(object = X)
	X@meta.data[, "batch"] <- batchname
	X@meta.data[,"labels"] <- label
	return(X)
}

