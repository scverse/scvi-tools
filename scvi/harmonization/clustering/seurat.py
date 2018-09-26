import numpy as np
import rpy2.robjects as ro
import warnings
from rpy2.rinterface import RRuntimeWarning
import rpy2.robjects.numpy2ri as numpy2ri
from scipy.io import mmwrite
class SEURAT():
    def __init__(self):
        warnings.filterwarnings("ignore", category=RRuntimeWarning)
        numpy2ri.activate()
        r_source = ro.r['source']
        r_source("scvi/harmonization/clustering/R/Seurat.functions.R")
        ro.r["library"]("Matrix")
        ro.r["library"]("RcppCNPy")
        ro.r["library"]("reticulate")

    def create_seurat(self, dataset, batchname):
        genenames = dataset.gene_names
        genenames, uniq = np.unique(genenames,return_index=True)
        labels = [dataset.cell_types[int(i)] for i in np.concatenate(dataset.labels)]
        matrix = dataset.X[:,uniq]
        mmwrite('temp.mtx', matrix.T)
        ro.r('X <- readMM("temp.mtx")')
        ro.r.assign("batchname", batchname)
        ro.r.assign("genenames", ro.StrVector(genenames))
        ro.r.assign("labels", ro.StrVector(labels))
        ro.r('seurat'+str(batchname)+'<- SeuratPreproc(X,labels,batchname,genenames)')
        return 1
    def get_pcs(self):
        ro.r('seurat1 <- RunPCA(seurat1, pc.genes = seurat1@var.genes, do.print = FALSE)')
        ro.r('seurat2 <- RunPCA(seurat2, pc.genes = seurat2@var.genes, do.print = FALSE)')
        pc1 = ro.r('GetDimReduction(object = seurat1, reduction.type = "pca", slot = "cell.embeddings")[,c(1:10)]')
        pc2 = ro.r('GetDimReduction(object = seurat2, reduction.type = "pca", slot = "cell.embeddings")[,c(1:10)]')
        return pc1,pc2


    def get_cca(self):
        ro.r('combined <- hvg_CCA(list(seurat1,seurat2))')
        latent = ro.r('combined[[1]]')
        labels = ro.r('combined[[3]]')
        batch_indices = ro.r('combined[[2]]')
        cell_types,labels = np.unique(labels,return_inverse=True)
        return latent,batch_indices,labels,cell_types
