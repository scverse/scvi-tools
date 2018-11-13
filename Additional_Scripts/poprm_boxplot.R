library(reshape2)
library(ggplot2)

dataname = 'NoOverlap'
dataname = 'PopRemove'
data1 = read.csv('NoOverlap.res.txt',header=F,sep='\t',as.is=T)
data2 = read.csv('NoOverlap.acc.res.txt',header=F,sep='\t',as.is=T)

data1 = read.csv('PopRemove.res.txt',header=F,sep='\t',as.is=T)
data2 = read.csv('PopRemove.acc.res.txt',header=F,sep='\t',as.is=T)


acc_res = lapply(unique(as.character(data2[,2])),function(celltype){
	# celltype = 'CD4 T cells'
	temp = data2[data2[,2]==celltype,]
	pops = as.character(temp[1,c(12:20)])
	temp = temp[,c(3:11)]
	colnames(temp) = pops
	temp = melt(temp)
	temp = cbind(rep(c('Seurat','vae'),9),temp)
	colnames(temp)=c('method','type','value')
	temp$celltype = c(rep('removed',2),rep('other',16))
	temp$removed = rep(celltype,18)
	return(temp)
})

acc_res = do.call(rbind,acc_res)
acc_res = acc_res[acc_res$type!='Other',]
p <- ggplot(data=acc_res, aes(x=removed, y=value,fill=method)) +
geom_boxplot() + 
theme_minimal() + 
theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
ylab('Per Cell Type Classification Accuracy') + 
xlab('Cell Type')
ggsave(paste(dataname, ".PerClusterACC.pdf"), width = 6, height = 4)


BE_res = lapply(unique(as.character(data1[,2])),function(celltype){
	# celltype = 'CD4 T cells'
	temp = data1[data1[,2]==celltype,]
	pops = as.character(temp[1,c(16:23)])
	temp = temp[,c(3:10)]
	colnames(temp) = pops
	temp = melt(temp)
	temp = cbind(rep(c('Seurat','vae'),8),temp)
	colnames(temp)=c('method','type','value')
	temp$celltype = c(rep('removed',2),rep('other',14))
	temp$removed = rep(celltype,16)
	return(temp)
})
BE_res = do.call(rbind,BE_res)
BE_res = BE_res[BE_res$type!='Other',]

if (dataname=='PopRemove'){
p <- ggplot(data=BE_res[BE_res$celltype=='other',], aes(x=removed, y=value,fill=method)) +
geom_boxplot() + 
theme_minimal() + 
theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
ylab('Per Cell Type Batch Entropy Mixing') + 
xlab('Cell Type')
ggsave(paste(dataname, ".PerClusterBE.pdf"), width = 6, height = 4)
removed_BE = BE_res[BE_res$celltype=='removed',]
p <- ggplot(data=removed_BE, aes(x=method, y=value,fill=method)) +
geom_boxplot() + 
theme_minimal() + 
theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
ylab('Batch Entropy Mixing') + 
xlab('Method')
ggsave(paste(dataname, ".removedBE.pdf"), width = 3, height = 4)

}else{
p <- ggplot(data=BE_res, aes(x=removed, y=value,fill=method)) +
geom_boxplot() + 
theme_minimal() + 
theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
ylab('Per Cell Type Batch Entropy Mixing') + 
xlab('Cell Type')
ggsave(paste(dataname, ".PerClusterBE.pdf"), width = 6, height = 4)

}



data1 	
    res = [BE1] + BE2 + [res_knn[x] for x in list(res_knn.keys())[:5]]
    g.write('vae' + '\t' + rmCellTypes + ("\t%.4f" * 13 + "\t%s"*8 + "\n") % tuple(res+cell_type))

data2
    f.write('vae' + '\t' + rmCellTypes + ("\t%.4f" * 9 + "\t%s" * 9 + "\n") % tuple(acc + list(cell_type)))
