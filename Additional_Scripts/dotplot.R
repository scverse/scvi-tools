library(reshape2)
library(ggplot2)

prop = read.table('celltypeprop.txt',sep='\t',row.names=NULL,as.is=T)
ncelltypes = length(prop[1,])-1
colnames(prop) = c(colnames(prop)[1:ncelltypes+1])
prop = prop[,1:ncelltypes]

others = read.csv('others.percluster.res.txt',sep='\t',row.names=NULL,as.is=T)
colnames(others) = c(colnames(others)[2:length(colnames(others))],'NA')
others = others[,1:length(others[1,])-1]

scvi = read.csv('scvi.percluster.res.txt',sep='\t',row.names=NULL,as.is=T)
colnames(scvi) = c(colnames(scvi)[2:length(colnames(scvi))],'NA')
scvi = scvi[,1:length(scvi[1,])-1]

Dotplot(others,scvi,'p1',c('SCMAP','SCANVI'),'percluster_scanvi_scmap1.pdf')
Dotplot(others,scvi,'p2',c('SCMAP','SCANVI'),'percluster_scanvi_scmap2.pdf')
Dotplot(others,scvi,'p1',c('SCMAP'),'percluster_scmap1.pdf')
Dotplot(others,scvi,'p2',c('SCMAP'),'percluster_scmap2.pdf')
Dotplot(others,scvi,'all',c('CCA','SCANVI'),'percluster_scanvi_cca.pdf')
Dotplot(others,scvi,'all',c('CCA'),'percluster_cca.pdf')


Dotplot <- function(others,scvi,ann,methods,plotname){
	scvi = scvi[scvi[,2]==ann,]
	others = others[others[,2]==ann,]
	if(ann=='p1'){prop_values = as.numeric(prop[2,])
		}else if(ann=='p2'){prop_values = as.numeric(prop[3,])
			}else{prop_values = as.numeric(prop[1,])}

	celltypes = colnames(prop)

	res = lapply(celltypes,function(celltype){
		temp =  c(others[,colnames(others)==celltype],
			scvi[,colnames(scvi)==celltype]
		)
		return(temp)
	})
	res = do.call(cbind,res)
	rownames(res) = c(others[,1],scvi[,1])

	celltypes = celltypes[order(prop_values)]
	res = res[,order(prop_values)]
	prop_values = prop_values[order(prop_values)]

	prop_values = prop_values[!is.na(colSums(res))]
	celltypes = celltypes[!is.na(colSums(res))]
	res = res[,!is.na(colSums(res))]

	prop_values = prop_values[celltypes!='nan']
	res = res[,celltypes!='nan']
	celltypes = celltypes[celltypes!='nan']

	df = data.frame(x=c(1:length(celltypes)), VAE =res[rownames(res)=='vae',],celltypes= celltypes)

	if('SCMAP' %in% methods){df$SCMAP = res[rownames(res)=='scmap',]}
	if('CCA' %in% methods){df$CCA = res[rownames(res)=='readSeurat',]}
	if('SCANVI' %in% methods){df$SCANVI = res[rownames(res)=='scanvi',]}

	df$prop = prop_values
	df = melt(df,id=c('x','celltypes','prop'))

	p = ggplot(df,aes(x,value)) + 
	geom_point(aes(size=prop,colour=variable),position=position_jitter(w=0.1,h=0.1)) + 
	scale_size_area(max_size = 20)+
	scale_x_continuous(breaks=df$x[df$variable=='VAE'],labels=celltypes) +
	theme(axis.text.x = element_text(angle = 90, hjust = 1))

	ggsave(plotname,p)

}


