library(reshape2)
library(ggplot2)

prop = read.table('celltypeprop.txt',sep='\t',row.names=NULL,as.is=T)
ncelltypes = length(prop[1,])-1
colnames(prop) = prop[1,]
prop = prop[2:4,1:ncelltypes]
colnames(prop)[colnames(prop)=='NaN']='nan'

others = read.csv('others.percluster.res.txt',sep='\t',row.names=NULL,as.is=T)
colnames(others) = c(colnames(others)[2:length(colnames(others))],'NA')
others = others[,1:length(others[1,])-1]

scvi = read.csv('scvi.percluster.res.txt',sep='\t',row.names=NULL,as.is=T)
colnames(scvi) = c(colnames(scvi)[2:length(colnames(scvi))],'NA')
scvi = scvi[,1:length(scvi[1,])-1]

scanvi = read.csv('scanvi_acc.txt',sep='\t',row.names=NULL,as.is=T)
# ncols = length(scanvi[1,])
# scanvi = scanvi[,c(1,c(3:(ncols-1)))]

Dotplot <- function(others,scvi,scanvi,ann,methods,plotname){
	scvi = scvi[scvi[,2]==ann & scvi[,1]=='vae',]
	others = others[others[,2]==ann,]
    celltypes = colnames(prop)
    if(ann=='p1'){
        scanvi = scanvi[scanvi[,1]=='scanvi1',]
    }else if(ann=='p2'){
        scanvi = scanvi[scanvi[,1]=='scanvi2',]
    }else if(ann=='p'){
        scanvi = scanvi[scanvi[,1]=='scanvi',]
    }
	if(ann=='p1'){prop_values = as.numeric(prop[2,])
		}else if(ann=='p2'){prop_values = as.numeric(prop[3,])
			}else{prop_values = as.numeric(prop[1,])}
	res = lapply(celltypes,function(celltype){
		temp =  c(others[,colnames(others)==celltype],
			scvi[,colnames(scvi)==celltype],
            scanvi[,colnames(scanvi)==celltype]
		)
		return(temp)
	})
	res = do.call(cbind,res)
	rownames(res) = c(others[,1],scvi[,1],'scanvi')
	celltypes = celltypes[order(prop_values)]
	res = res[,order(prop_values)]
	prop_values = prop_values[order(prop_values)]
    prop_values = prop_values[!is.na(colSums(res))]
	celltypes = celltypes[!is.na(colSums(res))]
	res = res[,!is.na(colSums(res))]
	prop_values = prop_values[celltypes!='nan']
	res = res[,celltypes!='nan']
	celltypes = celltypes[celltypes!='nan']
    df = data.frame(x=c(1:length(celltypes)), scVI=res[rownames(res)=='vae',],celltypes= celltypes)
	if('SCMAP' %in% methods){df$SCMAP = res[rownames(res)=='scmap',]}
	if('CCA' %in% methods){df$CCA = res[rownames(res)=='readSeurat',]}
	if('SCANVI' %in% methods){df$SCANVI = res[rownames(res)=='scanvi',]}
    if ('SCMAP' %in% methods){
            colors = c('red','darkgreen','orange','blue')
    }else if('CCA' %in% methods){
            colors = c('red','darkgreen','blue')
    }
	df$prop = prop_values
	df = melt(df,id=c('celltypes','prop','x'))
    df$variable = factor(df$variable, levels = c('scVI','SCANVI','SCMAP','CCA'))
	p = ggplot(df,aes(x,value)) +
	geom_point(aes(size=prop,colour=variable)) +
	scale_size_area(max_size = 10)+
	scale_x_continuous(breaks=df$x[df$variable=='scVI'],labels=celltypes) +
	theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_color_manual(values=colors)+ ylab("Cell Type Classification Accuracy")+xlab('Cell Type Labels')
	ggsave(plotname,p)
    p = ggplot(df,aes(x,value,fill=variable)) +
	scale_fill_manual(values=colors)+ geom_boxplot() +
    ylab("Cell Type Classification Accuracy")+xlab('Method')
	ggsave(paste('box_',plotname,sep=''),p)
}



Dotplot(others,scvi,scanvi,'p1',c('scVI','SCANVI','SCMAP','CCA'),'percluster_scmap_vae_scanvi_p1.pdf')
Dotplot(others,scvi,scanvi,'p2',c('scVI','SCANVI','SCMAP','CCA'),'percluster_scmap_vae_scanvi_p2.pdf')
Dotplot(others,scvi,scanvi,'p',c('scVI','SCANVI','CCA'),'percluster_cca_vae_scanvi.pdf')


# pdfjam Easy1/box_percluster_cca_vae_scanvi.pdf \
# Easy1/box_percluster_scmap_vae_scanvi_p1.pdf \
# Easy1/box_percluster_scmap_vae_scanvi_p2.pdf \
# Tech1/box_percluster_cca_vae_scanvi.pdf \
# Tech1/box_percluster_scmap_vae_scanvi_p1.pdf \
# Tech1/box_percluster_scmap_vae_scanvi_p2.pdf \
# Tech3/box_percluster_cca_vae_scanvi.pdf \
# Tech3/box_percluster_scmap_vae_scanvi_p1.pdf \
# Tech3/box_percluster_scmap_vae_scanvi_p2.pdf \
# Tech4/box_percluster_cca_vae_scanvi.pdf \
# Tech4/box_percluster_scmap_vae_scanvi_p1.pdf \
# Tech4/box_percluster_scmap_vae_scanvi_p2.pdf \
# --frame true --nup 3x4  --no-landscape \
# --scale 0.9 \
# --outfile acc_boxplot.combined.pdf

