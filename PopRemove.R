scvi = read.table('scvi.res.txt',header=TRUE)
others = read.table('others.res.txt',header=TRUE)
modelnames = c(as.character(scvi[,1]),as.character(others[,1]))
statnames = colnames(scvi)
scvi = scvi[,grep('res_jac',statnames)]
others = others[,grep('res_jac',statnames)]
data = rbind(scvi,others)
knearest = c(seq(10,100,10),seq(150,450,50))
colnames(data) = as.character(knearest)
rownames(data) = modelnames
library(reshape2)
library(ggplot2)
data = data[rownames(data)!='scmap',]
data$model = rownames(data)
df = melt(data)

ggplot(df, aes(x=as.integer(as.character(variable)),y=value,group=model))+
geom_point(aes(colour=model))+
geom_line(aes(colour=model),size=2) + 
theme_minimal()