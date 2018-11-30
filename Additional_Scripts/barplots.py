import pandas as pd
import numpy as np
import sys
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import matplotlib.pyplot as plt


dataname = str(sys.argv[1])
criteria = str(sys.argv[2])

scvi = pd.read_table('../' + dataname + '/scvi.res.txt',delim_whitespace=True)
others = pd.read_table('../' + dataname + '/others.res.txt',delim_whitespace=True)
stats = pd.concat([scvi,others])
model_types = stats['model_type']
stat_names = np.asarray(list(scvi.columns)[1:])
model_types = np.asarray(model_types)
#
# ['knn_asw', 'knn_nmi', 'knn_ari', 'knn_uca', 'knn_wuca',
#        'p_knn_asw', 'p_knn_nmi', 'p_knn_ari', 'p_knn_uca', 'p_knn_wuca',
#        'p1_knn_asw', 'p1_knn_nmi', 'p1_knn_ari', 'p1_knn_uca',
#        'p1_knn_wuca', 'p2_knn_asw', 'p2_knn_nmi', 'p2_knn_ari',
#        'p2_knn_uca', 'p2_knn_wuca', 'kmeans_asw', 'kmeans_nmi',
#        'kmeans_ari', 'kmeans_uca', 'kmeans_wuca', 'p_kmeans_asw',
#        'p_kmeans_nmi', 'p_kmeans_ari', 'p_kmeans_uca', 'p_kmeans_wuca',
#        'p1_kmeans_asw', 'p1_kmeans_nmi', 'p1_kmeans_ari', 'p1_kmeans_uca',
#        'p1_kmeans_wuca', 'p2_kmeans_asw', 'p2_kmeans_nmi',
#        'p2_kmeans_ari', 'p2_kmeans_uca', 'p2_kmeans_wuca',
#        'res_jaccard10', 'res_jaccard10.1', 'res_jaccard10.2',
#        'res_jaccard10.3', 'res_jaccard10.4', 'res_jaccard10.5',
#        'res_jaccard10.6', 'res_jaccard10.7', 'res_jaccard10.8',
#        'res_jaccard10.9', 'res_jaccard50', 'res_jaccard50.1',
#        'res_jaccard50.2', 'res_jaccard50.3', 'res_jaccard50.4',
#        'res_jaccard50.5', 'res_jaccard50.6', 'jaccard_score',
#        'likelihood', 'BE', 'classifier_acc']
res=[]
for x in np.unique(model_types):
    stat = np.mean(np.asarray(stats[model_types==x])[:, 1:], axis=0)
    res.append(stat)

model_types = np.unique(model_types)
res = np.asarray(res)

sorted_res=[]
model_order = ['vae', 'scanvi','scanvi1','scanvi2','readSeurat', 'MNN','Combat','PCA']
model_names = ['scVI','SCAN-VI','SCAN-VI1','SCAN-VI2','CCA','MNN','Combat','PCA']


# model_order = ['vae', 'readSeurat','MNN']
# model_names = ['scVI','CCA','MNN']

# model_order = ['vae']
# model_names = ['scVI','SCMAP']
# scmap = res[model_types=='scmap',:][0]

sorted_res = []
for x in model_order:
    sorted_res.append(res[model_types==x,:])

sorted_res = np.asarray(sorted_res)
sorted_res = sorted_res.squeeze()
if(len(model_order)==1):
    sorted_res = sorted_res.reshape(1,len(sorted_res))

# scmap1 = scmap[13]
# scmap2 =  scmap[18]


# criteria = 'BE'
filtered_res = np.asarray(sorted_res[:, stat_names==criteria]).ravel()
# filtered_res = np.append(filtered_res,scmap1)
#
# filtered_res = np.asarray(sorted_res[:, stat_names==criteria]).ravel()
#
plt.figure(figsize=(5, 3))
colors = ('r','g','mediumseagreen','yellowgreen','b','y','m','c')
plt.bar(np.arange(len(model_names)), filtered_res, color=colors)
plt.xticks(np.arange(len(model_names)), model_names,rotation=90)
plt.tight_layout()
plt.savefig("../%s/%s_%s_barplot.pdf" % (dataname,dataname,criteria))
#
# filtered_res = np.asarray(sorted_res[:, stat_names==criteria]).ravel()
# filtered_res = np.append(filtered_res,scmap2)
#
# filtered_res = np.asarray(sorted_res[:, stat_names==criteria]).ravel()
#
# plt.figure(figsize=(5, 5))
# plt.bar(np.arange(len(model_names)), filtered_res, color=colors[:len(filtered_res)])
# plt.xticks(np.arange(len(model_names)), model_names)
# plt.savefig("../%s/%s_%s_barplot_scmapvae.pdf" % (dataname,dataname,criteria))
#
