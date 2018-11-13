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
scmap = res[model_types=='scmap',:][0]

sorted_res = []
for x in model_order:
    sorted_res.append(res[model_types==x,:])

sorted_res = np.asarray(sorted_res)
sorted_res = sorted_res.squeeze()
if(len(model_order)==1):
    sorted_res = sorted_res.reshape(1,len(sorted_res))

scmap1 = scmap[13]
scmap2 =  scmap[18]


# criteria = 'BE'
filtered_res = np.asarray(sorted_res[:, stat_names==criteria]).ravel()
# filtered_res = np.append(filtered_res,scmap1)
#
# filtered_res = np.asarray(sorted_res[:, stat_names==criteria]).ravel()
#
plt.figure(figsize=(5, 3))
colors = ['red','green','blue','yellow','magenta','cyan']
plt.bar(np.arange(len(model_names)), filtered_res, color=colors[:len(filtered_res)])
plt.xticks(np.arange(len(model_names)), model_names)
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
