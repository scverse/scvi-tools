import pandas as pd
import numpy as np
import sys
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import matplotlib.pyplot as plt
from sklearn.preprocessing import scale,MinMaxScaler


dataname = str(sys.argv[1])

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
# model_order = ['vae', 'scanvi','readSeurat', 'MNN','Combat','PCA']
# model_names = ['scVI','SCAN-VI','CCA','MNN','Combat','PCA']
model_order = ['vae',  'readSeurat', 'MNN']
model_names = ['scVI','CCA','MNN']
for x in model_order:
    sorted_res.append(res[model_types==x,:])

sorted_res = np.asarray(sorted_res)
sorted_res = sorted_res.squeeze()

filtered_res = [np.asarray(sorted_res[:, i]) for i,x in enumerate(stat_names) if 'res_jaccard' in x]
filtered_res = np.asarray(filtered_res)

import matplotlib
KNeighbors = np.concatenate([np.arange(10, 100, 10), np.arange(100, 500, 50)])
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt

plt.figure(figsize=(5, 5))
colors = ('r','g','b','y','m','c')

for i,x in enumerate(model_names):
    plt.plot(KNeighbors, filtered_res[:,i],colors[i],label=x)

legend = plt.legend(loc='lower right', shadow=False)
plt.savefig("../%s/%s_compare4_KNN.pdf" % (dataname,dataname))
#
# plt.figure(figsize=(5, 5))
# colors = ('r','g','b','y','m','c')
# for i in [0,2,3]:
#     x = model_names[i]
#     plt.plot(KNeighbors, filtered_res[:,i],colors[i],label=x)
#
# legend = plt.legend(loc='lower right', shadow=False)
# plt.savefig("../%s/%s_compare3_KNN.pdf" % (dataname,dataname))
#
