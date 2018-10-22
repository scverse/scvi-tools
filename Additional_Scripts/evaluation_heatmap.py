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

# scmap2 = res[model_types=='scmap1',[16,17,18,19]]
# res[model_types=='scmap1',[16,17,18,19]] = -1
# temp = np.repeat(-1.0,len(res[0,:]))
# temp[[11,12,13,14]] = scmap2
# res = np.append(res, temp.reshape(1,len(temp)),axis=0)
# model_types = np.append(model_types,'scmap2')

sorted_res=[]
# model_order = ['vae','scanvi0', 'scanvi', 'scanvi1', 'scanvi2',
#        'scmap1', 'scmap2', 'readSeurat', 'Combat', 'MNN']
# model_order = ['vae','scanvi0', 'scanvi', 'scanvi1', 'scanvi2','readSeurat', 'Combat', 'MNN']
model_order = ['vae', 'scanvi', 'scanvi1', 'scanvi2',
       'scmap', 'readSeurat', 'Combat', 'MNN']
for x in model_order:
    sorted_res.append(res[model_types==x,:])

sorted_res = np.asarray(sorted_res)
sorted_res = sorted_res.squeeze()

# tabs = ['knn_asw','knn_uca','knn_wuca','kmeans_uca','kmeans_wuca','BE','jaccard_score']
# tabs = ['knn_asw','knn_uca','p_knn_uca','p1_knn_uca','p2_knn_uca','BE','jaccard_score','classifier_acc']
tabs = ['kmeans_asw','kmeans_uca','p_kmeans_uca','p1_kmeans_uca','p2_kmeans_uca','BE','jaccard_score','classifier_acc']
# for x in tabs:
#     print(x)
#     print(np.sum(stat_names==x))
filtered_res = [sorted_res[:, stat_names==x] for x in tabs]
filtered_res = np.concatenate(filtered_res,axis=1)

filtered_names = np.concatenate(np.asarray( [stat_names[[y in x for x in stat_names]] for y in tabs]))
filtered_names[0]='asw'
rm_values = (filtered_res==-1)


def impute(x):
    avg = np.mean(x[x != -1])
    x[x == -1] = avg
    return x


filtered_res = np.apply_along_axis(impute,0,filtered_res)



def Heatmap(value_matrix, cololor_matrix, rownames,colnames,title,filename):
    fig, ax = plt.subplots(figsize=(7,7))
    # We want to show all ticks...
    im = ax.imshow(cololor_matrix,aspect='auto')
    ax.set_xticks(np.arange(len(rownames)))
    ax.set_yticks(np.arange(len(colnames)))
    # ... and label them with the respective list entries
    ax.set_xticklabels(rownames)
    ax.set_yticklabels(colnames)
    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")
    # Loop over data dimensions and create text annotations.
    for i in range(len(colnames)):
        for j in range(len(rownames)):
            text = ax.text(j, i, "{:.2f}".format(value_matrix[i, j]),
                           ha="center", va="center", color="w")
    ax.set_title(title)
    fig.tight_layout()
    plt.savefig(filename, transparent=True)


filtered_res = np.asarray(filtered_res).astype('float')
scaler = MinMaxScaler()
scaled_res = scaler.fit_transform(filtered_res)
scaled_res[rm_values]=np.nan
filtered_res[rm_values]=np.nan
Heatmap(filtered_res,scaled_res,tabs,model_order,dataname,'../'+ dataname + '/'+dataname + '.heatmap.pdf')
