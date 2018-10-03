from sklearn.manifold import TSNE
from sklearn.neighbors import kneighbors_graph
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

import sys
method = str(sys.argv[1])

tusi_batches = np.load('5000genes_10d/Tusi.batchid.npy')
paul_labels = np.load('Paul.labels.npy')
tusi_labels = np.load('B.npy')
potential = np.load('V.npy')

# paul = np.load('5000genes_10d/Paul.latent.npy')
# tusi = np.load('5000genes_10d/Tusi.latent.npy')
# merged = np.load('5000genes_10d/Tusi_Paul.latent.npy')
if method == 'CCA':
    paul = np.genfromtxt('latent/Continuous.1.CCA.txt')
    tusi = np.genfromtxt('latent/Continuous.2.CCA.txt')
    merged = np.genfromtxt('latent/Continuous.CCA.txt')
elif method == 'vae':
    paul = np.load('latent/Paul.vae.latent.npy')
    tusi = np.load('latent/Tusi.vae.latent.npy')
    merged = np.load('latent/merged.vae.latent.npy')
elif method == 'scanvi':
    paul = np.load('latent/Paul.scanvi.latent.npy')
    tusi = np.load('latent/Tusi.scanvi.latent.npy')
    merged = np.load('latent/merged.scanvi1.latent.npy')


def show_viz(latent, algo=None, labels=None, clusters_cmap=7, cmap="tab10", return_layout=False):    
    if clusters_cmap == 0:
        cmap = plt.get_cmap(cmap)
    else:
        cmap = plt.get_cmap(cmap, clusters_cmap)        
    if labels is None:
        labels = c_train
    if algo == "tSNE":
        layout = TSNE().fit_transform(latent)
    elif algo == "kNN":
        ad = kneighbors_graph(latent, 5, include_self=False)
        graph = nx.from_scipy_sparse_matrix(ad)
        lag_node = 0.9
        pos = nx.spring_layout(graph, k=np.sqrt(1.0/(lag_node*ad.shape[0])), iterations=200)
        layout =  np.concatenate([x[:, np.newaxis] for x in pos.values()], axis=1).T   
    else:
        layout = latent
    plt.figure(figsize=(10, 10))
    cax = plt.scatter(layout[:, 0], layout[:, 1], c=labels, \
                                   cmap=cmap, edgecolors='none')
    plt.axis("off")
    if algo is None:
        st = "tSNE"
    else:
        st = algo
    plt.tight_layout()
    if return_layout:
        return layout

spring_paul = show_viz(paul,algo='kNN',labels=paul_labels,return_layout=True)
spring_tusi = show_viz(tusi,algo='kNN',labels = potential[tusi_batches==1],return_layout=True)
potential_merged = np.concatenate([np.repeat('NaN',2730),potential[tusi_batches==1]])
potential_merged = potential_merged.astype('float')
spring_merged = show_viz(merged,algo='kNN',labels = potential_merged,return_layout=True)


if method == 'CCA':
    np.savetxt('spring_paul.CCA.csv', spring_paul,delimiter=',')
    np.savetxt('spring_tusi.CCA.csv', spring_tusi,delimiter=',')
    np.savetxt('spring_merged.CCA.csv', spring_merged,delimiter=',')
elif method == 'vae':
    np.savetxt('spring_paul.vae.csv', spring_paul,delimiter=',')
    np.savetxt('spring_tusi.vae.csv', spring_tusi,delimiter=',')
    np.savetxt('spring_merged.vae.csv', spring_merged,delimiter=',')
elif method == 'scanvi':
    np.savetxt('spring_paul.scanvi.csv', spring_paul,delimiter=',')
    np.savetxt('spring_tusi.scanvi.csv', spring_tusi,delimiter=',')
    np.savetxt('spring_merged.scanvi.csv', spring_merged,delimiter=',')



# batchid = np.concatenate([np.repeat(0,2730),np.repeat(1,4016)])
# plt.figure(figsize=(10, 10))
# plt.scatter(spring_merged[:, 0], spring_merged[:, 1], c='gray', edgecolors='none')
# plt.scatter(spring_merged[batchid==1, 0], spring_merged[batchid==1, 1], c=potential[tusi_batches==1].argsort().argsort(), edgecolors='none',s=10)
# plt.show()


# plt.figure(figsize=(10, 10))
# plt.scatter(spring_tusi[:, 0], spring_tusi[:, 1], c=potential[tusi_batches==1].argsort().argsort(), edgecolors='none')
# plt.show()

# batchid = np.concatenate([np.repeat(0,2730),np.repeat(1,4016)])
# plt.figure(figsize=(10, 10))
# plt.scatter(spring_merged[:, 0], spring_merged[:, 1], c='gray', edgecolors='none')
# plt.scatter(spring_merged[batchid==0, 0], spring_merged[batchid==0, 1], c=paul_labels.astype('int').astype('str'), edgecolors='none',s=10,cmap=plt.get_cmap('tab20'))
# plt.show()


# plt.figure(figsize=(10, 10))
# plt.scatter(spring_paul[:, 0], spring_paul[:, 1], c=paul_labels, edgecolors='none',cmap=plt.get_cmap('tab20'))
# plt.show()



