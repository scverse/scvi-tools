from scvi.dataset.dataset import GeneExpressionDataset
from scvi.harmonization.utils_chenling import SubsetGenes,run_model
from scvi.harmonization.utils_chenling import KNNJaccardIndex
plotname = 'Tech4'
from copy import deepcopy
import numpy as np

from scvi.dataset.scanorama import DatasetSCANORAMA
dirs = ['data/pancreas/pancreas_inDrop', 'data/pancreas/pancreas_multi_celseq2_expression_matrix', 'data/pancreas/pancreas_multi_celseq_expression_matrix', 'data/pancreas/pancreas_multi_fluidigmc1_expression_matrix', 'data/pancreas/pancreas_multi_smartseq2_expression_matrix']

datasets = [DatasetSCANORAMA(d) for d in dirs]

labels = (open('/data/scanorama/data/cell_labels/pancreas_cluster.txt').read().rstrip().split())

all_dataset = GeneExpressionDataset.concat_datasets(*datasets)
batch_id = all_dataset.batch_indices.ravel()

all_dataset = GeneExpressionDataset.concat_datasets(datasets[0],datasets[1])
all_dataset.cell_types,all_dataset.labels = np.unique(np.asarray(labels)[(batch_id==0)+(batch_id==1)],return_inverse=True)
all_dataset.labels = all_dataset.labels.reshape(len(all_dataset.labels),1)
all_dataset.n_labels = len(np.unique(all_dataset.labels))

all_dataset.subsample_genes(all_dataset.nb_genes)
dataset1 = deepcopy(all_dataset)
dataset1.update_cells(dataset1.batch_indices.ravel()==0)
dataset2 = deepcopy(all_dataset)
dataset2.update_cells(dataset2.batch_indices.ravel()==1)

seurat1 = np.genfromtxt('../Seurat_data/' + plotname + '.1.CCA.txt')
seurat2 = np.genfromtxt('../Seurat_data/' + plotname + '.2.CCA.txt')
seurat, batch_indices, labels, keys, stats = run_model('readSeurat', all_dataset, dataset1, dataset2,
                                                       filename=plotname)

dataset1, dataset2, gene_dataset = SubsetGenes(dataset1, dataset2, all_dataset, plotname)
vae1, _, _, _, _ = run_model('vae', dataset1, 0, 0, filename=plotname, rep='vae1')
vae2, _, _, _, _ = run_model('vae', dataset2, 0, 0, filename=plotname, rep='vae2')
vae,  _, _, _, _ = run_model('vae', gene_dataset, dataset1, dataset2,filename=plotname, rep='0')

KNeighbors = np.concatenate([np.arange(10, 100, 10), np.arange(100, 500, 50)])
seurat_seurat = [KNNJaccardIndex(seurat1, seurat2, seurat, batch_indices, k)[0] for k in KNeighbors]
vae_seurat = [KNNJaccardIndex(vae1, vae2, seurat, batch_indices, k)[0] for k in KNeighbors]
vae_vae = [KNNJaccardIndex(vae1, vae2, vae, batch_indices, k)[0] for k in KNeighbors]
seurat_vae = [KNNJaccardIndex(seurat1, seurat2, vae, batch_indices, k)[0] for k in KNeighbors]

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
plt.figure(figsize=(10, 10))
plt.plot(KNeighbors, seurat_seurat,'r',label='Seurat_Seurat')
plt.plot(KNeighbors, vae_vae,'b',label='VAE_VAE')
plt.plot(KNeighbors, vae_seurat,'g',label='VAE_Seurat')
plt.plot(KNeighbors, seurat_vae,'y',label='Seurat_VAE')
legend = plt.legend(loc='lower right', shadow=False)
plt.savefig("../%s/%s.KNN.pdf" % (plotname,plotname))
