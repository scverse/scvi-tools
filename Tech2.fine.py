use_cuda = True

from scvi.dataset.BICCN import *
from scvi.dataset.dataset import GeneExpressionDataset
from scvi.harmonization.utils_chenling import run_model, eval_latent
import sys
print(model_type)

model_type = str(sys.argv[1])
plotname = 'Zeng.fine'
if model_type !='readSeurat':
    dataset1 = Zeng10X(coarse=False)
    dataset1.subsample_cells(6274)
    dataset2 = ZengSS2(coarse=False)
    dataset1.subsample_genes(dataset1.nb_genes)
    dataset2.subsample_genes(dataset2.nb_genes)
    gene_dataset = GeneExpressionDataset.concat_datasets(dataset1, dataset2)
    latent, batch_indices, labels, keys = run_model(model_type, gene_dataset, dataset1, dataset2, filename=plotname,ngenes=5000)
elif model_type == 'readSeurat':
    latent, batch_indices, labels, keys = run_model(model_type, 0, 0, 0, plotname)
    batch_indices = batch_indices-1


shared_labels = set(labels[batch_indices==0]).intersection(labels[batch_indices==1])
shared = np.asarray([x in shared_labels for x in labels])
latentx = latent[shared,:]
batch_indicesx = batch_indices[shared]
labelsx = labels[shared]
keyx = keys[np.sort([x for x in shared_labels])]
labels_dict = dict(zip(np.sort([x for x in shared_labels]), np.arange(len(shared_labels))))
labelsx = np.asarray([labels_dict[i] for i in labelsx])
eval_latent(batch_indicesx, labelsx, latentx,keyx , plotname + '.' + model_type)



labels = gene_dataset.labels.ravel()
batch_indices = gene_dataset.batch_indices.ravel()
np.unique(batch_indices[shared],return_counts=True)
np.unique(batch_indices,return_counts=True)
