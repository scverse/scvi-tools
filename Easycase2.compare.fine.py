from scvi.harmonization.utils_chenling import run_model, eval_latent
use_cuda = True
from scvi.dataset.BICCN import *
from scvi.dataset.dataset import GeneExpressionDataset

import sys
model_type = str(sys.argv[1])
plotname = 'Macosko_Regev.fine'
print(model_type)

if model_type!='readSeurat':
    dataset1 = MacoskoDataset(coarse=False)
    dataset2 = RegevDataset(coarse=False)
    dataset1.subsample_genes(dataset1.nb_genes)
    dataset2.subsample_genes(dataset2.nb_genes)
    gene_dataset = GeneExpressionDataset.concat_datasets(dataset1, dataset2)
    gene_dataset.subsample_genes(5000)
    if model_type=='writedata':
        latent, batch_indices, labels, keys = run_model(model_type, gene_dataset, dataset1, dataset2, plotname)
    else:
        latent, batch_indices, labels, keys = run_model(model_type, gene_dataset, dataset1, dataset2)
elif model_type=='readSeurat':
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

#vae
# asw 0.24857216
# nmi 0.816948223930474
# ari 0.7137941166798535
# uca 0.7684242984890534
# Entropy batch mixing : 0.5986809094007403
# asw 0.25624087
# nmi 0.8086121204710517
# ari 0.6934660066471109
# uca 0.7506167129201357
# Entropy batch mixing : 0.601227696792241

# scanvi
# asw 0.17032264
# nmi 0.7911944565358177
# ari 0.6668419534923542
# uca 0.7086678238229551
# Entropy batch mixing : 0.6204837658895879
# asw 0.30121857
# nmi 0.8543763034969332
# ari 0.7707112664847917
# uca 0.8364619790009545
# Entropy batch mixing : 0.6082580751080363

# asw 0.21227054989866825
# nmi 0.7967035834174173
# ari 0.70222569164588
# uca 0.7674889113687651

#Seurat
# asw 0.21227054989866825
# nmi 0.7967035834174173
# ari 0.70222569164588
# uca 0.7674889113687651
# Entropy batch mixing : 0.5940539831386287

#[1] "runtime for CCA = 2019.222" "runtime for CCA = 504.306"
# [3] "runtime for CCA = 2524.14"  "runtime for CCA = 0"
# [5] "runtime for CCA = 0"

