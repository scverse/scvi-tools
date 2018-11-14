use_cuda = True
from scvi.dataset.BICCN import *
from scvi.dataset.dataset import GeneExpressionDataset
from scvi.harmonization.utils_chenling import CompareModels
from scvi.harmonization.utils_chenling import SubsetGenes
import sys
models = str(sys.argv[1])
plotname = 'Zeng'


dataset1 = Zeng10X()
dataset2 = ZengSS2()
dataset1.subsample_genes(dataset1.nb_genes)
dataset2.subsample_genes(dataset2.nb_genes)
gene_dataset = GeneExpressionDataset.concat_datasets(dataset1, dataset2)

CompareModels(gene_dataset, dataset1, dataset2, plotname, models)

dataset1, dataset2, gene_dataset = SubsetGenes(dataset1, dataset2, gene_dataset, plotname)
from scvi.harmonization.utils_chenling import run_model
latent1, _, _, _, _ = run_model('vae', dataset1, 0, 0, filename=plotname, rep='vae1')
latent2, _, _, _, _ = run_model('vae', dataset2, 0, 0, filename=plotname, rep='vae2')

model_type='scanvi2'
latent, batch_indices, labels, keys, stats = run_model(model_type, gene_dataset, dataset1, dataset2,
                                                       filename=plotname, rep='0')

res_jaccard = [KNNJaccardIndex(latent1, latent2, latent, batch_indices, k)[0] for k in KNeighbors]
res_jaccard_score = np.sum(res_jaccard * K_int)
res_knn, res_knn_partial, res_kmeans, res_kmeans_partial = \
    eval_latent(batch_indices=batch_indices, labels=labels, latent=latent, keys=keys,
                labelled_idx=labelled_idx, unlabelled_idx=unlabelled_idx,
                plotname=plotname + '.' + model_type, plotting=True, partial_only=False)

_, res_knn_partial1, _, res_kmeans_partial1 = \
    eval_latent(batch_indices=batch_indices, labels=labels, latent=latent, keys=keys,
                labelled_idx=(batch_indices == 0), unlabelled_idx=(batch_indices == 1),
                plotname=plotname + '.' + model_type, plotting=False)

_, res_knn_partial2, _, res_kmeans_partial2 = \
    eval_latent(batch_indices=batch_indices, labels=labels, latent=latent, keys=keys,
                labelled_idx=(batch_indices == 1), unlabelled_idx=(batch_indices == 0),
                plotname=plotname + '.' + model_type, plotting=False)

res = [res_knn[x] for x in ['asw', 'nmi', 'ari', 'ca', 'weighted ca']] + \
      [res_knn_partial[x] for x in ['asw', 'nmi', 'ari', 'ca', 'weighted ca']] + \
      [res_knn_partial1[x] for x in ['asw', 'nmi', 'ari', 'ca', 'weighted ca']] + \
      [res_knn_partial2[x] for x in ['asw', 'nmi', 'ari', 'ca', 'weighted ca']] + \
      [res_kmeans[x] for x in ['asw', 'nmi', 'ari', 'uca', 'weighted uca']] + \
      [res_kmeans_partial[x] for x in ['asw', 'nmi', 'ari', 'uca', 'weighted uca']] + \
      [res_kmeans_partial1[x] for x in ['asw', 'nmi', 'ari', 'uca', 'weighted uca']] + \
      [res_kmeans_partial2[x] for x in ['asw', 'nmi', 'ari', 'uca', 'weighted uca']] + \
      res_jaccard + \
      [res_jaccard_score, stats[0], stats[1], stats[2]]
