from numpy import loadtxt
from scipy.io import mmread
import tables
import scipy.sparse as sparse

from scvi.inference import UnsupervisedTrainer, SemiSupervisedTrainer
from scvi.inference.posterior import *
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.manifold import TSNE
import sys
use_cuda = True
from scvi.metrics.clustering import select_indices_evenly, clustering_scores, entropy_batch_mixing

from scvi.models.scanvi import SCANVI
from scvi.models.vae import VAE
from sklearn.neighbors import NearestNeighbors

from scvi.dataset import GeneExpressionDataset
from copy import deepcopy

def assign_label(cellid, geneid, labels_map, count, cell_type, seurat):
    labels = seurat[1:, 4]
    labels = np.int64(np.asarray(labels))
    labels_new = deepcopy(labels)
    for i, j in enumerate(labels_map):
        labels_new[labels == i] = j
    temp = dict(zip(cellid, count))
    new_count = []
    for x in seurat[1:, 5]:
        new_count.append(temp[x])
    new_count = sparse.vstack(new_count)
    dataset = GeneExpressionDataset(*GeneExpressionDataset.get_attributes_from_matrix(new_count, labels=labels_new),
                                    gene_names=geneid, cell_types=cell_type)
    return dataset


def trainVAE(gene_dataset,nlayers=2):
    vae = VAE(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches, n_labels=gene_dataset.n_labels,
              n_hidden=128, n_latent=10, n_layers=nlayers, dispersion='gene')
    trainer = UnsupervisedTrainer(vae, gene_dataset, train_size=1.0)
    trainer.train(n_epochs=250)
    if gene_dataset.X.shape[0] < 20000:
        trainer.train(n_epochs=250)
    batch_entropy = trainer.train_set.entropy_batch_mixing()
    print("Entropy batch mixing :", batch_entropy)
    trainer.train_set.show_t_sne(color_by='batches',
                                 save_name='../' + filename + '.' + model_type + '.batch.png')
    full = trainer.create_posterior(vae, gene_dataset, indices=np.arange(len(gene_dataset)))
    return full

def SCANVI_pretrain(gene_dataset, nlayers=2):
    vae = VAE(gene_dataset.nb_genes, gene_dataset.n_batches, gene_dataset.n_labels, n_latent=10, n_layers=nlayers)
    trainer = UnsupervisedTrainer(vae, gene_dataset, train_size=1.0)
    trainer.train(n_epochs=250)
    scanvi = SCANVI(gene_dataset.nb_genes, gene_dataset.n_batches, gene_dataset.n_labels, n_layers=nlayers)
    scanvi.load_state_dict(vae.state_dict(), strict=False)
    trainer_scanvi = SemiSupervisedTrainer(scanvi, gene_dataset, classification_ratio=1,
                                           n_epochs_classifier=1, lr_classification=5 * 1e-3)
    return trainer_scanvi


def run_model(model_type, gene_dataset, dataset1, dataset2, filename='temp', nlayers=2):
    if model_type == 'vae':
        full = trainVAE(gene_dataset)
        batch_entropy = full.entropy_batch_mixing()
        print("Entropy batch mixing :", batch_entropy)
        full.show_t_sne(color_by='batches',save_name='../' + filename + '.' + model_type + '.batch.png')
        ll = full.ll(verbose=True)
        latent, batch_indices, labels = full.sequential().get_latent()
        batch_indices = batch_indices.ravel()
        labels = labels.ravel()
        batch_entropy = full.entropy_batch_mixing()
        stats = [ll, batch_entropy]
    elif model_type == 'scanvi1':
        trainer_scanvi = SCANVI_pretrain(gene_dataset)
        trainer_scanvi.labelled_set = trainer_scanvi.create_posterior(indices=(gene_dataset.batch_indices == 0))
        trainer_scanvi.unlabelled_set = trainer_scanvi.create_posterior(indices=(gene_dataset.batch_indices == 1))
        trainer_scanvi.train(n_epochs=50)
        acc = trainer_scanvi.unlabelled_set.accuracy()
        full = trainer_scanvi.create_posterior(trainer_scanvi.model, gene_dataset, indices=np.arange(len(gene_dataset)))
        ll = full.ll(verbose=True)
        batch_entropy = full.entropy_batch_mixing()
        stats = [ll, batch_entropy, acc]
        trainer_scanvi.full_dataset.show_t_sne(color_by='batches', save_name='../'+filename+'.'+model_type+'.batch.png')
        keys = gene_dataset.cell_types
        latent, batch_indices, labels = trainer_scanvi.unlabelled_set.get_latent()
        batch_indices = batch_indices.ravel()
    elif model_type == 'scanvi2':
        trainer_scanvi = SCANVI_pretrain(gene_dataset)
        trainer_scanvi.labelled_set = trainer_scanvi.create_posterior(indices=(gene_dataset.batch_indices == 1))
        trainer_scanvi.unlabelled_set = trainer_scanvi.create_posterior(indices=(gene_dataset.batch_indices == 0))
        trainer_scanvi.train(n_epochs=50)
        print('svaec acc =',  trainer_scanvi.unlabelled_set.accuracy())
        full = trainer_scanvi.create_posterior(scanvi, gene_dataset, indices=np.arange(len(gene_dataset)))
        full.ll(verbose=True)
        full.marginal_ll(verbose=True)
        batch_entropy = full.entropy_batch_mixing()
        print("Entropy batch mixing :", batch_entropy)
        trainer_scanvi.full_dataset.show_t_sne(color_by='batches', save_name='../'+filename+'.'+model_type+'.batch.png')
        keys = gene_dataset.cell_types
        latent, batch_indices, labels = trainer_scanvi.unlabelled_set.get_latent()
        batch_indices = batch_indices.ravel()
    elif model_type == 'scanvi':
        vae = VAE(gene_dataset.nb_genes, gene_dataset.n_batches, gene_dataset.n_labels, n_latent=10, n_layers=2)
        trainer = UnsupervisedTrainer(vae, gene_dataset, train_size=1.0)
        trainer.train(n_epochs=250)
        scanvi = SCANVI(gene_dataset.nb_genes, gene_dataset.n_batches, gene_dataset.n_labels, n_layers=2)
        scanvi.load_state_dict(vae.state_dict(), strict=False)
        trainer_scanvi = SemiSupervisedTrainer(scanvi, gene_dataset, classification_ratio=1,
                                               n_epochs_classifier=1, lr_classification=5 * 1e-3)
        trainer_scanvi.train(n_epochs=50)
        print('svaec acc =',  trainer_scanvi.unlabelled_set.accuracy())
        batch_entropy = trainer_scanvi.full_dataset.entropy_batch_mixing()
        print("Entropy batch mixing :", batch_entropy)
        trainer_scanvi.full_dataset.show_t_sne(color_by='batches', save_name='../'+filename+'.'+model_type+'.batch.png')
        trainer_scanvi.full_dataset.ll(verbose=True)
        trainer_scanvi.full_dataset.marginal_ll(verbose=True)
        keys = gene_dataset.cell_types
        latent, batch_indices, labels = trainer_scanvi.unlabelled_set.get_latent()
        batch_indices = batch_indices.ravel()
    elif model_type =='scanvi0':
        vae = VAE(gene_dataset.nb_genes, gene_dataset.n_batches, gene_dataset.n_labels, n_latent=10, n_layers=2)
        trainer = UnsupervisedTrainer(vae, gene_dataset, train_size=1.0)
        trainer.train(n_epochs=250)
        scanvi = SCANVI(gene_dataset.nb_genes, gene_dataset.n_batches, gene_dataset.n_labels, n_layers=2)
        scanvi.load_state_dict(vae.state_dict(), strict=False)
        trainer_scanvi = SemiSupervisedTrainer(scanvi, gene_dataset, classification_ratio=1,
                                               n_labelled_samples_per_class=0,n_epochs_classifier=1,
                                               lr_classification=5 * 1e-3)
        trainer_scanvi.train(n_epochs=50)
        print('svaec acc =', trainer_scanvi.unlabelled_set.accuracy())
        batch_entropy = trainer_scanvi.full_dataset.entropy_batch_mixing()
        print("Entropy batch mixing :", batch_entropy)
        keys = gene_dataset.cell_types
        latent, batch_indices, labels = trainer_scanvi.unlabelled_set.get_latent()
        batch_indices = batch_indices.ravel()
    elif model_type=='MNN':
        from scvi.harmonization.clustering.MNN import MNN
        from sklearn.decomposition import PCA
        mnn = MNN()
        out = mnn.fit_transform(gene_dataset.X.todense(), gene_dataset.batch_indices.ravel(), [0, 1])
        latent = PCA(n_components=10).fit_transform(out)
        batch_indices = gene_dataset.batch_indices.ravel()
        labels = gene_dataset.labels.ravel()
        keys = gene_dataset.cell_types
    elif model_type == 'Combat':
        from scvi.harmonization.clustering.combat import COMBAT
        combat = COMBAT()
        latent = combat.combat_pca(gene_dataset)
        latent = latent.T
        batch_indices = np.concatenate(gene_dataset.batch_indices)
        labels = np.concatenate(gene_dataset.labels)
        labels = labels.astype('int')
        keys = gene_dataset.cell_types
    elif model_type == 'Seurat':
        from scvi.harmonization.clustering.seurat import SEURAT
        seurat = SEURAT()
        seurat.create_seurat(dataset1, 1)
        seurat.create_seurat(dataset2, 2)
        latent, batch_indices, labels, keys = seurat.get_cca()
    elif model_type == 'SeuratPC':
        from scvi.harmonization.clustering.seurat import SEURAT
        seurat = SEURAT()
        seurat.create_seurat(dataset1, 1)
        seurat.create_seurat(dataset2, 2)
        latent, batch_indices = seurat.get_pcs()
        labels, keys, _, _ = seurat.get_cca()
        keys=0
    # elif model_type =='scmap':
    #     from scvi.harmonization.classification.scmap import SCMAP
    #     print("Starting scmap")
    #     scmap = SCMAP()
    #     count1 = np.asarray(gene_dataset.X[gene_dataset.batch_indices.ravel() == 0, :].todense())
    #     label1 = gene_dataset.labels[gene_dataset.batch_indices.ravel() == 0].ravel().astype('int')
    #     count2 = np.asarray(gene_dataset.X[gene_dataset.batch_indices.ravel() == 1, :].todense())
    #     label2 = gene_dataset.labels[gene_dataset.batch_indices.ravel() == 1].ravel().astype('int')
    #     for n_features in [100, 300, 500, 1000, 2000]:
    #         scmap.set_parameters(n_features=n_features)
    #         scmap.fit_scmap_cluster(count1, label1.astype(np.int))
    #         print("Score Dataset1->Dataset2:%.4f  [n_features = %d]\n" % (scmap.score(count2, label2), n_features))
    #         scmap.fit_scmap_cluster(count2, label2.astype(np.int))
    #         print("Score Dataset2->Dataset1:%.4f  [n_features = %d]\n" % (scmap.score(count1, label1), n_features))
    #     sys.exit()
    elif model_type =='writedata':
        from scipy.io import mmwrite
        count = gene_dataset.X
        count = count.astype('int')
        genenames = gene_dataset.gene_names.astype('str')
        labels = gene_dataset.labels.ravel().astype('int')
        cell_types = gene_dataset.cell_types.astype('str')
        batchid = gene_dataset.batch_indices.ravel().astype('int')
        mmwrite('../Seurat_data/' + filename + '.X.mtx', count)
        np.savetxt('../Seurat_data/'+ filename + '.celltypes.txt', cell_types, fmt='%s')
        np.save('../Seurat_data/' + '../Seurat_data/'+ filename + '.labels.npy', labels)
        np.save('../Seurat_data/' + filename + '.genenames.npy', genenames)
        np.save('../Seurat_data/' + filename + '.batch.npy', batchid)
        sys.exit()
    elif model_type =='readSeurat':
        latent = np.genfromtxt('../Seurat_data/' + filename + '.CCA.txt')
        labels = gene_dataset.labels.ravel()
        keys = gene_dataset.cell_types
        batch_indices = gene_dataset.batch_indices.ravel()
    return latent, batch_indices, labels, keys


def eval_latent(batch_indices, labels, latent, keys, plotname=None,plotting=False):
    res = clustering_scores(np.asarray(latent), labels, 'knn', len(np.unique(labels)))
    for x in res:
        print('KNN',x, res[x])
    res = clustering_scores(np.asarray(latent), labels, 'KMeans', len(np.unique(labels)))
    for x in res:
        print('KMenas' ,x, res[x])
    sample = select_indices_evenly(2000, batch_indices)
    batch_entropy = entropy_batch_mixing(latent[sample, :], batch_indices[sample])
    print("Entropy batch mixing :", batch_entropy)
    if plotting==True:
        sample = select_indices_evenly(2000, labels)
        if plotname is not None:
            colors = sns.color_palette('tab20')
            latent_s = latent[sample, :]
            label_s = labels[sample]
            batch_s = batch_indices[sample]
            if latent_s.shape[1] != 2:
                latent_s = TSNE().fit_transform(latent_s)
            fig, ax = plt.subplots(figsize=(13, 10))
            key_order = np.argsort(keys)
            for i,k in enumerate(key_order):
                ax.scatter(latent_s[label_s == k, 0], latent_s[label_s == k, 1], c=colors[i%20], label=keys[k],
                           edgecolors='none')
                ax.legend(bbox_to_anchor=(1.1, 0.5), borderaxespad=0, fontsize='x-large')
            fig.tight_layout()
            plt.savefig('../'+plotname+'.labels.png')
            plt.figure(figsize=(10, 10))
            plt.scatter(latent_s[:, 0], latent_s[:, 1], c=batch_s, edgecolors='none')
            plt.axis("off")
            plt.tight_layout()
            plt.savefig('../' + plotname + '.batchid.png')


def JaccardIndex(x1,x2):
    intersection = np.sum(x1*x2)
    union = np.sum((x1+x2)>0)
    return intersection/union


def KNNJaccardIndex(latent1, latent2,latent,batchid,nn):
    knn = NearestNeighbors(n_neighbors=nn, algorithm='auto')
    nbrs1 = knn.fit(latent1)
    nbrs1 = nbrs1.kneighbors_graph(latent1).toarray()
    np.fill_diagonal(nbrs1,0)
    nbrs2 = knn.fit(latent2)
    nbrs2 = nbrs2.kneighbors_graph(latent2).toarray()
    np.fill_diagonal(nbrs2,0)
    nbrs_1 = knn.fit(latent[batchid==0,:])
    nbrs_1 = nbrs_1.kneighbors_graph(latent[batchid==0,:]).toarray()
    np.fill_diagonal(nbrs_1,0)
    nbrs_2 = knn.fit(latent[batchid==1,:])
    nbrs_2 = nbrs_2.kneighbors_graph(latent[batchid==1,:]).toarray()
    np.fill_diagonal(nbrs_2,0)
    JI1 = [JaccardIndex(x1, x2) for x1, x2 in zip(nbrs1, nbrs_1)]
    JI2 = [JaccardIndex(x1, x2) for x1, x2 in zip(nbrs2, nbrs_2)]
    return [(np.mean(JI1)+np.mean(JI2))/2]


def get_matrix_from_dir(dirname):
    geneid = loadtxt('../'+ dirname +'/genes.tsv',dtype='str',delimiter="\t")
    cellid = loadtxt('../'+ dirname + '/barcodes.tsv',dtype='str',delimiter="\t")
    count = mmread('../'+ dirname +'/matrix.mtx')
    return count, geneid, cellid



def get_matrix_from_h5(filename, genome):
    with tables.open_file(filename, 'r') as f:
        try:
            group = f.get_node(f.root, genome)
        except tables.NoSuchNodeError:
            print("That genome does not exist in this file.")
            return None
        gene_names = getattr(group, 'gene_names').read()
        barcodes = getattr(group, 'barcodes').read()
        data = getattr(group, 'data').read()
        indices = getattr(group, 'indices').read()
        indptr = getattr(group, 'indptr').read()
        shape = getattr(group, 'shape').read()
        matrix = sparse.csc_matrix((data, indices, indptr), shape=shape)
        gene_names = gene_names.astype('<U18')
        barcodes = barcodes.astype('<U18')
        return gene_names, barcodes, matrix


def TryFindCells(dict, cellid, count):
    """

    :param dict: mapping from cell id to cluster id
    :param cellid: cell id
    :param count: count matrix in the same order as cell id (filter cells out if cell id has no cluster mapping)
    :return:
    """
    res = []
    new_count = []
    for i,key in enumerate(cellid):
        try:
            res.append(dict[key])
            new_count.append(count[i])
        except KeyError:
            continue
    new_count = sparse.vstack(new_count)
    return(new_count, np.asarray(res))

