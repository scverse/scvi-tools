from numpy import loadtxt
from scipy.io import mmread
import tables
import scipy.sparse as sparse

from scvi.inference import UnsupervisedTrainer, SemiSupervisedTrainer
from scvi.inference.posterior import *
use_cuda = True
from scvi.metrics.clustering import select_indices_evenly, clustering_scores,clustering_accuracy
from scvi.inference.posterior import entropy_batch_mixing
from scvi.models.scanvi import SCANVI
from scvi.models.vae import VAE
from sklearn.neighbors import NearestNeighbors

from scvi.dataset import GeneExpressionDataset
from copy import deepcopy
import os
import torch
from scvi.dataset.dataset import SubsetGenes


import seaborn as sns
from umap import UMAP
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt


def scmap_eval(labels_unlabelled,labels_pred,labels):
    res = {
                'nmi': NMI(labels_unlabelled, labels_pred),
                'ari': ARI(labels_unlabelled, labels_pred),
                'ca': clustering_accuracy(labels_unlabelled, labels_pred),
                'weighted ca': clustering_accuracy(labels_unlabelled, labels_pred, True),
                'clusteracc':[np.mean(labels_pred[labels_unlabelled == i] == i) for i in np.unique(labels)]

    }
    return res

def JaccardIndex(x1,x2):
    intersection = np.sum(x1*x2)
    union = np.sum((x1+x2)>0)
    return intersection/union



def KNNJaccardIndex(latent1, latent2,latent,batchid,nn,subsample=False,max_number=30000):
    if subsample == True:
        n_samples = len(latent)
        keep_idx = np.random.choice(np.arange(n_samples), size=min(len(latent), max_number), replace=False)
        batch0size = np.sum(batchid == 0)
        keep_idx1 = keep_idx[keep_idx < batch0size]
        keep_idx2 = keep_idx[keep_idx >= batch0size] - batch0size
        latent1 = latent1[keep_idx1]
        latent2 = latent2[keep_idx2]
        latent, batchid = latent[keep_idx], batchid[keep_idx]
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


def trainVAE(gene_dataset, filename, rep, nlayers=2,n_hidden=128):
    vae = VAE(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches, n_labels=gene_dataset.n_labels,
          n_hidden=n_hidden, n_latent=10, n_layers=nlayers, dispersion='gene')
    trainer = UnsupervisedTrainer(vae, gene_dataset, train_size=1.0)
    if os.path.isfile('../' + filename + '/' + 'vae' + '.rep'+str(rep)+'.pkl'):
        trainer.model = torch.load('../' + filename + '/' + 'vae' + '.rep'+str(rep)+'.pkl')
    else:
        trainer.train(n_epochs=250)
        torch.save(trainer.model,'../' + filename + '/' + 'vae' + '.rep'+str(rep)+'.pkl')
    full = trainer.create_posterior(trainer.model, gene_dataset, indices=np.arange(len(gene_dataset)))
    return full

def VAEstats(full):
    ll = full.ll(verbose=True)
    latent, batch_indices, labels = full.sequential().get_latent()
    batch_indices = batch_indices.ravel()
    if len(np.unique(batch_indices))==2:
        sample = select_indices_evenly(np.min(np.unique(batch_indices,return_counts=True)[1]), batch_indices)
        batch_entropy = entropy_batch_mixing(latent[sample, :], batch_indices[sample])
    else:
        batch_entropy = -1
    labels = labels.ravel()
    stats = [ll,batch_entropy,-1,-1,np.arange(0,len(labels))]
    return latent, batch_indices, labels, stats

def SCANVIstats(trainer_scanvi,gene_dataset):
    full = trainer_scanvi.create_posterior(trainer_scanvi.model, gene_dataset, indices=np.arange(len(gene_dataset)))
    ll = full.ll(verbose=True)
    latent, batch_indices, labels = full.sequential().get_latent()
    batch_indices = batch_indices.ravel()
    if len(np.unique(batch_indices))==2:
        sample = select_indices_evenly(np.min(np.unique(batch_indices,return_counts=True)[1]), batch_indices)
        batch_entropy = entropy_batch_mixing(latent[sample, :], batch_indices[sample])
    else:
        batch_entropy = -1
    labelled_idx = trainer_scanvi.labelled_set.indices
    unlabelled_idx = trainer_scanvi.unlabelled_set.indices
    trainer_scanvi.unlabelled_set = trainer_scanvi.create_posterior(trainer_scanvi.model, gene_dataset,
                                                                    indices=unlabelled_idx)
    acc = trainer_scanvi.unlabelled_set.accuracy()
    stats = [ll, batch_entropy, acc,labelled_idx, unlabelled_idx]
    return latent, batch_indices, labels, stats

def trainSCANVI(gene_dataset,model_type,filename,rep, nlayers=2):
    vae_posterior = trainVAE(gene_dataset,filename,rep)
    scanvi = SCANVI(gene_dataset.nb_genes, gene_dataset.n_batches, gene_dataset.n_labels, n_layers=nlayers)
    scanvi.load_state_dict(vae_posterior.model.state_dict(), strict=False)
    trainer_scanvi = SemiSupervisedTrainer(scanvi, gene_dataset, classification_ratio=50,
                                           n_epochs_classifier=1, lr_classification=5 * 1e-3)
    if model_type=='scanvi1':
        trainer_scanvi.labelled_set = trainer_scanvi.create_posterior(indices=(gene_dataset.batch_indices == 0))
        trainer_scanvi.unlabelled_set = trainer_scanvi.create_posterior(indices=(gene_dataset.batch_indices == 1))
    elif model_type=='scanvi2':
        trainer_scanvi.labelled_set = trainer_scanvi.create_posterior(indices=(gene_dataset.batch_indices == 1))
        trainer_scanvi.unlabelled_set = trainer_scanvi.create_posterior(indices=(gene_dataset.batch_indices == 0))
    elif model_type=='scanvi0':
        trainer_scanvi.labelled_set = trainer_scanvi.create_posterior(indices=(gene_dataset.batch_indices <0))
        trainer_scanvi.unlabelled_set = trainer_scanvi.create_posterior(indices=(gene_dataset.batch_indices >= 0))

    if os.path.isfile('../' + filename + '/' + model_type + '.rep'+str(rep)+'.pkl'):
        trainer_scanvi.model = torch.load('../' + filename + '/' + model_type + '.rep'+str(rep)+'.pkl')
    else:
        # if gene_dataset.X.shape[0]>40000:
        #     trainer_scanvi.train(n_epochs=10)
        # elif gene_dataset.X.shape[0] < 20000:
        #     trainer_scanvi.train(n_epochs=50)
        # else:
        #     trainer_scanvi.train(n_epochs=25)
        trainer_scanvi.train(n_epochs=50)
        torch.save(trainer_scanvi.model,'../' + filename + '/' + model_type + '.rep'+str(rep)+'.pkl')
    return trainer_scanvi


def run_model(model_type, gene_dataset, dataset1, dataset2, filename='temp',rep='0', plotting=False):
    keys = gene_dataset.cell_types
    batch_indices = gene_dataset.batch_indices.ravel().astype('int')
    labels = gene_dataset.labels.ravel().astype('int')
    stats = []

    if model_type == 'vae':
        full = trainVAE(gene_dataset, filename, rep)
        latent, batch_indices, labels, stats = VAEstats(full)

    elif model_type == 'scanvi1':
        trainer_scanvi = trainSCANVI(gene_dataset,model_type,filename,rep)
        if plotting==True and (os.path.isfile('../' + filename + '/' + model_type + '.batch.png') is False):
            trainer_scanvi.full_dataset.show_t_sne(color_by='batches',save_name='../' + filename + '/' + model_type + '.batch.png')
        latent, batch_indices, labels, stats = SCANVIstats(trainer_scanvi,gene_dataset)

    elif model_type == 'scanvi2':
        trainer_scanvi = trainSCANVI(gene_dataset,model_type,filename,rep)
        if plotting==True and (os.path.isfile('../' + filename + '/' + model_type + '.batch.png') is False):
            trainer_scanvi.full_dataset.show_t_sne(color_by='batches',save_name='../' + filename + '.' + model_type + '.batch.png')
        latent, batch_indices, labels, stats = SCANVIstats(trainer_scanvi,gene_dataset)

    elif model_type == 'scanvi':
        trainer_scanvi = trainSCANVI(gene_dataset,model_type,filename,rep)
        if plotting==True and (os.path.isfile('../' + filename + '/' + model_type + '.batch.png') is False):
            trainer_scanvi.full_dataset.show_t_sne(color_by='batches',save_name='../' + filename + '.' + model_type + '.batch.png')
        latent, batch_indices, labels, stats = SCANVIstats(trainer_scanvi,gene_dataset)

    elif model_type =='scanvi0':
        trainer_scanvi = trainSCANVI(gene_dataset,model_type,filename,rep)
        if plotting==True and (os.path.isfile('../' + filename + '/' + model_type + '.batch.png') is False):
            trainer_scanvi.full_dataset.show_t_sne(color_by='batches',save_name='../' + filename + '.' + model_type + '.batch.png')
        latent, batch_indices, labels, stats = SCANVIstats(trainer_scanvi,gene_dataset)

    elif model_type=='MNN':
        if os.path.isfile('../' + filename + '/' + 'MNN'  + '.npy'):
            latent = np.load('../' + filename + '/' + 'MNN'  + '.npy')
        else:
            from scvi.harmonization.clustering.MNN import MNN
            from sklearn.decomposition import PCA
            mnn = MNN()
            out = mnn.fit_transform(gene_dataset.X.todense(), gene_dataset.batch_indices.ravel(), [0, 1])
            latent = PCA(n_components=10).fit_transform(out)
            np.save('../' + filename + '/' + 'MNN' + '.npy',latent)

    elif model_type=='PCA':
        if os.path.isfile('../' + filename + '/' + 'PCA'  + '.npy'):
            latent = np.load('../' + filename + '/' + 'PCA'  + '.npy')
        else:
            from sklearn.decomposition import PCA
            X = np.log(1 + gene_dataset.X.todense())
            latent = PCA(n_components=10).fit_transform(X)
            np.save('../' + filename + '/' + 'PCA' + '.npy', latent)

    elif model_type == 'Combat':
        if os.path.isfile('../' + filename + '/' + 'Combat'  + '.npy'):
            latent = np.load('../' + filename + '/' + 'Combat'  + '.npy')
        else:
            from scvi.harmonization.clustering.combat import COMBAT
            combat = COMBAT()
            latent = combat.combat_pca(gene_dataset)
            latent = latent.T
            np.save('../' + filename + '/' + 'Combat' + '.npy',latent)

    elif model_type == 'readSeurat':
        latent = np.genfromtxt('../Seurat_data/' + filename + '.CCA.txt')

    elif model_type == 'Seurat':
        if os.path.isfile('../' + filename + '/' + 'Seurat'+'.'+ rep +'.npy'):
            latent = np.load('../' + filename + '/' + 'Seurat'+'.'+ rep +'.npy')
        else:
            from scvi.harmonization.clustering.seurat import SEURAT
            seurat = SEURAT('../' + filename + '/' + rep )
            seurat.create_seurat(dataset1, 1)
            seurat.create_seurat(dataset2, 2)
            latent, _, genes, cells = seurat.get_cca(filter_genes=True)
            np.save('../' + filename + '/' + 'Seurat'+'.'+ rep +'.npy',latent)
            np.save('../' + filename + '/' + 'Seurat' + '.' + rep + '.genes.npy', genes)
            np.save('../' + filename + '/' + 'Seurat' + '.' + rep + '.cells.npy', cells)

    elif model_type == 'SeuratPC':
        from scvi.harmonization.clustering.seurat import SEURAT
        seurat = SEURAT()
        seurat.create_seurat(dataset1, 1)
        seurat.create_seurat(dataset2, 2)
        # latent, batch_indices = seurat.get_pcs()
        latent, batch_indices, genes, cells = seurat.get_cca()

    elif model_type =='scmap':
        # if os.path.isfile('../' + filename + '/' + 'scmap.1.'+ rep +'.npy'):
        #     pred1 = np.load('../' + filename + '/' + 'scmap.1.'+ rep +'.npy')
        #     pred2 = np.load('../' + filename + '/' + 'scmap.2.'+ rep +'.npy')
        # else:
        from scvi.harmonization.classification.scmap import SCMAP
        print("Starting scmap")
        n_features = 500
        scmap = SCMAP()
        scmap.set_parameters(n_features=n_features)
        scmap.fit_scmap_cluster(gene_dataset, filename, batch=0)
        pred1 = scmap.predict_scmap_cluster(gene_dataset, filename, batch=1)
        # pred1cell = scmap.predict_scmap_cell(count2, label2.astype(np.int))
        scmap = SCMAP()
        scmap.set_parameters(n_features=n_features)
        scmap.fit_scmap_cluster(gene_dataset, filename, batch=1)
        pred2 = scmap.predict_scmap_cluster(gene_dataset, filename, batch=0)
        # pred2cell = scmap.predict_scmap_cell(count2, label2.astype(np.int))
        # np.save('../' + filename + '/' + 'scmap.1.' + rep + '.npy', pred1)
        # np.save('../' + filename + '/' + 'scmap.2.' + rep + '.npy', pred2)

        latent = pred1
        stats = pred2

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
        latent = 0
        stats = 0

    return latent, batch_indices, labels, keys, stats


def eval_latent(batch_indices, labels, latent, keys, labelled_idx=None,unlabelled_idx=None, plotname=None,plotting=False,partial_only=True):
    res_knn_partial = clustering_scores(latent, labels,'knn',True,labelled_idx,unlabelled_idx)
    res_kmeans_partial = clustering_scores(latent, labels,'KMeans',True,labelled_idx,unlabelled_idx)
    if partial_only==False:
        res_knn = clustering_scores(np.asarray(latent), labels, 'knn')
        res_kmeans = clustering_scores(np.asarray(latent), labels, 'KMeans')
    # sample = select_indices_evenly(2000, batch_indices)
    # batch_entropy = entropy_batch_mixing(latent[sample, :], batch_indices[sample])
    # print("Entropy batch mixing :", batch_entropy)
    if plotting==True and (os.path.isfile('../'+plotname+'.labels.pdf') is False):
        sample = select_indices_evenly(2000, labels)
        if plotname is not None:
            colors = sns.color_palette('bright') +\
                  sns.color_palette('muted') + \
                  sns.color_palette('dark') + \
                  sns.color_palette('pastel') + \
                  sns.color_palette('colorblind')
            latent_s = latent[sample, :]
            label_s = labels[sample]
            batch_s = batch_indices[sample]
            if latent_s.shape[1] != 2:
                latent_s = UMAP(spread=2).fit_transform(latent_s)
            fig, ax = plt.subplots(figsize=(18, 18))
            key_order = np.argsort(keys)
            for i,k in enumerate(key_order):
                ax.scatter(latent_s[label_s == k, 0], latent_s[label_s == k, 1], c=colors[i%30], label=keys[k],
                           edgecolors='none')
                # ax.legend(bbox_to_anchor=(1.1, 0.5), borderaxespad=0, fontsize='x-large')
            fig.patch.set_visible(False)
            ax.axis('off')
            fig.tight_layout()
            plt.savefig('../'+plotname+'.labels.pdf')
            plt.figure(figsize=(18, 18))
            plt.scatter(latent_s[:, 0], latent_s[:, 1], c=batch_s, edgecolors='none')
            plt.axis("off")
            plt.tight_layout()
            plt.savefig('../' + plotname + '.batchid.pdf')
    if partial_only==False:
        return res_knn,res_knn_partial, res_kmeans, res_kmeans_partial
    else:
        return 0, res_knn_partial, 0, res_kmeans_partial


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


def CompareModels(gene_dataset, dataset1, dataset2, plotname, models):
    KNeighbors = np.concatenate([np.arange(10, 100, 10), np.arange(100, 500, 50)])
    K_int = np.concatenate([np.repeat(10, 10), np.repeat(50, 7)])
    f = open('../' + plotname +'/' + models + '.res.txt', "w+")
    f.write("model_type " + \
            "knn_asw knn_nmi knn_ari knn_uca knn_wuca " + \
            "p_knn_asw p_knn_nmi p_knn_ari p_knn_uca p_knn_wuca " + \
            "p1_knn_asw p1_knn_nmi p1_knn_ari p1_knn_uca p1_knn_wuca " + \
            "p2_knn_asw p2_knn_nmi p2_knn_ari p2_knn_uca p2_knn_wuca " + \
            "kmeans_asw kmeans_nmi kmeans_ari kmeans_uca kmeans_wuca " + \
            "p_kmeans_asw p_kmeans_nmi p_kmeans_ari p_kmeans_uca p_kmeans_wuca " + \
            "p1_kmeans_asw p1_kmeans_nmi p1_kmeans_ari p1_kmeans_uca p1_kmeans_wuca " + \
            "p2_kmeans_asw p2_kmeans_nmi p2_kmeans_ari p2_kmeans_uca p2_kmeans_wuca " + \
            " ".join(['res_jaccard'+ x for x in np.concatenate([np.repeat(10, 10), np.repeat(50,7)]).astype('str')])+" " + \
            'jaccard_score likelihood BE classifier_acc\n'
            )
    g = open('../' + plotname +'/' + models + '.percluster.res.txt', "w+")
    g.write("model_type\tannotation\t" + "\t".join(gene_dataset.cell_types)+"\n")

    scanvi = SCANVI(gene_dataset.nb_genes, gene_dataset.n_batches, gene_dataset.n_labels)
    trainer_scanvi = SemiSupervisedTrainer(scanvi, gene_dataset, classification_ratio=1,
                                           n_epochs_classifier=1, lr_classification=5 * 1e-3)
    labelled_idx = trainer_scanvi.labelled_set.indices
    unlabelled_idx = trainer_scanvi.unlabelled_set.indices

    if models =='others':
        latent1 = np.genfromtxt('../Seurat_data/' + plotname + '.1.CCA.txt')
        latent2 = np.genfromtxt('../Seurat_data/' + plotname + '.2.CCA.txt')
        # for model_type in ['scmap']:
        for model_type in ['scmap','readSeurat','Combat','MNN','PCA']:
        # for model_type in [ 'readSeurat', 'Combat', 'MNN', 'PCA']:
            print(model_type)
            if model_type == 'scmap':
                latent, batch_indices, labels, keys, stats = run_model('scmap', gene_dataset, dataset1, dataset2,
                                                                       filename=plotname)
                pred1 = latent
                pred2 = stats
                res1 = scmap_eval(pred1, labels[batch_indices == 1], labels)
                res2 = scmap_eval(pred2, labels[batch_indices == 0], labels)
                g.write("%s\t" % (model_type) + "p1\t"+
                        ("%.4f\t" * len(gene_dataset.cell_types) % tuple(res1['clusteracc']) +
                         "\n"))
                g.write("%s\t" % (model_type) + "p2\t" +
                        ("%.4f\t" * len(gene_dataset.cell_types) % tuple(res2['clusteracc']) +
                         "\n"))
                res = [-1]*10 + \
                      [-1]+[res1[x] for x in ['nmi', 'ari', 'ca', 'weighted ca']] + \
                      [-1] + [res2[x] for x in ['nmi', 'ari', 'ca', 'weighted ca']] + \
                      [-1]*41

                f.write(model_type + (" %.4f" * 61 + "\n") % tuple(res))
            else:
                if model_type=='readSeurat':
                    dataset1, dataset2, gene_dataset = SubsetGenes(dataset1, dataset2, gene_dataset, plotname)

                latent, batch_indices, labels, keys, stats = run_model(model_type, gene_dataset, dataset1, dataset2,
                                                                           filename=plotname)

                res_jaccard = [KNNJaccardIndex(latent1, latent2, latent, batch_indices, k)[0] for k in KNeighbors]
                res_jaccard_score = np.sum(res_jaccard * K_int)
                res_knn, res_knn_partial, res_kmeans, res_kmeans_partial = \
                    eval_latent(batch_indices, labels, latent, keys,
                                labelled_idx, unlabelled_idx,
                                plotname=plotname + '.' + model_type, plotting=True,partial_only=False)

                _, res_knn_partial1, _, res_kmeans_partial1 = \
                    eval_latent(batch_indices, labels, latent, keys,
                                batch_indices == 0, batch_indices == 1,
                                plotname=plotname + '.' + model_type, plotting=False)

                _, res_knn_partial2, _, res_kmeans_partial2 = \
                    eval_latent(batch_indices, labels, latent, keys,
                                batch_indices == 1, batch_indices == 0,
                                plotname=plotname + '.' + model_type, plotting=False)

                sample = select_indices_evenly(np.min(np.unique(batch_indices,return_counts=True)[1]), batch_indices)
                batch_entropy = entropy_batch_mixing(latent[sample, :], batch_indices[sample])

                res = [res_knn[x] for x in ['asw', 'nmi', 'ari', 'ca', 'weighted ca']] + \
                      [res_knn_partial[x] for x in ['asw', 'nmi', 'ari', 'ca', 'weighted ca']] + \
                      [res_knn_partial1[x] for x in ['asw', 'nmi', 'ari', 'ca', 'weighted ca']] + \
                      [res_knn_partial2[x] for x in ['asw', 'nmi', 'ari', 'ca', 'weighted ca']] + \
                      [res_kmeans[x] for x in ['asw', 'nmi', 'ari', 'uca', 'weighted uca']] + \
                      [res_kmeans_partial[x] for x in ['asw', 'nmi', 'ari', 'uca', 'weighted uca']] + \
                      [res_kmeans_partial1[x] for x in ['asw', 'nmi', 'ari', 'uca', 'weighted uca']] + \
                      [res_kmeans_partial2[x] for x in ['asw', 'nmi', 'ari', 'uca', 'weighted uca']] + \
                      res_jaccard + \
                      [res_jaccard_score, -1, batch_entropy, -1]

                f.write(model_type + (" %.4f"*61+"\n") % tuple(res))
                g.write("%s\t" % (model_type) +'all\t'+
                        ("%.4f\t" * len(gene_dataset.cell_types) % tuple(res_knn['clusteracc']) +
                         "\n"))
                g.write("%s\t" % (model_type) +'p\t'+
                        ("%.4f\t" * len(gene_dataset.cell_types) % tuple(res_knn_partial['clusteracc']) +
                         "\n"))
                g.write("%s\t" % (model_type) +'p1\t'+
                        ("%.4f\t" * len(gene_dataset.cell_types) % tuple(res_knn_partial1['clusteracc']) +
                         "\n"))
                g.write("%s\t" % (model_type) +'p2\t'+
                        ("%.4f\t" * len(gene_dataset.cell_types) % tuple(res_knn_partial2['clusteracc']) +
                         "\n"))

    elif models=='scvi':
        if plotname != 'Sim1' and plotname != 'Sim2' and plotname != 'Sim3':
            dataset1, dataset2, gene_dataset = SubsetGenes(dataset1, dataset2, gene_dataset, plotname)
        latent1, _, _, _, _ = run_model('vae', dataset1, 0, 0, filename=plotname,rep='vae1')
        latent2, _, _, _, _ = run_model('vae', dataset2, 0, 0, filename=plotname,rep='vae2')

        for model_type in ['vae', 'scanvi', 'scanvi1', 'scanvi2', 'scanvi0']:
            print(model_type)
            latent, batch_indices, labels, keys, stats = run_model(model_type, gene_dataset, dataset1, dataset2,
                                                                   filename=plotname, rep='0')

            res_jaccard = [KNNJaccardIndex(latent1, latent2, latent, batch_indices, k)[0] for k in KNeighbors]
            res_jaccard_score = np.sum(res_jaccard * K_int)
            res_knn, res_knn_partial, res_kmeans, res_kmeans_partial = \
                eval_latent(batch_indices=batch_indices,labels=labels,latent=latent,keys=keys,
                            labelled_idx=labelled_idx,unlabelled_idx=unlabelled_idx,
                            plotname=plotname+'.'+model_type,plotting=True, partial_only=False)

            _, res_knn_partial1, _, res_kmeans_partial1 = \
                eval_latent(batch_indices= batch_indices, labels=labels, latent=latent, keys=keys,
                            labelled_idx=(batch_indices==0), unlabelled_idx=(batch_indices==1),
                            plotname=plotname + '.' + model_type, plotting=False)

            _, res_knn_partial2, _, res_kmeans_partial2 = \
                eval_latent(batch_indices= batch_indices, labels=labels, latent=latent, keys=keys,
                            labelled_idx=(batch_indices==1), unlabelled_idx=(batch_indices==0),
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

            f.write(model_type + (" %.4f" * 61 + "\n") % tuple(res))
            g.write("%s\t" % (model_type) + 'all\t' +
                    ("%.4f\t" * len(gene_dataset.cell_types) % tuple(res_knn['clusteracc']) +
                     "\n"))
            g.write("%s\t" % (model_type) + 'p\t' +
                    ("%.4f\t" * len(gene_dataset.cell_types) % tuple(res_knn_partial['clusteracc']) +
                     "\n"))
            g.write("%s\t" % (model_type) + 'p1\t' +
                    ("%.4f\t" * len(gene_dataset.cell_types) % tuple(res_knn_partial1['clusteracc']) +
                     "\n"))
            g.write("%s\t" % (model_type) + 'p2\t' +
                    ("%.4f\t" * len(gene_dataset.cell_types) % tuple(res_knn_partial2['clusteracc']) +
                     "\n"))
            # for i in [1, 2, 3]:
            #     latent, batch_indices, labels, keys, stats = run_model(model_type, gene_dataset, dataset1, dataset2,
            #                                                            filename=plotname, rep=str(i))
            #     res_jaccard, res_jaccard_score = KNNJaccardIndex(latent1, latent2, latent, batch_indices)
            #
            #     res_knn, res_knn_partial, res_kmeans, res_kmeans_partial = \
            #         eval_latent(batch_indices=batch_indices, labels=labels, latent=latent, keys=keys,
            #                     labelled_idx=labelled_idx, unlabelled_idx=unlabelled_idx,
            #                     plotname=plotname + '.' + model_type, plotting=False,partial_only=False)
            #
            #     _, res_knn_partial1, _, res_kmeans_partial1 = \
            #         eval_latent(batch_indices=batch_indices, labels=labels, latent=latent, keys=keys,
            #                     labelled_idx=(batch_indices == 0), unlabelled_idx=(batch_indices == 1),
            #                     plotname=plotname + '.' + model_type, plotting=False)
            #
            #     _, res_knn_partial2, _, res_kmeans_partial2 = \
            #         eval_latent(batch_indices=batch_indices, labels=labels, latent=latent, keys=keys,
            #                     labelled_idx=(batch_indices == 1), unlabelled_idx=(batch_indices == 0),
            #                     plotname=plotname + '.' + model_type, plotting=False)
            #
            #     res = [res_knn[x] for x in res_knn] + \
            #           [res_knn_partial[x] for x in res_knn_partial] + \
            #           [res_knn_partial1[x] for x in res_knn_partial1] + \
            #           [res_knn_partial2[x] for x in res_knn_partial2] + \
            #           [res_kmeans[x] for x in res_kmeans] + \
            #           [res_kmeans_partial[x] for x in res_kmeans_partial] + \
            #           [res_kmeans_partial1[x] for x in res_kmeans_partial1] + \
            #           [res_kmeans_partial2[x] for x in res_kmeans_partial2] + \
            #           res_jaccard + \
            #           [res_jaccard_score,stats[0], stats[1], stats[2]]
            #     f.write(model_type + (" %.4f" * 61 + "\n") % tuple(res))

    elif models=='writedata':
        _, _, _, _,_ = run_model('writedata', gene_dataset, dataset1, dataset2, filename=plotname)
    f.close()
    g.close()


def subsetByGenenames(dataset, subsetnames):
    genenames = dataset.gene_names
    filter = np.asarray([x in subsetnames for x in genenames])
    dataset.X = dataset.X[:, filter]
    dataset.update_genes(np.arange(len(filter))[filter])
    return dataset

# def SubsetGenes(dataset1,dataset2,dataset3,gene_dataset,plotname,ngenes=1000):
def SubsetGenes(dataset1, dataset2, gene_dataset, plotname, ngenes=1000):
    import pandas as pd
    genes1 = pd.read_table('../Seurat_data/' + plotname + '.1.hvg_info.csv', delimiter=',')
    geneid1 = np.asarray([x.replace('gene_', '') for x in genes1[genes1.keys()[0]]]).astype('int')
    genenames1 = genes1['genename']
    genes2 = pd.read_table('../Seurat_data/' + plotname + '.2.hvg_info.csv', delimiter=',')
    geneid2 = np.asarray([x.replace('gene_', '') for x in genes2[genes2.keys()[0]]]).astype('int')
    genenames2 = genes2['genename']
    # genes3 = pd.read_table('../Seurat_data/' + plotname + '.3.hvg_info.csv', delimiter=',')
    # geneid3 = np.asarray([x.replace('gene_', '') for x in genes3[genes3.keys()[0]]]).astype('int')
    # genenames3 = genes3['genename']
    assert np.sum(np.asarray(genenames1) == gene_dataset.gene_names) == len(gene_dataset.gene_names)
    assert np.sum(np.asarray(genenames2) == gene_dataset.gene_names) == len(gene_dataset.gene_names)
    # assert np.sum(np.asarray(genenames3) == gene_dataset.gene_names) == len(gene_dataset.gene_names)
    # geneid = np.union1d(np.union1d(geneid1[:ngenes], geneid2[:ngenes]),geneid3[:ngenes]) - 1
    geneid = np.union1d(geneid1[:ngenes], geneid2[:ngenes]) - 1
    genes = gene_dataset.gene_names[geneid]
    dataset1 = subsetByGenenames(dataset1,genes)
    dataset2 = subsetByGenenames(dataset2,genes)
    # dataset3 = subsetByGenenames(dataset2,genes)
    gene_dataset = subsetByGenenames(gene_dataset,genes)
    return dataset1,dataset2,gene_dataset
    # return dataset1,dataset2,dataset3,gene_dataset
