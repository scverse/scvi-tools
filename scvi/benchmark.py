import os

import numpy as np
from sklearn.decomposition import PCA

from scvi.dataset import CortexDataset, GeneExpressionDataset, SemiSupervisedDataLoaders
from scvi.dataset.BICNN import labels_groups
from scvi.inference import VariationalInference, VariationalInferenceFish, adversarial_wrapper, \
    SemiSupervisedVariationalInference
from scvi.metrics.classification import compute_accuracy_nn
from scvi.metrics.imputation import proximity_imputation
from scvi.models import VAE, VAEF, SCANVI


def cortex_benchmark(n_epochs=250, use_cuda=True, unit_test=False):
    cortex_dataset = CortexDataset()
    vae = VAE(cortex_dataset.nb_genes)
    infer_cortex_vae = VariationalInference(vae, cortex_dataset, use_cuda=use_cuda)
    infer_cortex_vae.train(n_epochs=n_epochs)

    infer_cortex_vae.ll('test')  # assert ~ 1200
    infer_cortex_vae.differential_expression('test')
    infer_cortex_vae.imputation('test', rate=0.1)  # assert ~ 2.3
    n_samples = 1000 if not unit_test else 10
    infer_cortex_vae.show_t_sne('test', n_samples=n_samples)
    return infer_cortex_vae


def benchmark(dataset, n_epochs=250, use_cuda=True):
    vae = VAE(dataset.nb_genes, n_batch=dataset.n_batches)
    infer = VariationalInference(vae, dataset, use_cuda=use_cuda)
    infer.train(n_epochs=n_epochs)
    infer.ll('test')
    infer.marginal_ll('test')
    infer.imputation('test', rate=0.1)  # assert ~ 2.1
    return infer


def harmonization_benchmarks(n_epochs=1, use_cuda=True):
    # retina_benchmark(n_epochs=n_epochs)
    pass


def annotation_benchmarks(n_epochs=1, use_cuda=True):
    # some cortex annotation benchmark
    pass


def all_benchmarks(n_epochs=250, use_cuda=True, unit_test=False):
    cortex_benchmark(n_epochs=n_epochs, use_cuda=use_cuda, unit_test=unit_test)

    harmonization_benchmarks(n_epochs=n_epochs, use_cuda=use_cuda)
    annotation_benchmarks(n_epochs=n_epochs, use_cuda=use_cuda)


def benchamrk_fish_scrna(gene_dataset_seq, gene_dataset_fish):
    gene_names = gene_dataset_fish.gene_names
    indexes_to_keep = np.arange(len(gene_names))
    vae = VAEF(gene_dataset_seq.nb_genes, indexes_to_keep, n_layers_decoder=2, n_latent=6,
               n_layers=2, n_hidden=256, reconstruction_loss='nb', dropout_rate=0.3, n_labels=7, n_batch=2,
               model_library=False)
    infer = VariationalInferenceFish(vae, gene_dataset_seq, gene_dataset_fish, train_size=0.9, verbose=True,
                                     frequency=5, weight_decay=0.35, n_epochs_even=100, n_epochs_kl=1000,
                                     cl_ratio=0, n_epochs_cl=100)
    infer = adversarial_wrapper(infer, scale=50, mode="smFISH")
    infer.train(n_epochs=1, lr=0.0008)
    concatenated_matrix = np.concatenate(
        (gene_dataset_fish.X[:, vae.indexes_to_keep], gene_dataset_seq.X[:, vae.indexes_to_keep]))
    concatenated_matrix = np.log(1 + concatenated_matrix)
    pca = PCA(n_components=9)
    latent_pca = pca.fit_transform(concatenated_matrix)
    pca_latent_fish = latent_pca[:gene_dataset_fish.X.shape[0], :]
    pca_latent_seq = latent_pca[gene_dataset_fish.X.shape[0]:, :]
    pca_values_seq = gene_dataset_seq.X
    pca_labels_seq = gene_dataset_seq.labels
    pca_labels_fish = gene_dataset_fish.labels
    _ = proximity_imputation(pca_latent_seq, pca_values_seq[:, 0], pca_latent_fish, k=5)
    _, _, = compute_accuracy_nn(pca_latent_seq, pca_labels_seq.ravel(), pca_latent_fish,
                                pca_labels_fish.ravel())


def benchmark_scanvi(source_datasets, target_dataset, params=dict(), max_acc=True):
    '''
    Annotate the target dataset from one, or multiple source datasets.
    Uses a warmup scheme for the classification only training the regular scVI model for the first epochs (faster)
    :param source_datasets: one dataset or a tuple / list of datasets
    :param target_dataset: dataset for which to predict the labels
    :param params: hyperparameters
    :param max_acc:
    :return:
    '''
    default_params = {
        'n_epoch_train_vae': 1,
        'lr_train_vae': 1e-3,
        'batch_size': 128,
        'n_layers': 2,
        'n_hidden': 256,
        'nb_genes': 100,
        'classifier_parameters': dict(),
        'weight_decay': 1e-6,

        'n_epoch_train_scanvi': 1,
        'lr_train_scanvi': 1e-4,
        'classification_ratio': 10,
        'lr_classification': 1e-3,
        'save_t_sne_folder': None
    }
    default_params.update(params)
    params = default_params
    print(params)

    if type(source_datasets) not in [tuple, list]:
        source_datasets = (source_datasets,)

    title = "%s->%s" % (' '.join([s.__class__.__name__ for s in source_datasets]), target_dataset.__class__.__name__)
    print(title)
    results = dict()

    dataset = GeneExpressionDataset.concat_datasets(target_dataset, *source_datasets)
    dataset.subsample_genes(new_n_genes=params['nb_genes'])
    if not max_acc:  # Then it is the Macosko-Regev dataset
        for group in labels_groups:
            dataset.merge_cell_types(group, group[0])

    data_loaders = SemiSupervisedDataLoaders(dataset, batch_size=params['batch_size'])
    data_loaders['unlabelled'] = data_loaders(indices=(dataset.batch_indices == 0).ravel().astype(np.bool))
    data_loaders['labelled'] = data_loaders['train'] = \
        data_loaders(indices=(dataset.batch_indices > 0).ravel().astype(np.bool))
    if max_acc:
        (_, labels_train), = data_loaders.raw_data(data_loaders['labelled'])
        (_, labels_test), = data_loaders.raw_data(data_loaders['unlabelled'])
        max_acc = np.mean([1 if l in np.unique(labels_train) else 0 for l in labels_test])
        print("Maximum Accuracy : ", max_acc)
        results.update({'max_acc': max_acc})

    # ~ equivalent to a warm-up for the classification
    vae = VAE(dataset.nb_genes, dataset.n_batches, dataset.n_labels,
              n_layers=params['n_layers'], n_hidden=params['n_hidden'], dropout_rate=0.1)
    infer = VariationalInference(vae, dataset, weight_decay=params['weight_decay'], verbose=True, frequency=50)

    infer.data_loaders = data_loaders
    infer.data_loaders.loop = ['all']
    infer.data_loaders['train'] = data_loaders['all']
    infer.train(params['n_epoch_train_vae'], lr=params['lr_train_vae'])
    print(infer.nn_latentspace('sequential'))
    results.update({
        'vae_latent_space_acc': infer.nn_latentspace('sequential'),
        'vae_entropy_batch_mixing': infer.entropy_batch_mixing('sequential'),
        'vae_clustering_scores': infer.clustering_scores('unlabelled')
    })

    scanvi = SCANVI(dataset.nb_genes, dataset.n_batches, dataset.n_labels,
                    n_layers=params['n_layers'], n_hidden=params['n_hidden'], dropout_rate=0.1,
                    classifier_parameters=params['classifier_parameters'])
    scanvi.load_state_dict(vae.state_dict(), strict=False)
    infer_scanvi = SemiSupervisedVariationalInference(
        scanvi, dataset, frequency=10, verbose=False, classification_ratio=params['classification_ratio'],
        n_epochs_classifier=1, lr_classification=params['lr_classification']
    )

    infer_scanvi.data_loaders = data_loaders
    data_loaders.loop = ['all', 'labelled']
    infer_scanvi.classifier_inference.data_loaders['labelled'] = data_loaders['labelled']
    infer_scanvi.classifier_inference.data_loaders.loop = ['labelled']


    infer_scanvi.train(params['n_epoch_train_scanvi'], lr=params['lr_train_scanvi'])

    results.update({
        'acc': infer_scanvi.accuracy('unlabelled'),
        'latent_space_acc': infer_scanvi.nn_latentspace('sequential'),
        'entropy_batch_mixing': infer_scanvi.entropy_batch_mixing('sequential'),
        'scanvi_clustering_scores': infer_scanvi.clustering_scores('unlabelled')
    })

    if params['save_t_sne_folder'] is not None:
        if not os.path.exists(params['save_t_sne_folder']):
            os.makedirs(params['save_t_sne_folder'])
        infer_scanvi.show_t_sne('sequential', color_by='batches and labels',
                                save_name=params['save_t_sne_folder'] + title + '.svg')

    return results
