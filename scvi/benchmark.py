from scvi.dataset import CortexDataset
from scvi.inference import VariationalInference, VariationalInferenceFish, adversarial_wrapper
from scvi.models import VAE, VAEF
import numpy as np
from sklearn.decomposition import PCA
from scvi.metrics.imputation import proximity_imputation
from scvi.metrics.classification import compute_accuracy_nn


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
