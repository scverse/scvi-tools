import numpy as np
from sklearn.decomposition import PCA

from scvi.dataset import CortexDataset
from scvi.inference import UnsupervisedTrainer, TrainerFish
from scvi.inference.annotation import compute_accuracy_nn
from scvi.inference.posterior import proximity_imputation
from scvi.models import VAE, VAEF, LDVAE


def cortex_benchmark(n_epochs=250, use_cuda=True, save_path='data/', show_plot=True):
    cortex_dataset = CortexDataset(save_path=save_path, total_genes=558)
    vae = VAE(cortex_dataset.nb_genes)
    trainer_cortex_vae = UnsupervisedTrainer(vae, cortex_dataset, use_cuda=use_cuda)
    trainer_cortex_vae.train(n_epochs=n_epochs)
    couple_celltypes = (4, 5)  # the couple types on which to study DE
    cell_idx1 = cortex_dataset.labels.ravel() == couple_celltypes[0]
    cell_idx2 = cortex_dataset.labels.ravel() == couple_celltypes[1]
    trainer_cortex_vae.train_set.differential_expression_score(cell_idx1, cell_idx2,
                                                               genes=["THY1", "MBP"])

    trainer_cortex_vae.test_set.reconstruction_error()  # assert ~ 1200
    vae = VAE(cortex_dataset.nb_genes)
    trainer_cortex_vae = UnsupervisedTrainer(vae, cortex_dataset, use_cuda=use_cuda)
    trainer_cortex_vae.corrupt_posteriors()
    trainer_cortex_vae.train(n_epochs=n_epochs)
    trainer_cortex_vae.uncorrupt_posteriors()
    trainer_cortex_vae.train_set.imputation_benchmark(save_path=save_path, show_plot=show_plot)

    n_samples = 10 if n_epochs == 1 else None  # n_epochs == 1 is unit tests
    trainer_cortex_vae.train_set.show_t_sne(n_samples=n_samples)
    return trainer_cortex_vae


def benchmark(dataset, n_epochs=250, use_cuda=True):
    vae = VAE(dataset.nb_genes, n_batch=dataset.n_batches)
    trainer = UnsupervisedTrainer(vae, dataset, use_cuda=use_cuda)
    trainer.train(n_epochs=n_epochs)
    trainer.test_set.reconstruction_error()
    trainer.test_set.marginal_ll()
    return trainer


def harmonization_benchmarks(n_epochs=1, use_cuda=True, save_path='data/'):
    # retina_benchmark(n_epochs=n_epochs)
    pass


def annotation_benchmarks(n_epochs=1, use_cuda=True, save_path='data/'):
    # some cortex annotation benchmark
    pass


def ldvae_benchmark(dataset, n_epochs, use_cuda=True):
    ldvae = LDVAE(dataset.nb_genes, n_batch=dataset.n_batches)
    trainer = UnsupervisedTrainer(ldvae, dataset, use_cuda=use_cuda)
    trainer.train(n_epochs=n_epochs)
    trainer.test_set.reconstruction_error()
    trainer.test_set.marginal_ll()

    ldvae.get_loadings()

    return trainer


def all_benchmarks(n_epochs=250, use_cuda=True, save_path='data/', show_plot=True):
    cortex_benchmark(n_epochs=n_epochs, use_cuda=use_cuda, save_path=save_path, show_plot=show_plot)

    harmonization_benchmarks(n_epochs=n_epochs, use_cuda=use_cuda, save_path=save_path)
    annotation_benchmarks(n_epochs=n_epochs, use_cuda=use_cuda, save_path=save_path)


def benchmark_fish_scrna(gene_dataset_seq, gene_dataset_fish):
    gene_names = gene_dataset_fish.gene_names
    indexes_to_keep = np.arange(len(gene_names))
    vae = VAEF(gene_dataset_seq.nb_genes, indexes_to_keep, n_layers_decoder=2, n_latent=6,
               n_layers=2, n_hidden=256, reconstruction_loss='nb', dropout_rate=0.3, n_labels=7, n_batch=2,
               model_library=False)
    trainer = TrainerFish(vae, gene_dataset_seq, gene_dataset_fish, train_size=0.9,
                          frequency=5, weight_decay=0.35, n_epochs_even=100, n_epochs_kl=1000,
                          cl_ratio=0, n_epochs_cl=100)
    trainer.train(n_epochs=1, lr=0.0008)
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
