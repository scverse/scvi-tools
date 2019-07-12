from scvi.dataset import CortexDataset
from scvi.inference import UnsupervisedTrainer
from scvi.models import VAE, LDVAE


def cortex_benchmark(n_epochs=250, use_cuda=True, save_path="data/", show_plot=True):
    cortex_dataset = CortexDataset(save_path=save_path, total_genes=558)
    vae = VAE(cortex_dataset.nb_genes)
    trainer_cortex_vae = UnsupervisedTrainer(vae, cortex_dataset, use_cuda=use_cuda)
    trainer_cortex_vae.train(n_epochs=n_epochs)
    couple_celltypes = (4, 5)  # the couple types on which to study DE
    cell_idx1 = cortex_dataset.labels.ravel() == couple_celltypes[0]
    cell_idx2 = cortex_dataset.labels.ravel() == couple_celltypes[1]
    trainer_cortex_vae.train_set.differential_expression_score(
        cell_idx1, cell_idx2, genes=["THY1", "MBP"]
    )

    trainer_cortex_vae.test_set.reconstruction_error()  # assert ~ 1200
    vae = VAE(cortex_dataset.nb_genes)
    trainer_cortex_vae = UnsupervisedTrainer(vae, cortex_dataset, use_cuda=use_cuda)
    trainer_cortex_vae.corrupt_posteriors()
    trainer_cortex_vae.train(n_epochs=n_epochs)
    trainer_cortex_vae.uncorrupt_posteriors()
    trainer_cortex_vae.train_set.imputation_benchmark(
        save_path=save_path, show_plot=show_plot
    )

    n_samples = 10 if n_epochs == 1 else None  # n_epochs == 1 is unit tests
    trainer_cortex_vae.train_set.show_t_sne(n_samples=n_samples)
    return trainer_cortex_vae


def ldvae_benchmark(dataset, n_epochs, use_cuda=True):
    ldvae = LDVAE(dataset.nb_genes, n_batch=dataset.n_batches)
    trainer = UnsupervisedTrainer(ldvae, dataset, use_cuda=use_cuda)
    trainer.train(n_epochs=n_epochs)
    trainer.test_set.reconstruction_error()
    trainer.test_set.marginal_ll()

    ldvae.get_loadings()

    return trainer


def all_benchmarks(n_epochs=250, use_cuda=True, save_path="data/", show_plot=True):
    cortex_benchmark(
        n_epochs=n_epochs, use_cuda=use_cuda, save_path=save_path, show_plot=show_plot
    )
