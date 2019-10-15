import numpy as np

from scvi.dataset import CortexDataset, SyntheticDataset
from scvi.inference import (
    JointSemiSupervisedTrainer,
    AlternateSemiSupervisedTrainer,
    ClassifierTrainer,
    UnsupervisedTrainer,
    AdapterTrainer,
    TotalTrainer,
)
from scvi.inference.annotation import compute_accuracy_rf, compute_accuracy_svc
from scvi.models import VAE, SCANVI, VAEC, LDVAE, TOTALVI, AutoZIVAE
from scvi.models.classifier import Classifier

use_cuda = True


def test_cortex(save_path):
    cortex_dataset = CortexDataset(save_path=save_path)
    vae = VAE(cortex_dataset.nb_genes, cortex_dataset.n_batches)
    trainer_cortex_vae = UnsupervisedTrainer(
        vae, cortex_dataset, train_size=0.5, use_cuda=use_cuda
    )
    trainer_cortex_vae.train(n_epochs=1)
    trainer_cortex_vae.train_set.reconstruction_error()
    trainer_cortex_vae.train_set.differential_expression_stats()

    trainer_cortex_vae.corrupt_posteriors(corruption="binomial")
    trainer_cortex_vae.corrupt_posteriors()
    trainer_cortex_vae.train(n_epochs=1)
    trainer_cortex_vae.uncorrupt_posteriors()

    trainer_cortex_vae.train_set.imputation_benchmark(
        n_samples=1, show_plot=False, title_plot="imputation", save_path=save_path
    )
    full = trainer_cortex_vae.create_posterior(
        vae, cortex_dataset, indices=np.arange(len(cortex_dataset))
    )
    x_new, x_old = full.generate(n_samples=10)
    assert x_new.shape == (cortex_dataset.nb_cells, cortex_dataset.nb_genes, 10)
    assert x_old.shape == (cortex_dataset.nb_cells, cortex_dataset.nb_genes)

    trainer_cortex_vae.train_set.imputation_benchmark(
        n_samples=1, show_plot=False, title_plot="imputation", save_path=save_path
    )

    svaec = SCANVI(
        cortex_dataset.nb_genes, cortex_dataset.n_batches, cortex_dataset.n_labels
    )
    trainer_cortex_svaec = JointSemiSupervisedTrainer(
        svaec, cortex_dataset, n_labelled_samples_per_class=3, use_cuda=use_cuda
    )
    trainer_cortex_svaec.train(n_epochs=1)
    trainer_cortex_svaec.labelled_set.accuracy()
    trainer_cortex_svaec.full_dataset.reconstruction_error()

    svaec = SCANVI(
        cortex_dataset.nb_genes, cortex_dataset.n_batches, cortex_dataset.n_labels
    )
    trainer_cortex_svaec = AlternateSemiSupervisedTrainer(
        svaec, cortex_dataset, n_labelled_samples_per_class=3, use_cuda=use_cuda
    )
    trainer_cortex_svaec.train(n_epochs=1, lr=1e-2)
    trainer_cortex_svaec.unlabelled_set.accuracy()
    data_train, labels_train = trainer_cortex_svaec.labelled_set.raw_data()
    data_test, labels_test = trainer_cortex_svaec.unlabelled_set.raw_data()
    compute_accuracy_svc(
        data_train,
        labels_train,
        data_test,
        labels_test,
        param_grid=[{"C": [1], "kernel": ["linear"]}],
    )
    compute_accuracy_rf(
        data_train,
        labels_train,
        data_test,
        labels_test,
        param_grid=[{"max_depth": [3], "n_estimators": [10]}],
    )

    cls = Classifier(cortex_dataset.nb_genes, n_labels=cortex_dataset.n_labels)
    cls_trainer = ClassifierTrainer(cls, cortex_dataset)
    cls_trainer.train(n_epochs=1)
    cls_trainer.train_set.accuracy()


def test_synthetic_1():
    synthetic_dataset = SyntheticDataset()
    synthetic_dataset.cell_types = np.array(["A", "B", "C"])
    svaec = SCANVI(
        synthetic_dataset.nb_genes,
        synthetic_dataset.n_batches,
        synthetic_dataset.n_labels,
    )
    trainer_synthetic_svaec = JointSemiSupervisedTrainer(
        svaec, synthetic_dataset, use_cuda=use_cuda
    )
    trainer_synthetic_svaec.train(n_epochs=1)
    trainer_synthetic_svaec.labelled_set.entropy_batch_mixing()
    trainer_synthetic_svaec.full_dataset.knn_purity()
    trainer_synthetic_svaec.labelled_set.show_t_sne(n_samples=5)
    trainer_synthetic_svaec.unlabelled_set.show_t_sne(n_samples=5, color_by="labels")
    trainer_synthetic_svaec.labelled_set.show_t_sne(
        n_samples=5, color_by="batches and labels"
    )
    trainer_synthetic_svaec.labelled_set.clustering_scores()
    trainer_synthetic_svaec.labelled_set.clustering_scores(prediction_algorithm="gmm")
    trainer_synthetic_svaec.unlabelled_set.unsupervised_classification_accuracy()
    trainer_synthetic_svaec.unlabelled_set.differential_expression_score(
        synthetic_dataset.labels.ravel() == 1,
        synthetic_dataset.labels.ravel() == 2,
        genes=["2", "4"],
        n_samples=2,
        M_permutation=10,
    )
    trainer_synthetic_svaec.unlabelled_set.one_vs_all_degenes(
        n_samples=2, M_permutation=10
    )


def test_synthetic_2():
    synthetic_dataset = SyntheticDataset()
    vaec = VAEC(
        synthetic_dataset.nb_genes,
        synthetic_dataset.n_batches,
        synthetic_dataset.n_labels,
    )
    trainer_synthetic_vaec = JointSemiSupervisedTrainer(
        vaec,
        synthetic_dataset,
        use_cuda=use_cuda,
        frequency=1,
        early_stopping_kwargs={
            "early_stopping_metric": "reconstruction_error",
            "on": "labelled_set",
            "save_best_state_metric": "reconstruction_error",
        },
    )
    trainer_synthetic_vaec.train(n_epochs=2)


def base_benchmark(gene_dataset):
    vae = VAE(gene_dataset.nb_genes, gene_dataset.n_batches, gene_dataset.n_labels)
    trainer = UnsupervisedTrainer(vae, gene_dataset, train_size=0.5, use_cuda=use_cuda)
    trainer.train(n_epochs=1)
    return trainer


def ldvae_benchmark(dataset, n_epochs, use_cuda=True):
    ldvae = LDVAE(dataset.nb_genes, n_batch=dataset.n_batches)
    trainer = UnsupervisedTrainer(ldvae, dataset, use_cuda=use_cuda)
    trainer.train(n_epochs=n_epochs)
    trainer.test_set.reconstruction_error()
    trainer.test_set.marginal_ll()

    ldvae.get_loadings()

    return trainer


def totalvi_benchmark(dataset, n_epochs, use_cuda=True):
    totalvae = TOTALVI(
        dataset.nb_genes, len(dataset.protein_names), n_batch=dataset.n_batches
    )
    trainer = TotalTrainer(totalvae, dataset, train_size=0.5, use_cuda=use_cuda)
    trainer.train(n_epochs=n_epochs)
    trainer.test_set.reconstruction_error()
    trainer.test_set.marginal_ll()

    trainer.test_set.get_protein_background_mean()
    trainer.test_set.get_latent()
    trainer.test_set.generate()
    trainer.test_set.get_sample_dropout()
    trainer.test_set.get_normalized_denoised_expression()
    trainer.test_set.imputation()

    return trainer


def test_synthetic_3():
    gene_dataset = SyntheticDataset()
    trainer = base_benchmark(gene_dataset)
    adapter_trainer = AdapterTrainer(
        trainer.model, gene_dataset, trainer.train_set, frequency=1
    )
    adapter_trainer.train(n_path=1, n_epochs=1)


def test_nb_not_zinb():
    synthetic_dataset = SyntheticDataset()
    svaec = SCANVI(
        synthetic_dataset.nb_genes,
        synthetic_dataset.n_batches,
        synthetic_dataset.n_labels,
        labels_groups=[0, 0, 1],
        reconstruction_loss="nb",
    )
    trainer_synthetic_svaec = JointSemiSupervisedTrainer(
        svaec, synthetic_dataset, use_cuda=use_cuda
    )
    trainer_synthetic_svaec.train(n_epochs=1)


def test_classifier_accuracy(save_path):
    cortex_dataset = CortexDataset(save_path=save_path)
    cls = Classifier(cortex_dataset.nb_genes, n_labels=cortex_dataset.n_labels)
    cls_trainer = ClassifierTrainer(
        cls,
        cortex_dataset,
        metrics_to_monitor=["accuracy"],
        frequency=1,
        early_stopping_kwargs={
            "early_stopping_metric": "accuracy",
            "save_best_state_metric": "accuracy",
        },
    )
    cls_trainer.train(n_epochs=2)
    cls_trainer.train_set.accuracy()


def test_LDVAE(save_path):
    synthetic_datset_one_batch = SyntheticDataset(n_batches=1)
    ldvae_benchmark(synthetic_datset_one_batch, n_epochs=1, use_cuda=False)
    synthetic_datset_two_batches = SyntheticDataset(n_batches=2)
    ldvae_benchmark(synthetic_datset_two_batches, n_epochs=1, use_cuda=False)


def test_sampling_zl(save_path):
    cortex_dataset = CortexDataset(save_path=save_path)
    cortex_vae = VAE(cortex_dataset.nb_genes, cortex_dataset.n_batches)
    trainer_cortex_vae = UnsupervisedTrainer(
        cortex_vae, cortex_dataset, train_size=0.5, use_cuda=use_cuda
    )
    trainer_cortex_vae.train(n_epochs=2)

    cortex_cls = Classifier((cortex_vae.n_latent + 1), n_labels=cortex_dataset.n_labels)
    trainer_cortex_cls = ClassifierTrainer(
        cortex_cls, cortex_dataset, sampling_model=cortex_vae, sampling_zl=True
    )
    trainer_cortex_cls.train(n_epochs=2)
    trainer_cortex_cls.test_set.accuracy()


def test_totalvi(save_path):
    synthetic_dataset_one_batch = SyntheticDataset(n_batches=1)
    totalvi_benchmark(synthetic_dataset_one_batch, n_epochs=1, use_cuda=use_cuda)
    synthetic_dataset_two_batches = SyntheticDataset(n_batches=2)
    totalvi_benchmark(synthetic_dataset_two_batches, n_epochs=1, use_cuda=use_cuda)


def test_autozi(save_path):
    data = SyntheticDataset(n_batches=1)

    for disp_zi in ["gene", "gene-label"]:
        autozivae = AutoZIVAE(
            n_input=data.nb_genes,
            dispersion=disp_zi,
            zero_inflation=disp_zi,
            n_labels=data.n_labels,
        )
        trainer_autozivae = UnsupervisedTrainer(
            model=autozivae, gene_dataset=data, train_size=0.5
        )
        trainer_autozivae.train(n_epochs=2, lr=1e-2)
        trainer_autozivae.test_set.elbo()
        trainer_autozivae.test_set.reconstruction_error()
        trainer_autozivae.test_set.marginal_ll()
