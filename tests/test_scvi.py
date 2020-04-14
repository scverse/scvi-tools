import numpy as np
import pandas as pd
import tempfile
import os
import pytest
import torch

from anndata import AnnData

from scvi.dataset import (
    AnnDatasetFromAnnData,
    CortexDataset,
    SyntheticDataset,
    GeneExpressionDataset,
    Dataset10X,
)
from scvi.inference import (
    JointSemiSupervisedTrainer,
    AlternateSemiSupervisedTrainer,
    ClassifierTrainer,
    UnsupervisedTrainer,
    AdapterTrainer,
    TotalTrainer,
    TotalPosterior,
)
from scvi.inference.posterior import unsupervised_clustering_accuracy
from scvi.inference.posterior_utils import load_posterior
from scvi.inference.annotation import compute_accuracy_rf, compute_accuracy_svc
from scvi.models import VAE, SCANVI, VAEC, LDVAE, TOTALVI, AutoZIVAE
from scvi.models.distributions import ZeroInflatedNegativeBinomial, NegativeBinomial
from scvi.models.classifier import Classifier
from scvi.models.log_likelihood import log_zinb_positive, log_nb_positive
from scvi import set_seed

set_seed(0)
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
    trainer_cortex_vae.train_set.generate_feature_correlation_matrix(
        n_samples=2, correlation_type="pearson"
    )
    trainer_cortex_vae.train_set.generate_feature_correlation_matrix(
        n_samples=2, correlation_type="spearman"
    )
    genes = cortex_dataset.gene_names[:3]
    sample_scale = trainer_cortex_vae.train_set.get_sample_scale(gene_list=genes)
    assert type(sample_scale) == pd.DataFrame
    assert np.array_equal(np.sort(np.array(sample_scale.columns)), np.sort(genes))

    trainer_cortex_vae.train_set.imputation(n_samples=1)
    trainer_cortex_vae.test_set.imputation(n_samples=5)

    trainer_cortex_vae.corrupt_posteriors(corruption="binomial")
    trainer_cortex_vae.corrupt_posteriors()
    trainer_cortex_vae.train(n_epochs=1)
    trainer_cortex_vae.uncorrupt_posteriors()

    trainer_cortex_vae.train_set.imputation_benchmark(
        n_samples=1, show_plot=False, title_plot="imputation", save_path=save_path
    )
    trainer_cortex_vae.train_set.generate_parameters()

    n_cells, n_genes = (
        len(trainer_cortex_vae.train_set.indices),
        cortex_dataset.nb_genes,
    )
    n_samples = 3
    (dropout, means, dispersions) = trainer_cortex_vae.train_set.generate_parameters()
    assert dropout.shape == (n_cells, n_genes) and means.shape == (n_cells, n_genes)
    assert dispersions.shape == (n_cells, n_genes)
    (dropout, means, dispersions) = trainer_cortex_vae.train_set.generate_parameters(
        n_samples=n_samples
    )
    assert dropout.shape == (n_samples, n_cells, n_genes)
    assert means.shape == (n_samples, n_cells, n_genes)
    (dropout, means, dispersions) = trainer_cortex_vae.train_set.generate_parameters(
        n_samples=n_samples, give_mean=True
    )
    assert dropout.shape == (n_cells, n_genes) and means.shape == (n_cells, n_genes)

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

    with tempfile.TemporaryDirectory() as temp_dir:
        posterior_save_path = os.path.join(temp_dir, "posterior_data")
        original_post = trainer_synthetic_svaec.labelled_set.sequential()
        original_post.save_posterior(posterior_save_path)
        new_svaec = SCANVI(
            synthetic_dataset.nb_genes,
            synthetic_dataset.n_batches,
            synthetic_dataset.n_labels,
        )
        new_post = load_posterior(posterior_save_path, model=new_svaec, use_cuda=False)
    assert np.array_equal(new_post.indices, original_post.indices)
    assert np.array_equal(new_post.gene_dataset.X, original_post.gene_dataset.X)
    assert np.array_equal(
        new_post.gene_dataset.labels, original_post.gene_dataset.labels
    )

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
    ldvae = LDVAE(
        dataset.nb_genes, n_batch=dataset.n_batches, latent_distribution="normal"
    )
    trainer = UnsupervisedTrainer(ldvae, dataset, use_cuda=use_cuda)
    trainer.train(n_epochs=n_epochs)
    trainer.test_set.reconstruction_error()
    trainer.test_set.marginal_ll()

    ldvae = LDVAE(dataset.nb_genes, n_batch=dataset.n_batches, latent_distribution="ln")
    trainer = UnsupervisedTrainer(ldvae, dataset, use_cuda=use_cuda)
    trainer.train(n_epochs=n_epochs)
    trainer.test_set.reconstruction_error()

    ldvae.get_loadings()

    return trainer


def totalvi_benchmark(dataset, n_epochs, use_cuda=True):
    totalvae = TOTALVI(
        dataset.nb_genes, len(dataset.protein_names), n_batch=dataset.n_batches
    )
    trainer = TotalTrainer(
        totalvae, dataset, train_size=0.5, use_cuda=use_cuda, early_stopping_kwargs=None
    )
    trainer.train(n_epochs=n_epochs)
    trainer.test_set.reconstruction_error()
    trainer.test_set.marginal_ll()

    trainer.test_set.get_protein_background_mean()
    trainer.test_set.get_latent()
    trainer.test_set.generate()
    trainer.test_set.get_sample_dropout()
    trainer.test_set.get_normalized_denoised_expression(transform_batch=0)
    trainer.test_set.get_normalized_denoised_expression(transform_batch=0)
    trainer.test_set.imputation()
    trainer.test_set.get_protein_mean()
    trainer.test_set.one_vs_all_degenes(n_samples=2, M_permutation=10)
    trainer.test_set.generate_feature_correlation_matrix(n_samples=2)
    trainer.test_set.generate_feature_correlation_matrix(n_samples=2, transform_batch=0)

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


def test_poisson_not_zinb():
    synthetic_dataset = SyntheticDataset()
    svaec = SCANVI(
        synthetic_dataset.nb_genes,
        synthetic_dataset.n_batches,
        synthetic_dataset.n_labels,
        labels_groups=[0, 0, 1],
        reconstruction_loss="poisson",
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


def test_annealing_procedures(save_path):
    cortex_dataset = CortexDataset(save_path=save_path)
    cortex_vae = VAE(cortex_dataset.nb_genes, cortex_dataset.n_batches)

    trainer_cortex_vae = UnsupervisedTrainer(
        cortex_vae,
        cortex_dataset,
        train_size=0.5,
        use_cuda=use_cuda,
        n_epochs_kl_warmup=1,
    )
    trainer_cortex_vae.train(n_epochs=2)
    assert trainer_cortex_vae.kl_weight >= 0.99, "Annealing should be over"

    trainer_cortex_vae = UnsupervisedTrainer(
        cortex_vae,
        cortex_dataset,
        train_size=0.5,
        use_cuda=use_cuda,
        n_epochs_kl_warmup=5,
    )
    trainer_cortex_vae.train(n_epochs=2)
    assert trainer_cortex_vae.kl_weight <= 0.99, "Annealing should be proceeding"

    # iter
    trainer_cortex_vae = UnsupervisedTrainer(
        cortex_vae,
        cortex_dataset,
        train_size=0.5,
        use_cuda=use_cuda,
        n_iter_kl_warmup=1,
        n_epochs_kl_warmup=None,
    )
    trainer_cortex_vae.train(n_epochs=2)
    assert trainer_cortex_vae.kl_weight >= 0.99, "Annealing should be over"


def test_differential_expression(save_path):
    dataset = CortexDataset(save_path=save_path)
    n_cells = len(dataset)
    all_indices = np.arange(n_cells)
    vae = VAE(dataset.nb_genes, dataset.n_batches)
    trainer = UnsupervisedTrainer(vae, dataset, train_size=0.5, use_cuda=use_cuda)
    trainer.train(n_epochs=2)
    post = trainer.create_posterior(vae, dataset, shuffle=False, indices=all_indices)

    with tempfile.TemporaryDirectory() as temp_dir:
        posterior_save_path = os.path.join(temp_dir, "posterior_data")
        post = post.sequential(batch_size=3)
        post.save_posterior(posterior_save_path)
        new_vae = VAE(dataset.nb_genes, dataset.n_batches)
        new_post = load_posterior(posterior_save_path, model=new_vae, use_cuda=False)
    assert new_post.data_loader.batch_size == 3
    assert np.array_equal(new_post.indices, post.indices)
    assert np.array_equal(new_post.gene_dataset.X, post.gene_dataset.X)

    # Sample scale example
    px_scales = post.scale_sampler(
        n_samples_per_cell=4, n_samples=None, selection=all_indices
    )["scale"]
    assert (
        px_scales.shape[1] == dataset.nb_genes
    ), "posterior scales should have shape (n_samples, n_genes)"

    # Differential expression different models
    idx_1 = [1, 2, 3]
    idx_2 = [4, 5, 6, 7]
    de_dataframe = post.differential_expression_score(
        idx1=idx_1,
        idx2=idx_2,
        n_samples=10,
        mode="vanilla",
        use_permutation=True,
        M_permutation=100,
    )

    de_dataframe = post.differential_expression_score(
        idx1=idx_1,
        idx2=idx_2,
        n_samples=10,
        mode="change",
        use_permutation=True,
        M_permutation=100,
        cred_interval_lvls=[0.5, 0.95],
    )
    print(de_dataframe.keys())
    assert (
        de_dataframe["lfc_confidence_interval_0.5_min"]
        <= de_dataframe["lfc_confidence_interval_0.5_max"]
    ).all()
    assert (
        de_dataframe["lfc_confidence_interval_0.95_min"]
        <= de_dataframe["lfc_confidence_interval_0.95_max"]
    ).all()

    # DE estimation example
    de_probabilities = de_dataframe.loc[:, "proba_de"]
    assert ((0.0 <= de_probabilities) & (de_probabilities <= 1.0)).all()

    # Test totalVI DE
    sp = os.path.join(save_path, "10X")
    dataset = Dataset10X(dataset_name="pbmc_10k_protein_v3", save_path=sp)
    n_cells = len(dataset)
    all_indices = np.arange(n_cells)
    vae = TOTALVI(
        dataset.nb_genes, len(dataset.protein_names), n_batch=dataset.n_batches
    )
    trainer = TotalTrainer(
        vae, dataset, train_size=0.5, use_cuda=use_cuda, early_stopping_kwargs=None
    )
    trainer.train(n_epochs=2)
    post = trainer.create_posterior(
        vae, dataset, shuffle=False, indices=all_indices, type_class=TotalPosterior
    )

    # Differential expression different models
    idx_1 = [1, 2, 3]
    idx_2 = [4, 5, 6, 7]
    de_dataframe = post.differential_expression_score(
        idx1=idx_1,
        idx2=idx_2,
        n_samples=10,
        mode="vanilla",
        use_permutation=True,
        M_permutation=100,
    )

    de_dataframe = post.differential_expression_score(
        idx1=idx_1,
        idx2=idx_2,
        n_samples=10,
        mode="change",
        use_permutation=True,
        M_permutation=100,
    )


def test_totalvi(save_path):
    synthetic_dataset_one_batch = SyntheticDataset(n_batches=1)
    totalvi_benchmark(synthetic_dataset_one_batch, n_epochs=1, use_cuda=use_cuda)
    synthetic_dataset_two_batches = SyntheticDataset(n_batches=2)
    totalvi_benchmark(synthetic_dataset_two_batches, n_epochs=1, use_cuda=use_cuda)

    # adversarial testing
    dataset = synthetic_dataset_two_batches
    totalvae = TOTALVI(
        dataset.nb_genes, len(dataset.protein_names), n_batch=dataset.n_batches
    )
    trainer = TotalTrainer(
        totalvae,
        dataset,
        train_size=0.5,
        use_cuda=use_cuda,
        early_stopping_kwargs=None,
        use_adversarial_loss=True,
    )
    trainer.train(n_epochs=1)

    with tempfile.TemporaryDirectory() as temp_dir:
        posterior_save_path = os.path.join(temp_dir, "posterior_data")
        original_post = trainer.create_posterior(
            totalvae,
            dataset,
            indices=np.arange(len(dataset)),
            type_class=TotalPosterior,
        )
        original_post.save_posterior(posterior_save_path)
        new_totalvae = TOTALVI(
            dataset.nb_genes, len(dataset.protein_names), n_batch=dataset.n_batches
        )
        new_post = load_posterior(
            posterior_save_path, model=new_totalvae, use_cuda=False
        )
        assert hasattr(new_post.gene_dataset, "protein_names")
        assert new_post.posterior_type == "TotalPosterior"
        assert np.array_equal(
            new_post.gene_dataset.protein_expression, dataset.protein_expression
        )


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


def test_multibatches_features():
    data = [
        np.random.randint(1, 5, size=(20, 10)),
        np.random.randint(1, 10, size=(20, 10)),
        np.random.randint(1, 10, size=(20, 10)),
        np.random.randint(1, 10, size=(30, 10)),
    ]
    dataset = GeneExpressionDataset()
    dataset.populate_from_per_batch_list(data)
    vae = VAE(dataset.nb_genes, dataset.n_batches)
    trainer = UnsupervisedTrainer(vae, dataset, train_size=0.5, use_cuda=use_cuda)
    trainer.train(n_epochs=2)
    trainer.test_set.imputation(n_samples=2, transform_batch=0)
    trainer.train_set.imputation(n_samples=2, transform_batch=[0, 1, 2])


def test_deprecated_munkres():
    y = np.array([0, 1, 0, 1, 0, 1, 1, 1])
    y_pred = np.array([0, 0, 0, 0, 1, 1, 1, 1])
    reward, assignment = unsupervised_clustering_accuracy(y, y_pred)
    assert reward == 0.625
    assert (assignment == np.array([[0, 0], [1, 1]])).all()

    y = np.array([1, 1, 2, 2, 0, 0, 3, 3])
    y_pred = np.array([1, 1, 2, 2, 3, 3, 0, 0])
    reward, assignment = unsupervised_clustering_accuracy(y, y_pred)
    assert reward == 1.0
    assert (assignment == np.array([[0, 3], [1, 1], [2, 2], [3, 0]])).all()


def test_zinb_distribution():
    theta = 100.0 + torch.rand(size=(2,))
    mu = 15.0 * torch.ones_like(theta)
    pi = torch.randn_like(theta)
    x = torch.randint_like(mu, high=20)
    log_p_ref = log_zinb_positive(x, mu, theta, pi)

    dist = ZeroInflatedNegativeBinomial(mu=mu, theta=theta, zi_logits=pi)
    log_p_zinb = dist.log_prob(x)
    assert (log_p_ref - log_p_zinb).abs().max().item() <= 1e-8

    torch.manual_seed(0)
    s1 = dist.sample((100,))
    assert s1.shape == (100, 2)
    s2 = dist.sample(sample_shape=(4, 3))
    assert s2.shape == (4, 3, 2)

    log_p_ref = log_nb_positive(x, mu, theta)
    dist = NegativeBinomial(mu=mu, theta=theta)
    log_p_nb = dist.log_prob(x)
    assert (log_p_ref - log_p_nb).abs().max().item() <= 1e-8

    s1 = dist.sample((1000,))
    assert s1.shape == (1000, 2)
    assert (s1.mean(0) - mu).abs().mean() <= 1e0
    assert (s1.std(0) - (mu + mu * mu / theta) ** 0.5).abs().mean() <= 1e0

    size = (50, 3)
    theta = 100.0 + torch.rand(size=size)
    mu = 15.0 * torch.ones_like(theta)
    pi = torch.randn_like(theta)
    x = torch.randint_like(mu, high=20)
    dist1 = ZeroInflatedNegativeBinomial(mu=mu, theta=theta, zi_logits=pi)
    dist2 = NegativeBinomial(mu=mu, theta=theta)
    assert dist1.log_prob(x).shape == size
    assert dist2.log_prob(x).shape == size

    with pytest.raises(ValueError):
        ZeroInflatedNegativeBinomial(mu=-mu, theta=theta, zi_logits=pi)
    with pytest.warns(UserWarning):
        dist1.log_prob(-x)  # ensures neg values raise warning
    with pytest.warns(UserWarning):
        dist2.log_prob(0.5 * x)  # ensures float values raise warning


def test_anndata_loader():
    x = np.random.randint(low=0, high=100, size=(15, 4))
    batch_ids = np.random.randint(low=0, high=2, size=(15,))
    n_batches = 2
    adata = AnnData(X=x, obs=dict(batch=batch_ids))
    _ = AnnDatasetFromAnnData(adata, batch_label="batch")
    dataset = AnnDatasetFromAnnData(adata, batch_label="batch")
    assert (
        dataset.n_batches == n_batches
    ), "AnnDatasetFromAnnData should not modify the anndata object"
