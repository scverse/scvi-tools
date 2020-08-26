# import numpy as np
# import pytest
# import torch

# import scvi
# from scvi.core.models.vae import VAE
# from scvi.core.trainers.inference import UnsupervisedTrainer
# from scvi.core.trainers.annotation import ClassifierTrainer
# from scvi.core.models.classifier import Classifier

# scvi.set_seed(0)
# use_cuda = True


# def test_cortex(save_path):
#     cortex_dataset = scvi.dataset.cortex(save_path=save_path)
#     scvi.dataset.setup_anndata(cortex_dataset, labels_key="labels")
#     stats = cortex_dataset.uns["scvi_summary_stats"]
#     vae = VAE(stats["n_genes"], stats["n_batch"])
#     trainer_cortex_vae = UnsupervisedTrainer(
#         vae, cortex_dataset, train_size=0.5, use_cuda=use_cuda
#     )
#     trainer_cortex_vae.train(n_epochs=1)
#     trainer_cortex_vae.train_set.reconstruction_error()
#     trainer_cortex_vae.train_set.differential_expression_stats()
#     trainer_cortex_vae.train_set.generate_feature_correlation_matrix(
#         n_samples=2, correlation_type="pearson"
#     )
#     trainer_cortex_vae.train_set.generate_feature_correlation_matrix(
#         n_samples=2, correlation_type="spearman"
#     )
#     genes = cortex_dataset.var_names[:3]
#     sample_scale = trainer_cortex_vae.train_set.get_sample_scale(gene_list=genes)
#     assert type(sample_scale) == pd.DataFrame
#     assert np.array_equal(np.sort(np.array(sample_scale.columns)), np.sort(genes))

#     trainer_cortex_vae.train_set.imputation(n_samples=1)
#     trainer_cortex_vae.test_set.imputation(n_samples=5)

#     trainer_cortex_vae.train_set.generate_parameters()

#     n_cells, n_genes = (
#         len(trainer_cortex_vae.train_set.indices),
#         cortex_dataset.uns["scvi_summary_stats"]["n_genes"],
#     )
#     n_samples = 3
#     (dropout, means, dispersions) = trainer_cortex_vae.train_set.generate_parameters()
#     assert dropout.shape == (n_cells, n_genes) and means.shape == (n_cells, n_genes)
#     assert dispersions.shape == (n_cells, n_genes)
#     (dropout, means, dispersions) = trainer_cortex_vae.train_set.generate_parameters(
#         n_samples=n_samples
#     )
#     assert dropout.shape == (n_samples, n_cells, n_genes)
#     assert means.shape == (n_samples, n_cells, n_genes)
#     (dropout, means, dispersions) = trainer_cortex_vae.train_set.generate_parameters(
#         n_samples=n_samples, give_mean=True
#     )
#     assert dropout.shape == (n_cells, n_genes) and means.shape == (n_cells, n_genes)

#     full = trainer_cortex_vae.create_posterior(
#         vae, cortex_dataset, indices=np.arange(len(cortex_dataset))
#     )
#     x_new, x_old = full.generate(n_samples=10)
#     assert x_new.shape == (
#         cortex_dataset.uns["scvi_summary_stats"]["n_cells"],
#         cortex_dataset.uns["scvi_summary_stats"]["n_genes"],
#         10,
#     )
#     assert x_old.shape == (
#         cortex_dataset.uns["scvi_summary_stats"]["n_cells"],
#         cortex_dataset.uns["scvi_summary_stats"]["n_genes"],
#     )

#     svaec = SCANVI(stats["n_genes"], stats["n_batch"], stats["n_labels"])
#     trainer_cortex_svaec = JointSemiSupervisedTrainer(
#         svaec, cortex_dataset, n_labelled_samples_per_class=3, use_cuda=use_cuda
#     )
#     trainer_cortex_svaec.train(n_epochs=1)
#     trainer_cortex_svaec.labelled_set.accuracy()
#     trainer_cortex_svaec.full_dataset.reconstruction_error()

#     svaec = SCANVI(stats["n_genes"], stats["n_batch"], stats["n_labels"])
#     trainer_cortex_svaec = AlternateSemiSupervisedTrainer(
#         svaec, cortex_dataset, n_labelled_samples_per_class=3, use_cuda=use_cuda
#     )
#     trainer_cortex_svaec.train(n_epochs=1, lr=1e-2)
#     trainer_cortex_svaec.unlabelled_set.accuracy()
#     data_train, labels_train = trainer_cortex_svaec.labelled_set.raw_data()
#     data_test, labels_test = trainer_cortex_svaec.unlabelled_set.raw_data()
#     compute_accuracy_svc(
#         data_train,
#         labels_train,
#         data_test,
#         labels_test,
#         param_grid=[{"C": [1], "kernel": ["linear"]}],
#     )
#     compute_accuracy_rf(
#         data_train,
#         labels_train,
#         data_test,
#         labels_test,
#         param_grid=[{"max_depth": [3], "n_estimators": [10]}],
#     )

#     cls = Classifier(stats["n_genes"], n_labels=stats["n_labels"])
#     cls_trainer = ClassifierTrainer(cls, cortex_dataset)
#     cls_trainer.train(n_epochs=1)
#     cls_trainer.train_set.accuracy()


# def test_synthetic_1():
#     synthetic_dataset = scvi.dataset.synthetic_iid()
#     scvi.dataset.setup_anndata(
#         synthetic_dataset, batch_key="batch", labels_key="labels"
#     )
#     stats = synthetic_dataset.uns["scvi_summary_stats"]

#     svaec = SCANVI(stats["n_genes"], stats["n_batch"], stats["n_labels"])
#     trainer_synthetic_svaec = JointSemiSupervisedTrainer(
#         svaec, synthetic_dataset, use_cuda=use_cuda
#     )
#     trainer_synthetic_svaec.train(n_epochs=1)
#     trainer_synthetic_svaec.labelled_set.entropy_batch_mixing()

#     with tempfile.TemporaryDirectory() as temp_dir:
#         posterior_save_path = os.path.join(temp_dir, "posterior_data")
#         original_post = trainer_synthetic_svaec.labelled_set.sequential()
#         original_post.save_posterior(posterior_save_path)
#         new_svaec = SCANVI(stats["n_genes"], stats["n_batch"], stats["n_labels"])
#         new_post = load_posterior(posterior_save_path, model=new_svaec, use_cuda=False)
#     assert np.array_equal(new_post.indices, original_post.indices)
#     assert np.array_equal(new_post.gene_dataset.X, original_post.gene_dataset.X)
#     assert np.array_equal(
#         new_post.gene_dataset.labels, original_post.gene_dataset.labels
#     )
#     trainer_synthetic_svaec.full_dataset.knn_purity()
#     trainer_synthetic_svaec.labelled_set.clustering_scores()
#     trainer_synthetic_svaec.labelled_set.clustering_scores(prediction_algorithm="gmm")
#     trainer_synthetic_svaec.unlabelled_set.unsupervised_classification_accuracy()
#     trainer_synthetic_svaec.unlabelled_set.differential_expression_score(
#         synthetic_dataset.obs["_scvi_labels"].ravel() == 1,
#         synthetic_dataset.obs["_scvi_labels"].ravel() == 2,
#         n_samples=2,
#         M_permutation=10,
#     )
#     trainer_synthetic_svaec.unlabelled_set.one_vs_all_degenes(
#         n_samples=2, M_permutation=10
#     )


# def test_synthetic_2():
#     synthetic_dataset = scvi.dataset.synthetic_iid()
#     scvi.dataset.setup_anndata(
#         synthetic_dataset, batch_key="batch", labels_key="labels"
#     )
#     stats = synthetic_dataset.uns["scvi_summary_stats"]

#     vaec = VAEC(stats["n_genes"], stats["n_batch"], stats["n_labels"])
#     trainer_synthetic_vaec = JointSemiSupervisedTrainer(
#         vaec,
#         synthetic_dataset,
#         use_cuda=use_cuda,
#         frequency=1,
#         early_stopping_kwargs={
#             "early_stopping_metric": "reconstruction_error",
#             "on": "labelled_set",
#             "save_best_state_metric": "reconstruction_error",
#         },
#     )
#     trainer_synthetic_vaec.train(n_epochs=2)


# def base_benchmark(adata):
#     stats = adata.uns["scvi_summary_stats"]
#     vae = VAE(stats["n_genes"], stats["n_batch"], stats["n_labels"])
#     trainer = UnsupervisedTrainer(vae, adata, train_size=0.5, use_cuda=use_cuda)
#     trainer.train(n_epochs=1)
#     return trainer


# def ldvae_benchmark(dataset, n_epochs, use_cuda=True):
#     ldvae = LDVAE(
#         dataset.uns["scvi_summary_stats"]["n_genes"],
#         n_batch=dataset.uns["scvi_summary_stats"]["n_batch"],
#         latent_distribution="normal",
#     )
#     trainer = UnsupervisedTrainer(ldvae, dataset, use_cuda=use_cuda)
#     trainer.train(n_epochs=n_epochs)
#     trainer.test_set.reconstruction_error()
#     trainer.test_set.marginal_ll()

#     ldvae = LDVAE(
#         dataset.uns["scvi_summary_stats"]["n_genes"],
#         n_batch=dataset.uns["scvi_summary_stats"]["n_batch"],
#         latent_distribution="ln",
#     )
#     trainer = UnsupervisedTrainer(ldvae, dataset, use_cuda=use_cuda)
#     trainer.train(n_epochs=n_epochs)
#     trainer.test_set.reconstruction_error()

#     ldvae.get_loadings()

#     return trainer


# def totalvi_benchmark(dataset, n_epochs, use_cuda=True):
#     totalvae = TOTALVI(
#         dataset.uns["scvi_summary_stats"]["n_genes"],
#         dataset.uns["scvi_summary_stats"]["n_proteins"],
#         n_batch=dataset.uns["scvi_summary_stats"]["n_batch"],
#     )
#     trainer = TotalTrainer(
#         totalvae, dataset, train_size=0.5, use_cuda=use_cuda, early_stopping_kwargs=None
#     )
#     trainer.train(n_epochs=n_epochs)
#     trainer.test_set.reconstruction_error()
#     trainer.test_set.marginal_ll()

#     trainer.test_set.get_protein_background_mean()
#     trainer.test_set.get_latent()
#     trainer.test_set.generate()
#     trainer.test_set.get_sample_dropout()
#     trainer.test_set.get_normalized_denoised_expression(transform_batch=0)
#     trainer.test_set.get_normalized_denoised_expression(transform_batch=0)
#     trainer.test_set.imputation()
#     trainer.test_set.get_protein_mean()
#     trainer.test_set.one_vs_all_degenes(n_samples=2, M_permutation=10)
#     trainer.test_set.generate_feature_correlation_matrix(n_samples=2)
#     trainer.test_set.generate_feature_correlation_matrix(n_samples=2, transform_batch=0)

#     return trainer


# def test_synthetic_3():
#     adata = scvi.dataset.synthetic_iid()
#     scvi.dataset.setup_anndata(adata, batch_key="batch", labels_key="labels")
#     trainer = base_benchmark(adata)
#     adapter_trainer = AdapterTrainer(
#         trainer.model, adata, trainer.train_set, frequency=1
#     )
#     adapter_trainer.train(n_path=1, n_epochs=1)


# def test_nb_not_zinb():
#     synthetic_dataset = scvi.dataset.synthetic_iid()
#     scvi.dataset.setup_anndata(
#         synthetic_dataset, batch_key="batch", labels_key="labels"
#     )
#     svaec = SCANVI(
#         synthetic_dataset.uns["scvi_summary_stats"]["n_genes"],
#         synthetic_dataset.uns["scvi_summary_stats"]["n_batch"],
#         synthetic_dataset.uns["scvi_summary_stats"]["n_labels"],
#         labels_groups=[0, 0, 1],
#         reconstruction_loss="nb",
#     )
#     trainer_synthetic_svaec = JointSemiSupervisedTrainer(
#         svaec, synthetic_dataset, use_cuda=use_cuda
#     )
#     trainer_synthetic_svaec.train(n_epochs=1)


# def test_poisson_not_zinb():
#     synthetic_dataset = scvi.dataset.synthetic_iid()
#     scvi.dataset.setup_anndata(
#         synthetic_dataset, batch_key="batch", labels_key="labels"
#     )
#     svaec = SCANVI(
#         synthetic_dataset.uns["scvi_summary_stats"]["n_genes"],
#         synthetic_dataset.uns["scvi_summary_stats"]["n_batch"],
#         synthetic_dataset.uns["scvi_summary_stats"]["n_labels"],
#         labels_groups=[0, 0, 1],
#         reconstruction_loss="poisson",
#     )
#     trainer_synthetic_svaec = JointSemiSupervisedTrainer(
#         svaec, synthetic_dataset, use_cuda=use_cuda
#     )
#     trainer_synthetic_svaec.train(n_epochs=1)


# def test_classifier_accuracy(save_path):
#     cortex_dataset = scvi.dataset.cortex(save_path=save_path)
#     scvi.dataset.setup_anndata(cortex_dataset, labels_key="labels")
#     cls = Classifier(
#         cortex_dataset.uns["scvi_summary_stats"]["n_genes"],
#         n_labels=cortex_dataset.uns["scvi_summary_stats"]["n_labels"],
#     )
#     cls_trainer = ClassifierTrainer(
#         cls,
#         cortex_dataset,
#         metrics_to_monitor=["accuracy"],
#         frequency=1,
#         early_stopping_kwargs={
#             "early_stopping_metric": "accuracy",
#             "save_best_state_metric": "accuracy",
#         },
#     )
#     cls_trainer.train(n_epochs=2)
#     cls_trainer.train_set.accuracy()


# def test_LDVAE(save_path):
#     synthetic_datset_one_batch = scvi.dataset.synthetic_iid(n_batches=1)
#     scvi.dataset.setup_anndata(synthetic_datset_one_batch, batch_key="batch")
#     ldvae_benchmark(synthetic_datset_one_batch, n_epochs=1, use_cuda=False)
#     synthetic_datset_two_batches = scvi.dataset.synthetic_iid(n_batches=2)
#     scvi.dataset.setup_anndata(synthetic_datset_two_batches, batch_key="batch")
#     ldvae_benchmark(synthetic_datset_two_batches, n_epochs=1, use_cuda=False)


# def test_differential_expression(save_path):
#     dataset = scvi.dataset.cortex(save_path=save_path)
#     scvi.dataset.setup_anndata(dataset, labels_key="cell_type")
#     n_cells = len(dataset)
#     all_indices = np.arange(n_cells)
#     vae = VAE(
#         dataset.uns["scvi_summary_stats"]["n_genes"],
#         dataset.uns["scvi_summary_stats"]["n_batch"],
#     )
#     trainer = UnsupervisedTrainer(vae, dataset, train_size=0.5, use_cuda=use_cuda)
#     trainer.train(n_epochs=2)
#     post = trainer.create_posterior(vae, dataset, shuffle=False, indices=all_indices)

#     with tempfile.TemporaryDirectory() as temp_dir:
#         posterior_save_path = os.path.join(temp_dir, "posterior_data")
#         post = post.sequential(batch_size=3)
#         post.save_posterior(posterior_save_path)
#         new_vae = VAE(
#             dataset.uns["scvi_summary_stats"]["n_genes"],
#             dataset.uns["scvi_summary_stats"]["n_batch"],
#         )
#         new_post = load_posterior(posterior_save_path, model=new_vae, use_cuda=False)
#     assert new_post.data_loader.sampler.batch_size == 3
#     assert np.array_equal(new_post.indices, post.indices)
#     assert np.array_equal(new_post.gene_dataset.X, post.gene_dataset.X)

#     # Sample scale example
#     px_scales = post.scale_sampler(
#         n_samples_per_cell=4, n_samples=None, selection=all_indices
#     )["scale"]
#     assert (
#         px_scales.shape[1] == dataset.uns["scvi_summary_stats"]["n_genes"]
#     ), "posterior scales should have shape (n_samples, n_genes)"

#     # Differential expression different models
#     idx_1 = [1, 2, 3]
#     idx_2 = [4, 5, 6, 7]
#     de_dataframe = post.differential_expression_score(
#         idx1=idx_1,
#         idx2=idx_2,
#         n_samples=10,
#         mode="vanilla",
#         use_permutation=True,
#         M_permutation=100,
#     )

#     de_dataframe = post.differential_expression_score(
#         idx1=idx_1,
#         idx2=idx_2,
#         n_samples=10,
#         mode="change",
#         use_permutation=True,
#         M_permutation=100,
#         cred_interval_lvls=[0.5, 0.95],
#     )
#     print(de_dataframe.keys())
#     assert (
#         de_dataframe["lfc_confidence_interval_0.5_min"]
#         <= de_dataframe["lfc_confidence_interval_0.5_max"]
#     ).all()
#     assert (
#         de_dataframe["lfc_confidence_interval_0.95_min"]
#         <= de_dataframe["lfc_confidence_interval_0.95_max"]
#     ).all()

#     # DE estimation example
#     de_probabilities = de_dataframe.loc[:, "proba_de"]
#     assert ((0.0 <= de_probabilities) & (de_probabilities <= 1.0)).all()

#     # Test totalVI DE
#     sp = os.path.join(save_path, "10X")
#     dataset = scvi.dataset.dataset10X(
#         dataset_name="pbmc_10k_protein_v3", save_path=sp, gex_only=False
#     )
#     scvi.dataset.organize_cite_seq_10x(dataset)
#     setup_anndata(dataset, protein_expression_obsm_key="protein_expression")

#     n_cells = len(dataset)
#     all_indices = np.arange(n_cells)
#     vae = TOTALVI(
#         dataset.uns["scvi_summary_stats"]["n_genes"],
#         dataset.uns["scvi_summary_stats"]["n_proteins"],
#         n_batch=dataset.uns["scvi_summary_stats"]["n_batch"],
#     )
#     trainer = TotalTrainer(
#         vae, dataset, train_size=0.5, use_cuda=use_cuda, early_stopping_kwargs=None
#     )
#     trainer.train(n_epochs=2)
#     post = trainer.create_posterior(
#         vae, dataset, shuffle=False, indices=all_indices, type_class=TotalPosterior
#     )

#     # Differential expression different models
#     idx_1 = [1, 2, 3]
#     idx_2 = [4, 5, 6, 7]
#     de_dataframe = post.differential_expression_score(
#         idx1=idx_1,
#         idx2=idx_2,
#         n_samples=10,
#         mode="vanilla",
#         use_permutation=True,
#         M_permutation=100,
#     )

#     de_dataframe = post.differential_expression_score(
#         idx1=idx_1,
#         idx2=idx_2,
#         n_samples=10,
#         mode="change",
#         use_permutation=True,
#         M_permutation=100,
#     )


# def test_totalvi(save_path):
#     synthetic_dataset_one_batch = scvi.dataset.synthetic_iid(n_batches=1)
#     scvi.dataset.setup_anndata(
#         synthetic_dataset_one_batch,
#         protein_expression_obsm_key="protein_expression",
#         protein_names_uns_key="protein_names",
#         labels_key="labels",
#     )
#     totalvi_benchmark(synthetic_dataset_one_batch, n_epochs=1, use_cuda=use_cuda)
#     synthetic_dataset_two_batches = scvi.dataset.synthetic_iid(n_batches=2)
#     scvi.dataset.setup_anndata(
#         synthetic_dataset_two_batches,
#         batch_key="batch",
#         protein_expression_obsm_key="protein_expression",
#         protein_names_uns_key="protein_names",
#         labels_key="labels",
#     )
#     totalvi_benchmark(synthetic_dataset_two_batches, n_epochs=1, use_cuda=use_cuda)

#     # adversarial testing
#     dataset = synthetic_dataset_two_batches
#     totalvae = TOTALVI(
#         dataset.uns["scvi_summary_stats"]["n_genes"],
#         dataset.uns["scvi_summary_stats"]["n_proteins"],
#         n_batch=dataset.uns["scvi_summary_stats"]["n_batch"],
#     )
#     trainer = TotalTrainer(
#         totalvae,
#         dataset,
#         train_size=0.5,
#         use_cuda=use_cuda,
#         early_stopping_kwargs=None,
#         use_adversarial_loss=True,
#     )
#     trainer.train(n_epochs=1)

#     with tempfile.TemporaryDirectory() as temp_dir:
#         posterior_save_path = os.path.join(temp_dir, "posterior_data")
#         original_post = trainer.create_posterior(
#             totalvae,
#             dataset,
#             indices=np.arange(len(dataset)),
#             type_class=TotalPosterior,
#         )
#         original_post.save_posterior(posterior_save_path)
#         new_totalvae = TOTALVI(
#             dataset.uns["scvi_summary_stats"]["n_genes"],
#             dataset.uns["scvi_summary_stats"]["n_proteins"],
#             n_batch=dataset.uns["scvi_summary_stats"]["n_batch"],
#         )
#         new_post = load_posterior(
#             posterior_save_path, model=new_totalvae, use_cuda=False
#         )
#         assert hasattr(new_post.gene_dataset, "protein_names")
#         assert new_post.posterior_type == "TotalPosterior"
#         assert np.array_equal(
#             new_post.gene_dataset.protein_expression, dataset.obsm["protein_expression"]
#         )


# def test_autozi(save_path):
#     data = scvi.dataset.synthetic_iid(n_batches=1)
#     scvi.dataset.setup_anndata(data, batch_key="batch", labels_key="labels")

#     for disp_zi in ["gene", "gene-label"]:
#         autozivae = AutoZIVAE(
#             n_input=data.uns["scvi_summary_stats"]["n_genes"],
#             dispersion=disp_zi,
#             zero_inflation=disp_zi,
#             n_labels=data.uns["scvi_summary_stats"]["n_labels"],
#         )
#         trainer_autozivae = UnsupervisedTrainer(
#             model=autozivae, adata=data, train_size=0.5
#         )
#         trainer_autozivae.train(n_epochs=2, lr=1e-2)
#         trainer_autozivae.test_set.elbo()
#         trainer_autozivae.test_set.reconstruction_error()
#         trainer_autozivae.test_set.marginal_ll()


# def test_multibatches_features():
#     dataset = scvi.dataset.synthetic_iid(n_batches=3)
#     scvi.dataset.setup_anndata(dataset, batch_key="batch", labels_key="labels")

#     vae = VAE(
#         dataset.uns["scvi_summary_stats"]["n_genes"],
#         dataset.uns["scvi_summary_stats"]["n_batch"],
#     )
#     trainer = UnsupervisedTrainer(vae, dataset, train_size=0.5, use_cuda=use_cuda)
#     trainer.train(n_epochs=2)
#     trainer.test_set.imputation(n_samples=2, transform_batch=0)
#     trainer.train_set.imputation(n_samples=2, transform_batch=[0, 1, 2])
