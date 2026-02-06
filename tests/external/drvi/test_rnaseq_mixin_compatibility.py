import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from scvi.external.drvi import DRVI


class TestRNASeqMixinCompatibility:
    """Test DRVI compatibility with scvi-tools RNASeqMixin methods."""

    n = 100  # Number of cells
    g = 50  # Number of genes
    c = 3  # Number of cell types
    b = 4  # Number of batches

    def _make_test_adata(self, is_sparse=True):
        """Create test AnnData object with RNA-seq like data."""
        N, G, C, B = self.n, self.g, self.c, self.b

        # Generate cell type assignments
        ct_list = np.sort(np.random.choice(range(C), N))[:, np.newaxis]
        g_exp_list = np.sort(np.random.choice(range(C), G))[:, np.newaxis]
        ct_array = (np.indices((N, C))[1] == ct_list) + 0.0
        g_exp_array = (np.indices((G, C))[1] == g_exp_list) + 0.0

        # Generate batch assignments
        batch_list = np.random.choice(range(B), N)[:, np.newaxis]
        exp_indicator = ct_array @ g_exp_array.T
        g_mean_list = np.exp((np.random.random(G) + 0.5) * 2)[:, np.newaxis]

        # Generate expression matrix (negative binomial-like)
        exp_matrix = np.random.negative_binomial(
            n=5, p=1 / (1 + exp_indicator * g_mean_list.T)
        ).astype(np.float32)

        adata = ad.AnnData(
            X=sparse.csr_matrix(exp_matrix) if is_sparse else exp_matrix,
            obs=pd.DataFrame(
                {
                    "cell_type": [f"ct_{ct}" for ct in ct_list[:, 0]],
                    "batch": [f"batch_{bid}" for bid in batch_list[:, 0]],
                },
                index=[f"cell_{i}" for i in range(N)],
            ),
            var=pd.DataFrame(
                {
                    "gene_mean": g_mean_list[:, 0],
                    "gene_active_signature": np.apply_along_axis(
                        lambda x: "".join(x), axis=1, arr=g_exp_array.astype(str)
                    ),
                },
                index=[f"gene_{i}" for i in range(G)],
            ),
        )

        adata.obs["total_counts"] = adata.X.sum(axis=1)
        adata.layers["counts"] = adata.X.copy()
        adata.layers["lognorm"] = np.log1p(adata.X)

        if not is_sparse:
            adata.X = adata.X.A
            for l in ["counts", "lognorm"]:
                adata.layers[l] = adata.layers[l].A

        return adata

    def _setup_and_train_model(
        self, adata, n_latent=8, max_epochs=10, gene_likelihood="nb", model_kwargs=None
    ):
        """Setup and train a DRVI model for testing."""
        # Setup model
        DRVI.setup_anndata(
            adata,
            categorical_covariate_keys=["batch"],
            layer="lognorm" if gene_likelihood.startswith("normal") else "counts",
            is_count_data=False if gene_likelihood.startswith("normal") else True,
        )

        _model_kwargs = dict(  # noqa: C408
            n_latent=n_latent,
            encoder_dims=[64],
            decoder_dims=[64],
            gene_likelihood=gene_likelihood,
        )
        if model_kwargs is not None:
            _model_kwargs.update(model_kwargs)
        model = DRVI(adata, **_model_kwargs)

        # Train briefly
        model.train(accelerator="cpu", max_epochs=max_epochs)

        return model

    def test_rnaseq_mixin_inheritance(self):
        """Test that DRVI properly inherits from RNASeqMixin."""
        from scvi.model.base import RNASeqMixin

        assert issubclass(DRVI, RNASeqMixin), "DRVI should inherit from RNASeqMixin"

        # Check that RNASeqMixin methods are available
        required_methods = [
            "get_normalized_expression",
            "differential_expression",
            "posterior_predictive_sample",
            "get_likelihood_parameters",
            "get_latent_library_size",
        ]

        for method in required_methods:
            assert hasattr(DRVI, method), f"DRVI missing method: {method}"

    def test_get_normalized_expression(self):
        """Test get_normalized_expression method."""
        adata = self._make_test_adata()
        model = self._setup_and_train_model(adata, max_epochs=2)

        # Test get_normalized_expression
        normalized_expr = model.get_normalized_expression(
            adata=adata, n_samples=1, return_mean=True, return_numpy=True
        )

        assert normalized_expr.shape == (adata.n_obs, adata.n_vars), (
            f"Expected shape {(adata.n_obs, adata.n_vars)}, got {normalized_expr.shape}"
        )
        assert np.all(normalized_expr >= 0), "Normalized expression should be non-negative"

        # Test with gene subset
        gene_subset = adata.var_names[:10]
        normalized_expr_subset = model.get_normalized_expression(
            adata=adata, gene_list=gene_subset, n_samples=1, return_mean=True, return_numpy=True
        )

        assert normalized_expr_subset.shape == (adata.n_obs, len(gene_subset)), (
            f"Expected shape {(adata.n_obs, len(gene_subset))}, got {normalized_expr_subset.shape}"
        )

    def test_differential_expression(self):
        """Test differential_expression method."""
        adata = self._make_test_adata()
        model = self._setup_and_train_model(adata, max_epochs=2)

        # Test differential expression
        de_results = model.differential_expression(
            adata=adata, groupby="cell_type", group1="ct_0", group2="ct_1", batch_size=32
        )

        print(de_results.columns)

        assert isinstance(de_results, pd.DataFrame), "DE results should be a DataFrame"
        assert len(de_results) == adata.n_vars, "DE results should have one row per gene"
        assert "bayes_factor" in de_results.columns

    def test_posterior_predictive_sample(self):
        """Test posterior_predictive_sample method."""
        adata = self._make_test_adata()
        model = self._setup_and_train_model(adata, max_epochs=2)

        # Test posterior predictive sampling
        samples = model.posterior_predictive_sample(adata=adata, n_samples=3, batch_size=32)

        assert samples.shape == (adata.n_obs, adata.n_vars, 3), (
            f"Expected shape {(adata.n_obs, adata.n_vars, 3)}, got {samples.shape}"
        )
        assert np.all(samples >= 0), "Samples should be non-negative"

        # Test with gene subset
        gene_subset = adata.var_names[:10]
        samples_subset = model.posterior_predictive_sample(
            adata=adata, gene_list=gene_subset, n_samples=2, batch_size=32
        )

        assert samples_subset.shape == (adata.n_obs, len(gene_subset), 2), (
            f"Expected shape {(adata.n_obs, len(gene_subset), 2)}, got {samples_subset.shape}"
        )

    def test_get_likelihood_parameters(self):
        """Test get_likelihood_parameters method."""
        adata = self._make_test_adata()
        model = self._setup_and_train_model(adata, max_epochs=2)

        for n_samples in [1, 2]:
            # Test likelihood parameters
            likelihood_params = model.get_likelihood_parameters(
                adata=adata, n_samples=n_samples, give_mean=True, batch_size=32
            )

            assert isinstance(likelihood_params, dict), "Likelihood parameters should be a dict"
            assert "mean" in likelihood_params, "Should contain mean parameter"
            assert "dispersions" in likelihood_params, "Should contain dispersions parameter"

            mean = likelihood_params["mean"]
            dispersions = likelihood_params["dispersions"]

            assert mean.shape == (adata.n_obs, adata.n_vars), (
                f"Mean shape should be {(adata.n_obs, adata.n_vars)}, got {mean.shape}"
            )
            assert dispersions.shape == (adata.n_obs, adata.n_vars), (
                f"Dispersions shape should be {(adata.n_obs, adata.n_vars)}, got {dispersions.shape}"
            )

            assert np.all(mean >= 0), "Mean should be non-negative"
            assert np.all(dispersions > 0), "Dispersions should be positive"

    def test_get_latent_library_size(self):
        """Test get_latent_library_size method."""
        adata = self._make_test_adata()
        model = self._setup_and_train_model(adata, max_epochs=2)

        # Test latent library size
        # Note: DRVI doesn't compute posterior distribution for library size (ql),
        # so we use give_mean=False to use the observed library size instead
        library_size = model.get_latent_library_size(adata=adata, give_mean=False, batch_size=32)

        assert library_size.shape == (adata.n_obs,), (
            f"Expected shape {(adata.n_obs,)}, got {library_size.shape}"
        )
        assert np.all(library_size > 0), "Library size should be positive"

    def test_get_feature_correlation_matrix(self):
        """Test get_feature_correlation_matrix method."""
        adata = self._make_test_adata()
        model = self._setup_and_train_model(adata, max_epochs=2)

        # Test feature correlation matrix
        corr_matrix = model.get_feature_correlation_matrix(
            adata=adata, n_samples=5, batch_size=32, correlation_type="pearson"
        )

        assert isinstance(corr_matrix, pd.DataFrame), "Correlation matrix should be a DataFrame"
        assert corr_matrix.shape == (adata.n_vars, adata.n_vars), (
            f"Expected shape {(adata.n_vars, adata.n_vars)}, got {corr_matrix.shape}"
        )
        assert np.allclose(corr_matrix.values, corr_matrix.values.T), (
            "Correlation matrix should be symmetric"
        )
        assert np.allclose(np.diag(corr_matrix.values), 1.0), "Diagonal should be 1.0"

    def _test_module_methods(self):
        """Test that DRVI module has required methods."""
        from scvi.external.drvi.module import DRVIModule

        # Test required methods
        required_methods = ["forward", "inference", "_get_inference_input", "generative", "sample"]
        for method in required_methods:
            assert hasattr(DRVIModule, method), f"DRVIModule missing method: {method}"

        # Test generative method signature
        import inspect

        gen_sig = inspect.signature(DRVIModule.generative)
        assert "transform_batch" in gen_sig.parameters, (
            "generative method missing transform_batch parameter"
        )

    def _test_module_properties(self):
        """Test that DRVI module has required properties."""
        from scvi.external.drvi.module import DRVIModule

        # Test instance properties
        module = DRVIModule(n_input=50, n_latent=8, gene_likelihood="nb")
        assert hasattr(module, "gene_likelihood"), (
            "Module instance missing gene_likelihood property"
        )
        assert module.gene_likelihood == "nb", (
            f"Expected gene_likelihood='nb', got '{module.gene_likelihood}'"
        )

    def test_module_compatibility(self):
        """Test that DRVI module has required methods and properties."""
        self._test_module_methods()
        self._test_module_properties()

    def test_end_to_end_rnaseq_workflow(self):
        """Test a complete RNA-seq analysis workflow using RNASeqMixin methods."""
        adata = self._make_test_adata()
        model = self._setup_and_train_model(adata, max_epochs=2)

        # Get latent representation
        latent = model.get_latent_representation(adata)
        assert latent.shape == (adata.n_obs, 8), (
            f"Expected latent shape {(adata.n_obs, 8)}, got {latent.shape}"
        )

        # Get normalized expression
        normalized_expr = model.get_normalized_expression(
            adata=adata, n_samples=1, return_mean=True, return_numpy=True
        )

        # Perform differential expression analysis
        de_results = model.differential_expression(
            adata=adata,
            groupby="cell_type",
            group1="ct_0",
            group2="ct_1",
        )

        # Get likelihood parameters
        likelihood_params = model.get_likelihood_parameters(
            adata=adata, n_samples=1, give_mean=True
        )

        # Get latent library size
        # Note: DRVI doesn't compute posterior distribution for library size (ql),
        # so we use give_mean=False to use the observed library size instead
        library_size = model.get_latent_library_size(adata=adata, give_mean=False)

        # Verify all results are consistent
        assert latent.shape[0] == normalized_expr.shape[0] == library_size.shape[0], (
            "All results should have same number of cells"
        )
        assert normalized_expr.shape[1] == likelihood_params["mean"].shape[1] == adata.n_vars, (
            "All results should have same number of genes"
        )
        assert len(de_results) == adata.n_vars, "DE results should have one row per gene"

        print("✓ End-to-end RNA-seq workflow completed successfully!")
        print(f"  - Latent representation: {latent.shape}")
        print(f"  - Normalized expression: {normalized_expr.shape}")
        print(f"  - DE results: {len(de_results)} genes")
        print(f"  - Library size: {library_size.shape}")
        print(f"  - Likelihood parameters: {list(likelihood_params.keys())}")

    def _end_to_end_rnaseq_workflow(self, **kwargs):
        """Test a complete RNA-seq analysis workflow using RNASeqMixin methods."""
        adata = self._make_test_adata()
        model = self._setup_and_train_model(adata, **kwargs)

        # Get latent representation
        latent = model.get_latent_representation(adata)
        assert latent.shape == (adata.n_obs, 8), (
            f"Expected latent shape {(adata.n_obs, 8)}, got {latent.shape}"
        )

        # Get normalized expression
        normalized_expr = model.get_normalized_expression(
            adata=adata, n_samples=1, return_mean=True, return_numpy=True
        )

        # Perform differential expression analysis
        de_results = model.differential_expression(
            adata=adata,
            groupby="cell_type",
            group1="ct_0",
            group2="ct_1",
        )

        # Get likelihood parameters
        likelihood_params = model.get_likelihood_parameters(
            adata=adata, n_samples=1, give_mean=True
        )

        # Get latent library size
        # Note: DRVI doesn't compute posterior distribution for library size (ql),
        # so we use give_mean=False to use the observed library size instead
        library_size = model.get_latent_library_size(adata=adata, give_mean=False)

        # Verify all results are consistent
        assert latent.shape[0] == normalized_expr.shape[0] == library_size.shape[0], (
            "All results should have same number of cells"
        )
        assert normalized_expr.shape[1] == likelihood_params["mean"].shape[1] == adata.n_vars, (
            "All results should have same number of genes"
        )
        assert len(de_results) == adata.n_vars, "DE results should have one row per gene"

        print("✓ End-to-end RNA-seq workflow completed successfully!")
        print(f"  - Latent representation: {latent.shape}")
        print(f"  - Normalized expression: {normalized_expr.shape}")
        print(f"  - DE results: {len(de_results)} genes")
        print(f"  - Library size: {library_size.shape}")
        print(f"  - Likelihood parameters: {list(likelihood_params.keys())}")

    def test_end_to_end_rnaseq_workflow_for_different_gene_likelihoods(self):
        """Test a complete RNA-seq analysis workflow using RNASeqMixin methods."""
        for gene_likelihood in ["nb", "poisson", "pnb", "normal_sv"]:
            print(f"Testing gene likelihood: {gene_likelihood}")
            self._end_to_end_rnaseq_workflow(max_epochs=2, gene_likelihood=gene_likelihood)
