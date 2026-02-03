import anndata as ad
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import sparse

from scvi.external.drvi import DRVI
from scvi.external.drvi.utils.plotting import (
    plot_latent_dimension_stats,
    plot_latent_dims_in_heatmap,
    show_top_differential_vars,
)
from scvi.external.drvi.utils.tools import (
    calculate_differential_vars,
    iterate_on_top_differential_vars,
    set_latent_dimension_stats,
    traverse_latent,
)


class TestInterpretabilityPipeline:
    n = 50
    g = 20
    c = 2
    b = 5

    def make_test_adata(self, is_sparse=True):
        N, G, C, B = self.n, self.g, self.c, self.b

        ct_list = np.sort(np.random.choice(range(C), N))[:, np.newaxis]
        g_exp_list = np.sort(np.random.choice(range(C), G))[:, np.newaxis]
        ct_array = (np.indices((N, C))[1] == ct_list) + 0.0
        g_exp_array = (np.indices((G, C))[1] == g_exp_list) + 0.0

        batch_list = np.random.choice(range(B), N)[:, np.newaxis]
        batch_list_2 = np.random.choice(range(B), N)[:, np.newaxis]
        exp_indicator = ct_array @ g_exp_array.T
        g_mean_list = np.exp((np.random.random(G) + 0.5) * 2)[:, np.newaxis]

        exp_matrix = np.random.poisson(exp_indicator * g_mean_list.T).astype(np.float32)

        adata = ad.AnnData(
            X=sparse.csr_matrix(exp_matrix),
            obs=pd.DataFrame(
                {
                    "cell_type": [f"ct_{ct}" for ct in ct_list[:, 0]],
                    "batch": [f"batch_{bid}" for bid in batch_list[:, 0]],
                    "batch_2": [f"batch_{bid}" for bid in batch_list_2[:, 0]],
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

    def _whole_integration_and_interpretability_pipeline(self, model_kwargs=None):
        adata = self.make_test_adata()

        DRVI.setup_anndata(
            adata,
            batch_key="batch",
            layer="counts",
            is_count_data=True,
        )
        _model_kwargs = dict(  # noqa: C408
            n_latent=8,
            encoder_dims=[128],
            decoder_dims=[128],
            gene_likelihood="pnb",
            decoder_reuse_weights="everywhere",
        )
        if model_kwargs is not None:
            _model_kwargs.update(model_kwargs)
        model = DRVI(adata, **_model_kwargs)
        model.train(accelerator="cpu", max_epochs=100)
        embed = ad.AnnData(model.get_latent_representation(), obs=adata.obs)

        set_latent_dimension_stats(model, embed)

        traverse_adata = traverse_latent(model, embed, n_samples=2, max_noise_std=0.0)
        calculate_differential_vars(traverse_adata)

        return {
            "adata": adata,
            "model": model,
            "embed": embed,
            "traverse_adata": traverse_adata,
        }

    def test_whole_integration_and_interpretability_pipeline(self):
        """Test the complete interpretability pipeline."""
        self._whole_integration_and_interpretability_pipeline()

    def test_set_latent_dimension_stats(self):
        """Test set_latent_dimension_stats function."""
        adata = self.make_test_adata()
        DRVI.setup_anndata(
            adata,
            batch_key="batch",
            layer="counts",
            is_count_data=True,
        )
        model = DRVI(
            adata,
            n_latent=8,
            encoder_dims=[128],
            decoder_dims=[128],
            gene_likelihood="pnb",
            decoder_reuse_weights="everywhere",
        )
        model.train(accelerator="cpu", max_epochs=100)
        embed = ad.AnnData(model.get_latent_representation(), obs=adata.obs)

        # Test set_latent_dimension_stats
        set_latent_dimension_stats(model, embed)

        # Check that required columns are added
        assert "original_dim_id" in embed.var
        assert "reconstruction_effect" in embed.var
        assert "order" in embed.var
        assert "max_value" in embed.var
        assert "mean" in embed.var
        assert "min" in embed.var
        assert "max" in embed.var
        assert "std" in embed.var
        assert "title" in embed.var
        assert "vanished" in embed.var

        # Check that order is correct
        assert embed.var["order"].min() == 0
        assert embed.var["order"].max() == embed.n_vars - 1

    def test_traverse_latent(self):
        """Test traverse_latent function."""
        train_results = self._whole_integration_and_interpretability_pipeline()
        model = train_results["model"]
        embed = train_results["embed"]

        # Test traverse_latent
        traverse_adata = traverse_latent(model, embed, n_samples=2, n_steps=4, max_noise_std=0.0)

        # Check that traverse_adata has the expected structure
        assert "control" in traverse_adata.layers
        assert "effect" in traverse_adata.layers
        assert "control_latent" in traverse_adata.obsm
        assert "effect_latent" in traverse_adata.obsm
        assert "dim_id" in traverse_adata.obs
        assert "sample_id" in traverse_adata.obs
        assert "step_id" in traverse_adata.obs
        assert "span_value" in traverse_adata.obs

    def test_calculate_differential_vars(self):
        """Test calculate_differential_vars function."""
        train_results = self._whole_integration_and_interpretability_pipeline()
        model = train_results["model"]
        embed = train_results["embed"]

        traverse_adata = traverse_latent(model, embed, n_samples=2, n_steps=4, max_noise_std=0.0)

        # Test calculate_differential_vars
        calculate_differential_vars(traverse_adata)

        # Check that required keys are added
        assert "max_possible_traverse_effect_stepwise" in traverse_adata.uns
        assert "min_possible_traverse_effect_stepwise" in traverse_adata.uns
        assert "combined_score_traverse_effect_stepwise" in traverse_adata.uns
        assert "max_possible_traverse_effect_pos" in traverse_adata.varm
        assert "max_possible_traverse_effect_neg" in traverse_adata.varm
        assert "min_possible_traverse_effect_pos" in traverse_adata.varm
        assert "min_possible_traverse_effect_neg" in traverse_adata.varm
        assert "combined_score_traverse_effect_pos" in traverse_adata.varm
        assert "combined_score_traverse_effect_neg" in traverse_adata.varm

    def test_iterate_on_top_differential_vars(self):
        """Test iterate_on_top_differential_vars function."""
        train_results = self._whole_integration_and_interpretability_pipeline()
        model = train_results["model"]
        embed = train_results["embed"]
        traverse_adata = train_results["traverse_adata"]

        # Test iterate_on_top_differential_vars
        dimensions_interpretability = iterate_on_top_differential_vars(
            traverse_adata, key="combined_score", score_threshold=0.0
        )

        # Check that we get a list of tuples
        assert isinstance(dimensions_interpretability, list)
        assert len(dimensions_interpretability) > 0

        # Check structure of each element
        for dim_title, gene_scores in dimensions_interpretability:
            assert isinstance(dim_title, str)
            assert isinstance(gene_scores, pd.Series)
            assert len(gene_scores) > 0

    def test_plotting_functions(self):
        """Test plotting functions."""
        train_results = self._whole_integration_and_interpretability_pipeline()
        embed = train_results["embed"]
        traverse_adata = train_results["traverse_adata"]

        # Test plot_latent_dimension_stats
        plot_latent_dimension_stats(embed, ncols=2, show=False)
        plt.close()

        # Test plot_latent_dims_in_heatmap
        plot_latent_dims_in_heatmap(embed, "cell_type", title_col="title", show=False)
        plt.close()

        # Test show_top_differential_vars
        show_top_differential_vars(traverse_adata, key="combined_score", score_threshold=0.0, show=False)
        plt.close()

        # Test iterate_on_top_differential_vars (already tested above, but also test with plotting)
        dimensions_interpretability = iterate_on_top_differential_vars(
            traverse_adata, key="combined_score", score_threshold=0.0
        )

        # If we have dimensions, test plotting with a subset
        if len(dimensions_interpretability) > 0:
            sample_dim = dimensions_interpretability[0][0]
            show_top_differential_vars(
                traverse_adata,
                key="combined_score",
                dim_subset=[sample_dim],
                score_threshold=0.0,
                show=False,
            )
            plt.close()

