import logging
from typing import List, Optional, Sequence, Union

import jax
import jax.numpy as jnp
import numpy as np
import numpyro.distributions as dist
import optax
import tqdm
from anndata import AnnData
from flax.training import train_state
from jax import random

from scvi import REGISTRY_KEYS
from scvi._compat import Literal
from scvi.data.anndata import AnnDataManager
from scvi.data.anndata.fields import (
    CategoricalJointObsField,
    CategoricalObsField,
    LayerField,
    NumericalJointObsField,
    NumericalObsField,
)
from scvi.dataloaders import DataSplitter
from scvi.model._utils import _init_library_size
from scvi.module import JaxVAE
from scvi.utils import setup_anndata_dsp

from .base import BaseModelClass

logger = logging.getLogger(__name__)


class JaxSCVI(BaseModelClass):
    """
    single-cell Variational Inference [Lopez18]_, but with a Jax backend.

    Parameters
    ----------
    adata
        AnnData object that has been registered via :meth:`~scvi.model.SCVI.setup_anndata`.
    n_hidden
        Number of nodes per hidden layer.
    n_latent
        Dimensionality of the latent space.
    n_layers
        Number of hidden layers used for encoder and decoder NNs.
    dropout_rate
        Dropout rate for neural networks.
    dispersion
        One of the following:

        * ``'gene'`` - dispersion parameter of NB is constant per gene across cells
        * ``'gene-batch'`` - dispersion can differ between different batches
        * ``'gene-label'`` - dispersion can differ between different labels
        * ``'gene-cell'`` - dispersion can differ for every gene in every cell
    gene_likelihood
        One of:

        * ``'nb'`` - Negative binomial distribution
        * ``'zinb'`` - Zero-inflated negative binomial distribution
        * ``'poisson'`` - Poisson distribution
    latent_distribution
        One of:

        * ``'normal'`` - Normal distribution
        * ``'ln'`` - Logistic normal distribution (Normal(0, I) transformed by softmax)
    **model_kwargs
        Keyword args for :class:`~scvi.module.JaxVAE`

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.model.JaxSCVI.setup_anndata(adata, batch_key="batch")
    >>> vae = scvi.model.SCVI(adata)
    >>> vae.train()
    >>> adata.obsm["X_scVI"] = vae.get_latent_representation()
    >>> adata.obsm["X_normalized_scVI"] = vae.get_normalized_expression()
    """

    def __init__(
        self,
        adata: AnnData,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell"] = "gene",
        gene_likelihood: Literal["zinb", "nb", "poisson"] = "zinb",
        latent_distribution: Literal["normal", "ln"] = "normal",
        **model_kwargs,
    ):
        super().__init__(adata)

        n_cats_per_cov = (
            self.adata_manager.get_state_registry(
                REGISTRY_KEYS.CAT_COVS_KEY
            ).n_cats_per_key
            if REGISTRY_KEYS.CAT_COVS_KEY in self.adata_manager.data_registry
            else None
        )
        n_batch = self.summary_stats.n_batch
        use_size_factor_key = (
            REGISTRY_KEYS.SIZE_FACTOR_KEY in self.adata_manager.data_registry
        )
        library_log_means, library_log_vars = None, None
        if not use_size_factor_key:
            library_log_means, library_log_vars = _init_library_size(
                self.adata_manager, n_batch
            )

        self.module_kwargs = dict(
            n_input=self.summary_stats.n_vars,
            n_batch=n_batch,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            is_training=False,
        )
        self.module_kwargs.update(model_kwargs)
        self._module = None

        self._model_summary_string = (
            "SCVI Model with the following params: \nn_hidden: {}, n_latent: {}, n_layers: {}, dropout_rate: "
            "{}, dispersion: {}, gene_likelihood: {}, latent_distribution: {}"
        ).format(
            n_hidden,
            n_latent,
            n_layers,
            dropout_rate,
            dispersion,
            gene_likelihood,
            latent_distribution,
        )
        self.init_params_ = self._get_init_params(locals())

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        layer: Optional[str] = None,
        batch_key: Optional[str] = None,
        labels_key: Optional[str] = None,
        size_factor_key: Optional[str] = None,
        categorical_covariate_keys: Optional[List[str]] = None,
        continuous_covariate_keys: Optional[List[str]] = None,
        **kwargs,
    ):
        """
        %(summary)s.

        Parameters
        ----------
        %(param_layer)s
        %(param_batch_key)s
        %(param_labels_key)s
        %(param_size_factor_key)s
        %(param_cat_cov_keys)s
        %(param_cont_cov_keys)s
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            CategoricalObsField(REGISTRY_KEYS.LABELS_KEY, labels_key),
            NumericalObsField(
                REGISTRY_KEYS.SIZE_FACTOR_KEY, size_factor_key, required=False
            ),
            CategoricalJointObsField(
                REGISTRY_KEYS.CAT_COVS_KEY, categorical_covariate_keys
            ),
            NumericalJointObsField(
                REGISTRY_KEYS.CONT_COVS_KEY, continuous_covariate_keys
            ),
        ]
        adata_manager = AnnDataManager(
            fields=anndata_fields, setup_method_args=setup_method_args
        )
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)

    def _get_module(self, kwargs=None):
        if kwargs is None:
            kwargs = self.module_kwargs
        return JaxVAE(**kwargs)

    @property
    def module(self):
        if self._module is None:
            self._module = self._get_module()
        return self._module

    def train(
        self,
        max_epochs: Optional[int] = None,
        use_gpu: Optional[Union[str, int, bool]] = None,
        train_size: float = 0.9,
        validation_size: Optional[float] = None,
        batch_size: int = 128,
        **trainer_kwargs,
    ):
        """
        Train the model.

        Parameters
        ----------
        max_epochs
            Number of passes through the dataset. If `None`, defaults to
            `np.min([round((20000 / n_cells) * 400), 400])`
        use_gpu
            Use default GPU if available (if None or True), or index of GPU to use (if int),
            or name of GPU (if str, e.g., `'cuda:0'`), or use CPU (if False).
        train_size
            Size of training set in the range [0.0, 1.0].
        validation_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + validation_size < 1`, the remaining cells belong to a test set.
        batch_size
            Minibatch size to use during training.
        **trainer_kwargs
            Other keyword args for :class:`~scvi.train.Trainer`.
        """
        if max_epochs is None:
            n_cells = self.adata.n_obs
            max_epochs = np.min([round((20000 / n_cells) * 400), 400])

        data_splitter = DataSplitter(
            self.adata_manager,
            train_size=train_size,
            validation_size=validation_size,
            batch_size=batch_size,
            use_gpu=use_gpu,
            iter_ndarray=True,
        )
        data_splitter.setup()
        train_loader = data_splitter.train_dataloader()

        module_kwargs = self.module_kwargs.copy()
        module_kwargs.update(dict(is_training=True))
        module = self._get_module(module_kwargs)

        self.rngs = {
            "params": random.PRNGKey(0),
            "dropout": random.PRNGKey(1),
            "z": random.PRNGKey(2),
        }
        params = module.init(self.rngs, next(iter(train_loader)))["params"]

        state = train_state.TrainState.create(
            apply_fn=module.apply,
            params=params,
            tx=optax.adam(1e-3),
        )

        @jax.jit
        def train_step(state, array_dict, rngs, kl_weight=1) -> jnp.ndarray:
            x = array_dict[REGISTRY_KEYS.X_KEY]
            rngs = {k: random.split(v)[1] for k, v in rngs.items()}

            def loss_fn(params):
                outputs = module.apply({"params": params}, array_dict, rngs=rngs)
                log_likelihood = outputs.nb.log_prob(x).sum(-1)
                mean = outputs.mean
                scale = outputs.stddev
                prior = dist.Normal(jnp.zeros_like(mean), jnp.ones_like(scale))
                posterior = dist.Normal(mean, scale)
                kl = dist.kl_divergence(posterior, prior).sum(-1)
                elbo = log_likelihood - kl_weight * kl
                return -jnp.mean(elbo)

            value, grads = jax.value_and_grad(loss_fn)(state.params)
            state = state.apply_gradients(grads=grads)
            return state, value, rngs

        history = []
        epoch = 0
        with tqdm.trange(1, max_epochs + 1) as t:
            for i in t:
                epoch += 1
                epoch_loss = 0
                counter = 0
                for data in train_loader:
                    kl_weight = min(1.0, epoch / 400.0)
                    # gets new key for each epoch
                    state, value, self.rngs = train_step(
                        state, data, self.rngs, kl_weight=kl_weight
                    )
                    epoch_loss += value
                    counter += 1
                history += [jax.device_get(epoch_loss) / counter]
                t.set_postfix_str(
                    f"Epoch {i}, Loss: {epoch_loss / counter}, KL weight: {kl_weight}"
                )

        self.train_state = state
        self.params = state.params
        self.history_ = history
        self.is_trained_ = True

        self.module_kwargs.update(dict(is_training=False))
        self._module = None

    def get_latent_representation(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        give_mean: bool = True,
        mc_samples: int = 5000,
        batch_size: Optional[int] = None,
    ) -> np.ndarray:
        r"""
        Return the latent representation for each cell.

        This is denoted as :math:`z_n` in our manuscripts.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        give_mean
            Give mean of distribution or sample from it.
        mc_samples
            For distributions with no closed-form mean (e.g., `logistic normal`), how many Monte Carlo
            samples to take for computing mean.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.

        Returns
        -------
        latent_representation : np.ndarray
            Low-dimensional representation for each cell
        """
        self._check_if_trained(warn=False)

        adata = self._validate_anndata(adata)
        scdl = self._make_data_loader(
            adata=adata, indices=indices, batch_size=batch_size, iter_ndarray=True
        )

        @jax.jit
        def _get_val(array_dict, rngs):
            rngs = {k: random.split(v)[1] for k, v in rngs.items()}
            out = self.module.apply({"params": self.params}, array_dict, rngs=rngs)
            return out.mean, rngs

        latent = []
        for array_dict in scdl:
            mean, self.rngs = _get_val(array_dict, self.rngs)
            latent.append(mean)
        latent = jnp.concatenate(latent)

        return np.array(jax.device_get(latent))
