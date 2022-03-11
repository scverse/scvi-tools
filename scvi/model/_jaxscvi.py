import logging
from typing import Any, Optional, Sequence, Union

import jax
import jax.numpy as jnp
import numpy as np
import optax
import pandas as pd
import tqdm
from anndata import AnnData
from flax.core import FrozenDict
from flax.training import train_state
from jax import random

from scvi import REGISTRY_KEYS
from scvi._compat import Literal
from scvi.data import AnnDataManager
from scvi.data.fields import CategoricalObsField, LayerField
from scvi.dataloaders import DataSplitter
from scvi.module import JaxVAE
from scvi.utils import setup_anndata_dsp

from .base import BaseModelClass

logger = logging.getLogger(__name__)


class TrainState(train_state.TrainState):
    batch_stats: FrozenDict[str, Any]


class JaxSCVI(BaseModelClass):
    """
    EXPERIMENTAL single-cell Variational Inference [Lopez18]_, but with a Jax backend.

    This implementation is in a very experimental state. API is completely subject to change.

    Parameters
    ----------
    adata
        AnnData object that has been registered via :meth:`~scvi.model.JaxSCVI.setup_anndata`.
    n_hidden
        Number of nodes per hidden layer.
    n_latent
        Dimensionality of the latent space.
    dropout_rate
        Dropout rate for neural networks.
    gene_likelihood
        One of:

        * ``'nb'`` - Negative binomial distribution
        * ``'poisson'`` - Poisson distribution
    **model_kwargs
        Keyword args for :class:`~scvi.module.JaxVAE`

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.model.JaxSCVI.setup_anndata(adata, batch_key="batch")
    >>> vae = scvi.model.SCVI(adata)
    >>> vae.train()
    >>> adata.obsm["X_scVI"] = vae.get_latent_representation()
    """

    def __init__(
        self,
        adata: AnnData,
        n_hidden: int = 128,
        n_latent: int = 10,
        dropout_rate: float = 0.1,
        gene_likelihood: Literal["nb", "poisson"] = "nb",
        **model_kwargs,
    ):
        super().__init__(adata)

        n_batch = self.summary_stats.n_batch

        self.module_kwargs = dict(
            n_input=self.summary_stats.n_vars,
            n_batch=n_batch,
            n_hidden=n_hidden,
            n_latent=n_latent,
            dropout_rate=dropout_rate,
            is_training=False,
            gene_likelihood=gene_likelihood,
        )
        self.module_kwargs.update(model_kwargs)
        self._module = None

        self._model_summary_string = ""
        self.init_params_ = self._get_init_params(locals())

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        layer: Optional[str] = None,
        batch_key: Optional[str] = None,
        **kwargs,
    ):
        """
        %(summary)s.

        Parameters
        ----------
        %(param_layer)s
        %(param_batch_key)s
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
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
        check_val_every_n_epoch: Optional[int] = None,
        use_gpu: Optional[Union[str, int, bool]] = None,
        train_size: float = 0.9,
        validation_size: Optional[float] = None,
        batch_size: int = 128,
        lr: float = 1e-3,
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
        """
        if max_epochs is None:
            n_cells = self.adata.n_obs
            max_epochs = np.min([round((20000 / n_cells) * 400), 400])

        data_splitter = DataSplitter(
            self.adata_manager,
            train_size=train_size,
            validation_size=validation_size,
            batch_size=batch_size,
            # for pinning memory only
            use_gpu=False,
            iter_ndarray=True,
        )
        data_splitter.setup()
        train_loader = data_splitter.train_dataloader()
        val_dataloader = data_splitter.val_dataloader()

        if val_dataloader is None:
            raise UserWarning(
                "No observations in the validation loop. No validation metrics will be recorded."
            )

        module_kwargs = self.module_kwargs.copy()
        module_kwargs.update(dict(is_training=True))
        train_module = self._get_module(module_kwargs)

        # if key is generated on CPU, model params will be on CPU
        # we have to pay the price of a JIT compilation though
        if use_gpu is False:
            key = jax.jit(lambda i: random.PRNGKey(i), backend="cpu")
        else:
            # dummy function
            def key(i: int):
                return random.PRNGKey(i)

        self.rngs = {
            "params": key(0),
            "dropout": key(1),
            "z": key(2),
        }
        module_init = train_module.init(self.rngs, next(iter(train_loader)))
        params = module_init["params"]
        batch_stats = module_init["batch_stats"]

        state = TrainState.create(
            apply_fn=train_module.apply,
            params=params,
            tx=optax.adamw(lr, eps=0.01, weight_decay=1e-6),
            batch_stats=batch_stats,
        )

        @jax.jit
        def train_step(state, array_dict, rngs, **kwargs):
            rngs = {k: random.split(v)[1] for k, v in rngs.items()}

            # batch stats can't be passed here
            def loss_fn(params):
                vars_in = {"params": params, "batch_stats": state.batch_stats}
                outputs, new_model_state = state.apply_fn(
                    vars_in, array_dict, rngs=rngs, mutable=["batch_stats"], **kwargs
                )
                loss_recorder = outputs[2]
                loss = loss_recorder.loss
                elbo = jnp.mean(
                    loss_recorder.reconstruction_loss + loss_recorder.kl_local
                )
                return loss, (elbo, new_model_state)

            (loss, (elbo, new_model_state)), grads = jax.value_and_grad(
                loss_fn, has_aux=True
            )(state.params)
            new_state = state.apply_gradients(
                grads=grads, batch_stats=new_model_state["batch_stats"]
            )
            return new_state, loss, elbo, rngs

        @jax.jit
        def validation_step(state, array_dict, rngs, **kwargs):
            # note that self.module has is_training = False
            module = self.module
            rngs = {k: random.split(v)[1] for k, v in rngs.items()}
            vars_in = {"params": state.params, "batch_stats": state.batch_stats}
            outputs = module.apply(vars_in, array_dict, rngs=rngs, **kwargs)
            loss_recorder = outputs[2]
            loss = loss_recorder.loss
            elbo = jnp.mean(loss_recorder.reconstruction_loss + loss_recorder.kl_local)

            return loss, elbo

        history = dict(
            elbo_train=[], loss_train=[], elbo_validation=[], loss_validation=[]
        )
        epoch = 0
        with tqdm.trange(1, max_epochs + 1) as t:
            try:
                for i in t:
                    epoch += 1
                    epoch_loss = 0
                    epoch_elbo = 0
                    counter = 0
                    for data in train_loader:
                        kl_weight = min(1.0, epoch / 400.0)
                        # gets new key for each epoch
                        state, loss, elbo, self.rngs = train_step(
                            state,
                            data,
                            self.rngs,
                            loss_kwargs=dict(kl_weight=kl_weight),
                        )
                        epoch_loss += loss
                        epoch_elbo += elbo
                        counter += 1
                    history["loss_train"] += [jax.device_get(epoch_loss) / counter]
                    history["elbo_train"] += [jax.device_get(epoch_elbo) / counter]
                    t.set_postfix_str(
                        f"Epoch {i}, Elbo: {epoch_elbo / counter}, KL weight: {kl_weight}"
                    )

                    # validation loop
                    if (
                        check_val_every_n_epoch is not None
                        and epoch % check_val_every_n_epoch == 0
                    ):
                        val_counter = 0
                        val_epoch_loss = 0
                        val_epoch_elbo = 0
                        for data in val_dataloader:
                            val_loss, val_elbo = validation_step(
                                state,
                                data,
                                self.rngs,
                                loss_kwargs=dict(kl_weight=kl_weight),
                            )
                            val_epoch_loss += val_loss
                            val_epoch_elbo += val_elbo
                            val_counter += 1
                        history["loss_validation"] += [
                            jax.device_get(val_epoch_loss) / val_counter
                        ]
                        history["elbo_validation"] += [
                            jax.device_get(val_epoch_elbo) / val_counter
                        ]

            except KeyboardInterrupt:
                logger.info(
                    "Keyboard interrupt detected. Attempting graceful shutdown."
                )

        self.train_state = state
        self.params = state.params
        self.batch_stats = state.batch_stats
        self.history_ = {k: pd.DataFrame(v, columns=[k]) for k, v in history.items()}
        self.is_trained_ = True

        self.module_kwargs.update(dict(is_training=False))
        self._module = None
        self.bound_module = self.module.bind(
            {"params": self.params, "batch_stats": self.batch_stats}, rngs=self.rngs
        )

    def get_latent_representation(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        give_mean: bool = True,
        mc_samples: int = 1,
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
        def _get_val(array_dict):
            inference_input = self.bound_module._get_inference_input(array_dict)
            out = self.bound_module.inference(**inference_input, n_samples=mc_samples)
            return out

        latent = []
        for array_dict in scdl:
            out = _get_val(array_dict)
            if give_mean:
                z = out["qz"].mean
            else:
                z = out["z"]
            latent.append(z)
        concat_axis = 0 if ((mc_samples == 1) or give_mean) else 1
        latent = jnp.concatenate(latent, axis=concat_axis)

        return np.array(jax.device_get(latent))

    def save(self):
        raise NotImplementedError

    def load(self):
        raise NotImplementedError

    def to_device(self):
        raise NotImplementedError

    @property
    def device(self):
        raise NotImplementedError
