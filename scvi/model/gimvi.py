import numpy as np
import os
import logging
import pickle
import torch
from anndata import AnnData

from typing import List, Optional
from scvi.core.modules import JVAE, Classifier
from scvi.core.trainers.jvae_trainer import JvaeDataLoader
from scvi.core.trainers import JVAETrainer
from scvi.core.models import VAEMixin, BaseModelClass
from scvi.model._utils import _get_var_names_from_setup_anndata
from scvi import _CONSTANTS

logger = logging.getLogger(__name__)


def _unpack_tensors(tensors):
    x = tensors[_CONSTANTS.X_KEY].squeeze_(0)
    local_l_mean = tensors[_CONSTANTS.LOCAL_L_MEAN_KEY].squeeze_(0)
    local_l_var = tensors[_CONSTANTS.LOCAL_L_VAR_KEY].squeeze_(0)
    batch_index = tensors[_CONSTANTS.BATCH_KEY].squeeze_(0)
    y = tensors[_CONSTANTS.LABELS_KEY].squeeze_(0)
    return x, local_l_mean, local_l_var, batch_index, y


class GIMVI(VAEMixin, BaseModelClass):
    """
    Joint VAE for imputing missing genes in spatial data [Lopez19]_.

    Parameters
    ----------
    adata_seq
        AnnData object that has been registered via :func:`~scvi.data.setup_anndata`
        and contains RNA-seq data.
    adata_spatial
        AnnData object that has been registered via :func:`~scvi.data.setup_anndata`
        and contains spatial data.
    n_hidden
        Number of nodes per hidden layer.
    generative_distributions
        List of generative distribution for adata_seq data and adata_spatial data.
    model_library_size
        List of bool of whether to model library size for adata_seq and adata_spatial.
    n_latent
        Dimensionality of the latent space.
    use_cuda
        Use the GPU or not.
    **model_kwargs
        Keyword args for :class:`~scvi.core.modules.JVAE`

    Examples
    --------
    >>> adata_seq = anndata.read_h5ad(path_to_anndata_seq)
    >>> adata_spatial = anndata.read_h5ad(path_to_anndata_spatial)
    >>> scvi.data.setup_anndata(adata_seq)
    >>> scvi.data.setup_anndata(adata_spatial)
    >>> vae = scvi.model.GIMVI(adata_seq, adata_spatial)
    >>> vae.train(n_epochs=400)
    """

    def __init__(
        self,
        adata_seq: AnnData,
        adata_spatial: AnnData,
        generative_distributions: List = ["zinb", "nb"],
        model_library_size: List = [True, False],
        n_latent: int = 10,
        use_cuda: bool = True,
        **model_kwargs,
    ):
        super(GIMVI, self).__init__(use_cuda=use_cuda)
        self.use_cuda = use_cuda and torch.cuda.is_available()
        self.adatas = [adata_seq, adata_spatial]
        seq_var_names = _get_var_names_from_setup_anndata(adata_seq)
        spatial_var_names = _get_var_names_from_setup_anndata(adata_spatial)
        spatial_gene_loc = [
            np.argwhere(seq_var_names == g)[0] for g in spatial_var_names
        ]
        spatial_gene_loc = np.concatenate(spatial_gene_loc)
        gene_mappings = [slice(None), spatial_gene_loc]
        sum_stats = [d.uns["_scvi"]["summary_stats"] for d in self.adatas]
        n_inputs = [s["n_genes"] for s in sum_stats]
        total_genes = adata_seq.uns["_scvi"]["summary_stats"]["n_genes"]
        n_batches = sum([s["n_batch"] for s in sum_stats])

        self.model = JVAE(
            n_inputs,
            total_genes,
            gene_mappings,
            generative_distributions,
            model_library_size,
            n_batch=n_batches,
            n_latent=n_latent,
            **model_kwargs,
        )

        self._model_summary_string = "gimVI model with params"
        self.init_params_ = self._get_init_params(locals())

    @property
    def _trainer_class(self):
        return JVAETrainer

    @property
    def _scvi_dl_class(self):
        return JvaeDataLoader

    def train(
        self,
        n_epochs: Optional[int] = 200,
        kappa: Optional[int] = 5,
        discriminator: Optional[Classifier] = None,
        train_size: float = 0.9,
        frequency: int = 1,
        n_epochs_kl_warmup: int = 400,
        train_fun_kwargs: dict = {},
        **kwargs,
    ):
        """
        Train the model.

        Parameters
        ----------
        n_epochs
            Number of passes through the dataset.
        kappa
            Scaling parameter for the discriminator loss.
        discriminator
            :class:`~scvi.core.modules.Classifier` object.
        train_size
            Size of training set in the range [0.0, 1.0].
        frequency
            Frequency with which metrics are computed on the data for train/test/val sets.
        n_epochs_kl_warmup
            Number of passes through dataset for scaling term on KL divergence to go from 0 to 1.
        train_fun_kwargs
            Keyword args for the train method of :class:`~scvi.core.trainers.trainer.Trainer`.
        **kwargs
            Other keyword args for :class:`~scvi.core.trainers.trainer.Trainer`.
        """
        train_fun_kwargs = dict(train_fun_kwargs)
        if discriminator is None:
            discriminator = Classifier(self.model.n_latent, 32, 2, 3, logits=True)
        self.trainer = JVAETrainer(
            self.model,
            discriminator,
            self.adatas,
            train_size,
            frequency=frequency,
            kappa=kappa,
            n_epochs_kl_warmup=n_epochs_kl_warmup,
        )
        self.trainer.train(n_epochs=n_epochs, **train_fun_kwargs)

        self.is_trained_ = True

    def _make_scvi_dls(self, adatas: List[AnnData] = None, batch_size=128):
        if adatas is None:
            adatas = self.adatas
        post_list = [
            self._make_scvi_dl(adata, mode=i) for i, adata in enumerate(adatas)
        ]

        return post_list

    def get_latent_representation(
        self,
        adatas: List[AnnData] = None,
        deterministic: bool = True,
        batch_size: int = 128,
    ) -> List[np.ndarray]:
        """
        Return the latent space embedding for each dataset.

        Parameters
        ----------
        adatas
            List of adata seq and adata spatial.
        deterministic
            If true, use the mean of the encoder instead of a Gaussian sample.
        batch_size
            Minibatch size for data loading into model.
        """
        if adatas is None:
            adatas = self.adatas
        scdls = self._make_scvi_dls(adatas, batch_size=batch_size)
        self.model.eval()
        latents = []
        for mode, scdl in enumerate(scdls):
            latent = []
            for tensors in scdl:
                (
                    sample_batch,
                    local_l_mean,
                    local_l_var,
                    batch_index,
                    label,
                    *_,
                ) = _unpack_tensors(tensors)
                latent.append(
                    self.model.sample_from_posterior_z(
                        sample_batch, mode, deterministic=deterministic
                    )
                )

            latent = torch.cat(latent).cpu().detach().numpy()
            latents.append(latent)

        return latents

    def get_imputed_values(
        self,
        adatas: List[AnnData] = None,
        deterministic: bool = True,
        normalized: bool = True,
        decode_mode: Optional[int] = None,
        batch_size: int = 128,
    ) -> List[np.ndarray]:
        """
        Return imputed values for all genes for each dataset.

        Parameters
        ----------
        adatas
            List of adata seq and adata spatial
        deterministic
            If true, use the mean of the encoder instead of a Gaussian sample for the latent vector.
        normalized
            Return imputed normalized values or not.
        decode_mode
            If a `decode_mode` is given, use the encoder specific to each dataset as usual but use
            the decoder of the dataset of id `decode_mode` to impute values.
        batch_size
            Minibatch size for data loading into model.
        """
        self.model.eval()

        if adatas is None:
            adatas = self.adatas
        scdls = self._make_scvi_dls(adatas, batch_size=batch_size)

        imputed_values = []
        for mode, scdl in enumerate(scdls):
            imputed_value = []
            for tensors in scdl:
                (
                    sample_batch,
                    local_l_mean,
                    local_l_var,
                    batch_index,
                    label,
                    *_,
                ) = _unpack_tensors(tensors)
                if normalized:
                    imputed_value.append(
                        self.model.sample_scale(
                            sample_batch,
                            mode,
                            batch_index,
                            label,
                            deterministic=deterministic,
                            decode_mode=decode_mode,
                        )
                    )
                else:
                    imputed_value.append(
                        self.model.sample_rate(
                            sample_batch,
                            mode,
                            batch_index,
                            label,
                            deterministic=deterministic,
                            decode_mode=decode_mode,
                        )
                    )

            imputed_value = torch.cat(imputed_value).cpu().detach().numpy()
            imputed_values.append(imputed_value)

        return imputed_values

    @classmethod
    def load(
        cls,
        adata_seq: AnnData,
        adata_spatial: AnnData,
        dir_path: str,
        use_cuda: bool = False,
    ):
        """
        Instantiate a model from the saved output.

        Parameters
        ----------
        adata_seq
            AnnData organized in the same way as data used to train model.
            AnnData must be registered via :func:`~scvi.data.setup_anndata`.
        adata_spatial
            AnnData organized in the same way as data used to train model.
            AnnData must be registered via :func:`~scvi.data.setup_anndata`.
        dir_path
            Path to saved outputs.
        use_cuda
            Whether to load model on GPU.

        Returns
        -------
        Model with loaded state dictionaries.

        Examples
        --------
        >>> vae = GIMVI.load(adata_seq, adata_spatial, save_path)
        >>> vae.get_latent_representation()
        """
        model_path = os.path.join(dir_path, "model_params.pt")
        # optimizer_path = os.path.join(dir_path, "optimizer_params.pt")
        setup_dict_path = os.path.join(dir_path, "attr.pkl")
        with open(setup_dict_path, "rb") as handle:
            attr_dict = pickle.load(handle)
        # get the parameters for the class init signiture
        init_params = attr_dict.pop("init_params_")
        # grab all the parameters execept for kwargs (is a dict)
        non_kwargs = {k: v for k, v in init_params.items() if not isinstance(v, dict)}
        # expand out kwargs
        kwargs = {k: v for k, v in init_params.items() if isinstance(v, dict)}
        kwargs = {k: v for (i, j) in kwargs.items() for (k, v) in j.items()}
        model = cls(adata_seq, adata_spatial, **non_kwargs, **kwargs)
        for attr, val in attr_dict.items():
            setattr(model, attr, val)
        use_cuda = use_cuda and torch.cuda.is_available()

        if use_cuda:
            model.model.load_state_dict(torch.load(model_path))
            model.model.cuda()
        else:
            device = torch.device("cpu")
            model.model.load_state_dict(torch.load(model_path, map_location=device))
        model.model.eval()
        return model
