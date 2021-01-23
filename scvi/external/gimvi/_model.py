import logging
import os
import pickle
from itertools import cycle
from typing import List, Optional

import numpy as np
import torch
from anndata import AnnData, read
from torch.utils.data import DataLoader

from scvi import _CONSTANTS, settings
from scvi.data import transfer_anndata_setup
from scvi.dataloaders import AnnDataLoader
from scvi.lightning import Trainer
from scvi.model._utils import _get_var_names_from_setup_anndata
from scvi.model.base import BaseModelClass, VAEMixin

from ._module import JVAE
from ._task import GIMVITrainingPlan

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
    use_gpu
        Use the GPU or not.
    **model_kwargs
        Keyword args for :class:`~scvi.modules.JVAE`

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
        use_gpu: bool = True,
        **model_kwargs,
    ):
        super(GIMVI, self).__init__(use_gpu=use_gpu)
        self.use_gpu = use_gpu and torch.cuda.is_available()
        self.adatas = [adata_seq, adata_spatial]
        self.scvi_setup_dicts_ = {
            "seq": adata_seq.uns["_scvi"],
            "spatial": adata_spatial.uns["_scvi"],
        }

        seq_var_names = _get_var_names_from_setup_anndata(adata_seq)
        spatial_var_names = _get_var_names_from_setup_anndata(adata_spatial)

        if not set(spatial_var_names) <= set(seq_var_names):
            raise ValueError("spatial genes needs to be subset of seq genes")

        spatial_gene_loc = [
            np.argwhere(seq_var_names == g)[0] for g in spatial_var_names
        ]
        spatial_gene_loc = np.concatenate(spatial_gene_loc)
        gene_mappings = [slice(None), spatial_gene_loc]
        sum_stats = [d.uns["_scvi"]["summary_stats"] for d in self.adatas]
        n_inputs = [s["n_vars"] for s in sum_stats]

        total_genes = adata_seq.uns["_scvi"]["summary_stats"]["n_vars"]

        # since we are combining datasets, we need to increment the batch_idx
        # of one of the datasets
        adata_seq_n_batches = adata_seq.uns["_scvi"]["summary_stats"]["n_batch"]
        adata_spatial.obs["_scvi_batch"] += adata_seq_n_batches

        n_batches = sum([s["n_batch"] for s in sum_stats])

        self.module = JVAE(
            n_inputs,
            total_genes,
            gene_mappings,
            generative_distributions,
            model_library_size,
            n_batch=n_batches,
            n_latent=n_latent,
            **model_kwargs,
        )

        self._model_summary_string = (
            "GimVI Model with the following params: \nn_latent: {}, n_inputs: {}, n_genes: {}, "
            + "n_batch: {}, generative distributions: {}"
        ).format(n_latent, n_inputs, total_genes, n_batches, generative_distributions)
        self.init_params_ = self._get_init_params(locals())

    def train(
        self,
        max_epochs: int = 200,
        use_gpu: Optional[bool] = None,
        kappa: int = 5,
        train_size: float = 0.9,
        validation_size: Optional[float] = None,
        batch_size: int = 128,
        plan_kwargs: Optional[dict] = None,
        plan_class: Optional[None] = None,
        **kwargs,
    ):
        """
        Train the model.

        Parameters
        ----------
        max_epochs
            Number of passes through the dataset. If `None`, defaults to
            `np.min([round((20000 / n_cells) * 400), 400])`
        use_gpu
            If `True`, use the GPU if available. Will override the use_gpu option when initializing model
        kappa
            Scaling parameter for the discriminator loss.
        train_size
            Size of training set in the range [0.0, 1.0].
        validation_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + validation_size < 1`, the remaining cells belong to a test set.
        batch_size
            Minibatch size to use during training.
        plan_kwargs
            Keyword args for model-specific Pytorch Lightning task. Keyword arguments passed to
            `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        **kwargs
            Other keyword args for :class:`~scvi.lightning.Trainer`.
        """
        if use_gpu is None:
            use_gpu = self.use_gpu
        else:
            use_gpu = use_gpu and torch.cuda.is_available()
        gpus = 1 if use_gpu else None
        pin_memory = (
            True if (settings.dl_pin_memory_gpu_training and use_gpu) else False
        )

        self.trainer = Trainer(
            max_epochs=max_epochs,
            gpus=gpus,
            **kwargs,
        )
        self.train_indices_, self.test_indices_, self.validation_indices_ = [], [], []
        train_dls, test_dls, val_dls = [], [], []
        for i, ad in enumerate(self.adatas):
            train, val, test = self._train_test_val_split(
                ad,
                train_size=train_size,
                validation_size=validation_size,
                pin_memory=pin_memory,
                batch_size=batch_size,
            )
            train_dls.append(train)
            test_dls.append(test)
            val.mode = i
            val_dls.append(val)
            self.train_indices_.append(train.indices)
            self.test_indices_.append(test.indices)
            self.validation_indices_.append(val.indices)
        train_dl = TrainDL(train_dls)

        plan_kwargs = plan_kwargs if isinstance(plan_kwargs, dict) else dict()
        self._pl_task = self._plan_class(
            self.module,
            len(self.train_indices_),
            adversarial_classifier=True,
            scale_adversarial_loss=kappa,
            **plan_kwargs,
        )

        if train_size == 1.0:
            # circumvent the empty data loader problem if all dataset used for training
            self.trainer.fit(self._pl_task, train_dl)
        else:
            # accepts list of val dataloaders
            self.trainer.fit(self._pl_task, train_dl, val_dls)
        try:
            self.history_ = self.trainer.logger.history
        except AttributeError:
            self.history_ = None
        self.module.eval()
        if use_gpu:
            self.module.cuda()
        self.is_trained_ = True

    def _make_scvi_dls(self, adatas: List[AnnData] = None, batch_size=128):
        if adatas is None:
            adatas = self.adatas
        post_list = [self._make_scvi_dl(ad) for ad in adatas]
        for i, dl in enumerate(post_list):
            dl.mode = i

        return post_list

    @torch.no_grad()
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
        self.module.eval()
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
                    self.module.sample_from_posterior_z(
                        sample_batch, mode, deterministic=deterministic
                    )
                )

            latent = torch.cat(latent).cpu().detach().numpy()
            latents.append(latent)

        return latents

    @torch.no_grad()
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
        self.module.eval()

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
                        self.module.sample_scale(
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
                        self.module.sample_rate(
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

    def save(
        self,
        dir_path: str,
        overwrite: bool = False,
        save_anndata: bool = False,
        **anndata_write_kwargs,
    ):
        """
        Save the state of the model.

        Neither the trainer optimizer state nor the trainer history are saved.
        Model files are not expected to be reproducibly saved and loaded across versions
        until we reach version 1.0.

        Parameters
        ----------
        dir_path
            Path to a directory.
        overwrite
            Overwrite existing data or not. If `False` and directory
            already exists at `dir_path`, error will be raised.
        save_anndata
            If True, also saves the anndata
        anndata_write_kwargs
            Kwargs for anndata write function
        """
        # get all the user attributes
        user_attributes = self._get_user_attributes()
        # only save the public attributes with _ at the very end
        user_attributes = {a[0]: a[1] for a in user_attributes if a[0][-1] == "_"}
        # save the model state dict and the trainer state dict only
        if not os.path.exists(dir_path) or overwrite:
            os.makedirs(dir_path, exist_ok=overwrite)
        else:
            raise ValueError(
                "{} already exists. Please provide an unexisting directory for saving.".format(
                    dir_path
                )
            )
        if save_anndata:
            dataset_names = ["seq", "spatial"]
            for i in range(len(self.adatas)):
                save_path = os.path.join(
                    dir_path, "adata_{}.h5ad".format(dataset_names[i])
                )
                self.adatas[i].write(save_path)
                varnames_save_path = os.path.join(
                    dir_path, "var_names_{}.csv".format(dataset_names[i])
                )

                var_names = self.adatas[i].var_names.astype(str)
                var_names = var_names.to_numpy()
                np.savetxt(varnames_save_path, var_names, fmt="%s")

        model_save_path = os.path.join(dir_path, "model_params.pt")
        attr_save_path = os.path.join(dir_path, "attr.pkl")

        torch.save(self.module.state_dict(), model_save_path)
        with open(attr_save_path, "wb") as f:
            pickle.dump(user_attributes, f)

    @classmethod
    def load(
        cls,
        dir_path: str,
        adata_seq: Optional[AnnData] = None,
        adata_spatial: Optional[AnnData] = None,
        use_gpu: Optional[bool] = None,
    ):
        """
        Instantiate a model from the saved output.

        Parameters
        ----------
        adata_seq
            AnnData organized in the same way as data used to train model.
            It is not necessary to run :func:`~scvi.data.setup_anndata`,
            as AnnData is validated against the saved `scvi` setup dictionary.
            AnnData must be registered via :func:`~scvi.data.setup_anndata`.
        adata_spatial
            AnnData organized in the same way as data used to train model.
            If None, will check for and load anndata saved with the model.
        dir_path
            Path to saved outputs.
        use_gpu
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
        setup_dict_path = os.path.join(dir_path, "attr.pkl")
        seq_data_path = os.path.join(dir_path, "adata_seq.h5ad")
        spatial_data_path = os.path.join(dir_path, "adata_spatial.h5ad")
        seq_var_names_path = os.path.join(dir_path, "var_names_seq.csv")
        spatial_var_names_path = os.path.join(dir_path, "var_names_spatial.csv")

        if adata_seq is None and os.path.exists(seq_data_path):
            adata_seq = read(seq_data_path)
        elif adata_seq is None and not os.path.exists(seq_data_path):
            raise ValueError(
                "Save path contains no saved anndata and no adata was passed."
            )
        if adata_spatial is None and os.path.exists(spatial_data_path):
            adata_spatial = read(spatial_data_path)
        elif adata_spatial is None and not os.path.exists(spatial_data_path):
            raise ValueError(
                "Save path contains no saved anndata and no adata was passed."
            )
        adatas = [adata_seq, adata_spatial]

        seq_var_names = np.genfromtxt(seq_var_names_path, delimiter=",", dtype=str)
        spatial_var_names = np.genfromtxt(
            spatial_var_names_path, delimiter=",", dtype=str
        )
        var_names = [seq_var_names, spatial_var_names]

        for i, adata in enumerate(adatas):
            saved_var_names = var_names[i]
            user_var_names = adata.var_names.astype(str)
            if not np.array_equal(saved_var_names, user_var_names):
                logger.warning(
                    "var_names for adata passed in does not match var_names of "
                    "adata used to train the model. For valid results, the vars "
                    "need to be the same and in the same order as the adata used to train the model."
                )

        with open(setup_dict_path, "rb") as handle:
            attr_dict = pickle.load(handle)

        scvi_setup_dicts = attr_dict.pop("scvi_setup_dicts_")
        transfer_anndata_setup(scvi_setup_dicts["seq"], adata_seq)
        transfer_anndata_setup(scvi_setup_dicts["spatial"], adata_spatial)

        # get the parameters for the class init signiture
        init_params = attr_dict.pop("init_params_")

        # update use_gpu from the saved model
        if use_gpu is None:
            use_gpu = torch.cuda.is_available()

        init_params["use_gpu"] = use_gpu

        # new saving and loading, enable backwards compatibility
        if "non_kwargs" in init_params.keys():
            # grab all the parameters execept for kwargs (is a dict)
            non_kwargs = init_params["non_kwargs"]
            kwargs = init_params["kwargs"]

            # update use_gpu from the saved model
            # we assume use_gpu is exposed and not a kwarg
            non_kwargs["use_gpu"] = use_gpu

            # expand out kwargs
            kwargs = {k: v for (i, j) in kwargs.items() for (k, v) in j.items()}
        else:
            init_params["use_gpu"] = use_gpu

            # grab all the parameters execept for kwargs (is a dict)
            non_kwargs = {
                k: v for k, v in init_params.items() if not isinstance(v, dict)
            }
            kwargs = {k: v for k, v in init_params.items() if isinstance(v, dict)}
            kwargs = {k: v for (i, j) in kwargs.items() for (k, v) in j.items()}
        model = cls(adata_seq, adata_spatial, **non_kwargs, **kwargs)

        for attr, val in attr_dict.items():
            setattr(model, attr, val)

        if use_gpu:
            model.module.load_state_dict(torch.load(model_path))
            model.module.cuda()
        else:
            device = torch.device("cpu")
            model.module.load_state_dict(torch.load(model_path, map_location=device))
        model.module.eval()
        return model

    @property
    def _data_loader_cls(self):
        return AnnDataLoader

    @property
    def _plan_class(self):
        return GIMVITrainingPlan


class TrainDL(DataLoader):
    def __init__(self, data_loader_list, **kwargs):
        self.data_loader_list = data_loader_list
        self.largest_train_dl_idx = np.argmax(
            [len(dl.indices) for dl in data_loader_list]
        )
        self.largest_dl = self.data_loader_list[self.largest_train_dl_idx]
        super().__init__(self.largest_dl, **kwargs)

    def __len__(self):
        return len(self.largest_dl)

    def __iter__(self):
        train_dls = [
            dl if i == self.largest_train_dl_idx else cycle(dl)
            for i, dl in enumerate(self.data_loader_list)
        ]
        return zip(*train_dls)
