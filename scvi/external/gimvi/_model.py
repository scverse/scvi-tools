import logging
import os
import warnings
from itertools import cycle
from typing import List, Optional, Union

import numpy as np
import torch
from anndata import AnnData
from torch.utils.data import DataLoader

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data._compat import registry_from_setup_dict
from scvi.data._constants import _MODEL_NAME_KEY, _SETUP_ARGS_KEY
from scvi.data.fields import CategoricalObsField, LayerField
from scvi.dataloaders import DataSplitter
from scvi.model._utils import _init_library_size, parse_use_gpu_arg
from scvi.model.base import BaseModelClass, VAEMixin
from scvi.train import Trainer
from scvi.utils import setup_anndata_dsp

from ._module import JVAE
from ._task import GIMVITrainingPlan
from ._utils import _load_legacy_saved_gimvi_files, _load_saved_gimvi_files

logger = logging.getLogger(__name__)


def _unpack_tensors(tensors):
    x = tensors[REGISTRY_KEYS.X_KEY].squeeze_(0)
    batch_index = tensors[REGISTRY_KEYS.BATCH_KEY].squeeze_(0)
    y = tensors[REGISTRY_KEYS.LABELS_KEY].squeeze_(0)
    return x, batch_index, y


class GIMVI(VAEMixin, BaseModelClass):
    """
    Joint VAE for imputing missing genes in spatial data :cite:p:`Lopez19`.

    Parameters
    ----------
    adata_seq
        AnnData object that has been registered via :meth:`~scvi.external.GIMVI.setup_anndata`
        and contains RNA-seq data.
    adata_spatial
        AnnData object that has been registered via :meth:`~scvi.external.GIMVI.setup_anndata`
        and contains spatial data.
    n_hidden
        Number of nodes per hidden layer.
    generative_distributions
        List of generative distribution for adata_seq data and adata_spatial data. Defaults to ['zinb', 'nb'].
    model_library_size
        List of bool of whether to model library size for adata_seq and adata_spatial. Defaults to [True, False].
    n_latent
        Dimensionality of the latent space.
    **model_kwargs
        Keyword args for :class:`~scvi.external.gimvi.JVAE`

    Examples
    --------
    >>> adata_seq = anndata.read_h5ad(path_to_anndata_seq)
    >>> adata_spatial = anndata.read_h5ad(path_to_anndata_spatial)
    >>> scvi.external.GIMVI.setup_anndata(adata_seq)
    >>> scvi.external.GIMVI.setup_anndata(adata_spatial)
    >>> vae = scvi.model.GIMVI(adata_seq, adata_spatial)
    >>> vae.train(n_epochs=400)

    Notes
    -----
    See further usage examples in the following tutorials:

    1. :doc:`/user_guide/notebooks/gimvi_tutorial`
    """

    def __init__(
        self,
        adata_seq: AnnData,
        adata_spatial: AnnData,
        generative_distributions: Optional[List[str]] = None,
        model_library_size: Optional[List[bool]] = None,
        n_latent: int = 10,
        **model_kwargs,
    ):
        super().__init__()
        if adata_seq is adata_spatial:
            raise ValueError(
                "`adata_seq` and `adata_spatial` cannot point to the same object. "
                "If you would really like to do this, make a copy of the object and pass it in as `adata_spatial`."
            )
        model_library_size = model_library_size or [True, False]
        generative_distributions = generative_distributions or ["zinb", "nb"]
        self.adatas = [adata_seq, adata_spatial]
        self.adata_managers = {
            "seq": self._get_most_recent_anndata_manager(adata_seq, required=True),
            "spatial": self._get_most_recent_anndata_manager(
                adata_spatial, required=True
            ),
        }
        self.registries_ = []
        for adm in self.adata_managers.values():
            self._register_manager_for_instance(adm)
            self.registries_.append(adm.registry)

        seq_var_names = adata_seq.var_names
        spatial_var_names = adata_spatial.var_names

        if not set(spatial_var_names) <= set(seq_var_names):
            raise ValueError("spatial genes needs to be subset of seq genes")

        spatial_gene_loc = [
            np.argwhere(seq_var_names == g)[0] for g in spatial_var_names
        ]
        spatial_gene_loc = np.concatenate(spatial_gene_loc)
        gene_mappings = [slice(None), spatial_gene_loc]
        sum_stats = [adm.summary_stats for adm in self.adata_managers.values()]
        n_inputs = [s["n_vars"] for s in sum_stats]

        total_genes = n_inputs[0]

        # since we are combining datasets, we need to increment the batch_idx
        # of one of the datasets
        adata_seq_n_batches = sum_stats[0]["n_batch"]
        adata_spatial.obs[
            self.adata_managers["spatial"]
            .data_registry[REGISTRY_KEYS.BATCH_KEY]
            .attr_key
        ] += adata_seq_n_batches

        n_batches = sum(s["n_batch"] for s in sum_stats)

        library_log_means = []
        library_log_vars = []
        for adata_manager in self.adata_managers.values():
            adata_library_log_means, adata_library_log_vars = _init_library_size(
                adata_manager, n_batches
            )
            library_log_means.append(adata_library_log_means)
            library_log_vars.append(adata_library_log_vars)

        self.module = JVAE(
            n_inputs,
            total_genes,
            gene_mappings,
            generative_distributions,
            model_library_size,
            library_log_means,
            library_log_vars,
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
        use_gpu: Optional[Union[str, int, bool]] = None,
        kappa: int = 5,
        train_size: float = 0.9,
        validation_size: Optional[float] = None,
        batch_size: int = 128,
        plan_kwargs: Optional[dict] = None,
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
            Use default GPU if available (if None or True), or index of GPU to use (if int),
            or name of GPU (if str, e.g., `'cuda:0'`), or use CPU (if False).
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
            Keyword args for model-specific Pytorch Lightning task. Keyword arguments passed
            to `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        **kwargs
            Other keyword args for :class:`~scvi.train.Trainer`.
        """
        accelerator, lightning_devices, device = parse_use_gpu_arg(use_gpu)

        self.trainer = Trainer(
            max_epochs=max_epochs,
            accelerator=accelerator,
            devices=lightning_devices,
            **kwargs,
        )
        self.train_indices_, self.test_indices_, self.validation_indices_ = [], [], []
        train_dls, test_dls, val_dls = [], [], []
        for i, adm in enumerate(self.adata_managers.values()):
            ds = DataSplitter(
                adm,
                train_size=train_size,
                validation_size=validation_size,
                batch_size=batch_size,
                use_gpu=use_gpu,
            )
            ds.setup()
            train_dls.append(ds.train_dataloader())
            test_dls.append(ds.test_dataloader())
            val = ds.val_dataloader()
            val_dls.append(val)
            val.mode = i
            self.train_indices_.append(ds.train_idx)
            self.test_indices_.append(ds.test_idx)
            self.validation_indices_.append(ds.val_idx)
        train_dl = TrainDL(train_dls)

        plan_kwargs = plan_kwargs if isinstance(plan_kwargs, dict) else dict()
        self._training_plan = GIMVITrainingPlan(
            self.module,
            adversarial_classifier=True,
            scale_adversarial_loss=kappa,
            **plan_kwargs,
        )

        if train_size == 1.0:
            # circumvent the empty data loader problem if all dataset used for training
            self.trainer.fit(self._training_plan, train_dl)
        else:
            # accepts list of val dataloaders
            self.trainer.fit(self._training_plan, train_dl, val_dls)
        try:
            self.history_ = self.trainer.logger.history
        except AttributeError:
            self.history_ = None
        self.module.eval()

        self.to_device(device)
        self.is_trained_ = True

    def _make_scvi_dls(self, adatas: List[AnnData] = None, batch_size=128):
        if adatas is None:
            adatas = self.adatas
        post_list = [self._make_data_loader(ad) for ad in adatas]
        for i, dl in enumerate(post_list):
            dl.mode = i

        return post_list

    @torch.inference_mode()
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

    @torch.inference_mode()
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
        prefix: Optional[str] = None,
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
        prefix
            Prefix to prepend to saved file names.
        overwrite
            Overwrite existing data or not. If `False` and directory
            already exists at `dir_path`, error will be raised.
        save_anndata
            If True, also saves the anndata
        anndata_write_kwargs
            Kwargs for anndata write function
        """
        if not os.path.exists(dir_path) or overwrite:
            os.makedirs(dir_path, exist_ok=overwrite)
        else:
            raise ValueError(
                "{} already exists. Please provide an unexisting directory for saving.".format(
                    dir_path
                )
            )

        file_name_prefix = prefix or ""

        seq_adata = self.adatas[0]
        spatial_adata = self.adatas[1]
        if save_anndata:
            seq_save_path = os.path.join(dir_path, f"{file_name_prefix}adata_seq.h5ad")
            seq_adata.write(seq_save_path)

            spatial_save_path = os.path.join(
                dir_path, f"{file_name_prefix}adata_spatial.h5ad"
            )
            spatial_adata.write(spatial_save_path)

        # save the model state dict and the trainer state dict only
        model_state_dict = self.module.state_dict()

        seq_var_names = seq_adata.var_names.astype(str).to_numpy()
        spatial_var_names = spatial_adata.var_names.astype(str).to_numpy()

        # get all the user attributes
        user_attributes = self._get_user_attributes()
        # only save the public attributes with _ at the very end
        user_attributes = {a[0]: a[1] for a in user_attributes if a[0][-1] == "_"}

        model_save_path = os.path.join(dir_path, f"{file_name_prefix}model.pt")

        torch.save(
            dict(
                model_state_dict=model_state_dict,
                seq_var_names=seq_var_names,
                spatial_var_names=spatial_var_names,
                attr_dict=user_attributes,
            ),
            model_save_path,
        )

    @classmethod
    def load(
        cls,
        dir_path: str,
        adata_seq: Optional[AnnData] = None,
        adata_spatial: Optional[AnnData] = None,
        use_gpu: Optional[Union[str, int, bool]] = None,
        prefix: Optional[str] = None,
        backup_url: Optional[str] = None,
    ):
        """
        Instantiate a model from the saved output.

        Parameters
        ----------
        dir_path
            Path to saved outputs.
        adata_seq
            AnnData organized in the same way as data used to train model.
            It is not necessary to run :meth:`~scvi.external.GIMVI.setup_anndata`,
            as AnnData is validated against the saved `scvi` setup dictionary.
            AnnData must be registered via :meth:`~scvi.external.GIMVI.setup_anndata`.
        adata_spatial
            AnnData organized in the same way as data used to train model.
            If None, will check for and load anndata saved with the model.
        use_gpu
            Load model on default GPU if available (if None or True),
            or index of GPU to use (if int), or name of GPU (if str), or use CPU (if False).
        prefix
            Prefix of saved file names.
        backup_url
            URL to retrieve saved outputs from if not present on disk.

        Returns
        -------
        Model with loaded state dictionaries.

        Examples
        --------
        >>> vae = GIMVI.load(adata_seq, adata_spatial, save_path)
        >>> vae.get_latent_representation()
        """
        _, _, device = parse_use_gpu_arg(use_gpu)

        (
            attr_dict,
            seq_var_names,
            spatial_var_names,
            model_state_dict,
            loaded_adata_seq,
            loaded_adata_spatial,
        ) = _load_saved_gimvi_files(
            dir_path,
            adata_seq is None,
            adata_spatial is None,
            prefix=prefix,
            map_location=device,
            backup_url=backup_url,
        )
        adata_seq = loaded_adata_seq or adata_seq
        adata_spatial = loaded_adata_spatial or adata_spatial
        adatas = [adata_seq, adata_spatial]
        var_names = [seq_var_names, spatial_var_names]

        for i, adata in enumerate(adatas):
            saved_var_names = var_names[i]
            user_var_names = adata.var_names.astype(str)
            if not np.array_equal(saved_var_names, user_var_names):
                warnings.warn(
                    "var_names for adata passed in does not match var_names of "
                    "adata used to train the model. For valid results, the vars "
                    "need to be the same and in the same order as the adata used to train the model."
                )

        registries = attr_dict.pop("registries_")
        for adata, registry in zip(adatas, registries):
            if (
                _MODEL_NAME_KEY in registry
                and registry[_MODEL_NAME_KEY] != cls.__name__
            ):
                raise ValueError(
                    "It appears you are loading a model from a different class."
                )

            if _SETUP_ARGS_KEY not in registry:
                raise ValueError(
                    "Saved model does not contain original setup inputs. "
                    "Cannot load the original setup."
                )

            cls.setup_anndata(
                adata, source_registry=registry, **registry[_SETUP_ARGS_KEY]
            )

        # get the parameters for the class init signature
        init_params = attr_dict.pop("init_params_")

        # new saving and loading, enable backwards compatibility
        if "non_kwargs" in init_params.keys():
            # grab all the parameters except for kwargs (is a dict)
            non_kwargs = init_params["non_kwargs"]
            kwargs = init_params["kwargs"]

            # expand out kwargs
            kwargs = {k: v for (i, j) in kwargs.items() for (k, v) in j.items()}
        else:
            # grab all the parameters except for kwargs (is a dict)
            non_kwargs = {
                k: v for k, v in init_params.items() if not isinstance(v, dict)
            }
            kwargs = {k: v for k, v in init_params.items() if isinstance(v, dict)}
            kwargs = {k: v for (i, j) in kwargs.items() for (k, v) in j.items()}
        model = cls(adata_seq, adata_spatial, **non_kwargs, **kwargs)

        for attr, val in attr_dict.items():
            setattr(model, attr, val)

        model.module.load_state_dict(model_state_dict)
        model.module.eval()
        model.to_device(device)
        return model

    @classmethod
    def convert_legacy_save(
        cls,
        dir_path: str,
        output_dir_path: str,
        overwrite: bool = False,
        prefix: Optional[str] = None,
    ) -> None:
        """
        Converts a legacy saved GIMVI model (<v0.15.0) to the updated save format.

        Parameters
        ----------
        dir_path
            Path to directory where legacy model is saved.
        output_dir_path
            Path to save converted save files.
        overwrite
            Overwrite existing data or not. If ``False`` and directory
            already exists at ``output_dir_path``, error will be raised.
        prefix
            Prefix of saved file names.
        """
        if not os.path.exists(output_dir_path) or overwrite:
            os.makedirs(output_dir_path, exist_ok=overwrite)
        else:
            raise ValueError(
                "{} already exists. Please provide an unexisting directory for saving.".format(
                    dir_path
                )
            )

        file_name_prefix = prefix or ""
        (
            model_state_dict,
            seq_var_names,
            spatial_var_names,
            attr_dict,
            _,
            _2,
        ) = _load_legacy_saved_gimvi_files(
            dir_path,
            file_name_prefix,
            load_seq_adata=False,
            load_spatial_adata=False,
        )
        if "scvi_setup_dicts_" in attr_dict:
            scvi_setup_dicts = attr_dict.pop("scvi_setup_dicts_")
            registries = []
            for scvi_setup_dict in scvi_setup_dicts:
                registries.append(registry_from_setup_dict(cls, scvi_setup_dict))
            attr_dict["registries_"] = registries

        model_save_path = os.path.join(output_dir_path, f"{file_name_prefix}model.pt")
        torch.save(
            dict(
                model_state_dict=model_state_dict,
                seq_var_names=seq_var_names,
                spatial_var_names=spatial_var_names,
                attr_dict=attr_dict,
            ),
            model_save_path,
        )

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        batch_key: Optional[str] = None,
        labels_key: Optional[str] = None,
        layer: Optional[str] = None,
        **kwargs,
    ):
        """
        %(summary)s.

        Parameters
        ----------
        %(param_batch_key)s
        %(param_labels_key)s
        %(param_layer)s
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            CategoricalObsField(REGISTRY_KEYS.LABELS_KEY, labels_key),
        ]
        adata_manager = AnnDataManager(
            fields=anndata_fields, setup_method_args=setup_method_args
        )
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)


class TrainDL(DataLoader):
    """Train data loader."""

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
