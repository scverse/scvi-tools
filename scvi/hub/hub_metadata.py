import json
import os
from dataclasses import dataclass, field
from typing import List, Optional, Union

import torch
from huggingface_hub import ModelCard, ModelCardData

from scvi.data import AnnDataManager
from scvi.data._utils import _is_minified
from scvi.model.base._utils import _load_saved_files

from ._constants import _SCVI_HUB
from ._url import validate_url
from .model_card_template import template


@dataclass
class HubMetadata:
    """Encapsulates the required metadata for `scvi-tools` hub models.

    Parameters
    ----------
    scvi_version
        The version of `scvi-tools` that the model was trained with.
    anndata_version
        The version of anndata used during model training.
    model_cls_name
        The name of the model class.
    training_data_url
        Link to the training data used to train the model, if it is too large to be uploaded to the hub. This can be
        a cellxgene explorer session url. However it cannot be a self-hosted session -- it must be from the cellxgene
        portal (https://cellxgene.cziscience.com/).
    model_parent_module
        The parent module of the model class. Defaults to `scvi.model`. Change this if you are using a model
        class that is not in the `scvi.model` module, for example, if you are using a model class from a custom module.
    """

    scvi_version: str
    anndata_version: str
    model_cls_name: str
    training_data_url: Optional[str] = None
    model_parent_module: str = _SCVI_HUB.DEFAULT_PARENT_MODULE

    @classmethod
    def from_dir(
        cls,
        local_dir: str,
        anndata_version: str,
        map_location: Optional[Union[torch.device, str, dict]] = "cpu",
        **kwargs,
    ):
        """Create a `HubMetadata` object from a local directory.

        Parameters
        ----------
        local_dir
            The local directory containing the model files.
        anndata_version
            The version of anndata used during model training.
        map_location
            The device to map model tensors to, passed into :meth:`~torch.load`.
        kwargs
            Additional keyword arguments to pass to the HubMetadata initializer.
        """
        attr_dict, _, _, _ = _load_saved_files(
            local_dir, load_adata=False, map_location=map_location
        )
        scvi_version = attr_dict["registry_"]["scvi_version"]
        model_cls_name = attr_dict["registry_"]["model_name"]

        return cls(
            scvi_version,
            anndata_version,
            model_cls_name,
            **kwargs,
        )

    def __post_init__(self):
        if self.training_data_url is not None:
            validate_url(self.training_data_url, raise_error=True)


@dataclass
class HubModelCardHelper:
    """A helper for creating a `ModelCard` for `scvi-tools` hub models.

    Parameters
    ----------
    license_info
        The license information for the model.
    model_cls_name
        The name of the model class.
    model_init_params
        The model initialization parameters.
    model_setup_anndata_args
        The arguments used to call ``setup_anndata`` during model training.
    model_summary_stats
        The model summary statistics.
    model_data_registry
        The model data registry.
    scvi_version
        The version of `scvi-tools` that the model was trained with.
    anndata_version
        The version of anndata used during model training.
    data_modalities
        The modalities of the training data.
    tissues
        The tissues of the training data.
    data_is_annotated
        Whether the training data is annotated.
    data_is_minified
        Whether the training data uploaded with the model has been minified.
    training_data_url
        Link to the training data used to train the model, if it is too large to be uploaded to the hub. This can be
        a cellxgene explorer session url. However it cannot be a self-hosted session -- it must be from the cellxgene
        portal (https://cellxgene.cziscience.com/).
    training_code_url
        Link to the code used to train the model.
    model_parent_module
        The parent module of the model class. Defaults to `scvi.model`. Change this if you are using a model
        class that is not in the `scvi.model` module, for example, if you are using a model class from a custom module.
    description
        A description of the model.
    references_
        A list of references for the model.

    Notes
    -----
    It is not required to use this class to create a `ModelCard`. But this helps you do so in a way that is
    consistent with most other `scvi-tools` hub models. You can think of this as a template. The actual template
    string used can be found in ``scvi.template``. The resulting huggingface :class:`~huggingface_hub.ModelCard`
    can be accessed via the :meth:`~scvi.hub.HubModelCardHelper.model_card` property.
    """

    license_info: str
    model_cls_name: str
    model_init_params: dict
    model_setup_anndata_args: dict
    model_summary_stats: dict
    model_data_registry: dict
    scvi_version: str
    anndata_version: str
    data_modalities: List[str] = field(default_factory=list)
    tissues: List[str] = field(default_factory=list)
    data_is_annotated: Optional[bool] = None
    data_is_minified: Optional[bool] = None
    training_data_url: Optional[str] = None
    training_code_url: Optional[str] = None
    model_parent_module: str = _SCVI_HUB.DEFAULT_PARENT_MODULE
    description: str = _SCVI_HUB.DEFAULT_MISSING_FIELD
    references: str = _SCVI_HUB.DEFAULT_MISSING_FIELD

    def __post_init__(self):
        self.model_card = self._to_model_card()

        if self.training_data_url is not None:
            validate_url(self.training_data_url, raise_error=True)
        if self.training_code_url is not None:
            validate_url(self.training_code_url, raise_error=True)

    @classmethod
    def from_dir(
        cls,
        local_dir: str,
        license_info: str,
        anndata_version: str,
        data_is_minified: Optional[bool] = None,
        map_location: Optional[Union[torch.device, str, dict]] = "cpu",
        **kwargs,
    ):
        """Create a `HubModelCardHelper` object from a local directory.

        Parameters
        ----------
        local_dir
            The local directory containing the model files.
        license_info
            The license information for the model.
        anndata_version
            The version of anndata used during model training.
        data_is_minified
            Whether the training data uploaded with the model has been minified.
        map_location
            The device to map model tensors to, passed into :meth:`~torch.load`.
        kwargs
            Additional keyword arguments to pass to the HubModelCardHelper initializer.
        """
        attr_dict, _, _, _ = _load_saved_files(
            local_dir, load_adata=False, map_location=map_location
        )
        model_init_params = attr_dict["init_params_"]
        registry = attr_dict["registry_"]
        model_cls_name = registry["model_name"]
        scvi_version = registry["scvi_version"]
        model_setup_anndata_args = registry["setup_args"]
        model_summary_stats = dict(
            AnnDataManager._get_summary_stats_from_registry(registry)
        )
        model_data_registry = dict(
            AnnDataManager._get_data_registry_from_registry(registry)
        )

        # get `is_minified` from the param if it is given, else from adata if it on disk, else set it to None
        is_minified = data_is_minified
        if is_minified is None and os.path.isfile(f"{local_dir}/adata.h5ad"):
            is_minified = _is_minified(f"{local_dir}/adata.h5ad")

        return cls(
            license_info,
            model_cls_name,
            model_init_params,
            model_setup_anndata_args,
            model_summary_stats,
            model_data_registry,
            scvi_version,
            anndata_version,
            data_is_minified=is_minified,
            **kwargs,
        )

    def _to_model_card(self) -> ModelCard:
        # define tags
        tags = [
            "biology",
            "genomics",
            "single-cell",
            _SCVI_HUB.MODEL_CLS_NAME_TAG.format(self.model_cls_name),
            _SCVI_HUB.SCVI_VERSION_TAG.format(self.scvi_version),
            _SCVI_HUB.ANNDATA_VERSION_TAG.format(self.anndata_version),
        ]
        for m in self.data_modalities:
            tags.append(_SCVI_HUB.MODALITY_TAG.format(m))
        for t in self.tissues:
            tags.append(_SCVI_HUB.TISSUE_TAG.format(t))
        if self.data_is_annotated is not None:
            tags.append(_SCVI_HUB.ANNOTATED_TAG.format(self.data_is_annotated))

        # define the card data, which is the header
        card_data = ModelCardData(
            license=self.license_info,
            library_name=_SCVI_HUB.HF_LIBRARY_NAME,
            tags=tags,
        )

        # flatten the model_init_params into a single dict
        # for example {'kwargs': {'model_kwargs': {'foo': 'bar'}}, 'non_kwargs': {'n_hidden': 128, 'n_latent': 10}}
        # becomes {'n_hidden': 128, 'n_latent': 10, 'foo': 'bar'}
        if "non_kwargs" in self.model_init_params.keys():
            non_kwargs = self.model_init_params["non_kwargs"]
            kwargs = self.model_init_params["kwargs"]
        else:
            non_kwargs = {
                k: v
                for k, v in self.model_init_params.items()
                if not isinstance(v, dict)
            }
            kwargs = {
                k: v for k, v in self.model_init_params.items() if isinstance(v, dict)
            }
        kwargs = {k: v for (i, j) in kwargs.items() for (k, v) in j.items()}
        # kwargs and non_kwargs keys should be disjoint but if not, we'll just use the original model_init_params
        if len(set(kwargs.keys()).intersection(set(non_kwargs.keys()))) == 0:
            flattened_model_init_params = {**non_kwargs, **kwargs}
        else:
            flattened_model_init_params = self.model_init_params

        # create the content from the template
        content = template.format(
            card_data=card_data.to_yaml(),
            description=self.description,
            model_init_params=json.dumps(flattened_model_init_params, indent=4),
            model_setup_anndata_args=json.dumps(
                self.model_setup_anndata_args, indent=4
            ),
            model_summary_stats=AnnDataManager._view_summary_stats(
                self.model_summary_stats, as_markdown=True
            ),
            model_data_registry=AnnDataManager._view_data_registry(
                self.model_data_registry, as_markdown=True
            ),
            model_parent_module=self.model_parent_module,
            data_is_minified=_SCVI_HUB.DEFAULT_MISSING_FIELD
            if self.data_is_minified is None
            else self.data_is_minified,
            training_data_url=self.training_data_url or _SCVI_HUB.DEFAULT_NA_FIELD,
            training_code_url=self.training_code_url or _SCVI_HUB.DEFAULT_NA_FIELD,
            references=self.references,
        )

        # finally create and return the actual card
        return ModelCard(content)
