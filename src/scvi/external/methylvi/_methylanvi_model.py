from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from scvi.external.methylvi import METHYLVI_REGISTRY_KEYS

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence
    from typing import Literal

    from anndata import AnnData
    from mudata import MuData

import numpy as np
import pandas as pd
import torch

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager, fields
from scvi.data._constants import _SETUP_ARGS_KEY
from scvi.external.methylvi._base_components import BSSeqMixin
from scvi.external.methylvi._utils import _context_cov_key, _context_mc_key
from scvi.model.base import (
    ArchesMixin,
    BaseModelClass,
    SemisupervisedTrainingMixin,
    VAEMixin,
)
from scvi.train import SemiSupervisedTrainingPlan
from scvi.utils import setup_anndata_dsp

from ._methylanvi_module import METHYLANVAE

logger = logging.getLogger(__name__)


class METHYLANVI(VAEMixin, BSSeqMixin, SemisupervisedTrainingMixin, ArchesMixin, BaseModelClass):
    """Methylation annotation using variational inference :cite:p:`Weinberger23`.

    Inspired from M1 + M2 model, as described in (https://arxiv.org/pdf/1406.5298.pdf).

    Parameters
    ----------
    mdata
        MuData object registered via :meth:`~scvi.external.methylvi.METHYLVI.setup_mudata`.
    n_hidden
        Number of nodes per hidden layer.
    n_latent
        Dimensionality of the latent space.
    n_layers
        Number of hidden layers used for encoder and decoder NNs.
    dropout_rate
        Dropout rate for neural networks.
    likelihood
        One of
        * ``'betabinomial'`` - BetaBinomial distribution
        * ``'binomial'`` - Binomial distribution
    dispersion
        One of the following
        * ``'region'`` - dispersion parameter of BetaBinomial is constant per region across cells
        * ``'region-cell'`` - dispersion can differ for every region in every cell
    linear_classifier
        If ``True``, uses a single linear layer for classification instead of a
        multi-layer perceptron.
    **model_kwargs
        Keyword args for :class:`~scvi.module.SCANVAE`

    Examples
    --------
    >>> mdata = mudata.read_h5mu(path_to_mudata)
    >>> scvi.external.methylvi.METHYLANVI.setup_mudata(
    ...     mdata, labels_key="labels", unlabeled_category="Unknown"
    ... )
    >>> vae = scvi.external.methylvi.METHYLANVI(mdata)
    >>> vae.train()
    >>> mdata.obsm["X_scVI"] = vae.get_latent_representation()
    >>> mdata.obs["pred_label"] = vae.predict()

    """

    _module_cls = METHYLANVAE
    _training_plan_cls = SemiSupervisedTrainingPlan

    def __init__(
        self,
        mdata: MuData,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        likelihood: Literal["betabinomial", "binomial"] = "betabinomial",
        dispersion: Literal["region", "region-cell"] = "region",
        linear_classifier: bool = False,
        **model_kwargs,
    ):
        super().__init__(mdata)
        methylanvae_model_kwargs = dict(model_kwargs)

        self._set_indices_and_labels()

        # ignores unlabeled category
        n_labels = self.summary_stats.n_labels - 1
        n_cats_per_cov = (
            self.adata_manager.get_state_registry(REGISTRY_KEYS.CAT_COVS_KEY).n_cats_per_key
            if REGISTRY_KEYS.CAT_COVS_KEY in self.adata_manager.data_registry
            else None
        )

        n_batch = self.summary_stats.n_batch

        self.contexts = self.get_anndata_manager(mdata, required=True).registry[_SETUP_ARGS_KEY][
            "methylation_contexts"
        ]
        self.num_features_per_context = [mdata[context].shape[1] for context in self.contexts]

        n_input = np.sum(self.num_features_per_context)

        self.module = self._module_cls(
            n_input=n_input,
            n_batch=n_batch,
            n_cats_per_cov=n_cats_per_cov,
            n_labels=n_labels,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            dispersion=dispersion,
            likelihood=likelihood,
            linear_classifier=linear_classifier,
            contexts=self.contexts,
            num_features_per_context=self.num_features_per_context,
            **methylanvae_model_kwargs,
        )

        self.unsupervised_history_ = None
        self.semisupervised_history_ = None

        self._model_summary_string = (
            f"MethylANVI Model with the following params: \nunlabeled_category: "
            f"{self.unlabeled_category_}, n_hidden: {n_hidden}, n_latent: {n_latent}"
            f", n_layers: {n_layers}, dropout_rate: {dropout_rate}, dispersion: "
            f"{dispersion}, likelihood: {likelihood}"
        )
        self.init_params_ = self._get_init_params(locals())
        self.was_pretrained = False
        self.n_labels = n_labels

    def predict(
        self,
        mdata: MuData | None = None,
        indices: Sequence[int] | None = None,
        soft: bool = False,
        batch_size: int | None = None,
        use_posterior_mean: bool = True,
    ) -> np.ndarray | pd.DataFrame:
        """Return cell label predictions.

        Parameters
        ----------
        mdata
            MuData object registered via :meth:`~scvi.external.methylvi.MethylANVI.setup_mudata`.
        indices
            Return probabilities for each class label.
        soft
            If True, returns per class probabilities
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        use_posterior_mean
            If ``True``, uses the mean of the posterior distribution to predict celltype
            labels. Otherwise, uses a sample from the posterior distribution - this
            means that the predictions will be stochastic.
        """
        mdata = self._validate_anndata(mdata)

        if indices is None:
            indices = np.arange(mdata.n_obs)

        scdl = self._make_data_loader(
            adata=mdata,
            indices=indices,
            batch_size=batch_size,
        )
        y_pred = []
        for _, tensors in enumerate(scdl):
            inference_inputs = self.module._get_inference_input(tensors)  # (n_obs, n_vars)

            mc = inference_inputs[METHYLVI_REGISTRY_KEYS.MC_KEY]
            cov = inference_inputs[METHYLVI_REGISTRY_KEYS.COV_KEY]
            batch = tensors[REGISTRY_KEYS.BATCH_KEY]

            pred = self.module.classify(
                mc,
                cov,
                batch_index=batch,
                use_posterior_mean=use_posterior_mean,
            )
            if self.module.classifier.logits:
                pred = torch.nn.functional.softmax(pred, dim=-1)
            if not soft:
                pred = pred.argmax(dim=1)
            y_pred.append(pred.detach().cpu())

        y_pred = torch.cat(y_pred).numpy()
        if not soft:
            predictions = []
            for p in y_pred:
                predictions.append(self._code_to_label[p])

            return np.array(predictions)
        else:
            n_labels = len(pred[0])
            pred = pd.DataFrame(
                y_pred,
                columns=self._label_mapping[:n_labels],
                index=mdata.obs_names[indices],
            )
            return pred

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        **kwargs,
    ) -> AnnData | None:
        """
        %(summary)s.

        Parameters
        ----------
        %(param_adata)s

        Returns
        -------
        %(returns)s
        """
        raise NotImplementedError("METHYLANVI must be used with a MuData object.")

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_mudata(
        cls,
        mdata: MuData,
        mc_layer: str,
        cov_layer: str,
        labels_key: str,
        unlabeled_category: str,
        methylation_contexts: Iterable[str],
        batch_key: str | None = None,
        categorical_covariate_keys: Iterable[str] | None = None,
        modalities=None,
        **kwargs,
    ):
        """%(summary_mdata)s.

        Parameters
        ----------
        %(param_mdata)s
        mc_layer
            Layer containing methylated cytosine counts for each set of methylation features.
        cov_layer
            Layer containing total coverage counts for each set of methylation features.
        labels_key
            Obs field in `mdata` object containing cell type labels
        unlabeled_category
            Value of `mdata.obs[labels_key]` representing an unknown cell type label
        methylation_contexts
            List of modality fields in `mdata` object representing different methylation contexts.
            Each context must be equipped with a layer containing the number of methylated counts
            (specified by `mc_layer`) and total number of counts (specified by `cov_layer`) for
            each genomic region feature.
        %(param_batch_key)s
        %(param_categorical_covariate_keys)s
        %(param_modalities)s

        Examples
        --------
        METHYLANVI.setup_mudata(
            mdata,
            mc_layer="mc",
            cov_layer="cov",
            labels_key="CellType",
            unlabeled_category="Unknown",
            methylation_contexts=["mCG", "mCH"],
            categorical_covariate_keys=["Platform"],
            modalities={
                "categorical_covariate_keys": "mCG"
            },
        )

        """
        if modalities is None:
            modalities = {}
        setup_method_args = METHYLANVI._get_setup_method_args(**locals())

        if methylation_contexts is None:
            raise ValueError("Methylation contexts cannot be None.")

        modalities_ = cls._create_modalities_attr_dict(modalities, setup_method_args)

        batch_field = fields.MuDataCategoricalObsField(
            REGISTRY_KEYS.BATCH_KEY,
            batch_key,
            mod_key=modalities_.batch_key,
        )

        cat_cov_field = fields.MuDataCategoricalJointObsField(
            REGISTRY_KEYS.CAT_COVS_KEY,
            categorical_covariate_keys,
            mod_key=modalities_.categorical_covariate_keys,
        )

        cell_type_field = fields.MuDataLabelsWithUnlabeledObsField(
            REGISTRY_KEYS.LABELS_KEY,
            labels_key,
            unlabeled_category,
            mod_key=modalities_.labels_key,
        )

        mc_fields = []
        cov_fields = []

        for context in methylation_contexts:
            mc_fields.append(
                fields.MuDataLayerField(
                    _context_mc_key(context),
                    mc_layer,
                    mod_key=context,
                    is_count_data=True,
                    mod_required=True,
                )
            )

            cov_fields.append(
                fields.MuDataLayerField(
                    _context_cov_key(context),
                    cov_layer,
                    mod_key=context,
                    is_count_data=True,
                    mod_required=True,
                )
            )

        mudata_fields = (
            mc_fields + cov_fields + [batch_field] + [cat_cov_field] + [cell_type_field]
        )
        adata_manager = AnnDataManager(fields=mudata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(mdata, **kwargs)

        cls.register_manager(adata_manager)
