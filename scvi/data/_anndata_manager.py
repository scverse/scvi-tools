import logging
import warnings
from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd
from anndata import AnnData
from pandas.api.types import CategoricalDtype
from scipy.sparse import isspmatrix

import scvi
from scvi import _CONSTANTS

from ._utils import (
    _assert_key_in_obs,
    _assert_key_in_obsm,
    _check_nonnegative_integers,
    _get_batch_mask_protein_data,
)

logger = logging.getLogger(__name__)


class AnnDataField:
    """
    Represents a field within AnnData (obs/var/obsm/varm/etc.).

    Parameters
    ----------
    key
        The key of the field within the anndata attribute it belongs to. For example an obs field
        with name "batch" (adata.obs["batch"]) has key "batch".
    recorded_key
        The key we use in anndata to store scvi-internal information related to this field. For
        example if the field is an obs field with categorical data and ``recorded_key`` "_scvi_batch",
        we add its categorical codes to adata.obs["_scvi_batch"].
        If there is no extra scvi-internal information, ``recorded-key`` can be None.
    registry_key
        The key we use to represent this field in the scvi data registry, which is a dictionary that we store
        in adata.uns["scvi"]. For example an obs field with ``recorded_key`` "_scvi_batch" and ``registry_key``
        "obs_batch" is stored in the scvi data registry as ``'obs_batch' : {'attr_name': 'obs', 'attr_key': '_scvi_batch'}``.
    is_categorical
        ``True`` if this is a categorical field, ``False`` otherwise
    categorical_dtype
        The type to cast the data into. Only applicable if ``is_categorical == True``
    check_data_format
        Verify (and correct if not) that the data is CSR if sparse, and C_CONTIGUOUS
    """

    def __init__(
        self,
        key: str,
        recorded_key: Optional[str],
        registry_key: str,
        is_categorical: bool = True,
        categorical_dtype: Optional[CategoricalDtype] = None,
        check_data_format: bool = True,
    ) -> None:
        self.key = key
        self.recorded_key = recorded_key
        self.registry_key = registry_key
        self.is_categorical = is_categorical
        self.categorical_dtype = categorical_dtype
        self.check_data_format = check_data_format


class AnnDataManager:
    """
    Manages the anndata object that we use in scvi models.
    This class is the interface for all other components in the package to
    operate on anndata. It helps with setting up anndata (i.e. registering
    model-specific components onto it), validating the setup of the anndata
    against a reference anndata, transferring setup metadata if needed, and
    other anndata-related tasks.

    Parameters
    ----------
    adata
        The AnnData object containing the data
    """

    def __init__(self, adata: AnnData):
        if adata.is_view:
            raise ValueError("Please run `adata = adata.copy()`.")
        self.adata = adata
        # Set up the setup dict and some of its keys
        self.adata.uns["_scvi"] = dict(
            data_registry={},
            summary_stats={},
            scvi_version=scvi.__version__,
            categorical_obsm_keys={},
            continuous_obsm_keys={},
            categorical_mappings={},
        )
        self._add_to_summary_stats("n_cells", self.adata.shape[0])
        self._add_to_summary_stats("n_vars", self.adata.shape[1])
        self._add_to_summary_stats("n_proteins", 0)
        self._add_to_summary_stats("n_categorical_covs", 0)
        self._add_to_summary_stats("n_continuous_covs", 0)

    @property
    def setup_dict(self):
        return self.adata.uns["_scvi"]

    def _add_to_data_registry(
        self, registry_key: str, adata_attr_name: str, adata_attr_key: str
    ):
        """
        Registers the tensor associated with ``adata_attr_name`` and ``adata_attr_key``
        in the scvi data registry by adding it to the registry dictionary.
        """
        self.setup_dict["data_registry"][registry_key] = dict(
            attr_name=adata_attr_name, attr_key=adata_attr_key
        )
        logger.debug("Registered key: {}".format(registry_key))

    def _add_to_summary_stats(self, key: str, val: Union[int, float]):
        self.setup_dict["summary_stats"][key] = val

    def setup_batch(
        self,
        batch_key: Optional[str],
    ):
        """
        Wrapper of ``setup_obs_field`` to setup batch (sample) data.

        Parameters
        ----------
        batch_key
            key in `adata.obs` for batch information. Categories will be automatically converted into integer
            categories and saved to `adata.obs["_scvi_batch"]`. If `None`, assigns the same batch to all the data.
        """
        recorded_key = "_scvi_batch"
        field = AnnDataField(batch_key, recorded_key, _CONSTANTS.BATCH_KEY)
        self.setup_obs_field(field)
        n_batch = len(
            np.unique(self.setup_dict["categorical_mappings"][recorded_key]["mapping"])
        )
        self._add_to_summary_stats("n_batch", n_batch)

    def setup_labels(
        self,
        labels_key: Optional[str],
    ):
        """
        Wrapper of ``setup_obs_field`` to setup label data.

        Parameters
        ----------
        labels_key
            key in `adata.obs` for label information. Categories will be automatically converted into integer
            categories and saved to `adata.obs['_scvi_labels']`. If `None`, assigns the same label to all the data.
        """
        recorded_key = "_scvi_labels"
        field = AnnDataField(labels_key, recorded_key, _CONSTANTS.LABELS_KEY)
        self.setup_obs_field(field)
        n_labels = len(
            np.unique(self.setup_dict["categorical_mappings"][recorded_key]["mapping"])
        )
        self.add_to_summary_stats("n_labels", n_labels)

    def setup_obs_field(
        self,
        field: AnnDataField,
    ):
        """
        Sets up a field in ``adata.obs``.

        Parameters
        ----------
        field
            ``AnnDataField`` to set up
        """
        if not field.is_categorical:
            raise NotImplementedError

        if field.key is None:
            logger.info(
                "No obs key inputted, assuming all cells are same category in {}".format(
                    field.recorded_key
                )
            )
            self.adata.obs[field.recorded_key] = np.zeros(
                self.adata.shape[0], dtype=np.int64
            )
            obs_key = field.recorded_key
        else:
            _assert_key_in_obs(self.adata, field.key)
            logger.info('Using data from adata.obs["{}"]'.format(field.key))
        self._make_obs_field_categorical(
            obs_key=obs_key,
            recorded_obs_key=field.recorded_key,
            categorical_dtype=field.categorical_dtype,
        )
        self._add_to_data_registry(field.registry_key, "obs", field.recorded_key)

    def _make_obs_field_categorical(
        adata, obs_key, recorded_obs_key, categorical_dtype=None
    ):
        """
        Makes the data in ``adata.obs[obs_key]`` all categorical.

        If ``adata.obs[column_key]`` is not categorical, will categorize
        and save to ``adata.obs[recorded_obs_key]``
        """
        categorical_obs = adata.obs[obs_key].astype(
            "category" if categorical_dtype is None else categorical_dtype
        )

        # put codes in .obs[recorded_obs_key]
        codes = categorical_obs.cat.codes
        mapping = categorical_obs.cat.categories.to_numpy(copy=True)
        if -1 in np.unique(codes):
            received_categories = adata.obs[obs_key].astype("category").cat.categories
            raise ValueError(
                'Making .obs["{}"] categorical failed. Expected categories: {}. '
                "Received categories: {}. ".format(
                    obs_key, mapping, received_categories
                )
            )
        adata.obs[recorded_obs_key] = codes

        # store categorical mappings in the scvi data registry
        registry_dict = {
            recorded_obs_key: {"original_key": obs_key, "mapping": mapping}
        }
        adata.uns["_scvi"]["categorical_mappings"].update(registry_dict)

        # make sure each category contains enough cells
        unique, counts = np.unique(adata.obs[recorded_obs_key], return_counts=True)
        if np.min(counts) < 3:
            category = unique[np.argmin(counts)]
            warnings.warn(
                "Category {} in adata.obs['{}'] has fewer than 3 cells. SCVI may not train properly.".format(
                    category, recorded_obs_key
                )
            )
        # possible check for continuous?
        if len(unique) > (adata.shape[0] / 3):
            warnings.warn(
                "Is adata.obs['{}'] continuous? SCVI doesn't support continuous obs yet."
            )

    def setup_protein_expression(
        self, protein_expression_obsm_key: str, batch_key: str
    ):
        prot_exp = self.adata.obsm[protein_expression_obsm_key]
        if _check_nonnegative_integers(prot_exp) is False:
            warnings.warn(
                "adata.obsm[{}] does not contain unnormalized count data. Are you sure this is what you want?".format(
                    protein_expression_obsm_key
                )
            )

        field = AnnDataField(
            protein_expression_obsm_key, None, _CONSTANTS.PROTEIN_EXP_KEY
        )
        self.setup_obsm_field(field)
        logger.info(
            "Using protein expression from adata.obsm['{}']".format(
                protein_expression_obsm_key
            )
        )

        # set up protein names
        if isinstance(self.adata.obsm[protein_expression_obsm_key], pd.DataFrame):
            logger.info(
                "Using protein names from columns of adata.obsm['{}']".format(
                    protein_expression_obsm_key
                )
            )
            protein_names = np.asarray(
                self.adata.obsm[protein_expression_obsm_key].columns.to_list()
            )
        else:
            logger.info("Generating sequential protein names")
            protein_names = np.arange(
                self.adata.obsm[protein_expression_obsm_key].shape[1]
            )

        self.setup_dict["protein_names"] = protein_names
        self._add_to_summary_stats("n_proteins", len(protein_names))

        # batch mask totalVI
        batch_mask = _get_batch_mask_protein_data(
            self.adata, protein_expression_obsm_key, batch_key
        )

        # check if it's actually needed
        if np.sum([~b[1] for b in batch_mask.items()]) > 0:
            logger.info("Found batches with missing protein expression")
            self.setup_dict["totalvi_batch_mask"] = batch_mask

    def setup_obsm_field(
        self,
        field: AnnDataField,
    ):
        """
        Sets up a field in ``adata.obsm``.

        Parameters
        ----------
        field
            ``AnnDataField`` to set up
        """
        _assert_key_in_obsm(self.adata, field.key)

        self._add_to_data_registry(field.registry_key, "obsm", field.key)
        if field.check_data_format:
            self._verify_and_correct_data_format(field.registry_key)

    def setup_x(self, layer: Optional[str]):
        adata = self.adata
        if layer is not None:
            if layer not in adata.layers.keys():
                raise KeyError("{} is not a valid key in adata.layers".format(layer))
            logger.info('Using data from adata.layers["{}"]'.format(layer))
            x_loc = "layers"
            x_key = layer
            x = adata.layers[x_key]
        else:
            logger.info("Using data from adata.X")
            x_loc = "X"
            x_key = "None"
            x = adata.X

        if _check_nonnegative_integers(x) is False:
            logger_data_loc = (
                "adata.X" if layer is None else "adata.layers[{}]".format(layer)
            )
            warnings.warn(
                "{} does not contain unnormalized count data. Are you sure this is what you want?".format(
                    logger_data_loc
                )
            )

        self._add_to_data_registry(_CONSTANTS.X_KEY, x_loc, x_key)
        self._verify_and_correct_data_format(_CONSTANTS.X_KEY)

    def record_extra_categorical_covariates(
        self,
        categorical_covariate_keys: List[str],
        recorded_key: str = "_scvi_extra_categoricals",
        category_dict: Dict[str, List[str]] = None,
        registry_key: str = _CONSTANTS.CAT_COVS_KEY,
    ):
        """
        Records extra categorical covariates (which are obs keys) in
        an ``adata.obsm`` dataframe.

        This helps set up a series of obs keys that will be loaded into
        models as one tensor.

        Parameters
        ----------
        categorical_covariate_keys
            List of keys in ``adata.obs`` with categorical data
        recorded_key
            Key in ``adata.obsm`` where we record the dataframe with copied data
            from ``adata.obs[categorical_covariate_keys]``
        category_dict
            Optional dictionary with keys being keys of categorical data in obs
            and values being precomputed categories for each obs vector
        registry_key
            Keys under which to register the corresponding obsm field in the scvi
            data registry
        """
        adata = self.adata
        for key in categorical_covariate_keys:
            _assert_key_in_obs(adata, key)

        categories = {}
        df = pd.DataFrame(index=adata.obs_names)
        for key in categorical_covariate_keys:
            if category_dict is None:
                categorical_obs = adata.obs[key].astype("category")
                mapping = categorical_obs.cat.categories.to_numpy(copy=True)
                categories[key] = mapping
            else:
                possible_cats = category_dict[key]
                categorical_obs = adata.obs[key].astype(
                    CategoricalDtype(categories=possible_cats)
                )
            codes = categorical_obs.cat.codes
            df[key] = codes

        adata.obsm[recorded_key] = df

        store_cats = categories if category_dict is None else category_dict
        info_store = {}
        info_store["mappings"] = store_cats
        # this preserves the order of the keys added to the df
        info_store["keys"] = categorical_covariate_keys

        # how many categories per key, in the preserved order
        n_cats_per_key = []
        for k in categorical_covariate_keys:
            n_cats_per_key.append(len(store_cats[k]))
        info_store["n_cats_per_key"] = n_cats_per_key

        self.setup_dict["extra_categoricals"] = info_store

        field = AnnDataField(
            recorded_key,
            None,
            registry_key,
            check_data_format=False,
        )
        self.setup_obsm_field(field)

        self._add_to_summary_stats(
            "n_categorical_covs", len(categorical_covariate_keys)
        )

    def record_extra_continuous_covariates(
        self,
        continuous_covariate_keys: List[str],
        recorded_key: str = "_scvi_extra_continuous",
        registry_key: str = _CONSTANTS.CAT_COVS_KEY,
    ):
        """
        Records extra continuous covariates (which are obs keys) in
        an ``adata.obsm`` dataframe.

        This helps set up a series of obs keys that will be loaded into
        models as one tensor.

        Parameters
        ----------
        continuous_covariate_keys
            List of keys in ``adata.obs`` with continuous data
        recorded_key
            Key in ``adata.obsm`` where we record the dataframe with copied data
            from ``adata.obs[continuous_covariate_keys]``
        registry_key
            Keys under which to register the corresponding obsm field in the scvi
            data registry
        """
        adata = self.adata
        for key in continuous_covariate_keys:
            _assert_key_in_obs(adata, key)

        series = []
        for key in continuous_covariate_keys:
            s = adata.obs[key]
            series.append(s)

        adata.obsm[recorded_key] = pd.concat(series, axis=1)

        # add info to setup dict
        self.setup_dict["extra_continuous_keys"] = adata.obsm[
            recorded_key
        ].columns.to_numpy()

        field = AnnDataField(
            recorded_key,
            None,
            registry_key,
            check_data_format=False,
        )
        self.setup_obsm_field(field)

        self._add_to_summary_stats("n_continuous_covs", len(continuous_covariate_keys))

    def _verify_and_correct_data_format(self, registry_key: str):
        """
        Will make sure that the user's anndata is C_CONTIGUOUS and csr if it is dense numpy or sparse respectively.

        Parameters
        ----------
        registry_key
            Data registry key that represents the data we want to verify
        """
        adata = self.adata
        k = registry_key
        data = self.get_from_registry(adata, k)
        if isspmatrix(data) and (data.getformat() != "csr"):
            logger.warning(
                "Training will be faster when sparse matrix is formatted as CSR. It is safe to cast before model initialization."
            )
        elif isinstance(data, np.ndarray) and (data.flags["C_CONTIGUOUS"] is False):
            logger.debug(
                "{} is not C_CONTIGUOUS. Overwriting to C_CONTIGUOUS.".format(k)
            )
            data = np.asarray(data, order="C")
            self._set_data_in_anndata(adata, data, k)
        elif isinstance(data, pd.DataFrame) and (
            data.to_numpy().flags["C_CONTIGUOUS"] is False
        ):
            logger.debug(
                "{} is not C_CONTIGUOUS. Overwriting to C_CONTIGUOUS.".format(k)
            )
            index = data.index
            vals = data.to_numpy()
            columns = data.columns
            data = pd.DataFrame(
                np.ascontiguousarray(vals), index=index, columns=columns
            )
            self._set_data_in_anndata(adata, data, k)

    @staticmethod
    def get_from_registry(adata: AnnData, key: str) -> np.ndarray:
        """
        Returns the object in AnnData associated with the key in ``.uns['_scvi']['data_registry']``.

        Parameters
        ----------
        adata
            anndata object already setup with ``setup_anndata``
        key
            key of object to get from ``adata.uns['_scvi]['data_registry']``

        Returns
        -------
        The requested data

        Examples
        --------
        >>> import scvi
        >>> adata = scvi.data.cortex()
        >>> adata.uns['_scvi']['data_registry']
        {'X': ['_X', None],
        'batch_indices': ['obs', 'batch'],
        'labels': ['obs', 'labels']}
        >>> batch = get_from_registry(adata, "batch_indices")
        >>> batch
        array([[0],
            [0],
            [0],
            ...,
            [0],
            [0],
            [0]])
        """
        data_loc = adata.uns["_scvi"]["data_registry"][key]
        attr_name, attr_key = data_loc["attr_name"], data_loc["attr_key"]

        data = getattr(adata, attr_name)
        if attr_key != "None":
            if isinstance(data, pd.DataFrame):
                data = data.loc[:, attr_key]
            else:
                data = data[attr_key]
        if isinstance(data, pd.Series):
            data = data.to_numpy().reshape(-1, 1)
        return data

    def _set_data_in_anndata(self, data, key):
        """
        Sets the data associated with key in adata.uns['_scvi']['data_registry'].keys() to data.

        Note: This is a dangerous method and can change the underlying data of the user's anndata
        depdending on what attributes of adata you choose to change. Currently used to make the
        user's anndata C_CONTIGUOUS and csr if it is dense numpy or sparse respectively.

        Parameters
        ----------
        adata
            anndata object to change data of
        data
            data to change to
        key
            key in adata.uns['_scvi]['data_registry'].keys() associated with the data
        """
        adata = self.adata
        data_loc = adata.uns["_scvi"]["data_registry"][key]
        attr_name, attr_key = data_loc["attr_name"], data_loc["attr_key"]

        if attr_key == "None":
            setattr(adata, attr_name, data)
        elif attr_key != "None":
            attribute = getattr(adata, attr_name)
            if isinstance(attribute, pd.DataFrame):
                attribute.loc[:, attr_key] = data
            else:
                attribute[attr_key] = data
            setattr(adata, attr_name, attribute)
