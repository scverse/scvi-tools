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

from ._anndata import get_from_registry
from ._utils import _check_nonnegative_integers

logger = logging.getLogger(__name__)


class AnnDataRecorder:
    """
    Backend to setup AnnData.

    Constructs the scvi-tools setup dictionary to be placed in
    `adata.uns["_scvi"]`.
    """

    def __init__(self, adata: AnnData) -> None:

        self.adata = adata
        self.setup_dict = dict(
            data_registry={}, summary_stats={}, scvi_version=scvi.__version__
        )
        self.add_to_summary_stats("n_cells", self.adata.shape[0])
        self.add_to_summary_stats("n_vars", self.adata.shape[1])

    def get_setup_dict(self):
        return self.setup_dict

    def setup_batch(
        self,
        batch_key: str,
    ):
        """
        Wrapper of `setup_categorical_obs_key` to setup batch (sample) data.

        Parameters
        ----------
        batch_key
            key in `adata.obs` for batch information. Categories will automatically be converted into integer
            categories and saved to `adata.obs['_scvi_batch']`. If `None`, assigns the same batch to all the data.
        """
        recorded_key = "_scvi_batch"
        self.setup_categorical_obs_key(
            obs_key=batch_key,
            setup_param_name="batch_key",
            registry_key=_CONSTANTS.BATCH_KEY,
            recorded_key=recorded_key,
            default_same_cat=True,
        )
        n_batch = len(
            np.unique(self.setup_dict["categorical_mappings"][recorded_key]["mapping"])
        )
        self.add_to_summary_stats("n_batch", n_batch)

    def setup_labels(
        self,
        labels_key: str,
    ):
        """
        Wrapper of `setup_categorical_obs_key` to setup label data.

        Parameters
        ----------
        labels_key
            key in `adata.obs` for batch information. Categories will automatically be converted into integer
            categories and saved to `adata.obs['_scvi_batch']`. If `None`, assigns the same batch to all the data.
        """
        recorded_key = "_scvi_labels"
        self.setup_categorical_obs_key(
            obs_key=labels_key,
            setup_param_name="labels_key",
            registry_key=_CONSTANTS.LABELS_KEY,
            recorded_key=recorded_key,
            default_same_cat=True,
        )
        n_labels = len(
            np.unique(self.setup_dict["categorical_mappings"][recorded_key]["mapping"])
        )
        self.add_to_summary_stats("n_labels", n_labels)

    def setup_categorical_obs_key(
        self,
        obs_key: str,
        setup_param_name: str,
        registry_key: str,
        recorded_key: Optional[str] = None,
        default_same_cat: bool = False,
    ):
        """
        Sets up a key in obs with categorical data.

        Parameters
        ----------
        obs_key
            Key inputted by user containing data
        setup_param_name
            Name of parameter in `setup_anndata()` where `obs_key` was passed
        registry_key
            Key to use to represent this tensor in setup dict data registry
        default_same_cat
            If `obs_key` is `None`, register it assuming data is all same category
        """
        recorded_key = (
            "_scvi_{}".format(setup_param_name.split("_")[0])
            if recorded_key is None
            else recorded_key
        )
        if setup_param_name is None and default_same_cat is True:
            logger.info(
                "No {} inputted, assuming all cells are same category".format(
                    setup_param_name
                )
            )
            self.adata.obs[recorded_key] = np.zeros(self.adata.shape[0], dtype=np.int64)
        else:
            _assert_key_in_obs(self.adata, obs_key)
            logger.info(
                'Using data from adata.obs["{}"] for {}'.format(
                    obs_key, setup_param_name
                )
            )
        recorded_key = self._make_obs_column_categorical(
            column_key=obs_key, alternate_column_key=recorded_key
        )
        self.add_to_data_registry(registry_key, "obs", recorded_key)

    def setup_continuous_obs_key(self):
        raise NotImplementedError

    def setup_categorical_obs_key_iterable(
        self,
        categorical_covariate_keys: List[str],
        recorded_key: str,
        category_dict: Dict[str, List[str]] = None,
        registry_key: Optional[str] = None,
    ):
        """
        Setup obsm df for extra categorical covariates.

        Parameters
        ----------
        categorical_covariate_keys
            List of keys in adata.obs with categorical data
        recorded_key
            New df added to `.obsm` with copied data from `categorical_covariate_keys`
        category_dict
            Optional dictionary with keys being keys of categorical data in obs
            and values being precomputed categories for each obs vector
        registry_key
            Key to add to registry
        """
        adata = self.adata
        for key in categorical_covariate_keys:
            _assert_key_in_obs(adata, key)

        info_store = {}

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
        info_store["mappings"] = store_cats
        # this preserves the order of the keys added to the df
        info_store["keys"] = categorical_covariate_keys

        # how many cats per key, in the preserved order
        n_cats_per_key = []
        for k in categorical_covariate_keys:
            n_cats_per_key.append(len(store_cats[k]))
        info_store["n_cats_per_key"] = n_cats_per_key

        self.setup_dict["categorical_obsm_keys"][recorded_key] = info_store

        self.add_to_data_registry(
            _CONSTANTS.CAT_COVS_KEY if registry_key is None else registry_key,
            "obsm",
            recorded_key,
        )

    def setup_continuous_obs_key_iterable(
        self,
        continuous_covariate_keys: List[str],
        recorded_key: str,
        registry_key: Optional[str] = None,
    ):
        """
        Setup obsm df for extra continuous covariates.

        Parameters
        ----------
        continuous_covariate_keys
            List of keys in adata.obs with continuous data
        recorded_key
            New df added to `.obsm` with copied data from `continuous_covariate_keys`
        registry_key
            Key to add to registry
        """
        adata = self.adata
        for key in continuous_covariate_keys:
            _assert_key_in_obs(adata, key)

        info_store = {}

        series = []
        for key in continuous_covariate_keys:
            s = adata.obs[key]
            series.append(s)

        adata.obsm[recorded_key] = pd.concat(series, axis=1)
        info_store["keys"] = adata.obsm[recorded_key].columns.to_numpy()

        # add info to setup dict
        self.setup_dict["continuous_obsm_keys"][recorded_key] = info_store
        self.add_to_data_registry(
            _CONSTANTS.CONT_COVS_KEY if registry_key is None else registry_key,
            "obsm",
            recorded_key,
        )

    def setup_x(self, layer: str):
        adata = self.adata
        if layer is not None:
            assert (
                layer in adata.layers.keys()
            ), "{} is not a valid key in adata.layers".format(layer)
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

        self.add_to_data_registry(_CONSTANTS.X_KEY, x_loc, x_key)
        self._verify_and_correct_data_format(_CONSTANTS.X_KEY)

    def setup_obsm_key(
        self, obsm_key: str, registry_key: str, check_data_format: bool = True
    ):
        """
        Setup key with data in `.obsm`.

        Parameters
        ----------
        obsm_key
            Key in `.obsm`
        registry_key
            Key to add to registry
        check_data_format
            Check and correct that data is CSR if sparse, and C_CONTIGUOUS
        """
        adata = self.adata
        assert (
            obsm_key in adata.obsm.keys()
        ), "{} is not a valid key in adata.obsm".format(obsm_key)

        self.add_to_data_registry(registry_key, "obsm", obsm_key)
        if check_data_format:
            self._verify_and_correct_data_format(registry_key)

    def add_to_data_registry(
        self, registry_key: str, adata_attr_name: str, adata_attr_key: str
    ):
        """Registers the tensor in the local setup dictionaries data registry."""
        self.setup_dict["data_registry"][registry_key] = dict(
            attr_name=adata_attr_name, attr_key=adata_attr_key
        )

    def add_to_summary_stats(self, key: str, val: Union[int, float]):
        self.setup_dict["summary_stats"][key] = val

    def _make_obs_column_categorical(
        self, column_key: str, alternate_column_key: str, categorical_dtype=None
    ):
        """
        Makes the data in column_key in obs all categorical.

        If adata.obs[column_key] is not categorical, will first categorize
        and then save to .obs[alternate_column_key].
        """
        if categorical_dtype is None:
            categorical_obs = self.adata.obs[column_key].astype("category")
        else:
            categorical_obs = self.adata.obs[column_key].astype(categorical_dtype)

        # put codes in .obs[alternate_column_key]
        codes = categorical_obs.cat.codes
        mapping = categorical_obs.cat.categories.to_numpy(copy=True)
        if -1 in np.unique(codes):
            received_categories = (
                self.adata.obs[column_key].astype("category").cat.categories
            )
            raise ValueError(
                'Making .obs["{}"] categorical failed. Expected categories: {}. '
                "Received categories: {}. ".format(
                    column_key, mapping, received_categories
                )
            )
        self.adata.obs[alternate_column_key] = codes

        # store categorical mappings
        store_dict = {
            alternate_column_key: {"original_key": column_key, "mapping": mapping}
        }
        if "categorical_mappings" not in self.setup_dict.keys():
            self.setup_dict.update({"categorical_mappings": store_dict})
        else:
            self.setup_dict["categorical_mappings"].update(store_dict)

        # make sure each category contains enough cells
        unique, counts = np.unique(
            self.adata.obs[alternate_column_key], return_counts=True
        )
        if np.min(counts) < 3:
            category = unique[np.argmin(counts)]
            warnings.warn(
                "Category {} in adata.obs['{}'] has fewer than 3 cells. SCVI may not train properly.".format(
                    category, alternate_column_key
                )
            )
        # possible check for continuous?
        if len(unique) > (self.adata.shape[0] / 3):
            warnings.warn(
                "Is adata.obs['{}'] continuous? Please provide categorical data."
            )
        return alternate_column_key

    def _verify_and_correct_data_format(self, registry_key: str):
        """
        Will make sure that the user's anndata is C_CONTIGUOUS and csr if it is dense numpy or sparse respectively.

        Parameters
        ----------
        registry_key
            Key in data registry with data "path"
        """
        adata = self.adata
        k = registry_key
        data = get_from_registry(adata, k)
        if isspmatrix(data) and (data.getformat() != "csr"):
            logger.warning(
                "Training will be faster when sparse matrix is formatted as CSR. It is safe to cast before model initialization."
            )
        elif isinstance(data, np.ndarray) and (data.flags["C_CONTIGUOUS"] is False):
            logger.debug(
                "{} is not C_CONTIGUOUS. Overwriting to C_CONTIGUOUS.".format(k)
            )
            data = np.asarray(data, order="C")
            _set_data_in_registry(adata, data, k)
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
            _set_data_in_registry(adata, data, k)


def _assert_key_in_obs(adata, key):
    assert key in adata.obs.keys(), "{} is not a valid key for in adata.obs".format(key)


def _set_data_in_registry(adata, data, key):
    """
    Sets the data associated with key in adata.uns['_scvi']['data_registry'].keys() to data.

    Note: This is a dangerous method and will change the underlying data of the user's anndata
    Currently used to make the user's anndata C_CONTIGUOUS and csr if it is dense numpy
    or sparse respectively.

    Parameters
    ----------
    adata
        anndata object to change data of
    data
        data to change to
    key
        key in adata.uns['_scvi]['data_registry'].keys() associated with the data
    """
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
