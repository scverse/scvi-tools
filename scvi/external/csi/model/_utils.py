from typing import Optional, Union

import pandas as pd


def prepare_metadata(
    meta_data: pd.DataFrame,
    cov_cat_keys: Optional[list] = None,
    cov_cat_embed_keys: Optional[list] = None,
    cov_cont_keys: Optional[list] = None,
    categ_orders: Optional[dict] = None,
    key_orders: Optional[dict] = None,
):
    """
    Prepare content of dataframe columns for model training (one hot encoding, encoding for embedding, ...)

    Parameters
    ----------
    meta_data
        Dataframe containing metadata columns.
    cov_cat_keys
        List of categorical covariates column names to be one-hot encoded.
    cov_cat_embed_keys
        List of categorical covariates column names to be embedded.
    cov_cont_keys
        List of continuous covariates column names.
    categ_orders
        Defined orders for categorical covariates. Dict with keys being
        categorical covariates keys and values being lists of categories. May contain more
        categories than data.
    key_orders
        Defines order of covariate columns. Dict with keys being 'categorical', 'categorical_embed', 'continuous'
        and values being lists of keys.

    Returns
    -------
    Tuple of: covariate data that does not require embedding,
    covariate data that requires embedding,
    dict with order of categories per covariate (as orders),
    dict with keys (categorical, categorical_embed, and continuous) specifying as values
    order of covariates used to construct the two covariate datas

    """

    def get_categories_order(values: pd.Series, categories: Union[list, None] = None):
        """
        Helper to get order of categories based on values and optional list of categories

        Parameters
        ----------
        values
            Categorical values
        categories
            Optional order of categories

        Returns
        -------
            Categories order
        """
        if categories is None:
            categories = pd.Categorical(values).categories.values
        else:
            missing = set(values.unique()) - set(categories)
            if len(missing) > 0:
                raise ValueError(
                    f"Some values of {values.name} are not in the specified categories order: {missing}"
                )
        return list(categories)

    if cov_cat_keys is None:
        cov_cat_keys = []
    if cov_cat_embed_keys is None:
        cov_cat_embed_keys = []
    if cov_cont_keys is None:
        cov_cont_keys = []

    # Check & set order of covariates and categories
    if key_orders is not None:
        assert set(key_orders["categorical"]) == set(cov_cat_keys)
        cov_cat_keys = key_orders["categorical"]
        assert set(key_orders["categorical_embed"]) == set(cov_cat_embed_keys)
        cov_cat_embed_keys = key_orders["categorical_embed"]
        assert set(key_orders["continuous"]) == set(cov_cont_keys)
        cov_cont_keys = key_orders["continuous"]
    cov_dict = {
        "categorical": cov_cat_keys,
        "categorical_embed": cov_cat_embed_keys,
        "continuous": cov_cont_keys,
    }

    if categ_orders is None:
        categ_orders = {}
        for cov_key in cov_cat_keys + cov_cat_embed_keys:
            categ_orders[cov_key] = get_categories_order(
                values=meta_data[cov_key], categories=None
            )

    def dummies_categories(values: pd.Series, categories: list):
        """
        Make dummies of categorical covariates. Use specified order of categories.

        Parameters
        ----------
        values
            Categorical vales for each observation.
        categories
            Order of categories to use

        Returns
        -------
        dummies - one-hot encoding of categories in same order as categories.
        """
        # Get dummies
        # Ensure ordering
        values = pd.Series(
            pd.Categorical(values=values, categories=categories, ordered=True),
            index=values.index,
            name=values.name,
        )
        # This is problematic if many covariates
        dummies = pd.get_dummies(values, prefix=values.name)

        return dummies

    # Covs that are not embedded: continuous and one-hot encoded categorical covariates
    if len(cov_cat_keys) > 0 or len(cov_cont_keys) > 0:
        cov_cat_data = []
        for cov_cat_key in cov_cat_keys:
            cat_dummies = dummies_categories(
                values=meta_data[cov_cat_key], categories=categ_orders[cov_cat_key]
            )
            cov_cat_data.append(cat_dummies)
        # Prepare single cov array for all covariates
        cov_data_parsed = pd.concat(cov_cat_data + [meta_data[cov_cont_keys]], axis=1)
    else:
        cov_data_parsed = None

    # Data of covariates to be embedded
    if len(cov_cat_embed_keys) > 0:
        cov_embed_data = []
        for cov_cat_embed_key in cov_cat_embed_keys:
            cat_order = categ_orders[cov_cat_embed_key]
            cat_map = dict(zip(cat_order, range(len(cat_order))))
            cov_embed_data.append(meta_data[cov_cat_embed_key].map(cat_map))
        cov_embed_data = pd.concat(cov_embed_data, axis=1)
    else:
        cov_embed_data = None

    return cov_data_parsed, cov_embed_data, categ_orders, cov_dict
