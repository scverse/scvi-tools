import random
from scipy.stats import pearsonr
import numpy as np
import tqdm
import pandas
import anndata as ad
import torch
import scvi

from typing import List, Dict

def create_anndata_pseudobulk(adata: ad.AnnData,
                              x: np.array) -> ad.AnnData:
    """Creates an anndata object from a pseudobulk sample.

    Parameters
    ----------
    adata: ad.AnnData
        AnnData aobject storing training set
    x: np.array
        pseudobulk sample

    Return
    ------
    ad.AnnData
        Anndata object storing the pseudobulk array
    """
    df_obs = pandas.DataFrame.from_dict([{col: adata.obs[col].value_counts().index[0] for col in adata.obs.columns}])
    adata_pseudobulk = ad.AnnData(X=x, obs=df_obs)
    adata_pseudobulk.layers["counts"] = np.copy(x)

    return adata_pseudobulk


def replace_inf(x: np.array,
                strategy: str = "mean") -> np.array:
    """Checks if array contains nan values and replaces them

    Parmaters
    ---------
    x: np.array
        Input numpy array
    stratrgy: str
        Imputaton strategy

    Returns
    -------
    np.array
        Array without nan values.
    """
    if np.isinf(x).any():
        # Calculate the mean of non-infinite values
        finite_values = x[np.isfinite(x)]
        mean_of_finite = np.mean(finite_values)
        # Replace infinite values with the mean of finite values
        x[np.isinf(x)] = mean_of_finite
    return x

@torch.inference_mode()
def sanity_checks_metrics_latent(model: scvi.model.SCVI,
                          adata: ad.AnnData,
                          batch_sizes: List[int],
                          n_repeats: int,
                          use_get_latent: bool=True,
                          impute_nans: bool=False,
                          ) -> Dict:
    """Computes sanity check metrics for a given scVI model.

    Parameters
    ----------
    model: scvi.model.SCVI
        fitted scvi model
    adata: ad.AnnData
        Anndata object sotring the training set
    batch_sizes: List[int]
        List of batch size (# on sampled cells)
    n_repeats: int
        Number of repeats
    use_get_latent: bool
        Whether to use get_latent scvi built in function
        to infer the latent states. Is set to `False`, will
        use the module encoder directly.
    impute_nans: bool
        Whether to impute NaNs

    Returns
    -------
    Dict
        Dictionnary storing latent space metrics
    """
    metrics = {"corr": [], "kl": []}
    errors = {"corr": [], "kl": []}
    for n in tqdm.tqdm(batch_sizes):
        current_metrics = {"corr": [], "kl": []}
        for i in range(n_repeats):
            # Calculate the number of cells to sample from each cell type proportionally
            sampled_cells_per_type = adata.obs['cell_type'].value_counts(normalize=True) * n
            sampled_cells_per_type = sampled_cells_per_type.astype(int)

            # Perform stratified sampling for each cell type
            sampled_cells = []
            for cell_type, num_cells in sampled_cells_per_type.items():
                seed = random.seed()
                sampled_cells.extend(adata.obs[adata.obs['cell_type'] == cell_type].sample(n=num_cells,
                                                                                        random_state=seed).index)

            # Select the sampled cells from the DataFrame
            adata_sampled = adata[sampled_cells]

            # embeddings of single-cells
            if use_get_latent:
                latent_sampled = model.get_latent_representation(adata_sampled)
                mean_sampled_z = latent_sampled.mean(axis=0)
            else:
                sampled_x = torch.from_numpy(adata_sampled.layers["counts"].toarray()).to("cuda:0")
                outputs = model.module._regular_inference(x=sampled_x,
                                                          batch_index=0)
                latent_sampled = outputs["z"].detach().cpu().numpy()
                mean_sampled_z = latent_sampled.mean(axis=0)
                if impute_nans:
                    mean_sampled_z = replace_inf(mean_sampled_z)

            # pseudo-bulk embedding
            pseudobulk_x = adata_sampled.layers["counts"].mean(axis=0) #.astype(int).astype(numpy.float32)
            if use_get_latent:
                adata_pseudobulk = create_anndata_pseudobulk(adata, pseudobulk_x)
                pseudobulk_z = model.get_latent_representation(adata_pseudobulk).flatten()
            else:
                outputs  = model.module._regular_inference(x=torch.from_numpy(pseudobulk_x).to("cuda:0"),
                                                           batch_index=0)
                pseudobulk_z = outputs["z"].detach().cpu().numpy().flatten()
                if impute_nans:
                    pseudobulk_z = replace_inf(pseudobulk_z)

            # Compute correlation
            pearson_corr = pearsonr(mean_sampled_z, pseudobulk_z)
            current_metrics["corr"].append(pearson_corr[0])

            # compute kl-divergence
            # kl = kl_divergence(dist_z,
            #                    dist_pseudobulk_z).sum(dim=-1)


        metrics["corr"].append(np.mean(current_metrics["corr"]))
        errors["corr"].append(np.std(current_metrics["corr"]))

    return metrics, errors


@torch.inference_mode()
def sanity_checks_metrics_feature(model: scvi.model.SCVI,
                          adata: ad.AnnData,
                          batch_sizes: List[int],
                          n_repeats: int,
                          use_get_normalized: bool=False,
                          ) -> Dict:
    """Computes feature space sanity check metrics for a given scVI model.

    Parameters
    ----------
    model: scvi.model.SCVI
        fitted scvi model
    adata: ad.AnnData
        Anndata object sotring the training set
    batch_sizes: List[int]
        List of batch size (# on sampled cells)
    n_repeats: int
        Number of repeats
    use_get_normalized: bool
        Whether to use get_latent scvi built in function
        to infer the latent states. Is set to `False`, will
        use the module encoder directly.
    impute_nans: bool
        Whether to impute NaNs

    Returns
    -------
    Dict
        Dictionnary storing feature space metrics
    """
    metrics = {"corr": []}
    errors = {"corr": []}
    for n in tqdm.tqdm(batch_sizes):
        current_metrics = {"corr": [], "kl": []}
        for i in range(n_repeats):
            # Calculate the number of cells to sample from each cell type proportionally
            sampled_cells_per_type = adata.obs['cell_type'].value_counts(normalize=True) * n
            sampled_cells_per_type = sampled_cells_per_type.astype(int)

            # Perform stratified sampling for each cell type
            sampled_cells = []
            for cell_type, num_cells in sampled_cells_per_type.items():
                seed = random.seed()
                sampled_cells.extend(adata.obs[adata.obs['cell_type'] == cell_type].sample(n=num_cells,
                                                                                        random_state=seed).index)

            # Select the sampled cells from the DataFrame
            adata_sampled = adata[sampled_cells]

            # reconstruction of single-cells
            if use_get_normalized:
                rec_sampled = model.get_normalized_expression(adata_sampled)
                mean_sampled_rec = rec_sampled.mean(axis=0).values
            else:
                sampled_x = torch.from_numpy(adata_sampled.layers["counts"].toarray()).to("cuda:0")
                sampled_rec_x = model.module.sample(tensors={
                                                    "X": sampled_x,
                                                    "batch": 0,
                                                    "extra_continuous_covs": None,
                                                    "extra_categorical_covs": None,
                                                    "labels": 0
                                                    })
                sampled_rec_x = sampled_rec_x.detach().cpu().numpy()
                mean_sampled_rec = sampled_rec_x.mean(axis=0)

            # reconstruction of pseudo-bulk
            pseudobulk_x = adata_sampled.layers["counts"].mean(axis=0) #.astype(int).astype(numpy.float32)
            if use_get_normalized:
                adata_pseudobulk = create_anndata_pseudobulk(adata, pseudobulk_x)
                pseudobulk_rec_x = model.get_normalized_expression(adata_pseudobulk).values.flatten()
            else:
                pseudobulk_x = torch.from_numpy(pseudobulk_x).to("cuda:0")
                pseudobulk_rec_x = model.module.sample(tensors={
                                                    "X": pseudobulk_x,
                                                    "batch": 0,
                                                    "extra_continuous_covs": None,
                                                    "extra_categorical_covs": None,
                                                    "labels": 0
                                                    })
                pseudobulk_rec_x = pseudobulk_rec_x.detach().cpu().numpy().flatten()

            # Compute correlation
            pearson_corr = pearsonr(mean_sampled_rec, pseudobulk_rec_x)
            current_metrics["corr"].append(pearson_corr[0])


        metrics["corr"].append(np.mean(current_metrics["corr"]))
        errors["corr"].append(np.std(current_metrics["corr"]))

    return metrics, errors
