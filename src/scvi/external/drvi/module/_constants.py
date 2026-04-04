from __future__ import annotations

from scvi.module._constants import _MODULE_KEYS as _SCVI_MODULE_KEYS


class _DRVI_MODULE_KEYS(_SCVI_MODULE_KEYS):
    # generative
    PX_PARAMS_KEY = "px_params"
    PX_UNAGGREGATED_PARAMS_KEY = "px_unaggregated_params"
    # Extra
    LIKELIHOOD_ADDITIONAL_PARAMS_KEY: str = "gene_likelihood_additional_info"
    X_MASK_KEY: str = "x_mask"
    N_SAMPLES_KEY: str = "n_samples"
    RECONSTRUCTION_INDICES: str = "reconstruction_indices"
    # Tensor IO structure
    CONT_COVS_TENSOR_KEY: str = "cont_full_tensor"
    CAT_COVS_TENSOR_KEY: str = "cat_full_tensor"
    # Loss
    MSE_LOSS_KEY: str = "mse"
    RECONSTRUCTION_LOSS_KEY: str = "reconstruction_loss"


MODULE_KEYS = _DRVI_MODULE_KEYS()
