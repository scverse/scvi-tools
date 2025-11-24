import pandas as pd
from sklearn.gaussian_process import GaussianProcessClassifier


def _gaussian_process_classifier(
    lfc_g1_g2: pd.Series,
    lfc_n1_g2: pd.Series,
    fdr_g1_n1: pd.Series,
    lenght_scale_init: float = 10,
    lenght_scale_bounds: tuple = (1e1, 1e2),
    alpha_init: float = 0.1,
    alpha_bounds: tuple = (1e-3, 1.0),
    n_restarts_optimizer: int = 20,
    save_data: bool = True,
    restrict_to_upregulated: bool = True,
) -> GaussianProcessClassifier:
    """Train a Gaussian Process Classifier on the log fold change values of two groups.

    Parameters
    ----------
    lfc_g1_g2
        Log fold change values of group 1 vs group 2.
    lfc_n1_g2
        Log fold change values of neighbors 1 vs group 2.
    fdr_g1_n1
        FDR values of group 1 vs neighbors 1 (boolean).
    lenght_scale_init
        Initial value for the length scale hyperparameter.
    lenght_scale_bounds
        Bounds for the length scale hyperparameter.
    alpha_init
        Initial value for the alpha hyperparameter.
    alpha_bounds
        Bounds for the alpha hyperparameter.
    n_restarts_optimizer
        Number of restarts for the optimizer.
    save_data
        If True, save the input data to the classifier.
    restrict_to_upregulated
        If True, restrict the classifier to upregulated genes,
        for the group 1 versus group 2 comparison.

    Returns
    -------
    GaussianProcessClassifier
        Trained Gaussian Process Classifier.
    """
    from sklearn.gaussian_process.kernels import (
        ConstantKernel as C,
    )
    from sklearn.gaussian_process.kernels import (
        RationalQuadratic,
    )

    kernel = C(1.0, (1e-3, 1e3)) * RationalQuadratic(
        length_scale=lenght_scale_init,
        length_scale_bounds=lenght_scale_bounds,
        alpha=alpha_init,
        alpha_bounds=alpha_bounds,
    )

    if restrict_to_upregulated:
        # Identify upregulated genes
        upregulated_g1_g2 = lfc_g1_g2[lfc_g1_g2 > 0]
        selected_genes = upregulated_g1_g2.index

        # Subset features and labels for training
        lfc_g1_g2_train = lfc_g1_g2.loc[selected_genes]
        lfc_n1_g2_train = lfc_n1_g2.loc[selected_genes]
        fdr_g1_n1_train = fdr_g1_n1.loc[selected_genes]
    else:
        selected_genes = lfc_g1_g2.index  # Use all genes if not restricting
        lfc_g1_g2_train = lfc_g1_g2
        lfc_n1_g2_train = lfc_n1_g2
        fdr_g1_n1_train = fdr_g1_n1

    # The design matrix X is concat of lfc_g1_g2 and lfc_n1_g2:
    X = pd.concat(
        [lfc_g1_g2, lfc_n1_g2],
        axis=1,
    )
    X_train = pd.concat(
        [lfc_g1_g2_train, lfc_n1_g2_train],
        axis=1,
    )

    gpc = GaussianProcessClassifier(
        kernel=kernel, n_restarts_optimizer=n_restarts_optimizer, random_state=0
    ).fit(X_train, fdr_g1_n1_train)

    gpc.train_score_ = gpc.score(X_train, fdr_g1_n1_train)
    gene_probas = pd.Series(0.0, index=lfc_g1_g2.index)
    gene_probas.loc[selected_genes] = gpc.predict_proba(X_train)[:, 1]
    gpc.gene_probas_ = gene_probas

    if save_data:
        gpc.X_ = X
        gpc.y_ = fdr_g1_n1

    return gpc
