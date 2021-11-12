import scanpy as sc
from sklearn.metrics import silhouette_score


def silhouette_metric(model, labels_key):
    model.is_trained_ = True
    latent = model.get_latent_representation()
    model.is_trained_ = False
    adata = model.adata
    adata.obsm["X_scvi"] = latent
    sc.pp.neighbors(adata, use_rep="X_scVI")
    return silhouette_score(
        adata.obsp["distances"],
        adata.obs[labels_key],
    )
