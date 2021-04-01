import scanpy as sc
from sklearn.metrics import silhouette_score


def silhouette_metric(model):
    model.is_trained_ = True
    latent = model.get_latent_representation()
    model.is_trained_ = False
    adata = model.adata
    adata.obsm["X_scVI"] = latent
    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.leiden(adata, key_added="leiden_scVI", resolution=0.5)
    return silhouette_score(
        adata.obsm["X_scVI"],
        adata.obs["leiden_scVI"],
    )
