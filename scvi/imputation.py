import numpy as np
from torch.utils.data import DataLoader


def imputation_error(X_pred, X, i, j, ix):
    """
    X_pred: imputed dataset
    X: original dataset
    X_zero: zeros dataset
    i, j, ix: indices of where dropout was applied
    ========
    returns:
    median L1 distance between datasets at indices given
    """
    all_index = i[ix], j[ix]
    x, y = X_pred[all_index], X[all_index]
    return np.median(np.abs(x - y))


def imputation(vae, gene_dataset):
    if gene_dataset.get_all().is_cuda:
        X = gene_dataset.get_all().cpu().numpy()
    else:
        X = gene_dataset.get_all().numpy()

    gene_dropout_dataset, i, j, ix = gene_dataset.dropout(rate=0.1)
    data_loader_dropout = DataLoader(gene_dropout_dataset, batch_size=len(gene_dropout_dataset),
                                     shuffle=True, num_workers=1)
    for sample_batch, local_l_mean, local_l_var, batch_index in data_loader_dropout:
        _, _, px_rate, _, _, _, _, _ = vae(sample_batch, batch_index)
    if px_rate.data.is_cuda:
        mae = imputation_error(px_rate.data.cpu().numpy(), X, i, j, ix)
    else:
        mae = imputation_error(px_rate.data.numpy(), X, i, j, ix)
    return mae
