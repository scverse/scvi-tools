from .dataset10X import Dataset10X
from .dataset import GeneExpressionDataset


class PbmcDataset(GeneExpressionDataset):

    def __init__(self, save_path='data/'):
        pbmc = GeneExpressionDataset.concat_datasets(
                    Dataset10X("pbmc8k", save_path=save_path),
                    Dataset10X("pbmc4k", save_path=save_path)
                )
        super(PbmcDataset, self).__init__(pbmc.X, pbmc.local_means, pbmc.local_vars,
                                          pbmc.batch_indices, pbmc.labels, pbmc.gene_names)
