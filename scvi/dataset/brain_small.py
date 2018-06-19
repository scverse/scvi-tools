from .dataset10X import Dataset10X


class BrainSmallDataset(Dataset10X):
    def __init__(self, save_path='data/'):
        super(BrainSmallDataset, self).__init__(name="neuron_9k", p_genes=3000, save_path=save_path)
