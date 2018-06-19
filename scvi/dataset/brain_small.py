from .dataset10X import Dataset10X


class BrainSmallDataset(Dataset10X):
    def __init__(self):
        super(BrainSmallDataset, self).__init__(name="neuron_9k", p_genes=3000)
