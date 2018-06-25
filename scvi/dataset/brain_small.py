from .dataset10X import Dataset10X


class BrainSmallDataset(Dataset10X):
    def __init__(self, save_path='data/'):
        super(BrainSmallDataset, self).__init__(filename="neuron_9k",
                                                save_path=save_path)
