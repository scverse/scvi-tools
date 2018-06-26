from .loom import LoomDataset


class RetinaDataset(LoomDataset):

    def __init__(self, save_path='data/'):
        super(RetinaDataset, self).__init__(filename='retina.loom',
                                            save_path=save_path,
                                            url='https://github.com/YosefLab/scVI-data/raw/master/retina.loom')
