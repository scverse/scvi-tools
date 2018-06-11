from .loom import LoomDataset
import os


class RetinaDataset(LoomDataset):
    url = 'https://github.com/YosefLab/scVI-data/raw/master/retina.loom'

    def __init__(self, unit_test=False):

        self.save_path = 'data/' if not unit_test else 'tests/data/'
        self.download_name = 'retina.loom'

        if not os.path.exists(self.save_path + self.download_name):
            self.download()

        print("Finished preprocessing Retina dataset")
        super(RetinaDataset, self).__init__(self.save_path + self.download_name)
