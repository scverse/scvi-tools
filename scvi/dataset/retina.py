from .loom import LoomDataset
import os


class RetinaDataset(LoomDataset):
    url = 'https://github.com/YosefLab/scVI-data/raw/master/retina.loom'

    def __init__(self, save_path='data/'):

        self.save_path = save_path
        self.download_name = 'retina.loom'

        if not os.path.exists(self.save_path + self.download_name):
            self.download()

        print("Finished preprocessing Retina dataset")
        super(RetinaDataset, self).__init__(filename=self.download_name, save_path=self.save_path)
