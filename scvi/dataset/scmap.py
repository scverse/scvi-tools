"""The pancreas datasets gathered in the scmap paper"""
from scvi.dataset import LoomDataset


class XinDataset(LoomDataset):
    def __init__(self, save_path='data/'):
        super(XinDataset, self).__init__(
            filename='xin.loom',
            save_path=save_path,
            url='https://github.com/YosefLab/scVI-data/raw/master/xin.loom'
        )


class SegerstolpeDataset(LoomDataset):
    def __init__(self, save_path='data/'):
        super(SegerstolpeDataset, self).__init__(
            filename='segerstolpe.loom',
            save_path=save_path,
            url='https://github.com/YosefLab/scVI-data/raw/master/segerstolpe.loom'
        )


class MuraroDataset(LoomDataset):
    def __init__(self, save_path='data/'):
        super(MuraroDataset, self).__init__(
            filename='muraro.loom',
            save_path=save_path,
            url='https://github.com/YosefLab/scVI-data/raw/master/muraro.loom'
        )


class BaronDataset(LoomDataset):
    def __init__(self, save_path='data/'):
        super(BaronDataset, self).__init__(
            filename='baron.loom',
            save_path=save_path,
            url='https://github.com/YosefLab/scVI-data/raw/master/baron.loom'
        )
