from unittest import TestCase

from scvi.dataset import BiomDataset
from .utils import unsupervised_training_one_epoch, unsupervised_nb_training_one_epoch


class TestBiomDataset(TestCase):
    # fix once NB is fixed
    # def test_populate(self):
    #     dataset = BiomDataset(filename="feature-table.biom",
    #                           save_path="tests/data")
    #     unsupervised_training_one_epoch(dataset)

    def test_populate_nb(self):
        dataset = BiomDataset(filename="feature-table.biom",
                              save_path="tests/data")
        unsupervised_nb_training_one_epoch(dataset)

    def test_ag_url(self):
        url = (r'https://github.com/biocore/American-Gut/blob/master/'
               r'data/AG/AG_100nt.biom?raw=true')
        dataset = BiomDataset(filename='AG_100nt.biom', url=url, save_path="tests/data")
        unsupervised_nb_training_one_epoch(dataset)

