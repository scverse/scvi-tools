from unittest import TestCase

from scvi.dataset import SeqFishPlusDataset
from .utils import unsupervised_training_one_epoch


class TestSeqFishPlusDataset(TestCase):
    def test_populate(self):
        for tissue_region in ["subventricular cortex", "olfactory bulb"]:
            dataset = SeqFishPlusDataset(
                tissue_region=tissue_region, save_path="tests/data"
            )
            unsupervised_training_one_epoch(dataset)
