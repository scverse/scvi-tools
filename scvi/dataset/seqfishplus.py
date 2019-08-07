import os
import zipfile

import pandas as pd

from scvi.dataset import DownloadableDataset, CellMeasurement


class SeqFishPlusDataset(DownloadableDataset):
    """seqFISH+ can image mRNAs for 10,000 genes in single cells—with high accuracy and
    sub-diffraction-limit resolution—in the cortex, subventricular zone
    and olfactory bulb of mouse brain

    :param tissue_region: Region of the mouse brain, Either "subventricular cortex" or "olfactory bulb"
    :param save_path: Location to use when saving/loading the SeqFish+ data.
    :param delayed_populating: Switch for delayed populating mechanism.
    """

    def __init__(
        self,
        tissue_region: str = "subventricular cortex",
        save_path: str = "data",
        delayed_populating: bool = False,
    ):

        self.tissue_region = tissue_region
        if tissue_region == "subventricular cortex":
            self.file_prefix = "cortex_svz"
        elif tissue_region == "olfactory bulb":
            self.file_prefix = "ob"
        else:
            raise ValueError(
                '`tissue_type` must be "subventricular cortex" or "olfactory bulb", but got {}'.format(
                    tissue_region
                )
            )
        super().__init__(
            urls="https://github.com/CaiGroup/seqFISH-PLUS/raw/master/sourcedata.zip",
            filenames="seqfishplus.zip",
            save_path=save_path,
            delayed_populating=delayed_populating,
        )

    def populate(self):
        counts_filename = "{}_counts.csv".format(self.file_prefix)
        coordinates_filename = "{}_cellcentroids.csv".format(self.file_prefix)
        data_path = os.path.join(self.save_path, "seqfishplus")
        if not os.path.exists(data_path):
            os.makedirs(data_path)
        with zipfile.ZipFile(os.path.join(self.save_path, self.filenames[0])) as f:
            f.extract(counts_filename, path=data_path)
            f.extract(coordinates_filename, path=data_path)
        df_counts = pd.read_csv(os.path.join(data_path, counts_filename))
        df_coordinates = pd.read_csv(os.path.join(data_path, coordinates_filename))
        coordinates = CellMeasurement(
            name="coords",
            data=df_coordinates[["X", "Y"]],
            columns_attr_name="axis",
            columns=["x", "y"],
        )
        cell_attributes_dict = {"Cell ID": None, "Field of View": None}
        if self.tissue_region == "subventricular cortex":
            cell_attributes_dict.update({"Region": None})
        for key in cell_attributes_dict.keys():
            cell_attributes_dict[key] = df_coordinates[key]
        self.populate_from_data(
            X=df_counts.values,
            gene_names=df_counts.columns,
            Ys=[coordinates],
            cell_attributes_dict=cell_attributes_dict,
        )
