import logging
import os
from collections import namedtuple

import numpy as np
import pandas as pd

from scvi.dataset.dataset import CellMeasurement, DownloadableDataset

logger = logging.getLogger(__name__)

available_datasets = {
    "cbmc": "CBMC_8K_13AB_10X",
    "pbmc": "PBMC_vs_flow_10X",
    "cd8": "CD8_merged",
}
CiteSeqFilenames = namedtuple(
    "CiteSeqFilenames", field_names=["rna", "adt", "adt_centered"]
)


class CiteSeqDataset(DownloadableDataset):
    """Allows to form 3 different CiteSeq datasets.

    Note that their centered log ratio transformation for ADT counts is different from
    the standard clr transformation: they explain they add pseudocounts (for 0 values),
    but do not explicit the actual transformation.
    It doesn't seem to be simply adding count 1 to all entries, or only 0 entries.

    :param name: Name of the CiteSeq dataset to load. Either "cbmc", "pbmc" or "cd8".
    :param save_path: Location to use when saving/loading the data.
    :param delayed_populating: Switch for delayed populating mechanism.
    """

    def __init__(
        self,
        name: str = "cbmc",
        save_path: str = "data/citeSeq/",
        delayed_populating: bool = False,
    ):
        s = available_datasets[name]
        filenames = CiteSeqFilenames(
            rna="%s_rna.csv.gz" % name,
            adt="%s_adt.csv.gz" % name,
            adt_centered="%s_adt_centered.csv.gz" % name,
        )
        super().__init__(
            urls=[
                "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100866/suppl/GSE100866_%s-RNA_umi.csv.gz"
                % s,
                "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100866/suppl/GSE100866_%s-ADT_umi.csv.gz"
                % s,
                "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100866/suppl/"
                "GSE100866_%s-ADT_clr-transformed.csv.gz" % s,
            ],
            filenames=filenames,
            save_path=os.path.join(save_path, name),
            delayed_populating=delayed_populating,
        )

    def populate(self):
        logger.info("Preprocessing data")
        self.expression = pd.read_csv(
            os.path.join(self.save_path, self.filenames.rna),
            index_col=0,
            compression="gzip",
        ).T

        # process protein measurements
        adt = pd.read_csv(
            os.path.join(self.save_path, self.filenames.adt),
            index_col=0,
            compression="gzip",
        )
        protein_names = np.asarray(adt.index).astype(np.str)
        protein_measurement = CellMeasurement(
            name="protein_expression",
            data=adt.T.values,
            columns_attr_name="protein_names",
            columns=protein_names
        )
        adt_centered = pd.read_csv(
            os.path.join(self.save_path, self.filenames.adt_centered),
            index_col=0,
            compression="gzip",
        )
        if not np.array_equal(np.asarray(adt_centered.index).astype(np.str), protein_names):
            raise ValueError("Protein names are not the same for raw and centered counts.")
        protein_measurement_centered = CellMeasurement(
            name="protein_expression_clr",
            data=adt_centered.T.values,
            columns_attr_name="protein_names_clr",
            columns=protein_names,
        )

        # keep only human genes (there are also mouse genes)
        gene_names = np.asarray(self.expression.columns, dtype=str)
        human_filter = np.asarray(
            [name.startswith("HUMAN") for name in gene_names], dtype=np.bool
        )
        logger.info(
            "Selecting only HUMAN genes ({} / {})".format(
                human_filter.sum(), len(human_filter)
            )
        )
        X = self.expression.values[:, human_filter]
        gene_names = gene_names[human_filter]
        gene_names = np.char.upper(
            np.asarray(
                [name.split("_")[-1] if "_" in name else name for name in gene_names],
                dtype=np.str,
            )
        )

        logger.info("Finish preprocessing data")

        self.populate_from_data(
            X=X,
            gene_names=gene_names,
            Ys=[protein_measurement, protein_measurement_centered]
        )

        self.filter_cells_by_count()


class CbmcDataset(CiteSeqDataset):
    """Loads cbmc dataset.

    This dataset that includes 8,617 cord blood mononuclear cells profiled using 10x along with for each cell 13
    well-characterized mononuclear antibodies. We kept the top 600 genes by variance.

    Args:
        :save_path: Save path of raw data file. Default: ``'data/'``.

    Examples:
        >>> gene_dataset = CbmcDataset()

    """

    def __init__(
        self, save_path: str = "data/citeSeq/", delayed_populating: bool = False
    ):
        super().__init__(
            name="cbmc", save_path=save_path, delayed_populating=delayed_populating
        )
