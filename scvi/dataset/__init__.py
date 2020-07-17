from ._anndata import setup_anndata
from ._preprocessing import (
    highly_variable_genes_seurat_v3,
    poisson_gene_selection,
    organize_cite_seq_10x,
)

from .cite_seq import pbmcs_10x_cite_seq
from .dataset10X import dataset10X
from .brain_large import brainlarge_dataset
from ._datasets import synthetic_iid
from .pbmc import pbmc_dataset
from .cortex import cortex
from .seqfish import seqfish
from .seqfishplus import seqfishplus
from .smfish import smfish
from .csv import breast_cancer_dataset, mouse_ob_dataset
from .loom import (
    retina,
    prefrontalcortex_starmap,
    frontalcortex_dropseq,
    annotation_simulation,
)
