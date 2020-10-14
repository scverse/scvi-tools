====
Data
====
.. currentmodule:: scvi

Data loading
~~~~~~~~~~~~

``scvi-tools`` now relies entirely on the AnnData_ format. For convenience, we have included data loaders from the AnnData_ API. Scanpy_ also has utilities_ to load data that are outputted by 10x's Cell Ranger software.

.. _AnnData: https://anndata.readthedocs.io/en/stable/
.. _Scanpy: https://scanpy.readthedocs.io/en/stable/index.html
.. _utilities: https://scanpy.readthedocs.io/en/stable/api/index.html#reading

.. autosummary::
   :toctree: reference/

   data.read_h5ad
   data.read_csv
   data.read_loom
   data.read_text

Basic preprocessing
~~~~~~~~~~~~~~~~~~~

For general single-cell preprocessing, we defer to our friends at Scanpy_, and specifically their preprocessing module (:mod:`scanpy.pp`).

All ``scvi-tools`` models require raw UMI count data. The count data can be safely stored in an AnnData layer as one of the first steps of a Scanpy single-cell workflow::

    adata.layers["counts"] = adata.X.copy()

Here we maintain a few package specific utilities for feature selection, etc.

.. _Scanpy: https://scanpy.readthedocs.io/en/stable/index.html

.. autosummary::
   :toctree: reference/

   data.poisson_gene_selection
   data.organize_cite_seq_10x


Data preparation
~~~~~~~~~~~~~~~~

Setting up an AnnData object is a prerequisite for running any ``scvi-tools`` model.

.. autosummary::
   :toctree: reference/

   data.setup_anndata
   data.transfer_anndata_setup
   data.register_tensor_from_anndata
   data.view_anndata_setup

Built in data
~~~~~~~~~~~~~

.. autosummary::
   :toctree: reference/

   data.pbmcs_10x_cite_seq
   data.spleen_lymph_cite_seq
   data.purified_pbmc_dataset
   data.dataset_10x
   data.brainlarge_dataset
   data.pbmc_dataset
   data.cortex
   data.seqfishplus
   data.seqfish
   data.smfish
   data.breast_cancer_dataset
   data.mouse_ob_dataset
   data.retina
   data.prefrontalcortex_starmap
   data.frontalcortex_dropseq

