====
User
====


Import scvi-tools as::

   import scvi


.. currentmodule:: scvi


Model
~~~~~

.. currentmodule:: scvi


.. autosummary::
   :toctree: reference/
   :nosignatures:

   model.AUTOZI
   model.LinearSCVI
   model.PEAKVI
   model.SCANVI
   model.SCVI
   model.TOTALVI


.. currentmodule:: scvi

External models
~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: reference/
   :nosignatures:

   external.CellAssign
   external.GIMVI
   external.RNAStereoscope
   external.SpatialStereoscope
   external.SOLO


Data loading
~~~~~~~~~~~~

``scvi-tools`` relies entirely on the AnnData_ format. For convenience, we have included data loaders from the AnnData_ API. Scanpy_ also has utilities_ to load data that are outputted by 10x's Cell Ranger software.

.. _AnnData: https://anndata.readthedocs.io/en/stable/
.. _Scanpy: https://scanpy.readthedocs.io/en/stable/index.html
.. _utilities: https://scanpy.readthedocs.io/en/stable/api/index.html#reading

.. autosummary::
   :toctree: reference/
   :nosignatures:

   data.read_h5ad
   data.read_csv
   data.read_loom
   data.read_text
   data.read_10x_atac


Data preparation
~~~~~~~~~~~~~~~~

Setting up an AnnData object is a prerequisite for running any ``scvi-tools`` model.

.. autosummary::
   :toctree: reference/
   :nosignatures:

   data.setup_anndata
   data.transfer_anndata_setup
   data.register_tensor_from_anndata
   data.view_anndata_setup


Basic preprocessing
~~~~~~~~~~~~~~~~~~~

For general single-cell preprocessing, we defer to our friends at Scanpy_, and specifically their preprocessing module (:mod:`scanpy.pp`).

All ``scvi-tools`` models require raw UMI count data. The count data can be safely stored in an AnnData layer as one of the first steps of a Scanpy single-cell workflow::

    adata.layers["counts"] = adata.X.copy()

Here we maintain a few package specific utilities for feature selection, etc.

.. _Scanpy: https://scanpy.readthedocs.io/en/stable/index.html

.. autosummary::
   :toctree: reference/
   :nosignatures:

   data.poisson_gene_selection
   data.organize_cite_seq_10x


Configuration
~~~~~~~~~~~~~

An instance of the :class:`~scvi._settings.ScviConfig` is available as ``scvi.settings`` and allows configuring scvi-tools.

.. autosummary::
   :toctree: reference/
   :nosignatures:

   _settings.ScviConfig
