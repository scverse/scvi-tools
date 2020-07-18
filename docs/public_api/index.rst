.. module:: scvi
.. automodule:: scvi
   :noindex:

API
===


Import scvi as::

   import scvi


Dataset
-------

.. module:: scvi.dataset
.. currentmodule:: scvi

Data Loading
~~~~~~~~~~~~

``scvi`` now relies entirely on the AnnData_ format. Please refer to the AnnData_ documentation for loading data from a variety of formats (e.g., `.csv`, `.loom`, `.h5ad`).

Scanpy_ also has utilities_ to load data that are outputted by 10X's Cell Ranger software.

.. _AnnData: https://anndata.readthedocs.io/en/stable/
.. _Scanpy: https://scanpy.readthedocs.io/en/stable/index.html
.. _utilities: https://scanpy.readthedocs.io/en/stable/api/index.html#reading

Data preparation
~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: .

   dataset.setup_anndata


Basic Preprocessing
~~~~~~~~~~~~~~~~~~~

For general single-cell preprocessing, we defer to our friends at Scanpy_, and specifically their preprocessing module (:mod:`scanpy.pp`).

Here we maintain a few package specific utilities for feature selection, etc.

.. _Scanpy: https://scanpy.readthedocs.io/en/stable/index.html

.. autosummary::
   :toctree: .

   dataset.highly_variable_genes_seurat_v3
   dataset.poisson_gene_selection



Inference
---------

.. module:: scvi.inference
.. currentmodule:: scvi


Posterior
~~~~~~~~~


.. autosummary::
   :toctree: .

   inference.Posterior
