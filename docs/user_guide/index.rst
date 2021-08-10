User guide
===========


Single-cell RNA-sequencing analysis
------------------------------------

scvi-tools hosts implementations of the following models:

.. list-table:: Overview of scRNA-seq models
   :widths: 35 100
   :header-rows: 1

   * - Model
     - Tasks
   * - scVI [Lopez18]_
     - Dimensionality reduction, removal of unwanted variation, integration across replicates, donors, and technologies, differential expression, imputation, normalization of other cell- and sample-level confounding factors
   * - CellAssign [Zhang19]_
     - Marker-based automated annotation

.. toctree::
    :maxdepth: 2

    models/cellassign


External models
----------------


