User guide
===========


scRNA-seq analysis
--------------------

scvi-tools hosts implementations of the following models:

.. list-table::
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

    models/scvi
    models/cellassign


CITE-seq analysis
-----------------

.. list-table::
   :widths: 35 100
   :header-rows: 1

   * - Model
     - Tasks
   * - totalVI [GayosoSteier21]_
     - Dimensionality reduction, removal of unwanted variation, integration across replicates, donors, and technologies, differential expression, protein imputation, imputation, normalization of other cell- and sample-level confounding factors

.. toctree::
    :maxdepth: 2

    models/totalvi


