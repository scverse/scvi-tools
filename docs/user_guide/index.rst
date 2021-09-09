User guide
===========


scRNA-seq analysis
--------------------

scvi-tools hosts implementations of the following models:

.. list-table::
   :widths: 15 15 15 100
   :header-rows: 1

   * - Model
     - Reference
     - Class
     - Tasks
   * - :doc:`/user_guide/models/scvi`
     - [Lopez18]_
     - :class:`~scvi.model.SCVI`
     - Dimensionality reduction, removal of unwanted variation, integration across replicates, donors, and technologies, differential expression, imputation, normalization of other cell- and sample-level confounding factors
   * - :doc:`/user_guide/models/linearscvi`
     - :class:`~scvi.model.LinearSCVI`
     - scVI tasks with linear decoder
   * - :doc:`/user_guide/models/cellassign`
     - [Zhang19]_
     - :class:`~scvi.external.CellAssign`
     - Marker-based automated annotation

.. toctree::
    :maxdepth: 2
    :hidden:

    models/scvi
    models/cellassign


CITE-seq analysis
-----------------

.. list-table::
   :widths: 15 15 15 100
   :header-rows: 1

   * - Model
     - Reference
     - Class
     - Tasks
   * - :doc:`/user_guide/models/totalvi`
     - [GayosoSteier21]_
     - :class:`~scvi.model.TOTALVI`
     - Dimensionality reduction, removal of unwanted variation, integration across replicates, donors, and technologies, differential expression, protein imputation, imputation, normalization of other cell- and sample-level confounding factors

.. toctree::
    :maxdepth: 2
    :hidden:

    models/totalvi



Background
-----------------

.. toctree::
    :maxdepth: 2

    background/variational_inference
    background/differential_expression
    background/counterfactual_prediction
    background/transfer_learning

