User guide
===========


scRNA-seq analysis
--------------------

.. list-table::
   :widths: 15 15 100
   :header-rows: 1

   * - Model
     - Reference
     - Tasks
   * - :doc:`/user_guide/models/scvi`
     - [Lopez18]_
     - Dimensionality reduction, removal of unwanted variation, integration across replicates, donors, and technologies, differential expression, imputation, normalization of other cell- and sample-level confounding factors
   * - :doc:`/user_guide/models/linearscvi`
     - [Svensson20]_
     - scVI tasks with linear decoder
   * - :doc:`/user_guide/models/cellassign`
     - [Zhang19]_
     - Marker-based automated annotation

.. toctree::
    :maxdepth: 2
    :hidden:

    models/scvi
    models/linearscvi
    models/cellassign


CITE-seq analysis
-----------------

.. list-table::
   :widths: 15 15 100
   :header-rows: 1

   * - Model
     - Reference
     - Tasks
   * - :doc:`/user_guide/models/totalvi`
     - [GayosoSteier21]_
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

