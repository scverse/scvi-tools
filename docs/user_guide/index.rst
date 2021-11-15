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
   * - :doc:`/user_guide/models/scanvi`
     - [Xu21]_
     - scVI tasks with cell type transfer from reference, seed labeling
   * - :doc:`/user_guide/models/linearscvi`
     - [Svensson20]_
     - scVI tasks with linear decoder
   * - :doc:`/user_guide/models/autozi`
     - [Clivio19]_
     -  for assessing gene-specific levels of zero-inflation in scRNA-seq data 
   * - :doc:`/user_guide/models/cellassign`
     - [Zhang19]_
     - Marker-based automated annotation
   * - :doc:`/user_guide/models/solo`
     - [Bernstein19]_
     - Doublet detection

.. toctree::
    :maxdepth: 2
    :hidden:

    models/scvi
    models/linearscvi
    models/cellassign


ATAC-seq analysis
-----------------

.. list-table::
   :widths: 15 15 100
   :header-rows: 1

   * - Model
     - Reference
     - Tasks
   * - :doc:`/user_guide/models/peakvi`
     - [Ashuach21]_
     - Dimensionality reduction, removal of unwanted variation, integration across replicates, donors, and technologies, differential expression, imputation, normalization of other cell- and sample-level confounding factors

.. toctree::
    :maxdepth: 2
    :hidden:

    models/peakvi


Multimodal analysis
--------------------

CITE-seq
^^^^^^^^^

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

Multiome
^^^^^^^^

.. list-table::
   :widths: 15 15 100
   :header-rows: 1

   * - Model
     - Reference
     - Tasks
   * - :doc:`/user_guide/models/multivi`
     - [AshuachGabitto21]_
     - Integration of paired/unpaired multiome data, missing modality imputation, normalization of other cell- and sample-level confounding factors

.. toctree::
    :maxdepth: 2
    :hidden:

    models/multivi


Spatial transcriptomics analysis
--------------------------------

.. list-table::
   :widths: 15 15 100
   :header-rows: 1

   * - Model
     - Reference
     - Tasks
   * - :doc:`/user_guide/models/destvi`
     - [Lopez21]_
     - Multi-resolution deconvolution, cell-type-specific gene expression imputation, comparative analysis
   * - :doc:`/user_guide/models/stereoscope`
     - [Andersson20]_
     - Deconvolution
   * - :doc:`/user_guide/models/gimvi`
     - [Lopez19]_
     - Imputation of missing spatial genes

.. toctree::
    :maxdepth: 2
    :hidden:

    models/destvi
    models/stereoscope
    model/gimvi


General purpose analysis
------------------------

.. list-table::
   :widths: 15 15 100
   :header-rows: 1

   * - Model
     - Reference
     - Tasks
   * - :doc:`/user_guide/models/amortizedlda`
     - [Blei03]_
     - Topic modeling 

.. toctree::
    :maxdepth: 2
    :hidden:

    models/amortizedlda


Background
-----------------

.. toctree::
    :maxdepth: 2

    background/variational_inference
    background/differential_expression
    background/counterfactual_prediction
    background/transfer_learning
    background/codebase_overview


