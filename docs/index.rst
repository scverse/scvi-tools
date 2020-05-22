.. include:: ../README.rst
   :end-line: 22
.. include:: _authors.rst

scVI is a package for end-to-end analysis of single-cell omics data. The package is composed of several deep generative models for omics data analysis, namely:

* scVI for analysis of single-cell RNA-seq data [Lopez18]_
* scANVI for cell annotation of scRNA-seq data using semi-labeled examples [Xu19]_
* totalVI for analysis of CITE-seq data [GayosoSteier20]_
* gimVI for imputation of missing genes in spatial transcriptomics from scRNA-seq data [Lopez19]_
* AutoZI for assessing gene-specific levels of zero-inflation in scRNA-seq data [Clivio19]_
* LDVAE for an interpretable linear factor model version of scVI [Svensson20]_

These models are able to simultaneously perform many downstream tasks such as learning low-dimensional cell representations, harmonizing datasets from different experiments, and identifying differential expressed features [Boyeau19]_. By levaraging advances in stochastic optimization, these models scale to millions of cells. We invite you to explore these models in our :doc:`tutorials <tutorials/index>`.

* If you find a model useful for your research, please consider citing the corresponding publication.

.. toctree::
   :maxdepth: 1
   :hidden:

   installation
   tutorials/index
   contributed_tutorials/index
   contributing
   history
   references

.. toctree::
   :maxdepth: 1
   :caption: Package Reference
   :hidden:

   scvi.dataset
   scvi.inference
   scvi.models
   scvi.models-modules

Indices and tables
==================
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
