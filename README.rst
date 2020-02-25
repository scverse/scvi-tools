====
scVI
====

|PyPI| |bioconda| |Docs| |Build Status| |Coverage| |Code Style| |Downloads|

.. |PyPI| image:: https://img.shields.io/pypi/v/scVI.svg
    :target: https://pypi.org/project/scvi
.. |bioconda| image:: https://img.shields.io/badge/bioconda-blue.svg
    :target: http://bioconda.github.io/recipes/scvi/README.html
.. |Docs| image:: https://readthedocs.org/projects/scvi/badge/?version=latest
    :target: https://scvi.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status
.. |Build Status| image:: https://travis-ci.org/YosefLab/scVI.svg?branch=master
    :target: https://travis-ci.org/YosefLab/scVI
.. |Coverage| image:: https://codecov.io/gh/YosefLab/scVI/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/YosefLab/scVI
.. |Code Style| image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/python/black
.. |Downloads| image:: https://pepy.tech/badge/scvi
   :target: https://pepy.tech/project/scvi

Single-cell Variational Inference

* Free software: MIT license
* Documentation: https://scvi.readthedocs.io.

scVI is a package for end-to-end analysis of single-cell omics data. It simultaneously performs preprocessing, harmonization, model fitting and a certain number of downstream tasks. The package is composed of several deep generative models for omics data analysis, namely:

    * scVI [Lopez18]_ for analysis of single-cell RNA-seq data,
    * scANVI [Xu19]_ for cell annotation of scRNA-seq data using semi-labeled examples,
    * totalVI [Gayoso19]_ for analysis of CITE-seq data,
    * gimVI [Lopez19]_ for imputation of missing genes in spatial transcriptomics from scRNA-seq data,
    * AutoZI [Clivio19]_ that assesses gene-specific level of zero-inflation in scRNA-seq data.


Follow along with our Jupyter notebooks to quickly get familiar with scVI!

   a. Getting started:
       * `data loading`__
       * `basic usage (scVI)`__
   b. Analyzing several datasets:
       * `harmonization (scVI)`__
       * `annotation (scANVI)`__
   c. Advanced topics:
       * `interaction with scanpy`__
       * `linear decoder for gene interpretation (LDVAE)`__
       * `imputation of unobserved gene expression (gimVI)`__
       * `automated hyperparameter search`__
       * `joint model for CITE-seq data (totalVI)`__
       * `detection of zero-inflated genes (AutoZI)`__


.. __: https://nbviewer.jupyter.org/github/YosefLab/scVI/blob/master/tests/notebooks/data_loading.ipynb
.. __: https://nbviewer.jupyter.org/github/YosefLab/scVI/blob/master/tests/notebooks/basic_tutorial.ipynb
.. __: https://nbviewer.jupyter.org/github/YosefLab/scVI/blob/master/tests/notebooks/harmonization.ipynb
.. __: https://nbviewer.jupyter.org/github/YosefLab/scVI/blob/master/tests/notebooks/annotation.ipynb
.. __: https://nbviewer.jupyter.org/github/YosefLab/scVI/blob/master/tests/notebooks/scanpy_pbmc3k.ipynb
.. __: https://nbviewer.jupyter.org/github/YosefLab/scVI/blob/master/tests/notebooks/linear_decoder.ipynb
.. __: https://nbviewer.jupyter.org/github/YosefLab/scVI/blob/master/tests/notebooks/gimvi_tutorial.ipynb
.. __: https://nbviewer.jupyter.org/github/YosefLab/scVI/blob/master/tests/notebooks/autotune_advanced_notebook.ipynb
.. __: https://nbviewer.jupyter.org/github/YosefLab/scVI/blob/master/tests/notebooks/totalVI.ipynb
.. __: https://nbviewer.jupyter.org/github/YosefLab/scVI/blob/master/tests/notebooks/AutoZI_tutorial.ipynb

