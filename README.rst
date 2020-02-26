====
scVI -- Single cell Variational Inference
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

scVI is a package for end-to-end analysis of single-cell omics data. The package is composed of several deep generative models for omics data analysis, namely:

* scVI for analysis of single-cell RNA-seq data (`Nature Methods (2018)`_)
* scANVI for cell annotation of scRNA-seq data using semi-labeled examples (`bioRxiv`_)
* totalVI for analysis of CITE-seq data
* gimVI for imputation of missing genes in spatial transcriptomics from scRNA-seq data
* AutoZI for assessing gene-specific levels of zero-inflation in scRNA-seq data
* LDVAE for an interpretable linear factor model version of scVI

Tutorials and API reference are available in the documentation_.
If you'd like to contribute by opening an issue or creating a pull request,
please take a look at our `contributing guide`_.
If you find a model useful for your research, please consider citing the corresponding publication.

.. _Nature Methods (2018): https://www.nature.com/articles/s41592-018-0229-2.epdf?author_access_token=5sMbnZl1iBFitATlpKkddtRgN0jAjWel9jnR3ZoTv0P1-tTjoP-mBfrGiMqpQx63aBtxToJssRfpqQ482otMbBw2GIGGeinWV4cULBLPg4L4DpCg92dEtoMaB1crCRDG7DgtNrM_1j17VfvHfoy1cQ%3D%3D
.. _bioRxiv: https://www.biorxiv.org/content/biorxiv/early/2019/01/29/532895.full.pdf
.. _documentation: https://scvi.readthedocs.io
.. _contributing guide: CONTRIBUTING.rst




