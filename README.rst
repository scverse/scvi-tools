========================================
scVI - Single cell Variational Inference
========================================

|Stars| |PyPI| |BioConda| |Docs| |Build Status| |Coverage| |Code Style| |Downloads|

.. |Stars| image:: https://img.shields.io/github/stars/YosefLab/scVI?logo=GitHub&color=yellow
   :target: https://github.com/YosefLab/scVI/stargazers
.. |PyPI| image:: https://img.shields.io/pypi/v/scVI.svg
    :target: https://pypi.org/project/scvi
.. |BioConda| image:: https://img.shields.io/conda/vn/bioconda/scVI
   :target: https://bioconda.github.io/recipes/scvi/README.html
.. |Docs| image:: https://readthedocs.org/projects/scvi/badge/?version=latest
    :target: https://scvi.readthedocs.io/en/stable/?badge=stable
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

* scVI_ for analysis of single-cell RNA-seq data, as well as its improved differential expression framework_
* scANVI_ for cell annotation of scRNA-seq data using semi-labeled examples
* totalVI_ for analysis of CITE-seq data
* gimVI_ for imputation of missing genes in spatial transcriptomics from scRNA-seq data
* AutoZI_ for assessing gene-specific levels of zero-inflation in scRNA-seq data
* LDVAE_ for an interpretable linear factor model version of scVI

Tutorials and API reference are available in the documentation_.
Please use the issues here to discuss usage, or submit bug reports.
If you'd like to contribute, please check out our `contributing guide`_.
If you find a model useful for your research, please consider citing the corresponding publication (linked above).

.. _documentation: https://scvi.readthedocs.io
.. _`contributing guide`: https://scvi.readthedocs.io/en/update_docs/contributing.html
.. _scVI: https://rdcu.be/bdHYQ
.. _scANVI: https://www.biorxiv.org/content/biorxiv/early/2019/01/29/532895.full.pdf
.. _totalVI: https://www.biorxiv.org/content/10.1101/2020.05.08.083337v1.full.pdf
.. _AutoZI: https://www.biorxiv.org/content/biorxiv/early/2019/10/10/794875.full.pdf
.. _LDVAE: https://www.biorxiv.org/content/10.1101/737601v1.full.pdf
.. _gimVI: https://arxiv.org/pdf/1905.02269.pdf
.. _framework: https://www.biorxiv.org/content/biorxiv/early/2019/10/04/794289.full.pdf

