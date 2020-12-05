==========
scvi-tools
==========

|Stars| |PyPI| |BioConda| |Docs| |Build Status| |Coverage| |Code Style| |Downloads|

.. |Stars| image:: https://img.shields.io/github/stars/YosefLab/scvi-tools?logo=GitHub&color=yellow
   :target: https://github.com/YosefLab/scvi-tools/stargazers
.. |PyPI| image:: https://img.shields.io/pypi/v/scvi-tools.svg
    :target: https://pypi.org/project/scvi-tools
.. |BioConda| image:: https://img.shields.io/conda/vn/bioconda/scvi-tools
   :target: https://bioconda.github.io/recipes/scvi-tools/README.html
.. |Docs| image:: https://readthedocs.org/projects/scvi/badge/?version=latest
    :target: https://scvi.readthedocs.io/en/stable/?badge=stable
    :alt: Documentation Status
.. |Build Status| image:: https://github.com/YosefLab/scvi-tools/workflows/scvi-tools/badge.svg
.. |Coverage| image:: https://codecov.io/gh/YosefLab/scvi-tools/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/YosefLab/scvi-tools
.. |Code Style| image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/python/black
.. |Downloads| image:: https://pepy.tech/badge/scvi-tools
   :target: https://pepy.tech/project/scvi-tools

scvi-tools (single-cell variational inference tools) is a package for end-to-end analysis of single-cell omics data. The package is composed of several deep generative models for omics data analysis, namely:

* scVI_ for analysis of single-cell RNA-seq data, as well as its improved differential expression framework_
* scANVI_ for cell annotation of scRNA-seq data using semi-labeled examples
* totalVI_ for analysis of CITE-seq data
* gimVI_ for imputation of missing genes in spatial transcriptomics from scRNA-seq data
* AutoZI_ for assessing gene-specific levels of zero-inflation in scRNA-seq data
* LDVAE_ for an interpretable linear factor model version of scVI

Tutorials, API reference, and installation guides are available in the documentation_.
Please use the issues here to submit bug reports.
For discussion of usage, checkout out our `forum`_.
If you'd like to contribute, check out our `development guide`_.
If you find a model useful for your research, please consider citing the corresponding publication (linked above).

Package transition
------------------

``scvi`` is now ``scvi-tools``. If you're looking for ``scvi`` source code, please use the branch tags. ``scvi`` is still available via pypi and bioconda, but there will be no future releases past ``0.6.8``. Documentation for ``scvi`` is available at the same documentation website.

.. _documentation: https://scvi-tools.org/
.. _`development guide`: https://scvi-tools.org/en/stable/development.html
.. _scVI: https://rdcu.be/bdHYQ
.. _scANVI: https://www.biorxiv.org/content/biorxiv/early/2019/01/29/532895.full.pdf
.. _totalVI: https://www.biorxiv.org/content/10.1101/2020.05.08.083337v1.full.pdf
.. _AutoZI: https://www.biorxiv.org/content/biorxiv/early/2019/10/10/794875.full.pdf
.. _LDVAE: https://www.biorxiv.org/content/10.1101/737601v1.full.pdf
.. _gimVI: https://arxiv.org/pdf/1905.02269.pdf
.. _framework: https://www.biorxiv.org/content/biorxiv/early/2019/10/04/794289.full.pdf
.. _forum: https://discourse.scvi-tools.org
