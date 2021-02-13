.. raw:: html

    <img src="https://github.com/YosefLab/scvi-tools/blob/master/docs/_static/scvi-tools-horizontal.svg?raw=true" alt="scvi-tools" width="400px">

|Stars| |PyPI| |Docs| |Build Status| |Coverage| |Code Style| |Downloads| |Gitter|

.. |scvi-tools| image:: https://github.com/YosefLab/scvi-tools/blob/master/docs/_static/scvi-tools-horizontal.svg?raw=true
  :width: 50
  :alt: scvi-tools
.. |PyPI| image:: https://img.shields.io/pypi/v/scvi-tools.svg
    :target: https://pypi.org/project/scvi-tools
.. |BioConda| image:: https://img.shields.io/conda/vn/bioconda/scvi-tools
   :target: https://bioconda.github.io/recipes/scvi-tools/README.html
.. |Stars| image:: https://img.shields.io/github/stars/YosefLab/scvi-tools?logo=GitHub&color=yellow
   :target: https://github.com/YosefLab/scvi-tools/stargazers
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
.. |Gitter| image:: https://badges.gitter.im/scvi-tools/development.svg
   :alt: Join the chat at https://gitter.im/scvi-tools/development
   :target: https://gitter.im/scvi-tools/development?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge

scvi-tools_ (single-cell variational inference tools) is a package for probabilistic modeling of single-cell omics data, built on top of `PyTorch
<https://pytorch.org>`_ and `Anndata <https://anndata.readthedocs.io/en/latest/>`_.

Available implementations of single-cell omics models
-----------------------------------------------------

scvi-tools contains scalable implementations of several models that facilitate a broad number of tasks across many omics, including:

* scVI_ for analysis of single-cell RNA-seq data, as well as its improved differential expression framework_.
* scANVI_ for cell annotation of scRNA-seq data using semi-labeled examples.
* totalVI_ for analysis of CITE-seq data.
* gimVI_ for imputation of missing genes in spatial transcriptomics from scRNA-seq data.
* AutoZI_ for assessing gene-specific levels of zero-inflation in scRNA-seq data.
* LDVAE_ for an interpretable linear factor model version of scVI.
* Stereoscope_ for deconvolution of spatial transcriptomics data.
* peakVI for analysis of ATAC-seq data.
* scArches_ for transfer learning from one single-cell atlas to a query dataset (currently supports scVI, scANVI and TotalVI).

All these implementations have a high-level API that interacts with `scanpy <http://scanpy.readthedocs.io/>`_, standard save/load functions, and support GPU acceleration.

Fast prototyping of novel probabilistic models
----------------------------------------------

scvi-tools contains the building blocks to prototype novel probablistic models. These building blocks are powered by popular probabilistic and machine learning frameworks such as `PyTorch lightning <https://www.pytorchlightning.ai/>`_, and `Pyro <https://pyro.ai/>`_.

We recommend checking out the `skeleton repository <https://github.com/YosefLab/scvi-tools-skeleton>`_, as a starting point for developing new models into scvi-tools.

Resources
----------

* Tutorials, API reference, and installation guides are available in the documentation_.
* For discussion of usage, checkout out our `forum`_.
* Please use the issues here to submit bug reports.
* If you'd like to contribute, check out our `development guide`_.
* If you find a model useful for your research, please consider citing the corresponding publication (linked above).

.. _scvi-tools: https://scvi-tools.org/
.. _documentation: https://scvi-tools.org/
.. _`development guide`: https://scvi-tools.org/en/stable/development.html
.. _scVI: https://rdcu.be/bdHYQ
.. _scANVI: https://www.biorxiv.org/content/biorxiv/early/2019/01/29/532895.full.pdf
.. _totalVI: https://www.biorxiv.org/content/10.1101/2020.05.08.083337v1.full.pdf
.. _AutoZI: https://www.biorxiv.org/content/biorxiv/early/2019/10/10/794875.full.pdf
.. _LDVAE: https://www.biorxiv.org/content/10.1101/737601v1.full.pdf
.. _gimVI: https://arxiv.org/pdf/1905.02269.pdf
.. _Stereoscope: https://www.nature.com/articles/s42003-020-01247-y
.. _scArches: https://www.biorxiv.org/content/10.1101/2020.07.16.205997v1
.. _framework: https://www.biorxiv.org/content/biorxiv/early/2019/10/04/794289.full.pdf
.. _forum: https://discourse.scvi-tools.org
