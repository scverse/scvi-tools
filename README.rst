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


Quick Start
-----------

1. Install Python 3.7. We typically use the Miniconda_ Python distribution and Linux.

.. _Miniconda: https://conda.io/miniconda.html

2. Install PyTorch_. If you have an Nvidia GPU, be sure to install a version of PyTorch that supports it -- scVI runs much faster with a discrete GPU.

.. _PyTorch: http://pytorch.org

3. Install scVI through conda:

    ``conda install scvi -c bioconda -c conda-forge``

   Alternatively, you may try pip (``pip install scvi``), or you may clone this repository and run ``python setup.py install``.

4. If you wish to use multiple GPUs for hyperparameter tuning, install MongoDb_.

.. _MongoDb: https://docs.mongodb.com/manual/installation/

5. Follow along with our Jupyter notebooks to quickly get familiar with scVI!

   a. Getting started:
       * `data loading`__
       * `basic usage`__
   b. Analyzing several datasets:
       * `harmonization (scVI)`__
       * `annotation (scANVI)`__
   c. Advanced topics:
       * `interaction with scanpy`__
       * `linear decoder for gene interpretation`__
       * `imputation of unobserved gene expression (gimVI)`__
       * `automated hyperparameter search`__


.. __: https://nbviewer.jupyter.org/github/YosefLab/scVI/blob/master/tests/notebooks/data_loading.ipynb
.. __: https://nbviewer.jupyter.org/github/YosefLab/scVI/blob/master/tests/notebooks/basic_tutorial.ipynb
.. __: https://nbviewer.jupyter.org/github/YosefLab/scVI/blob/master/tests/notebooks/harmonization.ipynb
.. __: https://nbviewer.jupyter.org/github/YosefLab/scVI/blob/master/tests/notebooks/annotation.ipynb
.. __: https://nbviewer.jupyter.org/github/YosefLab/scVI/blob/master/tests/notebooks/scanpy_pbmc3k.ipynb
.. __: https://nbviewer.jupyter.org/github/YosefLab/scVI/blob/master/tests/notebooks/linear_decoder.ipynb
.. __: https://nbviewer.jupyter.org/github/YosefLab/scVI/blob/master/tests/notebooks/gimvi_tutorial.ipynb
.. __: https://nbviewer.jupyter.org/github/YosefLab/scVI/blob/master/tests/notebooks/autotune_advanced_notebook.ipynb

References
----------

Romain Lopez, Jeffrey Regier, Michael Cole, Michael I. Jordan, Nir Yosef.
**"Deep generative modeling for single-cell transcriptomics."**
Nature Methods, 2018. `[pdf]`__

.. __: https://rdcu.be/bdHYQ

Chenling Xu∗, Romain Lopez∗, Edouard Mehlman∗, Jeffrey Regier, Michael I. Jordan, Nir Yosef.
**"Harmonization and Annotation of Single-cell Transcriptomics data with Deep Generative Models."**
Submitted, 2019. `[pdf]`__

.. __: https://www.biorxiv.org/content/biorxiv/early/2019/01/29/532895.full.pdf

Romain Lopez∗, Achille Nazaret∗, Maxime Langevin*, Jules Samaran*, Jeffrey Regier*, Michael I. Jordan, Nir Yosef.
**"A joint model of unpaired data from scRNA-seq and spatial transcriptomics for imputing missing gene expression measurements."**
ICML Workshop on Computational Biology, 2019. `[pdf]`__

.. __: https://arxiv.org/pdf/1905.02269.pdf

