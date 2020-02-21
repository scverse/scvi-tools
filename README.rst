====
scVI
====

|PyPI| |bioconda| |Open In Colab| |Docs| |Build Status| |Coverage| |Code Style| |Downloads|

.. |PyPI| image:: https://img.shields.io/pypi/v/scVI.svg
    :target: https://pypi.org/project/scvi
.. |bioconda| image:: https://img.shields.io/badge/bioconda-blue.svg
    :target: http://bioconda.github.io/recipes/scvi/README.html
.. |Open In Colab| image:: https://colab.research.google.com/assets/colab-badge.svg
    :target: https://colab.research.google.com/github/yoseflab/scVI/blob/master
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

scVI is a deep generative model for single-cell RNA that simultaneously embed preprocessing, harmonization, data fitting, downstream application steps in a single framework.

Several projects have extended the scVI model to integrate additional data knowledge or to answer specific biological questions.
These models, included in this package are:
    * scANVI for cell annotation using labeled examples.
    * totalVI for CITE-seq data analysis.
    * gimVI, used to integrate spatial transcriptomics.
    * AutoZI, that quantify genes' levels of zero-inflation


All these models also benefit from integrated hyperparameter autotuning features provided in this package.

* Free software: MIT license
* Documentation: https://scvi.readthedocs.io.

YosefLab's blog can be found here_.

.. _here: https://yoseflab.github.io/

Installation
---------------------------------

1. Install Python 3.7. We typically use the Miniconda_ Python distribution and Linux.

.. _Miniconda: https://conda.io/miniconda.html

2. Install PyTorch_. If you have an Nvidia GPU, be sure to install a version of PyTorch that supports it -- scVI runs much faster with a discrete GPU.

.. _PyTorch: http://pytorch.org

3. Install scVI in one of the following ways:

    1. Through conda ``conda install scvi -c bioconda -c conda-forge``
    2. Through pip ``pip install scvi``
    3. Through pip with packages to run notebooks ``pip install scvi[notebooks]``
    4. Nightly version - clone this repo and run ``pip install .``
    5. For development - clone this repo and run ``pip install -e .[test,notebooks]``

4. If you wish to use multiple GPUs for hyperparameter tuning, install MongoDb_.

.. _MongoDb: https://docs.mongodb.com/manual/installation/


Notebooks
-------------------------------------------------------

Follow along with our Jupyter notebooks to quickly get familiar with scVI 
(also available on Colab_)!

.. _Colab: https://colab.research.google.com/github/yoseflab/scVI/blob/master

To run these notebooks, we recommend installing scVI extra dependencies using the command

``pip install scvi[notebooks]``

Getting started
""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    * `data loading`__: presents how to load and manipulate built-in and custom datasets in different formats.
    * `basic usage (scVI)`__: highlights scVI's main features on mouse cortex cells.

Analyzing several datasets
""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * `harmonization (scVI)`__: illustrates how scVI can harmonize multiple datasets for clustering and differential expression applications on two human PBMC datasets.
    * `annotation (scANVI)`__: shows how to transfer the cell state knowledge to label an unlabelled dataset.

Advanced topics
""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    * `interaction with scanpy`__: explores the Scanpy PBMC 3K dataset with scVI and Scanpy.
    * `linear decoder for gene interpretation (LDVAE)`__: introduces a flavor of scVI that employs a linear decoder, which allows a better interpretability of latent dimensions.
    * `imputation of unobserved gene expression (gimVI)`__: presents gimVI for missing genes imputation in Spatial Data from Sequencing Data. 
    * `automated hyperparameter search`__: introduces the hyperparameter tuning functionalities of scVI.
    * `joint model for CITE-seq data (totalVI)`__: introduces totalVI for end-to-end CITE-seq data analysis.
    * `detection of zero-inflated genes (AutoZI)`__: shows how AutoZI, derived from scVI, can provide gene-specific zero-inflation analysis. 


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

Adam Gayoso, Romain Lopez, Zoë Steier, Jeffrey Regier, Aaron Streets, Nir Yosef.
**"A joint model of RNA expression and surface protein abundance in single cells."**
Machine Learning in Computational Biology (MLCB), 2019. `[pdf]`__

.. __: https://www.biorxiv.org/content/biorxiv/early/2019/10/07/791947.full.pdf

Oscar Clivio, Romain Lopez, Jeffrey Regier, Adam Gayoso, Michael I. Jordan, Nir Yosef.
**"Detecting zero-inflated genes in single-cell transcriptomics data."**
Machine Learning in Computational Biology (MLCB), 2019. `[pdf]`__

.. __: https://www.biorxiv.org/content/biorxiv/early/2019/10/10/794875.full.pdf

Pierre Boyeau, Romain Lopez, Jeffrey Regier, Adam Gayoso, Michael I. Jordan, Nir Yosef.
**"Deep generative models for detecting differential expression in single cells."**
Machine Learning in Computational Biology (MLCB), 2019. `[pdf]`__

.. __: https://www.biorxiv.org/content/biorxiv/early/2019/10/04/794289.full.pdf

Valentine Svensson, Lior Pachter.
**"Interpretable factor models of single-cell RNA-seq via variational autoencoders."**
bioRxiv, 2019. `[pdf]`__

.. __: https://www.biorxiv.org/content/10.1101/737601v1.full.pdf
