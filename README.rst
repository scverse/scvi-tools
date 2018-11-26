====
scVI
====

.. image:: https://travis-ci.org/YosefLab/scVI.svg?branch=master
    :target: https://travis-ci.org/YosefLab/scVI

.. image:: https://codecov.io/gh/YosefLab/scVI/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/YosefLab/scVI

.. image:: https://readthedocs.org/projects/scvi/badge/?version=latest
        :target: https://scvi.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

Single-cell Variational Inference

* Free software: MIT license
* Documentation: https://scvi.readthedocs.io.


Quick Start
-----------

1. Install Python 3.6 or later. We typically use the Miniconda_ Python distribution.

.. _Miniconda: https://conda.io/miniconda.html

2. Install PyTorch_. If you have an Nvidia GPU, be sure to install a version of PyTorch that supports it -- scVI runs much faster with a discrete GPU.

.. _PyTorch: http://pytorch.org

3. Install scVI through conda (``conda install scvi -c bioconda``) or through pip (``pip install scvi``). Alternatively, you may download or clone this repository and run ``python setup.py install``.

4. Follow along with our Jupyter notebooks to quickly get familiar with scVI!

   a. `data loading`__
   b. `basic usage`__ 
   c. `reproducing results from the paper`__ 

.. __: https://github.com/YosefLab/scVI/tree/master/tests/notebooks/data_loading.ipynb
.. __: https://github.com/YosefLab/scVI/tree/master/tests/notebooks/basic_tutorial.ipynb
.. __: https://github.com/YosefLab/scVI/blob/master/tests/notebooks/scVI_reproducibility.ipynb



References
----------

Romain Lopez, Jeffrey Regier, Michael Cole, Michael I. Jordan, Nir Yosef.
**"Deep generative modeling for single-cell transcriptomics"**
Nature Methods, in press (accepted Oct 26, 2018). 
Preprint available at https://www.biorxiv.org/content/early/2018/03/30/292037
