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
--------

1. Install Python 3.6 or later. We typically use the Miniconda Python distribution: https://conda.io/miniconda.html

2. Install PyTorch: http://pytorch.org

3. Install scvi through either conda or pip, or by cloning this repository and manually installing all the requirements listed in setup.py_.

.. _setup.py: https://github.com/YosefLab/scVI/tree/master/setup.py

    Through conda: ``conda install scvi -c bioconda``

    Through pip: ``pip install scvi``

    By cloning this repository: ``git clone https://github.com/YosefLab/scVI.git``

4. Refer to the `this Jupyter notebook`__ to see how to import dataset into scvi.

.. __: https://github.com/YosefLab/scVI/tree/master/examples/scVI-data-loading.ipynb

5. Refer to `this Jupyter notebook`__ to see how to train the scVI model, impute missing data, detect differential expression, and more!

.. __: https://github.com/YosefLab/scVI/tree/master/examples/scVI-dev.ipynb


Benchmarks
--------

To recreate the results appearing in the paper referenced below, run

.. code-block::

    python ./run_benchmarks.py --dataset=cortex 

References
--------

Romain Lopez, Jeffrey Regier, Michael B Cole, Michael Jordan, Nir Yosef.
**"Bayesian Inference for a Generative Model of Transcriptome Profiles from Single-cell RNA Sequencing."**
In submission. Preprint available at https://www.biorxiv.org/content/early/2018/03/30/292037
