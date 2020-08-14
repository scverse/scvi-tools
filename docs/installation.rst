Installation
------------

Prerequisites
~~~~~~~~~~~~~~
1. Install Python 3.7. We typically use the Miniconda_ Python distribution and Linux.

.. _Miniconda: https://conda.io/miniconda.html

2. Install PyTorch_. If you have an Nvidia GPU, be sure to install a version of PyTorch that supports it -- scvi runs much faster with a discrete GPU.

.. _PyTorch: http://pytorch.org


scvi installation
~~~~~~~~~~~~~~~~~~

Install scvi in one of the following ways:

Through conda::

    conda install scvi -c bioconda -c conda-forge

Through pip::

    pip install scvi

Through pip with packages to run notebooks. This installs scanpy, etc.::

    pip install scvi[notebooks]

Nightly version - clone this repo and run::

    pip install .

For development - clone this repo and run::

    pip install -e .[test,notebooks]

If you wish to use multiple GPUs for hyperparameter tuning, install MongoDb_.

.. _MongoDb: https://docs.mongodb.com/manual/installation/
