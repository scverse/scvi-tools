Installation
-----------

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
