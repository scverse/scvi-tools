=============
Release notes
=============

0.7.0 (2020-10-14)
-----------------
scvi is now scvi-tools. Version 0.7 introduces many breaking changes. The best way to learn how to use scvi-tools is with our documentation and tutorials.

* New high-level API and data loading, please see tutorials and examples for usage.
* `GeneExpressionDataset` and associated classes have been removed.
* Built-in datasets now return `AnnData` objects.
* `scvi-tools` now relies entirely on the AnnData_ format.
* `scvi.models` has been moved to `scvi.core.modules`.
* `Posterior` classes have been reduced to wrappers on `DataLoaders`
* `scvi.inference` has been split to `scvi.core.data_loaders` for `ScviDataLoader` classes and `scvi.core.trainers` for trainer classes.
* Usage of classes like `Trainer` and `ScviDataLoader` now require the `AnnData` data object as input.


0.6.7 (2020-8-05)
-----------------
* downgrade anndata>=0.7 and scanpy>=1.4.6 `@galen`_
* make loompy optional, raise sckmisc import error `@adam`_
* fix PBMCDataset download bug `@galen`_
* fix AnnDatasetFromAnnData _X in adata.obs bug `@galen`_

0.6.6 (2020-7-08)
-----------------
* add tqdm to within cluster DE genes `@adam`_
* restore tqdm to use simple bar instead of ipywidget `@adam`_
* move to numpydoc for doctstrings `@adam`_
* update issues templates `@adam`_
* Poisson variable gene selection `@valentine-svensson`_
* BrainSmallDataset set defualt save_path_10X `@gokcen-eraslan`_
* train_size must be float between 0.0 and 1.0 `@galen`_
* bump dependency versions `@galen`_
* remove reproducibility notebook `@galen`_
* fix scanVI dataloading `@pierre`_

0.6.5 (2020-5-10)
------------------
* updates to totalVI posterior functions and notebooks `@adam`_
* update seurat v3 HVG selection now using skmisc loess  `@adam`_

0.6.4 (2020-4-14)
------------------
* add back Python 3.6 support `@adam`_
* get_sample_scale() allows gene selection `@valentine-svensson`_
* bug fix to the dataset to anndata method with how cell measurements are stored `@adam`_
* fix requirements `@adam`_

0.6.3 (2020-4-01)
------------------
* bug in version for Louvian in setup.py `@adam`_

0.6.2 (2020-4-01)
------------------
* update highly variable gene selection to handle sparse matrices `@adam`_
* update DE docstrings `@pierre`_
* improve posterior save load to also handle subclasses `@pierre`_
* Create NB and ZINB distributions with torch and refactor code accordingly `@pierre`_
* typos in autozivae `@achille`_
* bug in csc sparse matrices in anndata data loader `@adam`_

0.6.1 (2020-3-13)
------------------
* handles gene and cell attributes with the same name `@han-yuan`_
* fixes anndata overwriting when loading `@adam`_, `@pierre`_
* formatting in basic tutorial `@adam`_

0.6.0 (2020-2-28)
------------------
* updates on TotalVI and LDVAE `@adam`_
* fix documentation, compatibility and diverse bugs `@adam`_, `@pierre`_ `@romain`_
* fix for external module on scanpy `@galen`_

0.5.0 (2019-10-17)
------------------
* do not automatically upper case genes `@adam`_
* AutoZI `@oscar`_
* Made the intro tutorial more user friendly `@adam`_
* Tests for LDVAE notebook `@adam`_
* black codebase `@achille`_ `@gabriel`_ `@adam`_
* fix compatibility issues with sklearn and numba `@romain`_
* fix Anndata `@francesco-brundu`_
* docstring, totalVI, totalVI notebook and CITE-seq data `@adam`_
* fix type `@eduardo-beltrame`_
* fixing installation guide `@jeff`_
* improved error message for dispersion `@stephen-flemming`_

0.4.1 (2019-08-03)
------------------

* docstring `@achille`_
* differential expression `@oscar`_ `@pierre`_

0.4.0 (2019-07-25)
------------------

* gimVI `@achille`_
* synthetic correlated datasets, fixed bug in marginal log likelihood `@oscar`_
* autotune, dataset enhancements `@gabriel`_
* documentation `@jeff`_
* more consistent posterior API, docstring, validation set `@adam`_
* fix anndataset `@michael-raevsky`_
* linearly decoded VAE `@valentine-svensson`_
* support for scanpy, fixed bugs, dataset enhancements `@achille`_
* fix filtering bug, synthetic correlated datasets, docstring, differential expression `@pierre`_
* better docstring `@jamie-morton`_
* classifier based on library size for doublet detection `@david-kelley`_

0.3.0 (2019-05-03)
------------------

* corrected notebook `@jules`_
* added UMAP and updated harmonization code `@chenling`_ `@romain`_
* support for batch indices in csvdataset `@primoz-godec`_
* speeding up likelihood computations `@william-yang`_
* better anndata interop `@casey-greene`_
* early stopping based on classifier accuracy `@david-kelley`_

0.2.4 (2018-12-20)
------------------

* updated to torch v1 `@jules`_
* added stress tests for harmonization `@chenling`_
* fixed autograd breaking `@romain`_
* make removal of empty cells more efficient `@john-reid`_
* switch to os.path.join `@casey-greene`_


0.2.2 (2018-11-08)
------------------

* added baselines and datasets for sMFISH imputation `@jules`_
* added harmonization content `@chenling`_
* fixing bugs on DE `@romain`_


0.2.0 (2018-09-04)
------------------

* annotation notebook `@eddie`_
* Memory footprint management `@jeff`_
* updated early stopping `@max`_
* docstring `@james-webber`_

0.1.6 (2018-08-08)
------------------

* MMD and adversarial inference wrapper `@eddie`_
* Documentation `@jeff`_
* smFISH data imputation `@max`_

0.1.5 (2018-07-24)
------------------

* Dataset additions `@eddie`_
* Documentation `@yining`_
* updated early stopping `@max`_


0.1.3 (2018-06-22)
------------------

* Notebook enhancement `@yining`_
* Semi-supervision `@eddie`_

0.1.2 (2018-06-13)
------------------

* First release on PyPi
* Skeleton code & dependencies `@jeff`_
* Unit tests `@max`_
* PyTorch implementation of scVI `@eddie`_ `@max`_
* Dataset preprocessing `@eddie`_ `@max`_ `@yining`_

0.1.0 (2017-09-05)
------------------

* First scVI TensorFlow version `@romain`_

.. _`@romain`: https://github.com/romain-lopez
.. _`@adam`: https://github.com/adamgayoso
.. _`@eddie`: https://github.com/Edouard360
.. _`@jeff`: https://github.com/jeff-regier
.. _`@pierre`: https://github.com/PierreBoyeau
.. _`@max`: https://github.com/maxime1310
.. _`@yining`: https://github.com/imyiningliu
.. _`@gabriel`: https://github.com/gabmis
.. _`@achille`: https://github.com/ANazaret
.. _`@chenling`: https://github.com/chenlingantelope
.. _`@jules`: https://github.com/jules-samaran
.. _`@david-kelley`: https://github.com/davek44
.. _`@william-yang`: https://github.com/triyangle
.. _`@oscar`: https://github.com/oscarclivio
.. _`@casey-greene`: https://github.com/cgreene
.. _`@jamie-morton`: https://github.com/mortonjt
.. _`@valentine-svensson`: https://github.com/vals
.. _`@stephen-flemming`: https://github.com/sjfleming
.. _`@michael-raevsky`: https://github.com/raevskymichail
.. _`@james-webber`: https://github.com/jamestwebber
.. _`@galen`: https://github.com/galenxing
.. _`@francesco-brundu`: https://github.com/fbrundu
.. _`@primoz-godec`: https://github.com/PrimozGodec
.. _`@eduardo-beltrame`: https://github.com/Munfred
.. _`@john-reid`: https://github.com/JohnReid
.. _`@han-yuan`: https://github.com/hy395
.. _`@gokcen-eraslan`: https://github.com/gokceneraslan
.. _AnnData: https://anndata.readthedocs.io/en/stable/
