# scvi History

The scvi-tools package used to be scvi. This page commemorates all the hard work on the scvi package by our numerous contributors.

## Contributors

- [@romain]
- [@adam]
- [@eddie]
- [@jeff]
- [@pierre]
- [@max]
- [@yining]
- [@gabriel]
- [@achille]
- [@chenling]
- [@jules]
- [@david-kelley]
- [@william-yang]
- [@oscar]
- [@casey-greene]
- [@jamie-morton]
- [@valentine-svensson]
- [@stephen-flemming]
- [@michael-raevsky]
- [@james-webber]
- [@galen]
- [@francesco-brundu]
- [@primoz-godec]
- [@eduardo-beltrame]
- [@john-reid]
- [@han-yuan]
- [@gokcen-eraslan]

## 0.6.7 (2020-8-05)

- downgrade anndata>=0.7 and scanpy>=1.4.6 [@galen]
- make loompy optional, raise sckmisc import error [@adam]
- fix PBMCDataset download bug [@galen]
- fix AnnDatasetFromAnnData \_X in adata.obs bug [@galen]

## 0.6.6 (2020-7-08)

- add tqdm to within cluster DE genes [@adam]
- restore tqdm to use simple bar instead of ipywidget [@adam]
- move to numpydoc for doctstrings [@adam]
- update issues templates [@adam]
- Poisson variable gene selection [@valentine-svensson]
- BrainSmallDataset set defualt save_path_10X [@gokcen-eraslan]
- train_size must be float between 0.0 and 1.0 [@galen]
- bump dependency versions [@galen]
- remove reproducibility notebook [@galen]
- fix scanVI dataloading [@pierre]

## 0.6.5 (2020-5-10)

- updates to totalVI posterior functions and notebooks [@adam]
- update seurat v3 HVG selection now using skmisc loess  [@adam]

## 0.6.4 (2020-4-14)

- add back Python 3.6 support [@adam]
- get_sample_scale() allows gene selection [@valentine-svensson]
- bug fix to the dataset to anndata method with how cell measurements are stored [@adam]
- fix requirements [@adam]

## 0.6.3 (2020-4-01)

- bug in version for Louvian in setup.py [@adam]

## 0.6.2 (2020-4-01)

- update highly variable gene selection to handle sparse matrices [@adam]
- update DE docstrings [@pierre]
- improve posterior save load to also handle subclasses [@pierre]
- Create NB and ZINB distributions with torch and refactor code accordingly [@pierre]
- typos in autozivae [@achille]
- bug in csc sparse matrices in anndata data loader [@adam]

## 0.6.1 (2020-3-13)

- handles gene and cell attributes with the same name [@han-yuan]
- fixes anndata overwriting when loading [@adam], [@pierre]
- formatting in basic tutorial [@adam]

## 0.6.0 (2020-2-28)

- updates on TotalVI and LDVAE [@adam]
- fix documentation, compatibility and diverse bugs [@adam], [@pierre] [@romain]
- fix for external module on scanpy [@galen]

## 0.5.0 (2019-10-17)

- do not automatically upper case genes [@adam]
- AutoZI [@oscar]
- Made the intro tutorial more user friendly [@adam]
- Tests for LDVAE notebook [@adam]
- black codebase [@achille] [@gabriel] [@adam]
- fix compatibility issues with sklearn and numba [@romain]
- fix Anndata [@francesco-brundu]
- docstring, totalVI, totalVI notebook and CITE-seq data [@adam]
- fix type [@eduardo-beltrame]
- fixing installation guide [@jeff]
- improved error message for dispersion [@stephen-flemming]

## 0.4.1 (2019-08-03)

- docstring [@achille]
- differential expression [@oscar] [@pierre]

## 0.4.0 (2019-07-25)

- gimVI [@achille]
- synthetic correlated datasets, fixed bug in marginal log likelihood [@oscar]
- autotune, dataset enhancements [@gabriel]
- documentation [@jeff]
- more consistent posterior API, docstring, validation set [@adam]
- fix anndataset [@michael-raevsky]
- linearly decoded VAE [@valentine-svensson]
- support for scanpy, fixed bugs, dataset enhancements [@achille]
- fix filtering bug, synthetic correlated datasets, docstring, differential expression [@pierre]
- better docstring [@jamie-morton]
- classifier based on library size for doublet detection [@david-kelley]

## 0.3.0 (2019-05-03)

- corrected notebook [@jules]
- added UMAP and updated harmonization code [@chenling] [@romain]
- support for batch indices in csvdataset [@primoz-godec]
- speeding up likelihood computations [@william-yang]
- better anndata interop [@casey-greene]
- early stopping based on classifier accuracy [@david-kelley]

## 0.2.4 (2018-12-20)

- updated to torch v1 [@jules]
- added stress tests for harmonization [@chenling]
- fixed autograd breaking [@romain]
- make removal of empty cells more efficient [@john-reid]
- switch to os.path.join [@casey-greene]

## 0.2.2 (2018-11-08)

- added baselines and datasets for sMFISH imputation [@jules]
- added harmonization content [@chenling]
- fixing bugs on DE [@romain]

## 0.2.0 (2018-09-04)

- annotation notebook [@eddie]
- Memory footprint management [@jeff]
- updated early stopping [@max]
- docstring [@james-webber]

## 0.1.6 (2018-08-08)

- MMD and adversarial inference wrapper [@eddie]
- Documentation [@jeff]
- smFISH data imputation [@max]

## 0.1.5 (2018-07-24)

- Dataset additions [@eddie]
- Documentation [@yining]
- updated early stopping [@max]

## 0.1.3 (2018-06-22)

- Notebook enhancement [@yining]
- Semi-supervision [@eddie]

## 0.1.2 (2018-06-13)

- First release on PyPi
- Skeleton code & dependencies [@jeff]
- Unit tests [@max]
- PyTorch implementation of scVI [@eddie] [@max]
- Dataset preprocessing [@eddie] [@max] [@yining]

## 0.1.0 (2017-09-05)

- First scVI TensorFlow version [@romain]

[@achille]: https://github.com/ANazaret
[@adam]: https://github.com/adamgayoso
[@casey-greene]: https://github.com/cgreene
[@chenling]: https://github.com/chenlingantelope
[@david-kelley]: https://github.com/davek44
[@eddie]: https://github.com/Edouard360
[@eduardo-beltrame]: https://github.com/Munfred
[@francesco-brundu]: https://github.com/fbrundu
[@gabriel]: https://github.com/gabmis
[@galen]: https://github.com/galenxing
[@gokcen-eraslan]: https://github.com/gokceneraslan
[@han-yuan]: https://github.com/hy395
[@james-webber]: https://github.com/jamestwebber
[@jamie-morton]: https://github.com/mortonjt
[@jeff]: https://github.com/jeff-regier
[@john-reid]: https://github.com/JohnReid
[@jules]: https://github.com/jules-samaran
[@max]: https://github.com/maxime1310
[@michael-raevsky]: https://github.com/raevskymichail
[@oscar]: https://github.com/oscarclivio
[@pierre]: https://github.com/PierreBoyeau
[@primoz-godec]: https://github.com/PrimozGodec
[@romain]: https://github.com/romain-lopez
[@stephen-flemming]: https://github.com/sjfleming
[@valentine-svensson]: https://github.com/vals
[@william-yang]: https://github.com/triyangle
[@yining]: https://github.com/imyiningliu
