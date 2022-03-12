<img src="https://github.com/scverse/scvi-tools/blob/master/docs/_static/scvi-tools-horizontal.svg?raw=true" width="400" alt="scvi-tools">

[![Stars](https://img.shields.io/github/stars/scverse/scvi-tools?logo=GitHub&color=yellow)](https://github.com/YosefLab/scvi-tools/stargazers)
[![PyPI](https://img.shields.io/pypi/v/scvi-tools.svg)](https://pypi.org/project/scvi-tools)
[![Documentation Status](https://readthedocs.org/projects/scvi/badge/?version=latest)](https://scvi.readthedocs.io/en/stable/?badge=stable)
![Build
Status](https://github.com/scverse/scvi-tools/workflows/scvi-tools/badge.svg)
[![Coverage](https://codecov.io/gh/scverse/scvi-tools/branch/master/graph/badge.svg)](https://codecov.io/gh/YosefLab/scvi-tools)
[![Code
Style](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/python/black)
[![Downloads](https://pepy.tech/badge/scvi-tools)](https://pepy.tech/project/scvi-tools)
[![Join the chat at https://gitter.im/scvi-tools/development](https://badges.gitter.im/scvi-tools/development.svg)](https://gitter.im/scvi-tools/development?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

[scvi-tools](https://scvi-tools.org/) (single-cell variational inference
tools) is a package for probabilistic modeling and analysis of single-cell omics
data, built on top of [PyTorch](https://pytorch.org) and
[AnnData](https://anndata.readthedocs.io/en/latest/).

# Analysis of single-cell omics data

scvi-tools is composed of models that can perform one or more tasks in single-cell omics data analysis. scvi-tools currently hosts implementations of:

-   [scVI](https://rdcu.be/bdHYQ) for analysis of single-cell RNA-seq
    data, as well as its improved differential expression
    [framework](https://www.biorxiv.org/content/biorxiv/early/2019/10/04/794289.full.pdf).
-   [scANVI](https://www.biorxiv.org/content/biorxiv/early/2019/01/29/532895.full.pdf)
    for cell annotation of scRNA-seq data using semi-labeled examples.
-   [totalVI](https://www.biorxiv.org/content/10.1101/2020.05.08.083337v1.full.pdf)
    for analysis of CITE-seq data.
-   [gimVI](https://arxiv.org/pdf/1905.02269.pdf) for imputation of
    missing genes in spatial transcriptomics from scRNA-seq data.
-   [AutoZI](https://www.biorxiv.org/content/biorxiv/early/2019/10/10/794875.full.pdf)
    for assessing gene-specific levels of zero-inflation in scRNA-seq
    data.
-   [LDVAE](https://www.biorxiv.org/content/10.1101/737601v1.full.pdf)
    for an interpretable linear factor model version of scVI.
-   [Stereoscope](https://www.nature.com/articles/s42003-020-01247-y)
    for deconvolution of spatial transcriptomics data.
-   [DestVI](https://www.biorxiv.org/content/10.1101/2021.05.10.443517v1) for multi-resolution deconvolution
    of spatial transcriptomics data.
-   [peakVI](https://www.biorxiv.org/content/10.1101/2021.04.29.442020v1) for analysis of scATAC-seq data.
-   [scArches](https://www.biorxiv.org/content/10.1101/2020.07.16.205997v1)
    for transfer learning from one single-cell atlas to a query dataset
    (currently supports scVI, scANVI and TotalVI).
-   [CellAssign](https://www.nature.com/articles/s41592-019-0529-1) for
    reference-based annotation of scRNA-seq data.
-   [Solo](https://www.sciencedirect.com/science/article/pii/S2405471220301952)
    for doublet detection in scRNA-seq data.

All these implementations have a high-level API that interacts with
[scanpy](http://scanpy.readthedocs.io/), standard save/load functions,
and support GPU acceleration.

# Rapid development of novel probabilistic models

scvi-tools contains the building blocks to develop and deploy novel probablistic
models. These building blocks are powered by popular probabilistic and
machine learning frameworks such as [PyTorch
Lightning](https://www.pytorchlightning.ai/) and
[Pyro](https://pyro.ai/). For an overview of how the scvi-tools package
is structured, you may refer to [this](https://docs.scvi-tools.org/en/stable/user_guide/background/codebase_overview.html) page.

We recommend checking out the [skeleton
repository](https://github.com/YosefLab/scvi-tools-skeleton) as a
starting point for developing and deploying new models with scvi-tools.

# Basic installation

For conda,
```
conda install scvi-tools -c conda-forge
```
and for pip,
```
pip install scvi-tools
```
Please be sure to install a version of [PyTorch](https://pytorch.org/) that is compatible with your GPU (if applicable).

# Resources

-   Tutorials, API reference, and installation guides are available in
    the [documentation](https://docs.scvi-tools.org/).
-   For discussion of usage, check out our
    [forum](https://discourse.scvi-tools.org).
-   Please use the [issues](https://github.com/scverse/scvi-tools/issues) to submit bug reports.
-   If you\'d like to contribute, check out our [contributing
    guide](https://docs.scvi-tools.org/en/stable/contributing/index.html).
-   If you find a model useful for your research, please consider citing
    the corresponding publication (linked above).

# Reference

If you used scvi-tools in your research, please consider citing

```
@article{Gayoso2022,
         author={Gayoso, Adam and Lopez, Romain and Xing, Galen and Boyeau, Pierre and Valiollah Pour Amiri, Valeh and Hong, Justin and Wu, Katherine and Jayasuriya, Michael and   Mehlman, Edouard and Langevin, Maxime and Liu, Yining and Samaran, Jules and Misrachi, Gabriel and Nazaret, Achille and Clivio, Oscar and Xu, Chenling and Ashuach, Tal and Gabitto, Mariano and Lotfollahi, Mohammad and Svensson, Valentine and da Veiga Beltrame, Eduardo and Kleshchevnikov, Vitalii and Talavera-L{\'o}pez, Carlos and Pachter, Lior and Theis, Fabian J. and Streets, Aaron and Jordan, Michael I. and Regier, Jeffrey and Yosef, Nir},
         title={A Python library for probabilistic analysis of single-cell omics data},
         journal={Nature Biotechnology},
         year={2022},
         month={Feb},
         day={07},
         issn={1546-1696},
         doi={10.1038/s41587-021-01206-w},
         url={https://doi.org/10.1038/s41587-021-01206-w}
}
```
along with the publicaton describing the model used. 

