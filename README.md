<img src="https://github.com/scverse/scvi-tools/blob/main/docs/_static/scvi-tools-horizontal.svg?raw=true" width="400" alt="scvi-tools">

[![Stars](https://img.shields.io/github/stars/scverse/scvi-tools?logo=GitHub&color=yellow)](https://github.com/YosefLab/scvi-tools/stargazers)
[![PyPI](https://img.shields.io/pypi/v/scvi-tools.svg)](https://pypi.org/project/scvi-tools)
[![Documentation Status](https://readthedocs.org/projects/scvi/badge/?version=latest)](https://scvi.readthedocs.io/en/stable/?badge=stable)
![Build
Status](https://github.com/scverse/scvi-tools/workflows/scvi-tools/badge.svg)
[![Coverage](https://codecov.io/gh/scverse/scvi-tools/branch/master/graph/badge.svg)](https://codecov.io/gh/YosefLab/scvi-tools)
[![Code
Style](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/python/black)
[![Downloads](https://pepy.tech/badge/scvi-tools)](https://pepy.tech/project/scvi-tools)
[![Project chat](https://img.shields.io/badge/zulip-join_chat-brightgreen.svg)](https://scverse.zulipchat.com/)

[scvi-tools](https://scvi-tools.org/) (single-cell variational inference
tools) is a package for probabilistic modeling and analysis of single-cell omics
data, built on top of [PyTorch](https://pytorch.org) and
[AnnData](https://anndata.readthedocs.io/en/latest/).

# Analysis of single-cell omics data

scvi-tools is composed of models that perform many analysis tasks across single- or multi-omics:

-   Dimensionality reduction
-   Data integration
-   Automated annotation
-   Factor analysis
-   Doublet detection
-   Spatial deconvolution
-   and more!

In the [user guide](https://docs.scvi-tools.org/en/stable/user_guide/index.html), we provide an overview of each model.
All model implementations have a high-level API that interacts with
[scanpy](http://scanpy.readthedocs.io/) and includes standard save/load functions, GPU acceleration, etc.

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
