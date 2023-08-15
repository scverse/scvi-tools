<img src="https://github.com/scverse/scvi-tools/blob/main/docs/_static/scvi-tools-horizontal.svg?raw=true" width="400" alt="scvi-tools">

[![Stars](https://img.shields.io/github/stars/scverse/scvi-tools?logo=GitHub&color=yellow)](https://github.com/scverse/scvi-tools/stargazers)
[![PyPI](https://img.shields.io/pypi/v/scvi-tools.svg)](https://pypi.org/project/scvi-tools)
[![PyPIDownloads](https://static.pepy.tech/badge/scvi-tools)](https://pepy.tech/project/scvi-tools)
[![CondaDownloads](https://img.shields.io/conda/dn/conda-forge/scvi-tools?logo=Anaconda)](https://anaconda.org/conda-forge/scvi-tools)
[![Docs](https://readthedocs.org/projects/scvi/badge/?version=latest)](https://scvi.readthedocs.io/en/stable/?badge=stable)
[![Build](https://github.com/scverse/scvi-tools/actions/workflows/build.yml/badge.svg)](https://github.com/scverse/scvi-tools/actions/workflows/build.yml/)
[![Coverage](https://codecov.io/gh/scverse/scvi-tools/branch/main/graph/badge.svg)](https://codecov.io/gh/scverse/scvi-tools)
[![Discourse](https://img.shields.io/discourse/posts?color=yellow&logo=discourse&server=https%3A%2F%2Fdiscourse.scverse.org)](https://discourse.scverse.org/)
[![Chat](https://img.shields.io/badge/zulip-join_chat-brightgreen.svg)](https://scverse.zulipchat.com/)
[![Powered by NumFOCUS][badge-numfocus]][link-numfocus]

[scvi-tools](https://scvi-tools.org/) (single-cell variational inference
tools) is a package for probabilistic modeling and analysis of single-cell omics
data, built on top of [PyTorch](https://pytorch.org) and
[AnnData](https://anndata.readthedocs.io/en/latest/).

[//]: # "numfocus-fiscal-sponsor-attribution"

scvi-tools is part of the scverse project ([website](https://scverse.org), [governance](https://scverse.org/about/roles)) and is fiscally sponsored by [NumFOCUS](https://numfocus.org/).
Please consider making a tax-deductible [donation](https://numfocus.org/donate-to-scverse) to help the project pay for developer time, professional services, travel, workshops, and a variety of other needs.

<a href="https://numfocus.org/project/scverse">
  <img
    src="https://raw.githubusercontent.com/numfocus/templates/master/images/numfocus-logo.png"
    width="200"
  >
</a>

[badge-numfocus]: https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A
[link-numfocus]: http://numfocus.org

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
repository](https://github.com/scverse/simple-scvi) as a
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

If you use `scvi-tools` in your work, please cite

> **A Python library for probabilistic analysis of single-cell omics data**
>
> Adam Gayoso, Romain Lopez, Galen Xing, Pierre Boyeau, Valeh Valiollah Pour Amiri, Justin Hong, Katherine Wu, Michael Jayasuriya, Edouard Mehlman, Maxime Langevin, Yining Liu, Jules Samaran, Gabriel Misrachi, Achille Nazaret, Oscar Clivio, Chenling Xu, Tal Ashuach, Mariano Gabitto, Mohammad Lotfollahi, Valentine Svensson, Eduardo da Veiga Beltrame, Vitalii Kleshchevnikov, Carlos Talavera-López, Lior Pachter, Fabian J. Theis, Aaron Streets, Michael I. Jordan, Jeffrey Regier & Nir Yosef
>
> _Nature Biotechnology_ 2022 Feb 07. doi: [10.1038/s41587-021-01206-w](https://doi.org/10.1038/s41587-021-01206-w).

along with the publicaton describing the model used.

You can cite the scverse publication as follows:

> **The scverse project provides a computational ecosystem for single-cell omics data analysis**
>
> Isaac Virshup, Danila Bredikhin, Lukas Heumos, Giovanni Palla, Gregor Sturm, Adam Gayoso, Ilia Kats, Mikaela Koutrouli, Scverse Community, Bonnie Berger, Dana Pe’er, Aviv Regev, Sarah A. Teichmann, Francesca Finotello, F. Alexander Wolf, Nir Yosef, Oliver Stegle & Fabian J. Theis
>
> _Nature Biotechnology_ 2023 Apr 10. doi: [10.1038/s41587-023-01733-8](https://doi.org/10.1038/s41587-023-01733-8).
