<a href="https://scvi-tools.org/">
  <img
    src="https://github.com/scverse/scvi-tools/blob/main/docs/_static/scvi-tools-horizontal.svg?raw=true"
    width="400"
    alt="scvi-tools"
  >
</a>

[![Stars][gh-stars-badge]][gh-stars-link]
[![PyPI][pypi-badge]][pypi-link]
[![PyPIDownloads][pepy-badge]][pepy-link]
[![CondaDownloads][conda-badge]][conda-link]
[![Docs][docs-badge]][docs-link]
[![Build][build-badge]][build-link]
[![Coverage][coverage-badge]][coverage-link]

[scvi-tools] (single-cell variational inference tools) is a package for probabilistic modeling and
analysis of single-cell omics data, built on top of [PyTorch] and [AnnData].

# Analysis of single-cell omics data

scvi-tools is composed of models that perform many analysis tasks across single-cell, multi, and
spatial omics data:

- Dimensionality reduction
- Data integration
- Automated annotation
- Factor analysis
- Doublet detection
- Spatial deconvolution
- and more!

In the [user guide], we provide an overview of each model. All model implementations have a
high-level API that interacts with [Scanpy] and includes standard save/load functions, GPU
acceleration, etc.

# Rapid development of novel probabilistic models

scvi-tools contains the building blocks to develop and deploy novel probabilistic models. These
building blocks are powered by popular probabilistic and machine learning frameworks such as
[PyTorch Lightning] and [Pyro]. For an overview of how the scvi-tools package is structured, you
may refer to the [codebase overview] page.

We recommend checking out the [skeleton repository] as a starting point for developing and
deploying new models with scvi-tools.

# Basic installation

For conda,

```bash
conda install scvi-tools -c conda-forge
```

and for pip,

```bash
pip install scvi-tools
```

Please be sure to install a version of [PyTorch] that is compatible with your GPU (if applicable).

# Resources

- Tutorials, API reference, and installation guides are available in the [documentation].
- For discussion of usage, check out our [forum].
- Please use the [issues] to submit bug reports.
- If you'd like to contribute, check out our [contributing guide].
- If you find a model useful for your research, please consider citing the corresponding
    publication.

# Reference

If you use `scvi-tools` in your work, please cite

> **A Python library for probabilistic analysis of single-cell omics data**
>
> Adam Gayoso, Romain Lopez, Galen Xing, Pierre Boyeau, Valeh Valiollah Pour Amiri, Justin Hong,
> Katherine Wu, Michael Jayasuriya, Edouard Mehlman, Maxime Langevin, Yining Liu, Jules Samaran,
> Gabriel Misrachi, Achille Nazaret, Oscar Clivio, Chenling Xu, Tal Ashuach, Mariano Gabitto,
> Mohammad Lotfollahi, Valentine Svensson, Eduardo da Veiga Beltrame, Vitalii Kleshchevnikov,
> Carlos Talavera-López, Lior Pachter, Fabian J. Theis, Aaron Streets, Michael I. Jordan,
> Jeffrey Regier & Nir Yosef
>
> _Nature Biotechnology_ 2022 Feb 07. doi: [10.1038/s41587-021-01206-w](https://doi.org/10.1038/s41587-021-01206-w).

along with the publication describing the model used.

You can cite the scverse publication as follows:

> **The scverse project provides a computational ecosystem for single-cell omics data analysis**
>
> Isaac Virshup, Danila Bredikhin, Lukas Heumos, Giovanni Palla, Gregor Sturm, Adam Gayoso,
> Ilia Kats, Mikaela Koutrouli, Scverse Community, Bonnie Berger, Dana Pe’er, Aviv Regev,
> Sarah A. Teichmann, Francesca Finotello, F. Alexander Wolf, Nir Yosef, Oliver Stegle &
> Fabian J. Theis
>
> _Nature Biotechnology_ 2023 Apr 10. doi: [10.1038/s41587-023-01733-8](https://doi.org/10.1038/s41587-023-01733-8).

scvi-tools is part of the scverse® project ([website](https://scverse.org),
[governance](https://scverse.org/about/roles)) and is fiscally sponsored by [NumFOCUS](https://numfocus.org/).

If you like scverse® and want to support our mission, please consider making a tax-deductible
[donation](https://numfocus.org/donate-to-scverse) to help the project pay for developer time,
professional services, travel, workshops, and a variety of other needs.

<div align="center">
<a href="https://numfocus.org/project/scverse">
  <img
    src="https://raw.githubusercontent.com/numfocus/templates/master/images/numfocus-logo.png"
    width="200"
  >
</a>
</div>

[anndata]: https://anndata.readthedocs.io/en/latest/
[build-badge]: https://github.com/scverse/scvi-tools/actions/workflows/build.yml/badge.svg
[build-link]: https://github.com/scverse/scvi-tools/actions/workflows/build.yml/
[codebase overview]: https://docs.scvi-tools.org/en/stable/user_guide/background/codebase_overview.html
[conda-badge]: https://img.shields.io/conda/dn/conda-forge/scvi-tools?logo=Anaconda
[conda-link]: https://anaconda.org/conda-forge/scvi-tools
[contributing guide]: https://docs.scvi-tools.org/en/stable/contributing/index.html
[coverage-badge]: https://codecov.io/gh/scverse/scvi-tools/branch/main/graph/badge.svg
[coverage-link]: https://codecov.io/gh/scverse/scvi-tools
[docs-badge]: https://readthedocs.org/projects/scvi/badge/?version=latest
[docs-link]: https://scvi.readthedocs.io/en/stable/?badge=stable
[documentation]: https://docs.scvi-tools.org/
[forum]: https://discourse.scvi-tools.org
[gh-stars-badge]: https://img.shields.io/github/stars/scverse/scvi-tools?style=flat&logo=GitHub&color=blue
[gh-stars-link]: https://github.com/scverse/scvi-tools/stargazers
[issues]: https://github.com/scverse/scvi-tools/issues
[pepy-badge]: https://static.pepy.tech/badge/scvi-tools
[pepy-link]: https://pepy.tech/project/scvi-tools
[pypi-badge]: https://img.shields.io/pypi/v/scvi-tools.svg
[pypi-link]: https://pypi.org/project/scvi-tools
[pyro]: https://pyro.ai/
[pytorch]: https://pytorch.org
[pytorch lightning]: https://lightning.ai/docs/pytorch/stable/
[scanpy]: http://scanpy.readthedocs.io/
[scvi-tools]: https://scvi-tools.org/
[skeleton repository]: https://github.com/scverse/simple-scvi
[user guide]: https://docs.scvi-tools.org/en/stable/user_guide/index.html
