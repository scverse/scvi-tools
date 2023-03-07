# Documentation

`scvi-tools` (single-cell variational inference tools) is a package for end-to-end analysis of single-cell omics data primarily developed and maintained by the [Yosef Lab](https://yoseflab.github.io/) at UC Berkeley. `scvi-tools` has two components

-   Interface for easy use of a range of probabilistic models for single-cell omics (e.g., scVI, scANVI, totalVI).
-   Tools to build new probabilistic models, which are powered by PyTorch, PyTorch Lightning, and Pyro.

If you find a model useful for your research, please consider citing the [scvi-tools manuscript](http://dx.doi.org/10.1038/s41587-021-01206-w) as well as the publication describing the model, which can be found in the corresponding documentation.

::::{grid} 1 2 3 3
:gutter: 2

:::{grid-item-card} Installation {octicon}`plug;1em;`
:link: installation
:link-type: doc

New to _scvi-tools_? Check out the installation guide.
:::

:::{grid-item-card} User guide {octicon}`info;1em;`
:link: user_guide/index
:link-type: doc

The user guide provides distilled mathematical descriptions of
the models implemented in scvi-tools and connects the math
with the code.
:::

:::{grid-item-card} API reference {octicon}`book;1em;`
:link: api/index
:link-type: doc

The API reference contains a detailed description of
the scvi-tools API.
:::

:::{grid-item-card} Tutorials {octicon}`play;1em;`
:link: tutorials/index
:link-type: doc

The tutorials walk you through real-world applications of scvi-tools models.
Developer tutorials help you build new probabilistic models.
:::

:::{grid-item-card} Discussion {octicon}`megaphone;1em;`
:link: https://discourse.scvi-tools.org

Need help? Reach out on our forum to get your questions answered!
:::

:::{grid-item-card} GitHub {octicon}`mark-github;1em;`
:link: https://github.com/yoseflab/scvi-tools

Find a bug? Interested in improving scvi-tools? Checkout our GitHub for the latest developments.
:::
::::

```{toctree}
:hidden: true
:maxdepth: 3
:titlesonly: true

installation
tutorials/index
user_guide/index
api/index
release_notes/index
references
contributing/index
Discussion <https://discourse.scvi-tools.org>
GitHub <https://github.com/YosefLab/scvi-tools>
```
