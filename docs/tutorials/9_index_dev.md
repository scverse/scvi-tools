# Development

This section features tutorials to guide through the construction of a model with scvi-tools. For an example of how scvi-tools can be used in an independent package, see the simple-scvi [example].

:::{note}
For questions about developing with scvi-tools, please use the [forum] or [zulip].
:::

```{toctree}
:maxdepth: 1

notebooks/dev/data_tutorial
notebooks/dev/module_user_guide
notebooks/dev/model_user_guide
```

[forum]: https://discourse.scvi-tools.org/
[zulip]: https://scverse.zulipchat.com/
[example]: https://github.com/scverse/simple-scvi


```{customcard}
:path: notebooks/dev/data_tutorial
:tags:

Learn about how data is handled in scvi-tools
```

```{customcard}
:path: notebooks/dev/module_user_guide
:tags:

Implement a novel statistical method for single-cell omics data as a module
```

```{customcard}
:path: notebooks/dev/model_user_guide
:tags:

Implement an scvi-tools model class to provide a convenient interface for the lower-level module objects
```
