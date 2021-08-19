#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# scvi documentation build configuration file, created by
# sphinx-quickstart on Fri Jun  9 13:47:02 2017.
#
# This file is execfile()d with the current directory set to its
# containing dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

# If extensions (or modules to document with autodoc) are in another
# directory, add these directories to sys.path here. If the directory is
# relative to the documentation root, use os.path.abspath to make it
# absolute, like shown here.
#
import sys
from pathlib import Path

HERE = Path(__file__).parent
sys.path[:0] = [str(HERE.parent), str(HERE / "extensions")]

import scvi  # noqa

# -- General configuration ---------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
needs_sphinx = "3.4"  # Nicer param docs

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "nbsphinx",
    "nbsphinx_link",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx_autodoc_typehints",  # needs to be after napoleon
    "sphinx.ext.autosummary",
    "scanpydoc.elegant_typehints",
    "scanpydoc.definition_list_typed_field",
    "scanpydoc.autosummary_generate_imported",
    *[p.stem for p in (HERE / "extensions").glob("*.py")],
    "sphinx_copybutton",
    "sphinx_gallery.load_style",
    "sphinx_tabs.tabs",
]

# nbsphinx specific settings
exclude_patterns = ["_build", "**.ipynb_checkpoints"]
nbsphinx_execute = "never"

templates_path = ["_templates"]
source_suffix = ".rst"

# Generate the API documentation when building
autosummary_generate = True
autodoc_member_order = "bysource"
napoleon_google_docstring = True  # for pytorch lightning
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = True  # having a separate entry generally helps readability
napoleon_use_param = True
napoleon_custom_sections = [("Params", "Parameters")]
todo_include_todos = False
numpydoc_show_class_members = False
annotate_defaults = True  # scanpydoc option, look into why we need this

# The master toctree document.
master_doc = "index"

intersphinx_mapping = dict(
    anndata=("https://anndata.readthedocs.io/en/stable/", None),
    ipython=("https://ipython.readthedocs.io/en/stable/", None),
    matplotlib=("https://matplotlib.org/", None),
    numpy=("https://numpy.org/doc/stable/", None),
    pandas=("https://pandas.pydata.org/docs/", None),
    python=("https://docs.python.org/3", None),
    scipy=("https://docs.scipy.org/doc/scipy/reference/", None),
    sklearn=("https://scikit-learn.org/stable/", None),
    torch=("https://pytorch.org/docs/master/", None),
    scanpy=("https://scanpy.readthedocs.io/en/stable/", None),
    pytorch_lightning=("https://pytorch-lightning.readthedocs.io/en/stable/", None),
    pyro=("http://docs.pyro.ai/en/stable/", None),
)


# General information about the project.
project = u"scvi-tools"
copyright = u"2021, Yosef Lab, UC Berkeley"
author = u"Romain Lopez, Adam Gayoso, Pierre Boyeau, Galen Xing"

# The version info for the project you're documenting, acts as replacement
# for |version| and |release|, also used in various other places throughout
# the built documents.
#
# The short X.Y version.
version = scvi.__version__
# The full version, including alpha/beta/rc tags.
release = scvi.__version__

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = None

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This patterns also effect to html_static_path and html_extra_path
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "default"
pygments_dark_style = "default"

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False


# -- Options for HTML output -------------------------------------------

html_show_sourcelink = True
html_theme = "pydata_sphinx_theme"

html_context = dict(
    # display_github=True,  # Integrate GitHub
    github_user="YosefLab",  # Username
    github_repo="scvi-tools",  # Repo name
    github_version="master",  # Version
    doc_path="docs/",  # Path in the checkout to the docs root
)
# Set link name generated in the top bar.
html_title = "scvi-tools"
html_logo = "_static/logo.png"

html_theme_options = {
    "github_url": "https://github.com/YosefLab/scvi-tools",
    "twitter_url": "https://twitter.com/YosefLab",
}
# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
html_css_files = ["css/override.css", "css/sphinx_gallery.css"]
html_show_sphinx = False


nbsphinx_prolog = r"""
.. raw:: html

{{% set docname = env.doc2path(env.docname, base=None).split("/")[-1] %}}

.. raw:: html

    <style>
        p {{
            margin-bottom: 0.5rem;
        }}
    </style>

.. raw:: html

    <div class="admonition note">
        <p class="admonition-title">Note</p>
        <p>
        This page was generated from
        <a class="reference external" href="https://github.com/yoseflab/scvi-tutorials/tree/{version}/">{docname}</a>.
        Interactive online version:
        <span style="white-space: nowrap;"><a href="https://colab.research.google.com/github/yoseflab/scvi_tutorials/blob/{version}/{docname}"><img alt="Colab badge" src="https://colab.research.google.com/assets/colab-badge.svg" style="vertical-align:text-bottom"></a>.</span>
        </p>
    </div>
""".format(
    version=version, docname="{{ docname|e }}"
)
nbsphinx_thumbnails = {
    "user_guide/notebooks/data_loading": "_static/tutorials/anndata.svg",
    "user_guide/notebooks/api_overview": "_static/tutorials/overview.svg",
    "user_guide/notebooks/linear_decoder": "_static/tutorials/ldvae.svg",
    "user_guide/notebooks/scvi_in_R": "_static/tutorials/Rlogo.png",
    "user_guide/notebooks/harmonization": "_static/tutorials/scvi_batch.png",
    "user_guide/notebooks/totalVI": "_static/tutorials/totalvi_cell.svg",
    "user_guide/notebooks/AutoZI_tutorial": "_static/tutorials/history.png",
    "user_guide/notebooks/gimvi_tutorial": "_static/tutorials/gimvi.png",
    "user_guide/notebooks/scarches_scvi_tools": "_static/tutorials/scarches.png",
    "user_guide/notebooks/cite_scrna_integration_w_totalVI": "_static/tutorials/cite_scrna.png",
    "user_guide/notebooks/scVI_DE_worm": "_static/tutorials/worm.png",
    "user_guide/notebooks/stereoscope_heart_LV_tutorial": "_static/tutorials/stereoscope.png",
    "user_guide/notebooks/seed_labeling": "_static/tutorials/seed.png",
    "user_guide/notebooks/cellassign_tutorial": "_static/tutorials/cellassign.png",
    "user_guide/notebooks/DestVI_tutorial": "_static/tutorials/destvi.png",
    "user_guide/notebooks/PeakVI": "_static/tutorials/peakvi.png",
}


def setup(app):
    # https://github.com/pradyunsg/furo/issues/49
    app.config.pygments_dark_style = "default"
