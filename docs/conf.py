import importlib.util
import inspect
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

HERE = Path(__file__).parent
sys.path[:0] = [str(HERE.parent), str(HERE / "extensions")]

import scvi  # noqa

# -- General configuration ---------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
needs_sphinx = "4.3"  # Nicer param docs

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.linkcode",
    "nbsphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx_autodoc_typehints",  # needs to be after napoleon
    "sphinx.ext.autosummary",
    "sphinxcontrib.bibtex",
    *[p.stem for p in (HERE / "extensions").glob("*.py")],
    "sphinx_copybutton",
    "sphinx_design",
    "sphinxext.opengraph",
    "hoverxref.extension",
]
github_repo = "https://github.com/scverse/scvi-tools"


# for sharing urls with nice info
ogp_site_url = "https://docs.scvi-tools.org/"
ogp_image = "https://docs.scvi-tools.org/en/stable/_static/logo.png"

# nbsphinx specific settings
exclude_patterns = ["_build", "**.ipynb_checkpoints"]
nbsphinx_execute = "never"

bibtex_bibfiles = ["references.bib"]
templates_path = ["_templates"]

# Generate the API documentation when building
autosummary_generate = True
autodoc_member_order = "bysource"
bibtex_reference_style = "author_year"
napoleon_google_docstring = True  # for pytorch lightning
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = True  # having a separate entry generally helps readability
napoleon_use_param = True
napoleon_custom_sections = [("Params", "Parameters")]
typehints_defaults = "braces"
todo_include_todos = False
numpydoc_show_class_members = False
annotate_defaults = True  # scanpydoc option, look into why we need this
myst_enable_extensions = [
    "colon_fence",
    "dollarmath",
    "amsmath",
]

# The master toctree document.
master_doc = "index"

intersphinx_mapping = {
    "anndata": ("https://anndata.readthedocs.io/en/stable/", None),
    "ipython": ("https://ipython.readthedocs.io/en/stable/", None),
    "matplotlib": ("https://matplotlib.org/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "python": ("https://docs.python.org/3", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/reference/", None),
    "sklearn": ("https://scikit-learn.org/stable/", None),
    "torch": ("https://pytorch.org/docs/master/", None),
    "scanpy": ("https://scanpy.readthedocs.io/en/stable/", None),
    "pytorch_lightning": ("https://pytorch-lightning.readthedocs.io/en/stable/", None),
    "pyro": ("http://docs.pyro.ai/en/stable/", None),
    "pymde": ("https://pymde.org/", None),
    "flax": ("https://flax.readthedocs.io/en/latest/", None),
    "jax": ("https://jax.readthedocs.io/en/latest/", None),
    "ml_collections": ("https://ml-collections.readthedocs.io/en/latest/", None),
    "mudata": ("https://mudata.readthedocs.io/en/latest/", None),
    "ray": ("https://docs.ray.io/en/latest/", None),
    "huggingface_hub": ("https://huggingface.co/docs/huggingface_hub/main/en", None),
}


# General information about the project.
project = "scvi-tools"
copyright = "2022, Yosef Lab, UC Berkeley"
author = "The scvi-tools development team"

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
pygments_dark_style = "native"

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False


# -- Options for HTML output -------------------------------------------

# html_show_sourcelink = True
html_theme = "furo"

# Set link name generated in the top bar.
html_title = "scvi-tools"
html_logo = "_static/logo.png"

html_theme_options = {
    "sidebar_hide_name": True,
    "light_css_variables": {
        "color-brand-primary": "#003262",
        "color-brand-content": "#003262",
        "admonition-font-size": "var(--font-size-normal)",
        "admonition-title-font-size": "var(--font-size-normal)",
        "code-font-size": "var(--font-size--small)",
    },
}
# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
html_css_files = ["css/override.css"]
html_show_sphinx = False

# -- Nbsphinx prolog -------------------------------------------


nbsphinx_prolog = r"""
.. raw:: html

{{% set docname = env.doc2path(env.docname, base=None).split("/")[-1] %}}

.. raw:: html

    <style>
        p {{
            margin-bottom: 0.5rem;
        }}
        /* Main index page overview cards */
        /* https://github.com/spatialaudio/nbsphinx/pull/635/files */
        .jp-RenderedHTMLCommon table,
        div.rendered_html table {{
        border: none;
        border-collapse: collapse;
        border-spacing: 0;
        font-size: 12px;
        table-layout: fixed;
        color: inherit;
        }}

        body:not([data-theme=light]) .jp-RenderedHTMLCommon tbody tr:nth-child(odd),
        body:not([data-theme=light]) div.rendered_html tbody tr:nth-child(odd) {{
        background: rgba(255, 255, 255, .1);
        }}
    </style>

.. raw:: html

    <div class="admonition note">
        <p class="admonition-title">Note</p>
        <p>
        This page was generated from
        <a class="reference external" href="https://github.com/scverse/scvi-tutorials/tree/{version}/">{docname}</a>.
        Interactive online version:
        <span style="white-space: nowrap;"><a href="https://colab.research.google.com/github/scverse/scvi-tutorials/blob/{version}/{docname}"><img alt="Colab badge" src="https://colab.research.google.com/assets/colab-badge.svg" style="vertical-align:text-bottom"></a>.</span>
        Some tutorial content may look better in light mode.
        </p>
    </div>
""".format(
    version=version, docname="{{ docname|e }}"
)

# -- Config for linkcode -------------------------------------------


def git(*args):
    """Run git command and return output as string."""
    return subprocess.check_output(["git", *args]).strip().decode()


# https://github.com/DisnakeDev/disnake/blob/7853da70b13fcd2978c39c0b7efa59b34d298186/docs/conf.py#L192
# Current git reference. Uses branch/tag name if found, otherwise uses commit hash
git_ref = None
try:
    git_ref = git("name-rev", "--name-only", "--no-undefined", "HEAD")
    git_ref = re.sub(r"^(remotes/[^/]+|tags)/", "", git_ref)
except Exception:  # noqa: B902
    pass

# (if no name found or relative ref, use commit hash instead)
if not git_ref or re.search(r"[\^~]", git_ref):
    try:
        git_ref = git("rev-parse", "HEAD")
    except Exception:  # noqa: B902
        git_ref = "master"

# https://github.com/DisnakeDev/disnake/blob/7853da70b13fcd2978c39c0b7efa59b34d298186/docs/conf.py#L192
_scvi_tools_module_path = os.path.dirname(importlib.util.find_spec("scvi").origin)  # type: ignore


def linkcode_resolve(domain, info):
    """Determine the URL corresponding to Python object."""
    if domain != "py":
        return None

    try:
        obj: Any = sys.modules[info["module"]]
        for part in info["fullname"].split("."):
            obj = getattr(obj, part)
        obj = inspect.unwrap(obj)

        if isinstance(obj, property):
            obj = inspect.unwrap(obj.fget)  # type: ignore

        path = os.path.relpath(inspect.getsourcefile(obj), start=_scvi_tools_module_path)  # type: ignore
        src, lineno = inspect.getsourcelines(obj)
    except Exception:  # noqa: B902
        return None

    path = f"{path}#L{lineno}-L{lineno + len(src) - 1}"
    return f"{github_repo}/blob/{git_ref}/scvi/{path}"


hoverx_default_type = "tooltip"
hoverxref_domains = ["py"]
hoverxref_role_types = dict.fromkeys(
    ["ref", "class", "func", "meth", "attr", "exc", "data"],
    "tooltip",
)
hoverxref_intersphinx = [
    "python",
    "numpy",
    "scanpy",
    "anndata",
    "pytorch_lightning",
    "scipy",
    "pandas",
    "ml_collections",
]
# use proxied API endpoint on rtd to avoid CORS issues
if os.environ.get("READTHEDOCS"):
    hoverxref_api_host = "/_"
