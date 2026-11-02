import importlib
import sys
import warnings
from pathlib import Path
from importlib.metadata import metadata
from datetime import datetime

from scverse_doc.config import DEFAULTS
from sphinx.ext.autodoc.mock import MockFinder
from scverse_doc.registry import intersphinx

HERE = Path(__file__).parent
sys.path[:0] = [str(HERE.parent / "src"), str(HERE.parent), str(HERE / "extensions")]

# -- Project information -----------------------------------------------------

info = metadata("scvi-tools")
project = info["Name"]
author = info["Author"]
copyright = f"{datetime.now():%Y}, {author}."
version = info["Version"]
urls = dict(pu.split(", ") for pu in info.get_all("Project-URL"))
source_repository = urls["Source"]

# The full version, including alpha/beta/rc tags
release = info["Version"]

bibtex_bibfiles = ["references.bib"]
templates_path = ["_templates"]

# The spatial/multimodal tutorial notebooks superseded by scviva-tools (resolVI, scVIVA,
# DestVI, gimVI, Tangram, Stereoscope, Cell2location, DiagVI) are no longer linked from any
# toctree now that their index cards point to the scviva-tools-hosted versions instead, but
# they're kept in the build (not `exclude_patterns`) since docs/user_guide/models/*.md still
# cross-reference them with {doc}. Remove this once those model pages are updated (planned for
# v1.6) and the notebooks themselves are removed.
suppress_warnings = ["toc.not_included"]

# Bare (unqualified) names used in type annotations that autodoc/autodoc-typehints
# should resolve against their fully-qualified, intersphinx-resolvable targets.
autodoc_type_aliases = {
    "AnnData": "anndata.AnnData",
    "Tensor": "torch.Tensor",
}

# Targets that only exist in docstrings we inherit from third-party base classes, which we
# render but cannot edit.
nitpick_ignore = [
    ("py:class", "AttributeDict"),
    ("py:class", "huggingface_hub.repocard.ModelCard"),
    ("py:class", "ml_collections.config_dict.config_dict.FrozenConfigDict"),
    ("py:class", "numpy._typing.TypeAliasType"),
    ("py:class", "numpy._typing._array_like.TypeAliasType"),
    ("py:class", "pathlib._local.Path"),
    ("py:class", "Distribution"),
    ("py:class", "LongTensor"),
    ("py:class", "LossRecord"),
    ("py:class", "Tensor"),
    ("py:class", "callable"),
    ("py:class", "optional"),
    ("py:class", "return"),
    ("py:meth", "save_hyperparameters"),
]

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = [
    "scverse_doc",
    "sphinx.ext.linkcode",
    "sphinx_autodoc_typehints",  # needs to be after napoleon
    "sphinx.ext.extlinks",
    "sphinxcontrib.bibtex",
    *[p.stem for p in (HERE / "extensions").glob("*.py")],
]

# for sharing urls with nice info
ogp_site_url = "https://docs.scvi-tools.org/"
ogp_image = "https://docs.scvi-tools.org/en/stable/_static/logo.png"


# Generate the API documentation when building
autodoc_member_order = "bysource"
bibtex_reference_style = "author_year"
napoleon_google_docstring = True  # for pytorch lightning
napoleon_custom_sections = [("Params", "Parameters")]
todo_include_todos = False
myst_url_schemes = ("http", "https", "mailto")
# Auto-generate #slug anchors for headings so in-page links like
# [prerequisites](#prerequisites) resolve (up to h4, the deepest level used in docs).
myst_heading_anchors = 4
typehints_defaults = "braces"
always_use_bars_union = True


def _importable(name: str) -> bool:
    try:
        importlib.import_module(name)
    except ModuleNotFoundError:
        return False
    except Exception as e:
        warnings.warn(
            f"{name!r} failed to import ({e!r}); mocking it for the documentation build.",
            stacklevel=2,
        )
        return False
    return True


# Optional dependencies that are not installable, or not loadable, on the machine building the
# docs. Mocking them keeps their signatures documented.
autodoc_mock_imports = [
    name
    for name in (
        "hyperopt",
        "ray",
        "ray.tune",
        "scib_metrics",
        "muon",
        "mlx",
        "lamindb",
        "tiledbsoma",
    )
    if not _importable(name)
]
# sphinx-autodoc-typehints looks objects up after autodoc's import-time mocks are gone,
# and scvi.model raises on access to mlxSCVI when mlx is missing, so keep them for the whole build.
sys.meta_path.insert(0, MockFinder(autodoc_mock_imports))

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = [
    *DEFAULTS.exclude_patterns,
    "tutorials/notebooks/.github/**",
    "tutorials/notebooks/README.md",
]

# extlinks config
extlinks = {
    "issue": (f"{source_repository}/issues/%s", "#%s"),
    "pr": (f"{source_repository}/pull/%s", "#%s"),
    "ghuser": ("https://github.com/%s", "@%s"),
}

intersphinx_mapping = intersphinx() | {
    "ipython": ("https://ipython.readthedocs.io/en/stable/", None),
    "sklearn": ("https://scikit-learn.org/stable/", None),
    "torch": ("https://pytorch.org/docs/master/", None),
    # lightning.ai intermittently 404s on objects.inv from RTD builders; fall back to a local snapshot
    "lightning": (
        "https://lightning.ai/docs/pytorch/stable/",
        (None, "_inventories/lightning.inv"),
    ),
    "pyro": ("https://docs.pyro.ai/en/stable/", None),
    "torchmetrics": (
        "https://lightning.ai/docs/torchmetrics/stable/",
        (None, "_inventories/torchmetrics.inv"),
    ),
    "flax": ("https://flax.readthedocs.io/en/latest/", None),
    "jax": ("https://jax.readthedocs.io/en/latest/", None),
    "ml_collections": ("https://ml-collections.readthedocs.io/en/latest/", None),
    "ray": ("https://docs.ray.io/en/latest/", None),
    "huggingface_hub": ("https://huggingface.co/docs/huggingface_hub/main/en", None),
    "sparse": ("https://sparse.pydata.org/en/stable/", None),
    "annbatch": ("https://annbatch.readthedocs.io/en/latest/", None),
    "rich": ("https://rich.readthedocs.io/en/stable/", None),
    "xarray": ("https://docs.xarray.dev/en/stable/", None),
    "torch_geometric": ("https://pytorch-geometric.readthedocs.io/en/latest/", None),
}

# -- Options for HTML output -------------------------------------------

html_logo = "_static/logo.png"
# The user guide numbers and references equations with \tag and \ref, which math-core does not support.
html_math_renderer = "mathjax"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
html_css_files = ["css/custom.css"]
html_js_files = ["js/custom.js"]
html_show_sphinx = False


def setup(app):
    """App setup hook."""
    app.add_config_value(
        "recommonmark_config",
        {
            "auto_toc_tree_section": "Contents",
            "enable_auto_toc_tree": True,
            "enable_math": True,
            "enable_inline_math": False,
            "enable_eval_rst": True,
        },
        True,
    )
