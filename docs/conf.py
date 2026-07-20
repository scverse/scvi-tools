import importlib.util
import inspect
import os
import re
import subprocess
import sys
import types
from pathlib import Path
from importlib.metadata import metadata
from datetime import datetime
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Any

HERE = Path(__file__).parent
sys.path[:0] = [str(HERE.parent / "src"), str(HERE.parent), str(HERE / "extensions")]

# -- Project information -----------------------------------------------------

info = metadata("scvi-tools")
project_name = info["Name"]
author = info["Author"]
copyright = f"{datetime.now():%Y}, {author}."
version = info["Version"]
repository_url = f"https://github.com/scverse/{project_name}"

# The full version, including alpha/beta/rc tags
release = info["Version"]

bibtex_bibfiles = ["references.bib"]
templates_path = ["_templates"]
nitpicky = True  # Warn about broken links
needs_sphinx = "4.0"

# Bare (unqualified) names used in type annotations that autodoc/autodoc-typehints
# should resolve against their fully-qualified, intersphinx-resolvable targets.
autodoc_type_aliases = {
    "AnnData": "anndata.AnnData",
    "Tensor": "torch.Tensor",
}

# Reference targets that can never resolve, either because they leak a
# third-party library's private submodule path via `__module__` (with no public
# alias available), because the target is inherited third-party docstring
# content rendered in our own doc pages, or because of a long-standing
# py:class/py:data domain mismatch between sphinx-autodoc-typehints and the
# target project's own inventory.
nitpick_ignore = [
    # sphinx-autodoc-typehints emits typing.Union as py:data, but CPython's own
    # inventory registers it as py:class.
    ("py:data", "typing.Union"),
    # numpy.array is a function, not a class; used informally as a type name in
    # some docstrings.
    ("py:class", "numpy.array"),
    ("py:class", "numpy._typing.TypeAliasType"),
    # Legacy torch tensor aliases/internals not present in torch's inventory.
    ("py:class", "torch.FloatTensor"),
    ("py:class", "torch._VariableFunctionsClass.tensor"),
    # Internal typing helpers not meant to have their own documentation page.
    ("py:class", "scvi.train._config.KwargsConfig"),
    ("py:class", "scvi.external.contrastivevi._contrastive_dataloader.ContrastiveDataLoader"),
    # Inherited members from lightning.pytorch base classes: their first-line
    # docstrings reference other members by short name, which only resolves in
    # lightning's own docs, not ours.
    ("py:meth", "save_hyperparameters"),
    ("py:class", "AttributeDict"),
    # Inherited members from torch.distributions.Distribution whose docstrings
    # use short names that only resolve in torch's own docs.
    ("py:class", "Tensor"),
    ("py:class", "Distribution"),
    ("py:class", "return"),
    # napoleon/numpydoc occasionally mis-splits a "type, optional" parameter type
    # string and tries to cross-reference the word "optional" itself; this also
    # occurs in torch's own inherited docstrings (e.g. torch.nn.functional.one_hot),
    # so it can't be fixed purely on our side.
    ("py:class", "optional"),
    ("py:class", "LongTensor"),
]
nitpick_ignore_regex = [
    # mudata.MuData is aliased from a private path (mudata._core.mudata.MuData)
    # that intersphinx can't always disambiguate depending on inventory caching.
    (r"py:class", r"mudata\._core\.mudata\.MuData"),
]

html_context = {
    "display_github": True,  # Integrate GitHub
    "github_user": "scverse",  # Username
    "github_repo": project_name,  # Repo name
    "github_version": "main",  # Version
    "conf_py_path": "/docs/",  # Path in the checkout to the docs root
}

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = [
    "myst_nb",
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.linkcode",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx_autodoc_typehints",  # needs to be after napoleon
    "sphinx.ext.extlinks",
    "sphinx.ext.autosummary",
    "sphinxcontrib.bibtex",
    *[p.stem for p in (HERE / "extensions").glob("*.py")],
    "sphinx_copybutton",
    "sphinx_design",
    "sphinxext.opengraph",
    "hoverxref.extension",
]


# for sharing urls with nice info
ogp_site_url = "https://docs.scvi-tools.org/"
ogp_image = "https://docs.scvi-tools.org/en/stable/_static/logo.png"


# Generate the API documentation when building
autosummary_generate = True
autodoc_member_order = "bysource"
bibtex_reference_style = "author_year"
napoleon_google_docstring = True  # for pytorch lightning
napoleon_numpy_docstring = True  # use numpydoc style
napoleon_include_init_with_doc = False
napoleon_use_rtype = True  # having a separate entry generally helps readability
napoleon_use_param = True
napoleon_custom_sections = [("Params", "Parameters")]
todo_include_todos = False
myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_image",
    "html_admonition",
]
myst_url_schemes = ("http", "https", "mailto")
# Auto-generate #slug anchors for headings so in-page links like
# [prerequisites](#prerequisites) resolve (up to h4, the deepest level used in docs).
myst_heading_anchors = 4
nb_output_stderr = "remove"
nb_execution_mode = "off"
nb_merge_streams = True
typehints_defaults = "braces"
autodoc_mock_imports = []
if os.environ.get("READTHEDOCS") == "True":
    autodoc_mock_imports += ["hyperopt", "ray", "ray.tune", "scib_metrics", "muon", "mlx"]
    try:
        import scvi.autotune
    except Exception:
        import scvi

        _autotune_error = ModuleNotFoundError(
            "Autotune requires optional dependencies; install scvi-tools[autotune]."
        )
        autotune_stub = types.ModuleType("scvi.autotune")

        class AutotuneExperiment:
            """Autotune requires optional dependencies; install scvi-tools[autotune]."""

            def __init__(self, *args, **kwargs):
                raise _autotune_error

        class ScibTuneReportCheckpointCallback:
            """Autotune requires optional dependencies; install scvi-tools[autotune]."""

            def __init__(self, *args, **kwargs):
                raise _autotune_error

        def run_autotune(*args, **kwargs):
            """Autotune requires optional dependencies; install scvi-tools[autotune]."""
            raise _autotune_error

        autotune_stub.AutotuneExperiment = AutotuneExperiment
        autotune_stub.ScibTuneReportCheckpointCallback = ScibTuneReportCheckpointCallback
        autotune_stub.run_autotune = run_autotune
        autotune_stub.__all__ = [
            "AutotuneExperiment",
            "ScibTuneReportCheckpointCallback",
            "run_autotune",
        ]
        sys.modules["scvi.autotune"] = autotune_stub
        scvi.autotune = autotune_stub

source_suffix = {
    ".rst": "restructuredtext",
    ".ipynb": "myst-nb",
    ".myst": "myst-nb",
}

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints"]

# extlinks config
extlinks = {
    "issue": (f"{repository_url}/issues/%s", "#%s"),
    "pr": (f"{repository_url}/pull/%s", "#%s"),
    "ghuser": ("https://github.com/%s", "@%s"),
}

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
    "lightning": ("https://lightning.ai/docs/pytorch/stable/", None),
    "pyro": ("http://docs.pyro.ai/en/stable/", None),
    "flax": ("https://flax.readthedocs.io/en/latest/", None),
    "jax": ("https://jax.readthedocs.io/en/latest/", None),
    "ml_collections": ("https://ml-collections.readthedocs.io/en/latest/", None),
    "mudata": ("https://mudata.readthedocs.io/stable/", None),
    "ray": ("https://docs.ray.io/en/latest/", None),
    "huggingface_hub": ("https://huggingface.co/docs/huggingface_hub/main/en", None),
    "sparse": ("https://sparse.pydata.org/en/stable/", None),
    "annbatch": ("https://annbatch.readthedocs.io/en/latest/", None),
    "rich": ("https://rich.readthedocs.io/en/stable/", None),
    "xarray": ("https://docs.xarray.dev/en/stable/", None),
    "torch_geometric": ("https://pytorch-geometric.readthedocs.io/en/latest/", None),
}

# -- Options for HTML output -------------------------------------------

# html_show_sourcelink = True
html_theme = "sphinx_book_theme"
html_title = project_name

html_logo = "_static/logo.png"

html_theme_options = {
    "repository_url": repository_url,
    "use_repository_button": True,
    "logo_only": True,
    "show_toc_level": 1,
    "launch_buttons": {"colab_url": "https://colab.research.google.com"},
    "path_to_docs": "docs/",
    "repository_branch": version,
}

pygments_style = "default"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
html_css_files = ["css/override.css"]
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
except Exception:
    pass

# (if no name found or relative ref, use commit hash instead)
if not git_ref or re.search(r"[\^~]", git_ref):
    try:
        git_ref = git("rev-parse", "HEAD")
    except Exception:
        git_ref = "main"

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
    except Exception:
        return None

    path = f"{path}#L{lineno}-L{lineno + len(src) - 1}"
    return f"{repository_url}/blob/{git_ref}/src/scvi/{path}"


# -- Config for hoverxref -------------------------------------------

hoverx_default_type = "tooltip"
hoverxref_domains = ["py"]
hoverxref_role_types = dict.fromkeys(
    ["ref", "class", "func", "meth", "attr", "exc", "data", "mod"],
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
    "ray",
]
# use proxied API endpoint on rtd to avoid CORS issues
if os.environ.get("READTHEDOCS"):
    hoverxref_api_host = "/_"
