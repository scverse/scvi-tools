import io
import sys
import types
from IPython import get_ipython
from nbformat import read
from IPython.core.interactiveshell import InteractiveShell
import os
import matplotlib.pyplot as plt
import re

base_path = os.getcwd()


def find_notebook(fullname, path=None):
    """find a notebook, given its fully qualified name and an optional path

    This turns "foo.bar" into "foo/bar.ipynb"
    and tries turning "Foo_Bar" into "Foo Bar" if Foo_Bar
    does not exist.
    """
    name = fullname.rsplit(".", 1)[-1]
    if not path:
        path = [""]
    for d in path:
        nb_path = os.path.join(d, name + ".ipynb")
        if os.path.isfile(nb_path):
            return nb_path
        # let import Notebook_Name find "Notebook Name.ipynb"
        nb_path = nb_path.replace("_", " ")
        if os.path.isfile(nb_path):
            return nb_path


class NotebookLoader(object):
    """Module Loader for Jupyter Notebooks"""

    def __init__(self, path=None):
        self.shell = InteractiveShell.instance()
        self.path = path

    def load_module(self, fullname):
        """import a notebook as a module"""
        path = find_notebook(fullname, [os.path.join(base_path, "tests/notebooks")])

        print("importing Jupyter notebook from %s" % path)

        # load the notebook object
        with io.open(path, "r", encoding="utf-8") as f:
            nb = read(f, 4)

        # create the module and add it to sys.modules
        mod = types.ModuleType(fullname)
        mod.__file__ = path
        mod.__loader__ = self
        mod.__dict__["get_ipython"] = get_ipython
        sys.modules[fullname] = mod

        # extra work to ensure that magics that would affect the user_ns
        # actually affect the notebook module's ns
        save_user_ns = self.shell.user_ns
        self.shell.user_ns = mod.__dict__

        try:
            i = 0
            for cell in nb.cells:
                if cell.cell_type == "code":
                    if (
                        i != 0
                    ):  # the first code cell is cd ../.. which is not needed here and would generate errors
                        # transform the input to executable Python
                        code = self.shell.input_transformer_manager.transform_cell(
                            cell.source
                        )
                        # replace parameters with test parameters to run faster
                        code = re.sub(r"n_epochs_all = None", "n_epochs_all = 1", code)
                        code = re.sub(r"n_cl = \d+", "n_cl = 3", code)
                        code = re.sub(
                            r"M_permutation = \d+", "M_permutation = 10", code
                        )
                        code = re.sub(r"n_samples = \d+", "n_samples = 1", code)
                        code = re.sub(
                            r"n_samples_tsne = \d+", "n_samples_tsne = 10", code
                        )
                        code = re.sub(
                            r"n_samples_posterior_density = \d+",
                            "n_samples_posterior_density = 2",
                            code,
                        )
                        code = re.sub(
                            "save_path = 'data/'",
                            "save_path = '" + os.getcwd() + "'",
                            code,
                        )
                        code = re.sub(
                            'save_path = "data/"',
                            "save_path = '" + os.getcwd() + "'",
                            code,
                        )
                        code = re.sub("show_plot = True", "show_plot = False", code)
                        code = re.sub("test_mode = False", "test_mode = True", code)
                        # run the code in themodule
                        exec(code, mod.__dict__)
                        plt.close("all")
                    i += 1
        finally:
            self.shell.user_ns = save_user_ns
        return mod


class NotebookFinder(object):
    """Module finder that locates Jupyter Notebooks"""

    def __init__(self):
        self.loaders = {}

    def find_module(self, fullname, path=None):
        nb_path = find_notebook(fullname, path)
        if not nb_path:
            return

        key = path
        if path:
            # lists aren't hashable
            key = os.path.sep.join(path)

        if key not in self.loaders:
            self.loaders[key] = NotebookLoader(path)
        return self.loaders[key]


sys.meta_path.append(NotebookFinder())


def test_notebooks_annotation(save_path):
    try:
        os.chdir(save_path)
        import notebooks.annotation

        notebooks.annotation.allow_notebook_for_test()
        plt.close("all")
    except BaseException:
        raise
    finally:
        os.chdir(path=base_path)


def test_notebooks_dataloading(save_path):
    try:
        os.chdir(save_path)
        import notebooks.data_loading

        notebooks.data_loading.allow_notebook_for_test()
        plt.close("all")
    except BaseException:
        raise
    finally:
        os.chdir(path=base_path)


def test_notebooks_basictutorial(save_path):
    try:
        os.chdir(save_path)
        import notebooks.basic_tutorial

        notebooks.basic_tutorial.allow_notebook_for_test()
        plt.close("all")
    except BaseException:
        raise
    finally:
        os.chdir(path=base_path)


def test_notebooks_gimvitutorial(save_path):
    try:
        os.chdir(save_path)
        import notebooks.gimvi_tutorial

        notebooks.gimvi_tutorial.allow_notebook_for_test()
        plt.close("all")
    except BaseException:
        raise
    finally:
        os.chdir(path=base_path)


def test_notebooks_reproducibility(save_path):
    try:
        os.chdir(save_path)
        import notebooks.scVI_reproducibility

        notebooks.scVI_reproducibility.allow_notebook_for_test()
        plt.close("all")
    except BaseException:
        raise
    finally:
        os.chdir(path=base_path)


def test_notebooks_harmonization(save_path):
    try:
        os.chdir(save_path)
        import notebooks.harmonization

        notebooks.harmonization.allow_notebook_for_test()
        plt.close("all")
    except BaseException:
        raise
    finally:
        os.chdir(path=base_path)


def test_notebooks_scanpy_api(save_path):
    try:
        os.chdir(save_path)
        import notebooks.scanpy_pbmc3k

        print(save_path)
        notebooks.scanpy_pbmc3k.allow_notebook_for_test()
        plt.close("all")
    except BaseException:
        raise
    finally:
        os.chdir(path=base_path)


def test_notebooks_autotune(save_path):
    try:
        os.chdir(save_path)
        import notebooks.autotune_advanced_notebook

        print(save_path)
        notebooks.autotune_advanced_notebook.allow_notebook_for_test()
        plt.close("all")
    except BaseException:
        raise
    finally:
        os.chdir(path=base_path)


def test_notebooks_totalvi(save_path):
    try:
        os.chdir(save_path)
        import notebooks.totalVI

        print(save_path)
        notebooks.totalVI.allow_notebook_for_test()
        plt.close("all")
    except BaseException:
        raise
    finally:
        os.chdir(path=base_path)
