import io
import sys
import types
from IPython import get_ipython
from nbformat import read
from IPython.core.interactiveshell import InteractiveShell
import os
import matplotlib.pyplot as plt
import pickle
import re


def find_notebook(fullname, path=None):
    """find a notebook, given its fully qualified name and an optional path

    This turns "foo.bar" into "foo/bar.ipynb"
    and tries turning "Foo_Bar" into "Foo Bar" if Foo_Bar
    does not exist.
    """
    name = fullname.rsplit('.', 1)[-1]
    if not path:
        path = ['']
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
        path = find_notebook(fullname, self.path)

        print("importing Jupyter notebook from %s" % path)

        # load the notebook object
        with io.open(path, 'r', encoding='utf-8') as f:
            nb = read(f, 4)

        # create the module and add it to sys.modules
        # if name in sys.modules:
        #    return sys.modules[name]
        mod = types.ModuleType(fullname)
        mod.__file__ = path
        mod.__loader__ = self
        mod.__dict__['get_ipython'] = get_ipython
        sys.modules[fullname] = mod

        # extra work to ensure that magics that would affect the user_ns
        # actually affect the notebook module's ns
        save_user_ns = self.shell.user_ns
        self.shell.user_ns = mod.__dict__

        save_path = pickle.load(open('tests/data/path_test', 'rb'))

        try:
            i = 0
            for cell in nb.cells:
                if cell.cell_type == 'code':
                    if i != 0:  # the first code cell is cd ../.. which is not needed here and would generate errors
                        # transform the input to executable Python
                        code = self.shell.input_transformer_manager.transform_cell(cell.source)
                        # replace parameters with test parameters to run faster
                        code = re.sub("n_epochs_all = None", "n_epochs_all = 1", code)
                        code = re.sub("n_cl = \d+", "n_cl = 3", code)
                        code = re.sub("M_permutation = \d+", "M_permutation = 10", code)
                        code = re.sub("M_sampling = \d+", "M_sampling = 1", code)
                        code = re.sub("n_samples_tsne = \d+", "n_samples_tsne = 10", code)
                        code = re.sub("n_samples_posterior_density = \d+", "n_samples_posterior_density = 2", code)
                        code = re.sub("save_path = 'data/'", "save_path = '"+save_path+"'", code)
                        # run the code in themodule
                        exec(code, mod.__dict__)
                        plt.close('all')
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


def test_notebooks(save_path):
    pickle.dump(save_path, open('tests/data/path_test', 'wb'))

    try:
        import notebooks.annotation
        notebooks.annotation.allow_notebook_for_test()
        plt.close('all')

        import notebooks.scRNA_and_smFISH
        notebooks.scRNA_and_smFISH.allow_notebook_for_test()
        plt.close('all')

        import notebooks.data_loading
        notebooks.data_loading.allow_notebook_for_test()
        plt.close('all')

        import notebooks.basic_tutorial
        notebooks.basic_tutorial.allow_notebook_for_test()
        plt.close('all')

        import notebooks.scVI_reproducibility
        notebooks.scVI_reproducibility.allow_notebook_for_test()
        plt.close('all')

    except BaseException:
        raise

    finally:
        os.remove('tests/data/path_test')
