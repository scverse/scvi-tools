import io
import sys
import types
from IPython import get_ipython
from nbformat import read
from IPython.core.interactiveshell import InteractiveShell
import os
import shutil
import matplotlib.pyplot as plt


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

        try:
            i = 0
            for cell in nb.cells:
                if cell.cell_type == 'code':
                    if i != 0:  # the first code cell is cd ../.. which is not needed here and would generate errors
                        # transform the input to executable Python
                        code = self.shell.input_transformer_manager.transform_cell(cell.source)
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

path = 'docs/notebooks/'
test_path = 'tests/notebooks/'


def test_notebooks():
    prefix = 'annotation'
    # if exists, save and overwrite
    config_filename = prefix + '.config.json'  # By default, benchmark config, but overridden for tests
    config_filename_tmp = prefix + '.config.tmp.json'

    # Save content in memory in tmp file
    existing_config_to_restore = os.path.exists(path + config_filename)
    if existing_config_to_restore:
        shutil.copy(path + config_filename, path + config_filename_tmp)

    # Overwritting
    shutil.copy(test_path + config_filename, path + config_filename)
    try:
        import docs.notebooks.annotation
        docs.notebooks.annotation.allow_notebook_for_test()
        plt.close('all')
    except BaseException:
        raise
    finally:
        if existing_config_to_restore:
            shutil.move(path + config_filename_tmp, path + config_filename)

    prefix = 'scRNA_and_smFISH'
    # if exists, save and overwrite
    config_filename = prefix + '.config.json'  # By default, benchmark config, but overridden for tests
    config_filename_tmp = prefix + '.config.tmp.json'

    # Save content in memory in tmp file
    existing_config_to_restore = os.path.exists(path + config_filename)
    if existing_config_to_restore:
        shutil.copy(path + config_filename, path + config_filename_tmp)

    # Overwritting
    shutil.copy(test_path + config_filename, path + config_filename)
    try:
        import docs.notebooks.scRNA_and_smFISH
        docs.notebooks.scRNA_and_smFISH.allow_notebook_for_test()
        plt.close('all')
    except BaseException:
        raise
    finally:
        if existing_config_to_restore:
            shutil.move(path + config_filename_tmp, path + config_filename)

    prefix = 'data_loading'
    # if exists, save and overwrite
    config_filename = prefix + '.config.json'  # By default, benchmark config, but overridden for tests
    config_filename_tmp = prefix + '.config.tmp.json'

    # Save content in memory in tmp file
    existing_config_to_restore = os.path.exists(path + config_filename)
    if existing_config_to_restore:
        shutil.copy(path + config_filename, path + config_filename_tmp)

    # Overwritting
    shutil.copy(test_path + config_filename, path + config_filename)
    try:
        import docs.notebooks.data_loading
        docs.notebooks.data_loading.allow_notebook_for_test()
        plt.close('all')
    except BaseException:
        raise
    finally:
        if existing_config_to_restore:
            shutil.move(path + config_filename_tmp, path + config_filename)

    prefix = 'basic_tutorial'
    # if exists, save and overwrite
    config_filename = prefix + '.config.json'  # By default, benchmark config, but overridden for tests
    config_filename_tmp = prefix + '.config.tmp.json'

    # Save content in memory in tmp file
    existing_config_to_restore = os.path.exists(path + config_filename)
    if existing_config_to_restore:
        shutil.copy(path + config_filename, path + config_filename_tmp)

    # Overwritting
    shutil.copy(test_path + config_filename, path + config_filename)
    try:
        import docs.notebooks.basic_tutorial
        docs.notebooks.basic_tutorial.allow_notebook_for_test()
        plt.close('all')
    except BaseException:
        raise
    finally:
        if existing_config_to_restore:
            shutil.move(path + config_filename_tmp, path + config_filename)

    prefix = 'scVI_reproducibility'
    # if exists, save and overwrite
    config_filename = prefix + '.config.json'  # By default, benchmark config, but overridden for tests
    config_filename_tmp = prefix + '.config.tmp.json'

    # Save content in memory in tmp file
    existing_config_to_restore = os.path.exists(path + config_filename)
    if existing_config_to_restore:
        shutil.copy(path + config_filename, path + config_filename_tmp)

    # Overwritting
    shutil.copy(test_path + config_filename, path + config_filename)
    try:
        import docs.notebooks.scVI_reproducibility
        docs.notebooks.scVI_reproducibility.allow_notebook_for_test()
        plt.close('all')
    except BaseException:
        raise
    finally:
        if existing_config_to_restore:
            shutil.move(path + config_filename_tmp, path + config_filename)
