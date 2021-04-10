.. highlight:: shell

============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.


Get Started!
------------

Ready to contribute? Here's how to set up `scvi-tools` for local development.

1. Fork the `scvi-tools` repo on GitHub.
2. Clone your fork locally::

    # Clone your fork of the repository (substitute in your username)
    git clone https://github.com/{your-username}/scvi-tools.git
    # Enter the cloned repository
    cd scvi-tools
    # Add our repository as a remote
    git remote add upstream https://github.com/yoseflab/scvi-tools.git
    # git branch --set-upstream-to "upstream/master"

3. Install your local copy into a virtualenv (or conda environment)::

    # If you have pyenv-virtualenv
    pyenv virtualenv scvi-tools-dev
    pyenv activate scvi-tools-dev
    # If you have conda
    conda create -n scvi-tools-dev
    conda activate scvi-tools-dev
    # Enter the cloned repository
    cd scvi-tools
    pip install -e ".[dev,docs,tutorials]"

4. **[Advanced users]** Install your local copy into a virtualenv with Poetry. Our preferred local installation method consists of using `pyenv-virtualenv` to create a virtualenv, and using `poetry` to create an editable local installation. If using this approach, please be sure to install `poetry` the `recommended <https://python-poetry.org/docs/#installation>`_ way. Once `poetry` is installed::

    pyenv virtualenv scvi-tools-dev
    pyenv activate scvi-tools-dev
    cd scvi-tools
    poetry install --extras "dev docs tutorials"

5. **[Optional]** Install a version of PyTorch that supports your GPU. This will even be the case if you use Poetry.

6. Create an ipykernel so you can use your environment with a Jupyter notebook. This will make this developement environment available through Jupyter notebook/lab. Inside your virtualenv::

    python -m ipykernel install --user --name=scvi-tools-dev

7. Install pre-commit, which will enforce the scvi-tools code style (black, flake8) on each of your commits::

    $ pre-commit install

8. Create a branch for local development::

    $ git checkout -b {your-branch-name}

   Now you can make your changes locally.

9. Add tests to the `/tests` directory. These files start with `test_` and contain functions that start similarly with `test_`.

10. When you're done making changes, run the tests using pytest::

    $ pytest tests/models/test_my_new_feature.py
    $ pytest tests/models/test_my_new_feature.py::test_particular_function_in_file

11. Commit your changes and push your branch to GitHub::

    $ git add <file> ...
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

12. Submit a pull request through the GitHub website.


Coding Standards
----------------
1. Don't duplicate code. Certainly no blocks longer than a couple of lines. It's almost always better to refactor than to duplicate blocks of code.
2. Almost all code should at least be run by a unit test. No pull request should decrease unit test coverage by much.
3. Document each new method and each new class with a docstring.
4. Don't commit commented-out code. Just delete it or store it somewhere outside of the repo. You probably aren't going to need it. At worse, it's stored in previous commits, from before it was commented out.
5. A pull request (PR) will typically close at least one Github issue. For these pull requests, write the issue it closes in the description, e.g. ``closes #210``. The issue will be automatically closed when the PR is merged.
6. Don't commit data to the repository, except perhaps a few small (< 50 KB) files of test data.
7. Respect the scvi-tools code style, the easiest way is to install pre-commit as described above.
8. Your PR will be evaluated by automated software for following coding standards, so it's easier to start with good practices.


Documenting Code
----------------
This section is under construction, but we use the same docstring style as Scanpy. See their `tutorial <https://scanpy.readthedocs.io/en/stable/dev/documentation.html#building-the-docs>`_ for more info.


Tips
----

1. `GitKraken <https://www.gitkraken.com/>`_ can be a useful GUI for using git locally.
2. ``git commit -m "my message" --no-verify`` allows overriding `pre-commit`.
3. Reach out on `gitter <https://gitter.im/scvi-tools/development>`_ if you need help.


Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated.
3. The pull request should work for Python 3.6-3.8. Your PR will be tested
   on these versions with our continuous integration checks.


Deploying
---------

A reminder for the maintainers on how to deploy.
Make sure all your changes are committed (including a release note entry).
First, please install Poetry.

Also, make sure you've tested your code using pytest by running::

$ pytest

Then run::

$ poetry version preversion # possible: major / minor / patch
$ poetry build
$ poetry publish

This will upload `scvi-tools` to PyPi. Also be sure to add a tag corresponding to the new version number on the tutorials repo, as the tagged repo is used for the Colab links.


Instructions on Uploading to conda
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
`scvi-tools` is available on bioconda channel. Typically, a PR will be automatically created once a new PyPi release is made.
This automated PR might need changes if we've changed dependencies. In that case, follow the below steps to upload a new version to bioconda channel.

Create a fork of bioconda-recipes on GitHub. Then::

$ git clone https://github.com/<USERNAME>/bioconda-recipes.git
$ git remote add upstream https://github.com/bioconda/bioconda-recipes.git

Update repo::

$ git checkout master
$ git pull origin master

Write a recipe::

$ git checkout -b my-recipe

Get the package's hash::

$ pip hash dist/scvi-tools-<NEW_VERSION_TAG>.tar.gz

Push changes, wait for tests to pass, submit pull request::

$ git push -u origin my-recipe

For this, it's easier to look at old scvi-tools PR's.

Write a GitHub release
~~~~~~~~~~~~~~~~~~~~~~

On the GitHub page, draft a release. This is important for ReadTheDocs, which uses the last release tag from GitHub as the stable version.
