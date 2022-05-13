```{highlight} shell
```

# Contributing

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.

## Get Started!

Ready to contribute? Here's how to set up `scvi-tools` for local development.

01. Fork the `scvi-tools` repo on GitHub.

02. Clone your fork locally:

    ```
    # Clone your fork of the repository (substitute in your username)
    git clone https://github.com/{your-username}/scvi-tools.git
    # Enter the cloned repository
    cd scvi-tools
    # Add our repository as a remote
    git remote add upstream https://github.com/yoseflab/scvi-tools.git
    # git branch --set-upstream-to "upstream/master"
    ```

03. Install your local copy into a virtualenv (or conda environment):

    ```
    # If you have pyenv-virtualenv
    pyenv virtualenv scvi-tools-dev
    pyenv activate scvi-tools-dev
    # If you have conda (omit the python parameter if you already have the relevant python version installed)
    conda create -n scvi-tools-dev python=3.8.8 # or any python >3.7 that is available (conda search python)
    conda activate scvi-tools-dev
    # Enter the cloned repository and install the package in editable mode
    cd scvi-tools
    pip install -e ".[dev,docs,tutorials]"
    ```

    To confirm that scvi-tools was successfully installed:

    ```
    # This should find the package. Note that other metadata (such as Version, Summary, etc.) might be missing. This
    # is expected because we use poetry instead of setup-tools. On a non-editable install, these would be populated.
    pip show scvi-tools
    ```

04. **\[Advanced users\]** Install your local copy into a virtualenv with Poetry. Our preferred local installation method consists of using `pyenv-virtualenv` to create a virtualenv, and using `poetry` to create an editable local installation. If using this approach, please be sure to install `poetry` the [recommended](https://python-poetry.org/docs/#installation) way. Once `poetry` is installed:

    ```
    pyenv virtualenv scvi-tools-dev
    pyenv activate scvi-tools-dev
    cd scvi-tools
    poetry install --extras "dev docs tutorials"
    ```

    To confirm that scvi-tools was successfully installed, proceed in the same way as above. This time, `pip show scvi-tools` should show all other metadata as well (Version, Summary, etc.).

05. **\[Optional\]** Install a version of PyTorch that supports your GPU. This will be the case even if you use Poetry.

06. Create an ipykernel so you can use your environment with a Jupyter notebook. This will make this developement environment available through Jupyter notebook/lab. Inside your virtualenv:

    ```
    python -m ipykernel install --user --name=scvi-tools-dev
    ```

07. Install pre-commit, which will enforce the scvi-tools code style (black, flake8) on each of your commits:

    ```
    $ pre-commit install
    ```

08. Create a branch for local development:

    ```
    $ git checkout -b {your-branch-name}
    ```

    Now you can make your changes locally.

09. Add tests to the `/tests` directory. These files start with `test_` and contain functions that start similarly with `test_`.

10. When you're done making changes, run the tests using pytest:

    ```
    $ pytest tests/models/test_my_new_feature.py
    $ pytest tests/models/test_my_new_feature.py::test_particular_function_in_file
    ```

11. Commit your changes and push your branch to GitHub:

    ```
    $ git add <file> ...
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature
    ```

12. Submit a pull request through the GitHub website.

## Coding Standards

1. Don't duplicate code. Certainly no blocks longer than a couple of lines. It's almost always better to refactor than to duplicate blocks of code.
2. Almost all code should at least be run by a unit test. No pull request should decrease unit test coverage by much.
3. Document each new method and each new class with a docstring.
4. Don't commit commented-out code. Just delete it or store it somewhere outside of the repo. You probably aren't going to need it. At worse, it's stored in previous commits, from before it was commented out.
5. A pull request (PR) will typically close at least one Github issue. For these pull requests, write the issue it closes in the description, e.g. `closes #210`. The issue will be automatically closed when the PR is merged.
6. Don't commit data to the repository, except perhaps a few small (\< 50 KB) files of test data.
7. Respect the scvi-tools code style, the easiest way is to install pre-commit as described above.
8. Your PR will be evaluated by automated software for following coding standards, so it's easier to start with good practices.

## Documenting Code

This section is under construction, but we use the same docstring style as Scanpy. See their [tutorial](https://scanpy.readthedocs.io/en/stable/dev/documentation.html#building-the-docs) for more info.

For the actual documentation pages (e.g., the user guide), we use markdown files via the [myst-parser](https://myst-parser.readthedocs.io/en/latest/index.html).

## Tips

1. [GitKraken](https://www.gitkraken.com/) can be a useful GUI for using git locally.
2. `git commit -m "my message" --no-verify` allows overriding `pre-commit`.
3. Reach out on [gitter](https://gitter.im/scvi-tools/development) if you need help.

## Pull Request Guidelines

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated.
3. The pull request should work for Python 3.6-3.8. Your PR will be tested
   on these versions with our continuous integration checks.

## Deploying

First, please install Poetry.

A reminder for the maintainers on how to deploy. Make sure all your changes are committed (including a release note entry).

Additionally, make sure to commit a version bump to [pyproject.toml](https://github.com/YosefLab/scvi-tools/blob/master/pyproject.toml) which can be updated by running:

```
$ poetry version preversion # possible: major / minor / patch
```

Then, make sure you've tested your code using pytest by running:

```
$ pytest
```

Subsequently run:

```
$ poetry build
$ poetry publish
```

This will upload `scvi-tools` to PyPi. Also be sure to add a tag corresponding to the new version number on the tutorials repo, as the tagged repo is used for the Colab links.

### Instructions on Uploading to conda

`scvi-tools` is available on conda-forge channel. Typically, a PR will be automatically created once a new PyPI release is made.
This automated PR might need changes if we've changed dependencies. In that case, follow the below steps to upload a new version to conda-forge channel.

Create a fork of the scvi-tools feedstock [repo] on GitHub and follow instructions in the README there.

### Writing a GitHub release

On the GitHub page, draft a release. This is important for ReadTheDocs, which uses the last release tag from GitHub as the stable version.

## Backporting

This is a guide for the maintainers on how we backport patches.

The mainstream development branch is the master branch. We snap releases off of release branches created off of master.

We use the MeeseeksDev GitHub bot for automatic backporting. The way it works, in a nutshell, is that the bot listens to certain web events - for example commits containing “@meeseeksdev backport to \[BRANCHNAME\]” on a PR - and automatically opens a PR to that repo/branch. (Note: They open the PR sourced from a fork of the repo under the [MeeseeksMachine](https://github.com/meeseeksmachine) organization, into the repo/branch of interest. That’s why under MeeseeksMachine you see a collection of repo's that are forks of the repo's that use MeeseeksDev).

For each release, we create a branch \[MAJOR\].\[MINOR\].x where MAJOR and MINOR are the Major and Minor version numbers for that release, respectively, and x is the literal “x”. Every time a bug fix PR is merged into master, we evaluate whether it is worthy of being backported into the current release and if so use MeeseeksDev to do it for us if it can. How? Simply leave a comment on the PR that was merged into master that says: “@meeseeksdev backport to \[MAJOR\].\[MINOR\].x” (for example “@meeseeksdev backport to 0.14.x” if we are on a release from the 0.14 series.
Note: Auto backporting can also be triggered if you associate the PR with a Milestone or Label the description of which contains “on-merge: backport to \[BRANCHNAME\]”.

```{highlight} none
```

```
feature foo <- head of branch master, main development branch
|
bug fix
|
feature bar <- head of branch 0.14.x, release branch for the 0.14.x release series, also tagged as v0.14.0 (release)
\
  my hotfix <- backported from master
  |
  my other hotfix <- backported from master, also tagged as v0.14.1 (release)
|
feature baz
|
my hotfix
|
another bug fix
|
my other hotfix
```

```{highlight} shell
```

### Manually backporting a patch

If MeeseeksDev cannot automatically cherry-pick the PR (e.g. due to conflicts requiring manual resolution), it will let us know. In that case we need to cherry-pick the commit ourselves. [Here](https://github.com/search?q=label%3A%22Still+Needs+Manual+Backport%22+is%3Aopen&state=closed&type=Issues) are examples of such cases, and [here](https://github.com/pandas-dev/pandas/wiki/Backporting) is one resource explaining how to do it, but there are probably a lot more on the web.

[repo]: https://github.com/conda-forge/scvi-tools-feedstock
