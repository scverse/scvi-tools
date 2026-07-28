# Contributing code

This guide is intended for developers who wish to contribute to our codebase.
Here's how to set up a local development environment:

:::{tip} The *hatch* project manager

We highly recommend to familiarize yourself with [hatch].
Hatch is a Python project manager that

- manages virtual environments, separately for development, testing and building the documentation.
  Separating the environments is useful to avoid dependency conflicts.
- allows to run tests locally in different environments (e.g. different Python versions).
- allows to run tasks defined in `pyproject.toml`, e.g. to build the documentation.

All of our environments and test suites are declared in the `[tool.hatch]` table of `pyproject.toml`, which makes it the single point of truth for local development and continuous integration alike.
While the project is set up with `hatch` in mind, it is still possible to use different tools to manage dependencies, such as `uv` or `pip`.

:::

## Setting up a development environment

1. Fork our [repository] on GitHub

2. Locally clone your forked repository (replace `your-username` with your GitHub username):

    ```bash
    git clone https://github.com/{your-username}/scvi-tools.git
    # or with SSH
    git clone git@github.com:{your-username}/scvi-tools.git

    cd scvi-tools
    ```

3. Add the main repository as a remote:

    ```bash
    git remote add upstream https://github.com/scverse/scvi-tools.git
    ```

4. Install the development dependencies and the package in editable mode into a virtual environment with Python 3.12 - 3.14.

    :::::{tab-set}
    ::::{tab-item} Hatch
    :sync: hatch

    Hatch resolves the environments in the background the first time you use them:

    ```bash
    hatch test              # defined by [[tool.hatch.envs.hatch-test.matrix]] in pyproject.toml
    hatch run docs:build    # defined by [tool.hatch.envs.docs] in pyproject.toml
    ```

    To get a list of all environments, including the ones for the opt-in test suites, run

    ```bash
    hatch env show
    ```

    As the main development environment, we recommend `hatch-test` with the latest supported Python version, for example `hatch-test.py3.14-stable`.
    Create it and print its path with

    ```bash
    hatch env create hatch-test.py3.14-stable
    hatch env find hatch-test.py3.14-stable
    ```

    and point your IDE at the `python` binary inside that directory.

    ::::

    ::::{tab-item} uv
    :sync: uv

    ```bash
    uv sync --extra optional --group dev --group test
    ```

    The `.venv` directory is typically automatically discovered by IDEs such as VS Code.

    ::::

    ::::{tab-item} Pip
    :sync: pip

    ```bash
    python3 -m venv .venv
    source .venv/bin/activate
    pip install -e ".[optional]" --group dev --group test
    ```

    ::::
    :::::

    Don't know how to set up a virtual environment? Check out our [installation] guide!

5. (Optional) confirm that the installation was successful:

    ```bash
    pip show scvi-tools
    ```

6. (Optional) Set up the git hooks.
    We use [pre-commit]-style hooks to enforce a consistent code style and recommend running them with [prek], a fast drop-in replacement that reads the same `.pre-commit-config.yaml`:

    ```bash
    prek install
    ```

    This will run the checks before each commit, including code formatting and linting.
    Alternatively, you can run the checks manually with:

    ```bash
    prek  # check modified files
    # or
    prek run --all-files  # check all files
    ```

    Alternatively, you can rely on the [pre-commit.ci] service, which is enabled on GitHub and pushes fixes to your pull request automatically.

## Scoping changes

Before you start working on a new feature or bug fix, we recommend opening an [issue] (if one does not already exist) to discuss the proposed changes.
This will help ensure that your changes are aligned with the project's goals and that you are not duplicating work.

We don't guarantee that all changes will be accepted, but we will do our best to provide feedback and guidance on how to improve your contributions.

## Adding code changes

We only accept code changes that are made through pull requests. To contribute, follow these steps:

1. Create a new branch for your changes:

    ```bash
    git checkout -b my-change
    ```

2. Make your changes and commit them:

    ```bash
    git add .
    git commit -m "My change"
    ```

3. (Optional) If your changes add a new feature or address an existing bug, we require that you add tests to cover the new code, which should be added under the `tests` directory.

    To run the tests, use:

    ```bash
    hatch test                                          # run all tests with the latest supported Python version
    hatch test --all                                    # run all tests with every supported Python version
    hatch test tests/test_my_change.py                  # run tests in a specific file
    hatch test tests/test_my_change.py::test_my_change  # run a specific test
    ```

    Most of our test suites are opt-in because they need extra dependencies or hardware.
    Each of them has its own hatch environment and its own workflow in `.github/workflows`, and is triggered on a pull request by applying the matching label (for example `cuda tests`, or `all tests` for every suite).
    To run one locally, invoke its environment directly:

    ```bash
    hatch run test-cuda:run --accelerator cuda --devices auto  # requires a CUDA device
    hatch run hatch-test.py3.13-stable:run --optional          # the slower, optional-dependency tests
    ```

4. (Optional) If your changes add a new function or class to the public API, please include docstrings that describe the purpose, usage, and parameters of the new code, and update the API reference (`docs/api`) accordingly.

5. Include a description of your changes in the release notes (`CHANGELOG.md`).
    If you are unsure where to place your changes, please ask in the pull request.

6. Push your changes to your fork:

    ```bash
    git push origin my-change
    ```

7. Open a pull request on the main repository.
    Make sure to include a detailed description of your changes in the body and reference any related issues.

## Building the documentation

```bash
hatch run docs:build  # build the documentation into docs/_build
hatch run docs:open   # open the result in a browser
hatch run docs:clean  # remove all generated files
```

The build treats warnings as errors, which is also how it runs on Read the Docs, so a locally clean build is a good predictor of a green documentation build.

## Standards and conventions

- We use [ruff] for formatting and linting Python files, which closely mirrors the [Black] code style.
- We use [biome] for formatting JSON and JSONC files, and [pyproject-fmt] for `pyproject.toml`.
- We use [zizmor] to catch security issues in our GitHub Actions workflows.
- We use the [numpydoc] style for docstrings.
    All public functions and classes must have a docstring that describes their purpose, usage, and parameters.
- Although not all parts of our codebase are type-annotated yet, we recommend that all new code be annotated with type hints according to the [PEP 484] and [PEP 526] guidelines.
- We generally don't commit data files, except if they are small and necessary for testing.
    This will be assessed on a case-by-case basis.
- Starting from version 1.2, we format commits into the main branch according to [Conventional Commits].
    Don't worry if you're not familiar with this convention as a maintainer will format any commits before merging into the main branch.
    For more details, see our [maintenance guide].
- Starting with version 0.20, we follow the [Keep a Changelog] convention for `CHANGELOG.md`.
- We version our releases according to [Semantic Versioning].

## Repository layout

This repository follows the [scverse cookiecutter template][template], which is what most other scverse packages are based on.
The template is tracked in `.cruft.json`, and updates are proposed automatically as pull requests by the template repository.
Do not edit `.cruft.json` by hand.
Files that we intentionally keep different from the template are listed under `[tool.cruft]` in `pyproject.toml` so that template updates leave them alone.

[repository]: https://github.com/scverse/scvi-tools
[issue]: https://github.com/scverse/scvi-tools/issues
[hatch]: https://hatch.pypa.io/latest/
[pre-commit]: https://pre-commit.com/
[prek]: https://prek.j178.dev/
[pre-commit.ci]: https://pre-commit.ci/
[ruff]: https://github.com/astral-sh/ruff
[black]: https://github.com/psf/black
[biome]: https://biomejs.dev/
[pyproject-fmt]: https://github.com/tox-dev/pyproject-fmt
[zizmor]: https://docs.zizmor.sh/
[numpydoc]: https://numpydoc.readthedocs.io/en/latest/format.html
[PEP 484]: https://peps.python.org/pep-0484/
[PEP 526]: https://peps.python.org/pep-0526/
[Conventional Commits]: https://www.conventionalcommits.org/
[Keep a Changelog]: https://keepachangelog.com/
[Semantic Versioning]: https://semver.org/
[template]: https://github.com/scverse/cookiecutter-scverse
[installation]: ../installation.md
[maintenance guide]: ./maintenance.md
