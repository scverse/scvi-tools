# Contributing code

This guide is intended for developers who wish to contribute to our codebase. Here's how to set
up a local development environment:

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

4. Install the development dependencies and the package in editable mode into a virtual
    environment with Python 3.10 - 3.13:

    ```bash
    pip install -e ".[dev]"
    # or with uv
    uv pip install -e ".[dev]"
    ```

    Don't know how to set up a virtual environment? Check out our [installation] guide!

5. (Optional) confirm that the installation was successful:

    ```bash
    pip show scvi-tools
    ```

6. (Optional) Set up pre-commit git hooks:

    ```bash
    pre-commit install
    ```

    This will run pre-commit checks before each commit, including code formatting and linting.
    Alternatively, you can run the checks manually with:

    ```bash
    pre-commit  # check modified files
    # or
    pre-commit run --all  # check all files
    ```

## Scoping changes

Before you start working on a new feature or bug fix, we recommend opening an [issue] (if one does
not already exist) to discuss the proposed changes. This will help ensure that your changes are
aligned with the project's goals and that you are not duplicating work.

We don't guarantee that all changes will be accepted, but we will do our best to provide feedback
and guidance on how to improve your contributions.

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

3. (Optional) If your changes add a new feature or address an existing bug, we require that you add
    tests to cover the new code, which should be added under the `tests` directory.

    To run the tests, use:

    ```bash
    pytest  # run all tests
    # or
    pytest tests/test_my_change.py  # run tests in a specific file
    # or
    pytest tests/test_my_change.py::test_my_change  # run a specific test
    ```

4. (Optional) If your changes add a new function or class to the public API, please include
    docstrings that describe the purpose, usage, and parameters of the new code, and update the
    API reference (`docs/api`) accordingly.

5. Include a description of your changes in the release notes (`CHANGELOG.md`). If you are unsure
    where to place your changes, please ask in the pull request.

6. Push your changes to your fork:

    ```bash
    git push origin my-change
    ```

7. Open a pull request on the main repository. Make sure to include a detailed description of your
    changes in the body and reference any related issues.

## Standards and conventions

- We use [ruff] for formatting and linting Python files, which closely mirrors the [Black] code
    style.
- We use [Prettier] for formatting YAML files.
- We use [Mdformat] and [markdownlint] for formatting and linting Markdown files.
- We use the [numpydoc] style for docstrings. All public functions and classes must have a
    docstring that describes their purpose, usage, and parameters.
- Although not all parts of our codebase are type-annotated yet, we recommend that all new code be
    annotated with type hints according to the [PEP 484] and [PEP 526] guidelines.
- We generally don't commit data files, except if they are small and necessary for testing. This
    will be assessed on a case-by-case basis.
- Starting from version 1.2, we format commits into the main branch according to
    [Conventional Commits]. Don't worry if you're not familiar with this convention as a maintainer
    will format any commits before merging into the main branch. For more details, see our
    [maintenance guide].
- Starting with version 0.20, we follow the [Keep a Changelog] convention for `CHANGELOG.md`
- We version our releases according to [Semantic Versioning].

[repository]: https://github.com/scverse/scvi-tools
[issue]: https://github.com/scverse/scvi-tools/issues
[ruff]: https://github.com/astral-sh/ruff
[black]: https://github.com/psf/black
[numpydoc]: https://numpydoc.readthedocs.io/en/latest/format.html
[PEP 484]: https://peps.python.org/pep-0484/
[PEP 526]: https://peps.python.org/pep-0526/
[Conventional Commits]: https://www.conventionalcommits.org/
[Prettier]: https://prettier.io/
[Mdformat]: https://github.com/executablebooks/mdformat
[markdownlint]: https://github.com/igorshubovych/markdownlint-cli
[Keep a Changelog]: https://keepachangelog.com/
[Semantic Versioning]: https://semver.org/
[installation]: ../installation.md
[maintenance guide]: ./maintenance.md
