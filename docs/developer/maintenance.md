# Maintenance guide

This guide includes various sections that are applicable to maintainers of the project.

## Updating tutorials

The [tutorials] repository is included as a git submodule in the main repository under
`docs/tutorials/notebooks`. This allows us to lower the size of the main repository as the
notebooks can be quite large. This also means that changes committed to the tutorials repository
are not immediately reflected in the main repository. In order to update the tutorials in the
main repository (and subsequently the documentation), run the [update tutorials workflow], which
will create a PR updating the head reference of the submodule to the latest commit in the
tutorials repository.

## Releases

We follow [Semantic Versioning] for naming releases. In short, this means that each release is
tagged with a version number in the format `MAJOR.MINOR.PATCH`.

### Release branches

We create release branches for each increment in the major or minor version, formatted as
`MAJOR.MINOR.x`, _e.g._, `1.0.x`. This means that all patch releases (_e.g._ `1.0.1`, `1.0.2`,
etc.) will be made from the `1.0.x` branch.

This allows us to merge new features and/or breaking changes into the `main` branch without
affecting the release branches. We can then backport bug fixes and other changes to the release
branches as needed.

### Backporting

The development branch is the `main` branch, and we create release branches off it. Since pull
requests only come into the `main` branch, backporting is used to bring new changes into release
branches.

See the diagram below for an example of how this works:

```text
feature foo <- head of branch main, main development branch
|
bug fix
|
feature bar <- head of branch 0.14.x, release branch for the 0.14.x release series
\
  my hotfix <- backported from main
  |
  my other hotfix <- backported from main, tagged as v0.14.1 (release)
|
feature baz
|
my hotfix
|
another bug fix
|
my other hotfix
```

#### Automatic backporting with MeeseeksDev

MeeseeksDev is a GitHub bot that handles automatic backporting. It is configured to backport a pull
request if one of the following conditions is met:

1. The pull request contains a label with the title and body set to
    `on-merge: backport to MAJOR.MINOR.x`, where `MAJOR.MINOR.x` is an existing release branch.
    The PR must be labeled _before_ it is merged in order to trigger the backport.

2. The pull request contains a comment with the following text:

    ```text
    @meeseeksdev backport to MAJOR.MINOR.x
    ```

    where `MAJOR.MINOR.x` is an existing release branch. This comment must be added _after_ the PR
    is merged.

Once the bot is triggered, it will create a new pull request into the specified release branch with
the changes from the original pull request. This will be done from a fork of the repository, so it
must be approved and merged by a maintainer.

Note that this process can be repeated for multiple release branches, if, for example, a patch
change is needed for multiple versions.

#### Manual backporting

It is possible that an automatic backport fails due to merge conflicts or other issues. In this
case, the bot will attach a `Still Needs Manual Backport` label and will include a comment in the
original pull request with instructions on how to manually backport the changes. You can see an
example of this in [#2584].

Follow the instructions and remove the `Still Needs Manual Backport` label once the backport is
complete!

### Making a release

For convenience, we have a [release checklist] that outlines the steps to take when making a new
release. It is recommended to create an issue from this template to track the progress of the
release (see [#2327] for an example). This section provides an overview of the steps involved.

#### (Optional) creating a release branch

As mentioned above, if the release increments the major or minor version, a new release branch
should be created from `main`. This branch should be named according to the new version, with the
literal `x` in place of the patch version, _e.g._, `1.0.x`. Our GitHub rulesets will automatically
protect this branch from direct pushes and will require pull requests for changes.

#### Bumping the version

The next step is to bump the version in the `pyproject.toml` file. This should be done via a PR
into `main` and then an appropriate backport into the release branch.

#### (Optional) Re-run the tutorials

It is recommended to re-run all the [tutorials] for major and minor releases, and affected
tutorials for patch releases. This ensures that the tutorials are up-to-date with the latest
changes.

First, trigger a [Docker image build] targeting the release branch. This will build and upload
an image to the registry with the new changes.

Then, [run the tutorials] using the new image. This will create individual PRs for each tutorial
that has changed, which must be reviewed and merged. For convenience, there is a
[tutorial checklist] for tracking the progress of the tutorials (see [#210] for an example).

It is possible that there are new bugs or issues in the tutorials revealed when running these
automated updates. Any new bugs should be addressed, and the tutorials re-run until they pass
successfully.

#### Publish a release off the tutorials repository

Once all relevant tutorials have been updated and merged, create a new release on the tutorials
repository targeting `main`. This release should be named according to the new version, _e.g._,
`1.0.0`.

#### Updating the tutorials

Then, update the tutorials submodule in the main repository using the [update tutorials workflow].

#### Run the release workflow

Once the version is bumped in `pyproject.toml` and the tutorials have been updated, it is time to
run the [release workflow]. To run the workflow, you will need to provide two arguments: the
semantic versioning tag (_e.g._ 1.1.5) and the release branch name (_e.g._, `1.0.x`). The workflow
will fail if the version tag does not match the one in `pyproject.toml`.

The workflow will build the package, tag the version with Git, publish the release on GitHub with
a default body, and upload the package to PyPI.

At this point, check that the version updates correctly on [PyPI]. If necessary, follow the
instructions in the next section. Additionally, check that [Read the Docs] builds correctly and
tags with the `stable` version.

:::{important}
Once a particular version tag is published to PyPI, it **cannot** be republished, even if the
release is yanked or deleted. This is a constraint of PyPI, so it is important to ensure that the
version is correct before publishing.

If a mistake is made, the version should be incremented and a new release created and published.
Yanking should only be used in exceptional circumstances, but note that this will not remove the
version from the PyPI index.
:::

#### (Optional) Manual release

If the release workflow fails for any reason, it is possible to manually build and upload the
package to PyPI. This can be done by running the following commands:

```bash
hatch build
hatch publish
```

#### Updating the conda-forge feedstock

Our package is also available on the conda-forge channel. The [feedstock repository] contains
the recipe for building the package.

Typically a PR into the [feedstock] will be automatically created a couple of days after the
PyPI release (see [#32] for an example) and will be automatically merged.

If there are any issues with the recipe (_e.g._ due to dependency changes), it may need to be
updated manually. This can be done by following the instructions in the PR (see [#38] for an
example).

#### Update Docker images

Finally, build new Docker images with the correct branch tag using the [Docker image build].

## Continuous integration

We use a combination of GitHub-hosted and self-hosted runners for our continuous integration
(CI) workflows. The GitHub-hosted runners are used for most of our workflows, while the self-hosted
runners are used for workflows that require a GPU.

Each `.yaml` file under `.github/workflows/` corresponds to a specific workflow, where the
tags under `on:` specify the events that trigger the workflow. Additionally, our branch rulesets
(`Settings > Rules > Rulesets`) specify which workflows are required to pass before merging into
a branch.

- [`build`](https://github.com/scverse/scvi-tools/blob/main/.github/workflows/build.yml): Runs
  on pushes and pull requests into `main` and release branches. This workflow attempts to build the
  package.
- [`release`](https://github.com/scverse/scvi-tools/blob/main/.github/workflows/release.yml):
  Triggered by a new release on GitHub, this workflow builds the package and uploads it to PyPI
- [`test`](https://github.com/scverse/scvi-tools/blob/main/.github/workflows/test_linux.yml): Runs
  the test suite on Ubuntu using a GitHub-hosted runner. This workflow is triggered on pushes and
  pull requests into `main` and release branches, as well as based on a cron schedule and manual
  triggers.
- [`test (cuda)`](https://github.com/scverse/scvi-tools/blob/main/.github/workflows/test_linux_cuda.yml):
  Same as the `test` workflow, but runs on a self-hosted runner with a GPU. This workflow is only
  triggered on pull requests with the `cuda tests` or `all tests` labels, which can only be added
  by maintainers.
- [`test (private)`](https://github.com/scverse/scvi-tools/blob/main/.github/workflows/test_linux_private.yml):
  Runs on a GitHub-hosted runner for tests that require authentication (_e.g._, downloading and
  uploading to AWS S3). This workflow is triggered on pull requests with the `private tests` or
  `all tests` labels, which can only be added by maintainers.
- [`test (resolution)`](https://github.com/scverse/scvi-tools/blob/main/.github/workflows/test_linux_resolution.yml):
  Optionally triggered on pull requests with the `resolution tests` or `all tests` labels, this
  workflow runs for various package resolution settings (_e.g._, pre-releases, lowest suppported
  dependencies). Useful for catching dependency issues early if packages have pre-releases.
- [`test (macos)`](https://github.com/scverse/scvi-tools/blob/main/.github/workflows/test_macos.yml):
  Runs the test suite on a GitHub-hosted macOS runner. Useful for catching macOS-specific issues.
- [`test (windows)`](https://github.com/scverse/scvi-tools/blob/main/.github/workflows/test_windows.yml):
  Runs the test suite on a GitHub-hosted Windows runner. Useful for catching Windows-specific
  issues.

## Documentation

Documentation is built and hosted on [Read the Docs] (RTD), and the configuration can be found in
`.readthedocs.yaml`. While the documentation will be automatically built for pull requests into
`main` and release branches, it can be useful to build the documentation locally since RTD only
allows one build at a time.

Make sure to install the documentation dependencies first:

```bash
pip install -e ".[docsbuild]"
```

Run the following command to build the documentation locally from
the root of the repository:

```bash
python -m sphinx -b html docs docs/_build
```

The documentation can be viewed by opening `docs/_build/index.html` in a web browser.

## Conventional commits

Starting with version 1.2, we use [conventional commits] for commits into the main branch. We
support the following types (loosely based on the [Angular convention]):

- `build`: Changes that primarily affect depdenencies
- `ci`: Changes to our CI configuration files and scripts, typically in the `.github/` directory
- `docs`: Documentation changes, typically in the `docs/` directory
- `feat`: A new feature
- `fix`: A bug fix
- `misc`: Changes that don't fit into any other category
- `perf`: A code change that improves performance
- `refactor`: A code change that neither fixes a bug nor adds a feature
- `release`: Changes necessary for a release
- `style`: Changes that do not affect the meaning of the code (_e.g._ whitespace, formatting
  changes in pre-commit configurations)
- `test`: Changes that primarily add or modify tests

We support the following scopes, based on our top-level package structure:

- `autotune`: Changes primarly affecting `scvi.autotune`
- `criticism`: Changes primarily affecting `scvi.criticism`
- `data`: Changes primarily affecting `scvi.data`
- `dataloaders`: Changes primarily affecting `scvi.dataloaders`
- `distributions`: Changes primarily affecting `scvi.distributions`
- `external`: Changes primarily affecting `scvi.external`
- `hub`: Changes primarily affecting `scvi.hub`
- `model`: Changes primarily affecting `scvi.model`
- `module`: Changes primarily affecting `scvi.module`
- `nn`: Changes primarily affecting `scvi.nn`
- `train`: Changes primarily affecting `scvi.train`
- `utils`: Changes primarily affecting `scvi.utils`

We use the `BREAKING CHANGE` footer to indicate that a commit introduces a breaking change.

[#32]: https://github.com/conda-forge/scvi-tools-feedstock/pull/32
[#38]: https://github.com/conda-forge/scvi-tools-feedstock/pull/38
[#210]: https://github.com/scverse/scvi-tutorials/issues/210
[#2327]: https://github.com/scverse/scvi-tools/issues/2327
[#2584]: https://github.com/scverse/scvi-tools/pull/2584
[Semantic Versioning]: https://semver.org/
[release checklist]: https://github.com/scverse/scvi-tools/blob/main/.github/ISSUE_TEMPLATE/release_checklist.md
[tutorials]: https://github.com/scverse/scvi-tutorials
[Docker image build]: https://github.com/scverse/scvi-tools/actions/workflows/build_image_latest.yaml
[run the tutorials]: https://github.com/scverse/scvi-tutorials/actions/workflows/run_notebook_individual.yaml
[tutorial checklist]: https://github.com/scverse/scvi-tutorials/blob/main/.github/ISSUE_TEMPLATE/release_checklist.md
[release workflow]: https://github.com/scverse/scvi-tools/actions/workflows/release.yml
[PyPI]: https://pypi.org/project/scvi-tools/
[feedstock repository]: https://github.com/conda-forge/scvi-tools-feedstock
[Read the Docs]: https://readthedocs.org/projects/scvi/
[conventional commits]: https://www.conventionalcommits.org/
[Angular convention]: https://github.com/angular/angular/blob/22b96b9/CONTRIBUTING.md#-commit-message-guidelines
[update tutorials workflow]: https://github.com/scverse/scvi-tools/actions/workflows/tutorials.yaml
