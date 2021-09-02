.. highlight:: shell

============
Maintenance
============

This page is a guide for the maintainers of scvi-tools.


Deploying
---------

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

Writing a GitHub release
~~~~~~~~~~~~~~~~~~~~~~~~

On the GitHub page, draft a release. This is important for ReadTheDocs, which uses the last release tag from GitHub as the stable version.


Backporting
-----------

The mainstream development branch is the master branch. We snap releases off of release branches created off of master.

We use the MeeseeksDev GitHub bot for automatic backporting. The way it works, in a nutshell, is that the bot listens to certain web events - for example commits containing “@meeseeksdev backport to [BRANCHNAME]” on a PR - and automatically opens a PR to that repo/branch. (Note: They open the PR sourced from a fork of the repo under the `MeeseeksMachine <https://github.com/meeseeksmachine>`_ organization, into the repo/branch of interest. That’s why under MeeseeksMachine you see a collection of repo's that are forks of the repo's that use MeeseeksDev).

For each release, we create a branch [MAJOR].[MINOR].x where MAJOR and MINOR are the Major and Minor version numbers for that release, respectively, and x is the literal “x”. Every time a bug fix PR is merged into master, we evaluate whether it is worthy of being backported into the current release and if so use MeeseeksDev to do it for us if it can. How? Simply leave a comment on the PR that was merged into master that says: “@meeseeksdev backport to [MAJOR].[MINOR].x” (for example “@meeseeksdev backport to 0.14.x” if we are on a release from the 0.14 series.
The PR also needs to be associated with a Milestone the description of which contains “on-merge: backport to [BRANCHNAME]”.

.. highlight:: none

::

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

.. highlight:: shell

Manually backporting a patch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If MeeseeksDev cannot automatically cherry-pick the PR (e.g. due to conflicts requiring manual resolution), it will let us know. In that case we need to cherry-pick the commit ourselves. `Here <https://github.com/search?q=label%3A%22Still+Needs+Manual+Backport%22+is%3Aopen&state=closed&type=Issues>`_ are examples of such cases, and `here <https://github.com/pandas-dev/pandas/wiki/Backporting>`_ is one resource explaining how to do it, but there are probably a lot more on the web.