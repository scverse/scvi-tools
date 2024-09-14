---
name: Release checklist
about: Checklist for developers
title: ""
labels: releases
assignees: ""
---

- [ ] Bump version in `pyproject.toml`
- [ ] If patch release, backport version bump PR into the appropriate branch. Else, create a new
    branch off `main` with the appropriate rules
- [ ] Trigger a Docker image build in [`scvi-tools-docker`](https://github.com/YosefLab/scvi-tools-docker)
    targeting the release branch
- [ ] After image builds and pushes to the registry, run the [tutorials](https://github.com/scverse/scvi-tutorials)
    using the new image
- [ ] Publish a new release on the tutorials repo off `main` after all tutorials changes have been merged
- [ ] Create a new branch off `main` in the main repo and run `git submodule update --remote`, and
    then merge the PR, with an appropriate backport as needed
- [ ] Run the release workflow
- [ ] Check that the version updates correctly on [PyPI](https://pypi.org/project/scvi-tools/)
- [ ] Build new Docker images with the `stable` and semantic versioning tags
- [ ] Check that the [feedstock repo](https://github.com/conda-forge/scvi-tools-feedstock) updates correctly
- [ ] (Optional) Post threads on Discourse and Twitter
