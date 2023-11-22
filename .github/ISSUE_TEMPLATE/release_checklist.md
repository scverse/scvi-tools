---
name: Release checklist
about: Checklist for developers
title: ""
labels: releases
assignees: ""
---

-   [ ] Bump version in `pyproject.toml`.
-   [ ] If patch release, backport version bump PR into the appropriate branch.
        Else, create a new branch off main with the appropriate rules.
-   [ ] Create a new GitHub draft release targeting the release branch with the
        same body as the previous release (don't release yet).
-   [ ] Trigger a Docker image build in
        [`scvi-tools-docker`](https://github.com/YosefLab/scvi-tools-docker)
        targeting the release branch.
-   [ ] After image builds and pushes to the registry, run the
        [tutorials](https://github.com/scverse/scvi-tutorials) using the new
        image.
-   [ ] Publish a new release on the tutorials repo off main after all tutorial
        changes have been merged.
-   [ ] Double check everything and then publish the draft release on the main
        repo.
