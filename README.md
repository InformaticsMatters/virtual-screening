# Virtual screening tools

![build](https://github.com/InformaticsMatters/virtual-screening/workflows/build/badge.svg)
![publish-stable](https://github.com/InformaticsMatters/virtual-screening/workflows/publish-stable/badge.svg)
![publish-latest](https://github.com/InformaticsMatters/virtual-screening/workflows/publish-latest/badge.svg)

![GitHub tag (latest SemVer)](https://img.shields.io/github/v/tag/informaticsmatters/virtual-screening)

This repo contains a set of tools for running virtual screening operations.
The tools can be used *standalone*, but many are also packaged up as Squonk2 *jobs*.

More information can be found in the following docs:

* [User guide](USER_GUIDE.md)
* [Developer guide](https://informaticsmatters.gitlab.io/squonk2-data-manager/1-1/creating-new-jobs.html)

These tools are available under the [Apache2.0 license](LICENSE).

## Building
The image builds are accomplished using GitHib workflows in this repository.

The `rdkit-base` image on which a number of the other images depend on is
built manually and pushed to docker hub when the RDKit release needs to be updated.
At the moment we push the image `informaticsmatters/rdkit-base:latest`, and this is
controlled by the dockerfile `Dockerfile-rdkit-base`.

The other images are built automatically by the GitHub workflows, where: -

- Every change to the `main` branch results in a series of `:stable` images
- Every change to the `staging` branch results in a series of `:latest` images
