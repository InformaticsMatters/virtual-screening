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

## Source layout
All the Python modules, Nextflow workflows and shell scripts live under `src/`,
including the `moldb`, `im_mordred` and `dmpk` packages and the Nextflow
processes in `src/nf-processes`. The Dockerfiles copy the contents of `src/` into
`/code` in the images, so in a container the files are still at `/code/<file>`
(e.g. `/code/screen.py`, `python -m moldb.filter`) and the Job Definitions refer
to them by those paths. The Nextflow Jobs run their workflow from
`{{ CODE_DIRECTORY ~ '/src' if CODE_DIRECTORY is defined else '/code' }}`
because `jote` sets `CODE_DIRECTORY` to the repository root when it tests them.

To run the tools outside a container, put `src` on the `PYTHONPATH`: -

    $ export PYTHONPATH=$PWD/src
    $ python -m screen --help

## Building
The image builds are accomplished using GitHib workflows in this repository.

The `vs-rdkit-base` image on which a number of the other images depend on is
built manually and pushed to docker hub when the RDKit release needs to be updated.
It is controlled by the dockerfile `Dockerfile-rdkit-base`, and is pushed with
the release tag that the dependent images will use: -

    $ IMAGE_TAG=2.0.0 docker-compose -f docker-compose-manual.yaml build
    $ IMAGE_TAG=2.0.0 docker-compose -f docker-compose-manual.yaml push

Publish it *before* running the release workflow below - `Dockerfile-fns`,
`-moldb`, `-mordred`, `-oddt` and `-prep` are all `FROM` it and will fail to
build without it.

The other images are built automatically by the GitHub workflows, where: -

- The `publish-tag` workflow, run manually, results in a series of images
  carrying the tag you give it (e.g. `:2.0.0`). This is how releases are made:
  the Job Definitions pin these static tags, so **a tag must never be reused**
- Every change to the `main` branch results in a series of `:stable` images
- Every change to the `staging` branch results in a series of `:latest` images

`:stable` and `:latest` are dynamic tags, intended for development only. See
[versioning](https://github.com/InformaticsMatters/squonk2-jobs/blob/main/docs/versioning.md)
in the `squonk2-jobs` umbrella repository.
