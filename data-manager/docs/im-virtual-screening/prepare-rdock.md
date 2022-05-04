# Job: prepare-rdock

This describes how to run the `prepare-rdock` job from the `virtual screening` category in the `im-virtual-screening` collection.

## What the job does

This job can be used to generate the rDock configuration file (.prm file) and the cavity definition file (.as file) that
are needed as inputs to the [run-rdock](run-rdock.md) job.
Note that this job only generates a basic *default* version of these files and you might want to prepare your own
versions of these files if you have more specific needs.

## Implementation details

* Python module: [prepare_rdock.py](/prepare_rdock.py)
* Job definition: `jobs.prepare-rdock` in [im-virtual-screening.yaml](/data-manager/im-virtual-screening.yaml)

## How to run the job

### Inputs
* **MOL2 file for receptor**: the prepared receptor in MOL2 format.
* **Molfile file for ligand**: a candidate ligand that is used to defined the cavity in Molfile format.

### Options
* **Base name for output files**: that base name for the output files