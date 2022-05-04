# Job: pdb2pqr

This describes how to run the `pdb2pqr` job from the `virtual screening` category in the `im-virtual-screening` collection.

## What the job does

This job can be used to prepare a protein for docking, using the pdb2pqr tool which protonates the protein according to
the local environment. It can also fix minor problems in the input such as missing heavy atoms.

The input is a PDB file and the fixed protein is written in PDB and PQR formats.

## Implementation details

* pdb2pqr GitHub repo: https://github.com/Electrostatics/pdb2pqr
* Job definition: `jobs.pdb2pqr` in [im-virtual-screening.yaml](/data-manager/im-virtual-screening.yaml])

## How to run the job

### Inputs
* **Input PDB file**: PDB file with the input protein.

### Options
* **pH to protonate at**: the pH to protonate the protein at.
* **Base filename for output**: The base filename for the output files (do not include the extension).