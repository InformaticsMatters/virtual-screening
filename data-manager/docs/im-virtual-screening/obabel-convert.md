This describes how to run the `obabel-convert` job from the `virtual screening` category in the `im-virtual-screening` collection.

## What the job does

This job can be used to convert molecule formats, but is primarily used to prepare a protein for docking.
e.g. protonating the protein and/or converting from PDB to MOL2 format (or vice-versa).

Note that the file formats are typically auto-detected from the file extension, but this can be specified using command 
line options (specified as the *Additional options* option).

Also note that the [pdb2pqr](https://discourse.squonk.it/t/job-pdb2pqr/76) job provides a  more sophisticated approach to protonation of a protein.

## Implementation details

OpenBabel options: https://open-babel.readthedocs.io/en/latest/Command-line_tools/babel.html
Job definition: `jobs.obabel-convert` in [/data-manager/virtual-screening.yaml]()

## How to run the job

### Inputs

**Input  file**: Input file e.g. protein.pdb

### Options
**Filename for output**: ouput filename with extension e.g. protein.mol2
**Additional options**: additional commandline options e.g. -p 7.0