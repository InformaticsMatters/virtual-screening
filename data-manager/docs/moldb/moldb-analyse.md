# Job: moldb-extract-confs

This describes how to run the `moldb-analyse` job from the `ligand preparation` category in the `moldb` collection.

## What the job does

The analyses the contents of the MolDB database and generates a report telling you the status of the database.
From this report you can work out which molecues need property calculation, enumeration or conformer generation.

## Implementation details

* Job implementation: [filter.nf](/moldb/analyse.py)
* Job definition: `jobs.moldb-analyse` in [moldb.yaml](/data-manager/moldb.yaml)

## How to run the job

For general information on using this job and interpreting the output look at the *Finding out what needs doing*
section [here](https://discourse.squonk.it/t/about-moldb/138).

### Inputs

* **Filter specification file**: file containing the specification of molecular property filters

### Options

* **Output file name**: output file name (.txt)
* **Skip enumeration analysis**: do not perform analysis of the enumeration table
* **Skip conformer analysis**: do not perform analysis of the conformer table

## Related topics

* [About MOlDB](https://discourse.squonk.it/t/about-moldb/138)
* [moldb-enumerate-mols job](moldb-enumerate-mols.md)
* [moldb-gen-confs job](moldb-gen-confs.md)