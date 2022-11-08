# Job: moldb-extract-confs

This describes how to run the `moldb-extract-confs` job from the `ligand preparation` category in the `moldb` collection.

## What the job does

The extracts low energy 3D conformers from the MolDB database and saves them to a file 
suitable for virtual screening. Files can be saved in CXSMILES or SDF formats.

**NOTE**: Before running this job ensure you have run the
`moldb-gen-confs` jobs to completion to generate a full
set of molecules otherwise you will not get a complete set.

## Implementation details

* Job implementation: [filter.nf](/moldb/filter.py)
* Job definition: `jobs.moldb-extract-confs` in [moldb.yaml](/data-manager/moldb.yaml)

## How to run the job

For general information on using this job and interpreting the output look at the *Extracting 3D conformers*
section [here](https://discourse.squonk.it/t/about-moldb/138).

### Options

* **Output file name**: output file name (.sdf
* **Max number of molecules to extract**: the maximum number of molecules to extract.
* **Enumerated types**: at least one of Base, Microstate, Tautomer, Stereoisomer
* **Minimum/Maximum XYZ**: Min or Max values for the various molecular property filters

## Related topics

* [About MOlDB](https://discourse.squonk.it/t/about-moldb/138)
* [moldb-enumerate-mols job](moldb-enumerate-mols.md)