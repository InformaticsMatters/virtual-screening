# Job: moldb-extract-molecules

This describes how to run the `moldb-extract-molecules` job from the `ligand preparation` category in the `moldb` collection.

## What the job does

The extracts molecules from the MolDB database and saves them to a file 
suitable for further work.

**NOTE**: Before running this job ensure you have run the `moldb-calc-props` job to completion to generate a full
set of molecules otherwise you will not get a complete set.

## Implementation details

* Job implementation: [moldb-extract-molecules.py](/moldb/moldb-extract-molecules.py)
* Job definition: `jobs.moldb-extract-molecules` in [moldb.yaml](/data-manager/moldb.yaml)

## How to run the job

For general information on using this job and interpreting the output look at the *Extracting molecules*
section [here](https://discourse.squonk.it/t/about-moldb/138).

### Options

* **Output file name**: output file name (.sdf or .cxsmi)
* **Max number of molecules to extract**: the maximum number of molecules to extract.

## Related topics

* [About MOlDB](https://discourse.squonk.it/t/about-moldb/138)
* [moldb-extract-enums job](moldb-extract-enums.md)
* [moldb-extract-confs job](moldb-extract-confs.md)