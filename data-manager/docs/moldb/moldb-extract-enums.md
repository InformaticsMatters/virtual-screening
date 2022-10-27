# Job: moldb-extract-enums

This describes how to run the `moldb-extract-enums` job from the `ligand preparation` category in the `moldb` collection.

## What the job does

The extracts enumerated forms of molecules from the MolDB database and saves them to a file 
suitable for virtual screening. Each molecule has a single low energy conformer. Files can 
be saved in CXSMILES or SDF formats.

**NOTE**: Before running this job ensure you have run the `moldb-calc-props` and 
`moldb-enumerate-mols` jobs to completion to generate a full
set of molecules otherwise you will not get a complete set.

## Implementation details

* Job implementation: [filter.nf](/moldb/filter.py)
* Job definition: `jobs.moldb-extract-enums` in [moldb.yaml](/data-manager/moldb.yaml)

## How to run the job

### Options

* **Output file name**: output file name (.sdf or .cxsmi)
* **Max number of molecules to extract**: the maximum number of molecules to extract.
* **Enumerated types**: at least one of Base, Microstate, Tautomer, Stereoisomer
* **Minimum/Maximum XYZ**: Min or Max values for the various molecular property filters

## Related topics

* [moldb-enumerate-mols job](moldb-enumerate-mols.md)