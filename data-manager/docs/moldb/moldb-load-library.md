# Job: moldb-load-library

This describes how to run the `moldb-load-library` job from the `ligand preparation` category in the `moldb` collection.

## What the job does

The job takes molecules in smiles format as input and:
* standardizes the molecules
* loads them into the MolDB database

Following this you should calculate molecular properties using the `moldb-calc-props` job.

**NOTE**: After running this job to load a new library ensure you run the `moldb-calc-props` to 
completion to generate a full set of molecular properties otherwise you will not get a complete
set of molecules when you run the `moldb-enumerate-mols` and `moldb-gen-confs` jobs.

## Implementation details

This job is implemented as a [Nextflow](https://www.nextflow.io/) workflow.

* Job implementation: [load_library.nf](/moldb/load_library.nf)
* Job definition: `jobs.moldb-load-library` in [moldb.yaml](/data-manager/moldb.yaml)

## How to run the job

### Inputs
* **Molecules to load**: a file from a vendor in SMILES format. An example file with a small subset of ChemSpace can be found [here](https://github.com/InformaticsMatters/virtual-screening/blob/main/data/100000.smi).

### Options
* **Library name**: the name of the library e.g. `ChemSpace`
* **Index of ID field**: the zero based index of the column in the input that contains the vendor's code for that molecule e.g. `1`
* **Delimiter**: separator for the file. e.g. `tab`, `comma`, `space`, `pipe`
* **Skip header line**: file has a header line that should be skipped
* **Chunk size for splitting**: number of molecules per 'chunk' when standardizing in parallel


## Related topics

* [moldb-calc-props job](moldb-calc-props.md)