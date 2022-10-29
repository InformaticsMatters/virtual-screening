# Job: moldb-enumerate-mols

This describes how to run the `enumerate-candidates` job from the `virtual screening` category in the `im-virtual-screening` collection.

## What the job does
This job enumerates microstates, tautomers and undefined chiral centres of a molecules to generate a set of variants of the molecules suitable for virtual screening.

The microstates are enumerated using [Dimorphite-DL](https://durrantlab.pitt.edu/dimorphite-dl/).
Tautomers are enumerated using [RDKit's tautomer generator](http://rdkit.org/docs/source/rdkit.Chem.MolStandardize.rdMolStandardize.html?highlight=tautomerenumerator#rdkit.Chem.MolStandardize.rdMolStandardize.TautomerEnumerator).
Chiral centres are enumerated using Python code using RDKit.

Molecules needing enumeration are extracted from the MolDB database, the different forms enumerated and 
those forms then loaded into the database (the `enumeration` table).

Typically about 10 molecules are generated for each input, but this number varies considerably.
Which molecules to process is defined by a set of molecular property filters, allowing you to, 
for instance, restrict processing to molecules with a particular heavy atom count range and no
more than a certain number of rotatable bonds.

A single low energy conformer of the molecule is created. These conformers are usually suitable as input
to docking runs (e.g. `run-rdock` and `run-smina` jobs).

**NOTE**: Before running this job ensure you have run the `moldb-calc-props` to completion to generate a full
set of molecular properties otherwise you will not get a complete set of molecules.

## Implementation details

This job is implemented as a [Nextflow](https://www.nextflow.io/) workflow.

* Nextflow workflow: [enumerate_mols.nf](/moldb/enumerate_mols.nf)
* Python module for extracting molecules: [filter.py](/moldb/filter.py)
* Python module for enumeration: [enumerate.py](/enumerate.py)
* Python module for loading molecules: [db_load.py](/moldb/db_load.py)
* Job definition: `jobs.moldb-enumerate-mols` in [moldb.yaml](/data-manager/moldb.yaml)

## How to run the job

### Options

* **Max number of molecules to process**: the maximum number of molecules to extract.
* **Chunk size for splitting**: number of molecules per 'chunk' when calculating properties in parallel
* **Minimum/Maximum XYZ**: Min or Max values for the various molecular property filters

## Related topics

* [moldb-gen-confs job](moldb-gen-confs.md)