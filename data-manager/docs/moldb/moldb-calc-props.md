# Job: moldb-calc-props

This describes how to run the `moldb-calc-props` job from the `ligand preparation` category in the `moldb` collection.

## What the job does

The job reads a set of molecules from the MolDB database that have not yet had molecular
properties calculated, calculates the molecular properties and then loads those properties
into teh databaase.

It is typically run immediately after the `moldb-load-library` job.

You can choose how many molecules to extract. For large libraries you will need to run this
job several times to get all molecules calculated.

These properties are calculated:

- Heavy atom count
- Number of rotatable bonds
- Number of rings
- Number of aromatic rings
- Number of chiral centres
- Number of undefined chiral centres
- Number of SP3 hybridised carbon atoms
- Crippen LogP
- Topological polar surface area

## Implementation details

This job is implemented as a [Nextflow](https://www.nextflow.io/) workflow.

* Job implementation: [calc_molprops.nf](/moldb/calc_molprops.nf)
* Job definition: `jobs.moldb-load-library` in [moldb.yaml](/data-manager/moldb.yaml)

## How to run the job

For general information on using this job and interpreting the output look at the *Calculating molecular properties*
section [here](https://discourse.squonk.it/t/about-moldb/138).

### Options

* **Max number of molecules to process**: the maximum number of molecules to extract.
* **Chunk size for splitting**: number of molecules per 'chunk' when calculating properties in parallel

## Related topics

* [About MOlDB](https://discourse.squonk.it/t/about-moldb/138)
* [moldb-load-library job](moldb-load-library.md)
* [rdk-molprops job](../rdkit/rdk-molprops.md)