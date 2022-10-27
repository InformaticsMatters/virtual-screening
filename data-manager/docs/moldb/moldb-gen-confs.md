# Job: moldb-gen-confs

This describes how to run the `moldb-gen-confs` job from the `ligand preparation` category in the `moldb` collection.

## What the job does

This job generates a multiple energy minimised 3D conformers of molecules from the MolDB database.
It is typically used to generate input for shape similarity screening using USR and it's derivatives.

The conformers are generated using RDKit's ETKDG conformer generator and energy minimised using the MMFF94 force field.
The basics are described [here](http://rdkit.org/docs/GettingStartedInPython.html#working-with-3d-molecules)

The number of conformers generated is based on the rules defined in Ejeber *et al.* J. Chem. Inf. Model. 2012, 52, 1146-1158
[doi:abs/10.1021/ci2004658](https://pubs.acs.org/doi/abs/10.1021/ci2004658).
Following conformer generation the structures are energy minimised using the MMFF94 force field. Finally, conformers with
very low RMS differences are removed (keeping the lowest energy one).  Note that we use a different RMS threshold as the
default, 1.0 as opposed to 0.35 that is used in Ejeber *et al.*.

Enumerated forms needing conformer generation are extracted from the MolDB database, conformers are generated and
those forms then loaded into the database (the `conformer` table).

**NOTE**: Before running this job ensure you have run the `moldb-enumerate-mols` to completion using the
same set of molecular property filters otherwise you will not get a complete set of molecules.

## Implementation details

This job is implemented as a [Nextflow](https://www.nextflow.io/) workflow.

* Nextflow workflow: [generate_confs.nf](/moldb/generate_confs.nf)
* Python module for extracting molecules: [filter.py](/moldb/filter.py)
* Python module for conformer generation: [conformers.py](/moldb/conformers.py)
* Python module for loading molecules: [db_load.py](/moldb/db_load.py)
* Job definition: `jobs.moldb-gen-confs` in [moldb.yaml](/data-manager/moldb.yaml)

## How to run the job

### Options

* **Max number of molecules to process**: the maximum number of molecules to extract.
* **Chunk size for splitting**: number of molecules per 'chunk' when calculating properties in parallel
* **Minimum/Maximum XYZ**: Min or Max values for the various molecular property filters

## Related topics

* [moldb-enumerate-mols job](moldb-enumerate-mols.md)