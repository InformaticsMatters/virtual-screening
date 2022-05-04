This describes how to run the *minimize-mmff94* job from the *comp chem* category in the *rdkit* collection.

## What the job does

This job reads a set of molecules form a SD file and performs an energy minimisation using RDKit and the MMFF94 force field. The minimized molecules are aligned to the original molecule.

## Implementation details

* Job implementation: [/minimize.py]()
* Job definition: `jobs.minimize-mmff94` in [/data-manager/rdkit.yaml]()

## How to run the job

### Inputs

* **Molecules to minimize**: Input molecules (*.sdf).

### Options
* **Output file name**: The name for the output file  (*.sdf).
* **Number of cycles**: The number of minimisation cycles
* **Remove Hs from outputs**: Remove hydrogen atoms from the outputs
* **Molecule to write**: which molecules to write: *minimized* = the minimized molecule, *original* = the original molecule, *merged* = both merged into a single molecule

## Outputs

A SD-file is output containing the molecules specified by the *Molecule to write* option, all the original properties plus the following properties are written:
* **RMSD**: the RMSD for the alignment of the minimised molecule to the original
* **CONVERGED**: did the minimisation converge
* **ENERGY_MIN**: the energy of the minimized molecule
* **ENERGY_DELTA**: the energy difference between the original and minimized molecules